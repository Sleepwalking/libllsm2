/*
  libllsm2 - Low Level Speech Model (version 2)
  ===
  Copyright (c) 2017 Kanru Hua.

  libllsm2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libllsm2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libllsm. If not, see <http://www.gnu.org/licenses/>.
*/

#ifdef USE_PTHREAD
#include <pthread.h>
#endif

#include "llsmrt.h"
#include "buffer.h"
#include "dsputils.h"
#include "external/ciglet/ciglet.h"

// Implemented in layer0.c
FP_TYPE* llsm_synthesize_harmonic_frame_auto(llsm_soptions* options,
  FP_TYPE* ampl, FP_TYPE* phse, int nhar, FP_TYPE f0, int nx);

typedef struct {
  int nout;      // number of output samples left in buffer_out
  int nchannel;  // number of noise channels
  int ntemplate; // size of the noise template
  int ninternal; // capacity of the internal buffers
  llsm_container* conf; // a copy of the configurations
  llsm_soptions opt;    // a copy of the synthesis options

  FP_TYPE fs;    // sampling rate (in Hz)
  FP_TYPE thop;  // hop size (in seconds)
  FP_TYPE cycle; // current position relative to floored sample position
  int curr_nhop; // rounding-adjusted current hop size (in samples)
  int exc_cycle; // current position in the noise template
  int sin_pos;   // read position of the sinusoid buffer
  FP_TYPE* win;  // overlap-add window
  int nfft;      // FFT size for noise filtering
  llsm_nmframe* prev_nm; // the previous noise model frame

  llsm_ringbuffer*  buffer_out;       // output buffer of audio samples

  FP_TYPE** exc_template_comps;       // noise template for each channel
  llsm_ringbuffer** buffer_mod_comps; // noise modulation buffer for each
                                      //   channel
  llsm_ringbuffer*  buffer_exc_mix;   // buffer for the sum of noise excitation
  llsm_ringbuffer*  buffer_noise;     // buffer for the filtered noise
  llsm_ringbuffer*  buffer_sin;       // buffer for the sinusoidal component

  FP_TYPE* buffer_psd; // size: nspec
  FP_TYPE* buffer_fft; // size: nfft * 4
  FP_TYPE* buffer_rawexc; // size: ninternal
  FP_TYPE* buffer_rawmod; // size: ninternal
  FP_TYPE* warp_axis;  // size: nspec

# ifdef USE_PTHREAD
  pthread_mutex_t buffer_out_mtx;
  pthread_cond_t  buffer_out_cv;
# endif
} llsm_rtsynth_buffer_;

static FP_TYPE llsm_get_circular_noise(FP_TYPE* x, int nx, int i) {
  int overlap = 32;
  i = i % nx;
  if(i >= overlap) return x[i];

  FP_TYPE y = x[i];
  FP_TYPE r = (FP_TYPE)i / overlap;
  y *= 1.0 - r;
  y += x[nx - overlap + i] * r;
  y /= sqrt(2 * r * (r - 1) + 1);
  return y;
}

static void llsm_make_exc_template(llsm_rtsynth_buffer_* dst,
  FP_TYPE* chanfreq) {
  FP_TYPE fs = dst -> fs;
  for(int c = 0; c < dst -> nchannel; c ++) {
    FP_TYPE fmin = c == 0 ? 0 : chanfreq[c - 1];
    FP_TYPE fmax = c == dst -> nchannel - 1 ? fs / 2.0 : chanfreq[c];
    if(fmin >= fs / 2.0) break;
    FP_TYPE* x = llsm_generate_bandlimited_noise(dst -> ntemplate,
      fmin / fs * 2.0, fmax / fs * 2.0);
    for(int j = 0; j < dst -> ntemplate; j ++)
      dst -> exc_template_comps[c][j] = llsm_get_circular_noise(x,
        dst -> ntemplate, j);
    free(x);
  }
}

// Prepare for synthesizing the next frame.
static void llsm_update_cycle(llsm_rtsynth_buffer_* dst) {
  for(int c = 0; c < dst -> nchannel; c ++)
    llsm_ringbuffer_appendblank(dst -> buffer_mod_comps[c], dst -> curr_nhop);
  llsm_ringbuffer_appendblank(dst -> buffer_sin, dst -> curr_nhop);
  llsm_ringbuffer_appendblank(dst -> buffer_noise, dst -> curr_nhop);
  dst -> cycle += dst -> thop;
  dst -> curr_nhop = floor(dst -> cycle * dst -> fs);
  dst -> cycle -= (FP_TYPE)dst -> curr_nhop / dst -> fs;
  int nwin = dst -> curr_nhop * 2;
  free(dst -> win);
  dst -> win = hanning_2(nwin);
}

// Assuming the modulation component buffers all have been loaded with the
//   next chunk of samples, load noise from the band-wise noise templates
//   and mix them down to the exc_mix buffer.
static void llsm_run_excitation_buffers(llsm_rtsynth_buffer_* dst, int nx) {
  FP_TYPE* x = dst -> buffer_rawexc; // noise
  FP_TYPE* m = dst -> buffer_rawmod; // modulation
  memset(x, 0, nx * sizeof(FP_TYPE));
  for(int c = 0; c < dst -> nchannel; c ++) {
    llsm_ringbuffer_readchunk(dst -> buffer_mod_comps[c],
      -dst -> curr_nhop - nx, nx, m);
    for(int i = 0; i < nx; i ++)
      x[i] += sqrt(m[i]) * dst -> exc_template_comps[c]
        [(dst -> exc_cycle + i) % dst -> ntemplate];
  }
  llsm_ringbuffer_appendchunk(dst -> buffer_exc_mix, nx, x);
  dst -> exc_cycle = (dst -> exc_cycle + nx) % dst -> ntemplate;
}

static void llsm_fill_excitation_buffers(llsm_rtsynth_buffer_* dst) {
  for(int i = 0; i < dst -> ninternal - 1; i ++)
    for(int c = 0; c < dst -> nchannel; c ++)
      llsm_ringbuffer_append(dst -> buffer_mod_comps[c], 1e-5);
  for(int i = 0; i < 5; i ++)
    llsm_run_excitation_buffers(dst, dst -> ninternal / 5);
}

llsm_rtsynth_buffer* llsm_create_rtsynth_buffer(llsm_soptions* options,
  llsm_container* conf, int capacity_samples) {

  int* nchannel = llsm_container_get(conf, LLSM_CONF_NCHANNEL);
  FP_TYPE* thop = llsm_container_get(conf, LLSM_CONF_THOP);
  FP_TYPE* chanfreq = llsm_container_get(conf, LLSM_CONF_CHANFREQ);
  if(nchannel == NULL || thop == NULL || chanfreq == NULL) return NULL;

  llsm_rtsynth_buffer_* ret = malloc(sizeof(llsm_rtsynth_buffer_));
  ret -> nout = 0;
  ret -> nchannel = *nchannel;
  ret -> ntemplate = options -> fs;
  ret -> ninternal = options -> fs * 0.2; // 0.2 sec buffers
  ret -> conf = llsm_copy_container(conf);
  ret -> opt = *options;

  ret -> fs = options -> fs;
  ret -> thop = *thop;
  ret -> cycle = 0;
  ret -> curr_nhop = 0;
  ret -> exc_cycle = 0;
  ret -> win = NULL;
  ret -> nfft = pow(2, ceil(log2(*thop * ret -> fs * 2.2 + 32)));
  ret -> prev_nm = NULL;

  ret -> buffer_out = llsm_create_ringbuffer(capacity_samples);
  ret -> buffer_exc_mix = llsm_create_ringbuffer(ret -> ninternal);
  ret -> buffer_noise = llsm_create_ringbuffer(ret -> ninternal);
  ret -> buffer_sin   = llsm_create_ringbuffer(ret -> ninternal);

  ret -> exc_template_comps = malloc2d(*nchannel, ret -> ntemplate,
    sizeof(FP_TYPE));
  ret -> buffer_mod_comps = malloc(*nchannel * sizeof(llsm_ringbuffer*));
  for(int i = 0; i < *nchannel; i ++)
    ret -> buffer_mod_comps[i] = llsm_create_ringbuffer(ret -> ninternal);

  int npsd = *((int*)llsm_container_get(conf, LLSM_CONF_NPSD));
  FP_TYPE fnyq = *((FP_TYPE*)llsm_container_get(conf, LLSM_CONF_FNYQ));
  FP_TYPE noswarp = *((FP_TYPE*)llsm_container_get(conf, LLSM_CONF_NOSWARP));
  ret -> buffer_psd = calloc(ret -> nfft / 2 + 1, sizeof(FP_TYPE));
  ret -> buffer_fft = calloc(ret -> nfft * 4, sizeof(FP_TYPE));
  ret -> buffer_rawexc = calloc(ret -> ninternal, sizeof(FP_TYPE));
  ret -> buffer_rawmod = calloc(ret -> ninternal, sizeof(FP_TYPE));
  ret -> warp_axis = llsm_warp_frequency(0, fnyq, npsd, noswarp);

# ifdef USE_PTHREAD
  pthread_mutex_init(& ret -> buffer_out_mtx, NULL);
  pthread_cond_init(& ret -> buffer_out_cv, NULL);
# endif

  // fill in curr_nhop and reset the cycle counter
  ret -> curr_nhop = 1;
  llsm_update_cycle(ret);
  ret -> cycle = 0;
  ret -> sin_pos = -ret -> curr_nhop * 2 - ret -> nfft / 2;

  llsm_make_exc_template(ret, chanfreq);
  llsm_fill_excitation_buffers(ret);

  return (llsm_rtsynth_buffer*)ret;
}

void llsm_delete_rtsynth_buffer(llsm_rtsynth_buffer* dstptr) {
  if(dstptr == NULL) return;
  llsm_rtsynth_buffer_* dst = dstptr;
  llsm_delete_container(dst -> conf);
  llsm_delete_ringbuffer(dst -> buffer_out);
  llsm_delete_ringbuffer(dst -> buffer_exc_mix);
  llsm_delete_ringbuffer(dst -> buffer_noise);
  llsm_delete_ringbuffer(dst -> buffer_sin);
  for(int i = 0; i < dst -> nchannel; i ++)
    llsm_delete_ringbuffer(dst -> buffer_mod_comps[i]);
  llsm_delete_nmframe(dst -> prev_nm);
  free(dst -> win);
  free(dst -> buffer_mod_comps);
  free2d(dst -> exc_template_comps, dst -> nchannel);
  free(dst -> buffer_psd);
  free(dst -> buffer_fft);
  free(dst -> buffer_rawexc);
  free(dst -> buffer_rawmod);
  free(dst -> warp_axis);
# ifdef USE_PTHREAD
  pthread_mutex_destroy(& dst -> buffer_out_mtx);
  pthread_cond_destroy(& dst -> buffer_out_cv);
# endif
  free(dst);
}

// Synthesize noise modulation components.
static void llsm_rtsynth_buffer_feed_modcomps(llsm_rtsynth_buffer_* dst,
  llsm_nmframe* nm, FP_TYPE f0) {
  int nwin = dst -> curr_nhop * 2;
  llsm_hmframe* unvoiced_hm = llsm_create_hmframe(0);
  for(int c = 0; c < dst -> nchannel; c ++) {
    llsm_hmframe* hm = f0 > 0 ? nm -> eenv[c] : unvoiced_hm;
    FP_TYPE* x = llsm_synthesize_harmonic_frame_auto(& dst -> opt,
      hm -> ampl, hm -> phse, hm -> nhar, f0 / dst -> fs, nwin);
    FP_TYPE x_min = minfp(x, nwin);
    FP_TYPE offset = nm -> edc[c];
    if(x_min < - nm -> edc[c]) offset = -x_min;
    for(int i = 0; i < nwin; i ++) {
      x[i] = (x[i] + offset + 1e-8) * dst -> win[i];
    }
    llsm_ringbuffer_addchunk(dst -> buffer_mod_comps[c], -nwin, nwin, x);
    free(x);
  }
  llsm_delete_hmframe(unvoiced_hm);
}

// Synthesize sinusoidal component.
static void llsm_rtsynth_buffer_feed_sinusoids(llsm_rtsynth_buffer_* dst,
  llsm_container* frame) {
  FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
  llsm_hmframe* hm = llsm_container_get(frame, LLSM_FRAME_HM);
  if(hm != NULL && f0 != NULL && *f0 > 0) {
    FP_TYPE* x = llsm_synthesize_harmonic_frame_auto(& dst -> opt,
      hm -> ampl, hm -> phse, hm -> nhar, *f0 / dst -> fs,
      dst -> curr_nhop * 2);
    for(int i = 0; i < dst -> curr_nhop * 2; i ++)
      x[i] *= dst -> win[i];
    llsm_ringbuffer_addchunk(dst -> buffer_sin, -dst -> curr_nhop * 2,
      dst -> curr_nhop * 2, x);
    free(x);
  }
  llsm_nmframe* nm = llsm_container_get(frame, LLSM_FRAME_NM);
  if(nm != NULL)
    llsm_rtsynth_buffer_feed_modcomps(dst, nm, f0 == NULL ? 0 : *f0);
}

static void llsm_rtsynth_buffer_feed_filter(llsm_rtsynth_buffer_* dst) {
  const int nfade = 16;
  int nfft = dst -> nfft;
  int lag = nfft / 2;
  int nspec = lag + 1;
  int nwin = dst -> curr_nhop * 2;
  FP_TYPE wsqr = 0;
  for(int i = 0; i < nwin; i ++)
    wsqr += dst -> win[i] * dst -> win[i];

  FP_TYPE* psd = dst -> buffer_psd;
  FP_TYPE* fftbuffer = dst -> buffer_fft;
  FP_TYPE* x_re = fftbuffer;
  FP_TYPE* x_im = fftbuffer + nfft;
  memset(fftbuffer, 0, nfft * 4 * sizeof(FP_TYPE));

  int npsd = *((int*)llsm_container_get(dst -> conf, LLSM_CONF_NPSD));
  FP_TYPE* warp_axis = dst -> warp_axis;

  llsm_nmframe* nm = dst -> prev_nm;
  if(nm != NULL) {
    // STFT
    llsm_ringbuffer_readchunk(dst -> buffer_exc_mix, -lag - dst -> curr_nhop,
      nwin, x_re + nfft / 2 - nwin / 2);
    for(int i = 0; i < nwin; i ++)
      x_re[i - nwin / 2 + nfft / 2] *= dst -> win[i];
    fft(x_re, NULL, x_re, x_im, nfft, fftbuffer + nfft * 2);

    // PSD -> warp -> diff
    llsm_fft_to_psd(x_re, x_im, nfft, wsqr, psd);
    FP_TYPE* env = llsm_spectral_mean(psd, nspec, dst -> fs / 2.0,
      warp_axis, npsd);
    for(int i = 0; i < npsd; i ++)
      env[i] = exp_2(nm -> psd[i] / 20.0 * 2.3026) / sqrt(env[i] + 1e-8);

    // filter
    FP_TYPE* H = llsm_spectrum_from_envelope(warp_axis, env, npsd, nspec - 1,
      dst -> fs / 2.0);
    for(int i = 0; i < nspec - 1; i ++) {
      x_re[i] *= H[i]; x_im[i] *= H[i];
    }
    x_re[nspec - 1] = x_re[nspec - 2]; complete_symm (x_re, nfft);
    x_im[nspec - 1] = x_im[nspec - 2]; complete_asymm(x_im, nfft);

    // ISTFT
    ifft(x_re, x_im, x_re, NULL, nfft, fftbuffer + nfft * 2);
    for(int i = 0; i < nfade; i ++) {
      x_re[i] *= (FP_TYPE)i / nfade;
      x_re[nfft - i - 1] *= 1.0 - (FP_TYPE)i / nfade;
    }

    llsm_ringbuffer_addchunk(dst -> buffer_noise, -nfft, nfft, x_re);
    free(H); free(env);
  }
}

static void llsm_rtsynth_buffer_feed_mix(llsm_rtsynth_buffer_* dst) {
  FP_TYPE* x_nos = calloc(dst -> curr_nhop, sizeof(FP_TYPE));
  FP_TYPE* x_sin = calloc(dst -> curr_nhop, sizeof(FP_TYPE));
  llsm_ringbuffer_readchunk(dst -> buffer_noise,
    -dst -> nfft, dst -> curr_nhop, x_nos);
  llsm_ringbuffer_readchunk(dst -> buffer_sin,
    dst -> sin_pos, dst -> curr_nhop, x_sin);
  for(int i = 0; i < dst -> curr_nhop; i ++)
    x_nos[i] += x_sin[i];
# ifdef USE_PTHREAD
  pthread_mutex_lock(& dst -> buffer_out_mtx);
  while(dst -> nout > dst -> buffer_out -> capacity - dst -> curr_nhop)
    pthread_cond_wait(& dst -> buffer_out_cv, & dst -> buffer_out_mtx);
# endif
  llsm_ringbuffer_appendchunk(dst -> buffer_out, dst -> curr_nhop, x_nos);
  dst -> nout += dst -> curr_nhop;
# ifdef USE_PTHREAD
  pthread_cond_broadcast(& dst -> buffer_out_cv);
  pthread_mutex_unlock(& dst -> buffer_out_mtx);
# endif
  free(x_nos);
  free(x_sin);
}

void llsm_rtsynth_buffer_feed(llsm_rtsynth_buffer* ptr,
  llsm_container* frame) {
  llsm_rtsynth_buffer_* dst = ptr;
  llsm_update_cycle(dst);
  llsm_rtsynth_buffer_feed_sinusoids(dst, frame);
  llsm_run_excitation_buffers(dst, dst -> curr_nhop);
  llsm_rtsynth_buffer_feed_filter(dst);
  llsm_rtsynth_buffer_feed_mix(dst);
  if(dst -> prev_nm != NULL) llsm_delete_nmframe(dst -> prev_nm);
  dst -> prev_nm = llsm_container_get(frame, LLSM_FRAME_NM);
  if(dst -> prev_nm != NULL)
    dst -> prev_nm = llsm_copy_nmframe(dst -> prev_nm);
}

int llsm_rtsynth_buffer_fetch(llsm_rtsynth_buffer* ptr, FP_TYPE* dst) {
  llsm_rtsynth_buffer_* src = ptr;
# ifdef USE_PTHREAD
  pthread_mutex_lock(& src -> buffer_out_mtx);
# endif
  if(src -> nout > 0) {
    *dst = llsm_ringbuffer_read(src -> buffer_out, -src -> nout);
    src -> nout --;
#   ifdef USE_PTHREAD
    pthread_mutex_unlock(& src -> buffer_out_mtx);
    pthread_cond_broadcast(& src -> buffer_out_cv);
#   endif
    return 1;
  } else
#   ifdef USE_PTHREAD
    pthread_mutex_unlock(& src -> buffer_out_mtx);
#   endif
    return 0;
}

int llsm_rtsynth_buffer_getlatency(llsm_rtsynth_buffer* ptr) {
  llsm_rtsynth_buffer_* src = ptr;
  return -src -> sin_pos - src -> curr_nhop;
}

int llsm_rtsynth_buffer_numoutput(llsm_rtsynth_buffer* ptr) {
  llsm_rtsynth_buffer_* src = ptr;
  return src -> nout;
}

void llsm_rtsynth_buffer_clear(llsm_rtsynth_buffer* ptr) {
  llsm_rtsynth_buffer_* dst = ptr;
  int capacity_samples = dst -> buffer_out -> capacity;
  dst -> nout = 0;
  dst -> cycle = 0;
  dst -> exc_cycle = 0;
  llsm_delete_ringbuffer(dst -> buffer_out);
  llsm_delete_ringbuffer(dst -> buffer_exc_mix);
  llsm_delete_ringbuffer(dst -> buffer_noise);
  llsm_delete_ringbuffer(dst -> buffer_sin);
  dst -> buffer_out = llsm_create_ringbuffer(capacity_samples);
  dst -> buffer_exc_mix = llsm_create_ringbuffer(dst -> ninternal);
  dst -> buffer_noise = llsm_create_ringbuffer(dst -> ninternal);
  dst -> buffer_sin   = llsm_create_ringbuffer(dst -> ninternal);
  dst -> curr_nhop = 1;
  llsm_update_cycle(dst);
  dst -> sin_pos = -dst -> curr_nhop * 2 - dst -> nfft / 2;
}

/*
  libllsm2 - Low Level Speech Model (version 2)
  ===
  Copyright (c) 2017-2019 Kanru Hua.

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

#include <ciglet/ciglet.h>

#include "llsmrt.h"
#include "buffer.h"
#include "dsputils.h"
#include "llsmutils.h"
#include "constants.h"

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
                 //   (in seconds)
  FP_TYPE pulse; // the most recent pulse location relative to floored sample
                 //   position (in samples)
  int curr_nhop; // rounding-adjusted current hop size (in samples)
  int next_nhop; // rounding-adjusted next hop size (in samples)
  int exc_cycle; // current position in the noise template
  int sin_pos;   // read position of the sinusoid buffer
  FP_TYPE* win;  // overlap-add window
  int nfft;      // FFT size for noise filtering
  int pbp_offset; // sample offset for mixing pulses into buffer_sin
  int pbp_state;  // on = 1, off = 0
  llsm_nmframe* prev_nm; // the previous noise model frame

  llsm_ringbuffer*  buffer_out;       // output buffer of audio samples

  FP_TYPE** exc_template_comps;       // noise template for each channel
  llsm_ringbuffer** buffer_mod_comps; // noise modulation buffer for each
                                      //   channel
  llsm_ringbuffer*  buffer_exc_mix;   // buffer for the sum of noise excitation
  llsm_ringbuffer*  buffer_noise;     // buffer for the filtered noise
  llsm_ringbuffer*  buffer_sin;       // buffer for the sinusoidal component
  llsm_dualbuffer*  buffer_pulse;     // buffer for the sum of pulses (PBPSYN)

  FP_TYPE* buffer_psd; // size: nspec
  FP_TYPE* buffer_fft; // size: nfft * 4
  FP_TYPE* buffer_rawexc; // size: ninternal
  FP_TYPE* buffer_rawmod; // size: ninternal
  FP_TYPE* buffer_phase;  // size: nfft
  FP_TYPE* psd_axis;  // size: npsd

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
    FP_TYPE* x = llsm_generate_bandlimited_noise(
      dst -> ntemplate, fmin / fs, fmax / fs);
    for(int j = 0; j < dst -> ntemplate; j ++)
      dst -> exc_template_comps[c][j] = llsm_get_circular_noise(x,
        dst -> ntemplate, j);
    free(x);
  }
}

// Prepare for synthesizing the next frame.
static void llsm_update_cycle(llsm_rtsynth_buffer_* dst) {
  int prev_nhop = dst -> curr_nhop;
  dst -> cycle += dst -> thop;
  dst -> curr_nhop = floor(dst -> cycle * dst -> fs);

  dst -> cycle -= (FP_TYPE)prev_nhop / dst -> fs;
  dst -> pulse -= prev_nhop;
  if(dst -> pbp_state && dst -> pbp_offset > dst -> sin_pos + dst -> curr_nhop)
    dst -> pbp_offset -= prev_nhop;
  int nwin = dst -> curr_nhop * 2;
  free(dst -> win); dst -> win = hanning_2(nwin);

  dst -> next_nhop = floor((dst -> cycle + dst -> thop) * dst -> fs);

  for(int c = 0; c < dst -> nchannel; c ++)
    llsm_ringbuffer_appendblank(dst -> buffer_mod_comps[c], dst -> curr_nhop);
  llsm_ringbuffer_appendblank(dst -> buffer_sin, dst -> curr_nhop);
  llsm_ringbuffer_appendblank(dst -> buffer_noise, dst -> curr_nhop);
  llsm_dualbuffer_forward(dst -> buffer_pulse, dst -> curr_nhop);
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
  ret -> pulse = 0;
  ret -> curr_nhop = 0;
  ret -> next_nhop = 0;
  ret -> exc_cycle = 0;
  ret -> win = NULL;
  ret -> nfft = pow(2, ceil(log2(*thop * ret -> fs * 2.2 + 32)));
  ret -> pbp_offset = 0;
  ret -> pbp_state = 0;
  ret -> prev_nm = NULL;

  ret -> buffer_out = llsm_create_ringbuffer(capacity_samples);
  ret -> buffer_exc_mix = llsm_create_ringbuffer(ret -> ninternal);
  ret -> buffer_noise = llsm_create_ringbuffer(ret -> ninternal);
  ret -> buffer_sin   = llsm_create_ringbuffer(ret -> ninternal);
  ret -> buffer_pulse = llsm_create_dualbuffer(ret -> ninternal);

  ret -> exc_template_comps = malloc2d(*nchannel, ret -> ntemplate,
    sizeof(FP_TYPE));
  ret -> buffer_mod_comps = malloc(*nchannel * sizeof(llsm_ringbuffer*));
  for(int i = 0; i < *nchannel; i ++)
    ret -> buffer_mod_comps[i] = llsm_create_ringbuffer(ret -> ninternal);

  int npsd = *((int*)llsm_container_get(conf, LLSM_CONF_NPSD));
  FP_TYPE fnyq = *((FP_TYPE*)llsm_container_get(conf, LLSM_CONF_FNYQ));
  ret -> buffer_psd = calloc(ret -> nfft / 2 + 1, sizeof(FP_TYPE));
  ret -> buffer_fft = calloc(ret -> nfft * 4, sizeof(FP_TYPE));
  ret -> buffer_rawexc = calloc(ret -> ninternal, sizeof(FP_TYPE));
  ret -> buffer_rawmod = calloc(ret -> ninternal, sizeof(FP_TYPE));
  ret -> buffer_phase  = calloc(ret -> nfft, sizeof(FP_TYPE));
  ret -> psd_axis = linspace(0, fnyq, npsd);

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
  llsm_delete_dualbuffer(dst -> buffer_pulse);
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
  free(dst -> buffer_phase);
  free(dst -> psd_axis);
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
    FP_TYPE offset = nm -> edc[c];
    for(int i = 0; i < nwin; i ++)
      x[i] = (max(x[i] + offset, 1e-8)) * dst -> win[i];
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
  FP_TYPE* phase = dst -> buffer_phase;
  if(hm != NULL && f0 != NULL && *f0 > 0) {
    FP_TYPE phase_shift = dst -> cycle * 2 * M_PI * (*f0);
    int nhar = min(hm -> nhar, dst -> nfft);
    for(int k = 0; k < nhar; k ++)
      phase[k] = hm -> phse[k] - phase_shift * (k + 1.0);
    FP_TYPE* x = llsm_synthesize_harmonic_frame_auto(& dst -> opt,
      hm -> ampl, phase, nhar, *f0 / dst -> fs, dst -> curr_nhop * 2);
    for(int i = 0; i < dst -> curr_nhop * 2; i ++)
      x[i] *= dst -> win[i];
    llsm_ringbuffer_addchunk(dst -> buffer_sin, -dst -> curr_nhop * 2,
      dst -> curr_nhop * 2, x);
    free(x);
  }
}

// Synthesize deterministic component (semi-harmonic excitation and
//   noise envelope components).
static void llsm_rtsynth_buffer_feed_deterministic(llsm_rtsynth_buffer_* dst,
  llsm_container* frame) {
  FP_TYPE* f0 = llsm_container_get(frame, LLSM_FRAME_F0);
  llsm_nmframe* nm = llsm_container_get(frame, LLSM_FRAME_NM);
  if(nm != NULL)
    llsm_rtsynth_buffer_feed_modcomps(dst, nm, f0 == NULL ? 0 : *f0);
  if(! dst -> opt.use_l1) {
    llsm_rtsynth_buffer_feed_sinusoids(dst, frame);
    return;
  }
  int nhop = dst -> curr_nhop;

  FP_TYPE* vsphse = llsm_container_get(frame, LLSM_FRAME_VSPHSE);
  FP_TYPE* rd = llsm_container_get(frame, LLSM_FRAME_RD);
  int* pbpsyn = llsm_container_get(frame, LLSM_FRAME_PBPSYN);
  llsm_pbpeffect* pbpeff = llsm_container_get(frame, LLSM_FRAME_PBPEFF);
  if(vsphse == NULL || rd == NULL || *f0 == 0) return;
  int pbp_on = pbpsyn != NULL && pbpsyn[0] == 1;
  int nspec = *((int*)llsm_container_get(dst -> conf, LLSM_CONF_NSPEC));

  // update locations of pulses locked onto the first source harmonic
  FP_TYPE len_period = dst -> fs / f0[0];
  FP_TYPE t_period = 1.0 / f0[0];
  lfmodel source_model = lfmodel_from_rd(*rd, t_period, 1.0);
  FP_TYPE source_p0 = 0;
  free(lfmodel_spectrum(source_model, f0, 1, & source_p0));
  source_p0 -= 0.5 * M_PI; // integrate (flow derivative to flow velocity)
  FP_TYPE p0 = wrap(vsphse[0]);
  FP_TYPE p0_dist = phase_diff(source_p0, p0);
  if(p0_dist < 0) p0_dist += 2.0 * M_PI;
  // the next position where a glottal flow cycle begins (relative to
  //   -curr_nhop in the buffer)
  FP_TYPE pulse_projected = p0_dist / 2 / M_PI * len_period;
  // reset the pulse tracker after unvoiced part
  int len_reset = max(len_period, nhop) * 2;
  if(pulse_projected - dst -> pulse > len_reset)
    dst -> pulse = pulse_projected - len_reset;
  int num_periods = round((pulse_projected - dst -> pulse) / len_period);
  if(num_periods > 0)
    len_period = (pulse_projected - dst -> pulse) / num_periods;
  int pulse_size = pow(2, ceil(log2(max(len_period * 2, nspec))));

  FP_TYPE* fnyq = llsm_container_get(dst -> conf, LLSM_CONF_FNYQ);
  FP_TYPE* liprad = llsm_container_get(dst -> conf, LLSM_CONF_LIPRADIUS);

  int pbp_onset = 0;
  int pbp_termination = 0;
  if(pbp_on && ! dst -> pbp_state) {
    pbp_onset = 1;
    dst -> pbp_state = 1;
    dst -> pbp_offset = -nhop;
    if(llsm_container_get(frame, LLSM_FRAME_HM) == NULL)
      llsm_frame_tolayer0(frame, dst -> conf);
    llsm_rtsynth_buffer_feed_sinusoids(dst, frame);
  }
  if(! pbp_on && dst -> pbp_state) {
    pbp_termination = 1;
    dst -> pbp_state = 0;
    num_periods += ceil((-dst -> pbp_offset) / len_period);
  }

  // Pulse OLA
  if(dst -> pbp_state || pbp_termination) {
    int period_begin = pbp_onset ? -2 : 0;
    int period_end = num_periods;
    int num_pulses = period_end - period_begin;
    int pre_rotate = min(len_period, nhop * 2);
    if(num_pulses > 0) {
      FP_TYPE* offsets = calloc(num_pulses, sizeof(FP_TYPE));
      lfmodel* sources = calloc(num_pulses, sizeof(lfmodel));
      for(int i = 0; i < num_pulses; i ++) {
        FP_TYPE delta_t = 0;
        if(pbpeff != NULL) {
          llsm_gfm g = llsm_lfmodel_to_gfm(source_model);
          pbpeff -> modifier(& g, & delta_t, pbpeff -> info, frame);
          sources[i] = llsm_gfm_to_lfmodel(g);
        } else
          sources[i] = source_model;
        offsets[i] = dst -> pulse + (i + period_begin) * len_period
                   + delta_t * dst -> fs;
      }
      int pulse_base = offsets[0];
      for(int i = 0; i < num_pulses; i ++) offsets[i] -= pulse_base;
      FP_TYPE* y = llsm_make_filtered_pulse(frame, sources, offsets,
        num_pulses, pre_rotate, pulse_size, *fnyq, *liprad, dst -> fs);
      llsm_dualbuffer_addchunk(dst -> buffer_pulse,
        pulse_base - pre_rotate - nhop, pulse_size, y);
      free(y);
      free(offsets);
      free(sources);
    }
  }
  if(! dst -> pbp_state) {
    if(llsm_container_get(frame, LLSM_FRAME_HM) == NULL)
      llsm_frame_tolayer0(frame, dst -> conf);
    llsm_rtsynth_buffer_feed_sinusoids(dst, frame);
  }
  dst -> pulse = pulse_projected;

  if(dst -> pbp_state &&
     dst -> pbp_offset <= dst -> sin_pos + nhop) {
    FP_TYPE* x = calloc(nhop * 2, sizeof(FP_TYPE));
    llsm_dualbuffer_readchunk(
      dst -> buffer_pulse, dst -> pbp_offset, nhop * 2, x);
    for(int i = 0; i < nhop * 2; i ++)
      x[i] *= dst -> win[i];
    llsm_ringbuffer_addchunk(
      dst -> buffer_sin, dst -> pbp_offset, nhop * 2, x);
    free(x);
  }
  if(pbp_termination) {
    // Overlap-add PbP results using a trapezoid-like window to catch up
    //   with the harmonic model.
    int size = -nhop - dst -> pbp_offset;
    FP_TYPE* x = calloc(size, sizeof(FP_TYPE));
    llsm_dualbuffer_readchunk(
      dst -> buffer_pulse, dst -> pbp_offset, size, x);
    for(int i = 0; i < nhop; i ++) {
      x[i] *= dst -> win[i];
      x[size - nhop + i] *= dst -> win[i + nhop];
    }
    llsm_ringbuffer_addchunk(
      dst -> buffer_sin, dst -> pbp_offset, size, x);
    free(x);
  }
}

static void llsm_rtsynth_buffer_feed_filter(llsm_rtsynth_buffer_* dst) {
  const int nfade = 16;
  int nfft = dst -> nfft;
  int nspec = nfft / 2 + 1;
  int nhop = dst -> curr_nhop;
  int nwin = nhop * 2;
  FP_TYPE wsqr = 0;
  for(int i = 0; i < nwin; i ++)
    wsqr += dst -> win[i] * dst -> win[i];

  FP_TYPE* psd = dst -> buffer_psd;
  FP_TYPE* fftbuffer = dst -> buffer_fft;
  FP_TYPE* x_re = fftbuffer;
  FP_TYPE* x_im = fftbuffer + nfft;
  memset(fftbuffer, 0, nfft * 4 * sizeof(FP_TYPE));

  int npsd = *((int*)llsm_container_get(dst -> conf, LLSM_CONF_NPSD));
  FP_TYPE* psd_axis = dst -> psd_axis;

  llsm_nmframe* nm = dst -> prev_nm;
  if(nm != NULL) {
    FP_TYPE peak = maxfp(nm -> psd, npsd);
    if(peak < -100) return; // -100 dB noise floor

    // STFT
    llsm_ringbuffer_readchunk(dst -> buffer_exc_mix, -nhop * 2,
      nwin, x_re + nfft / 2 - nhop);
    for(int i = 0; i < nwin; i ++)
      x_re[i - nhop + nfft / 2] *= dst -> win[i];
    fft(x_re, NULL, x_re, x_im, nfft, fftbuffer + nfft * 2);

    // PSD -> diff
    llsm_fft_to_psd(x_re, x_im, nfft, wsqr, psd);
    FP_TYPE* env = moving_avg(psd, nspec, 3);
    FP_TYPE* H = llsm_spectrum_from_envelope(
      psd_axis, nm -> psd, npsd, nspec - 1, dst -> fs / 2.0);
    for(int j = 0; j < nspec - 1; j ++)
      H[j] = exp(DB2LOG(H[j])) / sqrt(env[j] * 44100 / dst -> fs + 1e-8);

    // filter
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
  FP_TYPE* x_nos = calloc(dst -> next_nhop, sizeof(FP_TYPE));
  FP_TYPE* x_sin = calloc(dst -> next_nhop, sizeof(FP_TYPE));
  llsm_ringbuffer_readchunk(dst -> buffer_noise,
    -dst -> nfft, dst -> next_nhop, x_nos);
  llsm_ringbuffer_readchunk(dst -> buffer_sin,
    dst -> sin_pos, dst -> next_nhop, x_sin);
  for(int i = 0; i < dst -> next_nhop; i ++)
    x_nos[i] += x_sin[i];
# ifdef USE_PTHREAD
  pthread_mutex_lock(& dst -> buffer_out_mtx);
  while(dst -> nout > dst -> buffer_out -> capacity - dst -> next_nhop)
    pthread_cond_wait(& dst -> buffer_out_cv, & dst -> buffer_out_mtx);
# endif
  llsm_ringbuffer_appendchunk(dst -> buffer_out, dst -> next_nhop, x_nos);
  dst -> nout += dst -> next_nhop;
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
  llsm_rtsynth_buffer_feed_deterministic(dst, frame);
  llsm_run_excitation_buffers(dst, dst -> curr_nhop);
  llsm_rtsynth_buffer_feed_filter(dst);
  llsm_rtsynth_buffer_feed_mix(dst);
  if(dst -> prev_nm != NULL) llsm_delete_nmframe(dst -> prev_nm);
  dst -> prev_nm = llsm_container_get(frame, LLSM_FRAME_NM);
  if(dst -> prev_nm != NULL)
    dst -> prev_nm = llsm_copy_nmframe(dst -> prev_nm);
  FP_TYPE* resvec = llsm_container_get(frame, LLSM_FRAME_PSDRES);
  if(resvec != NULL)
  for(int j = 0; j < dst -> prev_nm -> npsd; j ++)
    dst -> prev_nm -> psd[j] += resvec[j] - LOG2IN(LOGRESBIAS);
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
  llsm_delete_ringbuffer(dst -> buffer_out);
  llsm_delete_ringbuffer(dst -> buffer_exc_mix);
  llsm_delete_ringbuffer(dst -> buffer_noise);
  llsm_delete_ringbuffer(dst -> buffer_sin);
  llsm_delete_dualbuffer(dst -> buffer_pulse);
  dst -> buffer_out = llsm_create_ringbuffer(capacity_samples);
  dst -> buffer_exc_mix = llsm_create_ringbuffer(dst -> ninternal);
  dst -> buffer_noise = llsm_create_ringbuffer(dst -> ninternal);
  dst -> buffer_sin   = llsm_create_ringbuffer(dst -> ninternal);
  dst -> buffer_pulse = llsm_create_dualbuffer(dst -> ninternal);
  dst -> curr_nhop = 1;
  dst -> pbp_offset = 0;
  dst -> pbp_state = 0;
  llsm_update_cycle(dst);
  dst -> cycle = 0;
  dst -> pulse = 0;
  dst -> exc_cycle = 0;
  dst -> sin_pos = -dst -> curr_nhop * 2 - dst -> nfft / 2;
}

#include "llsm.h"
#include "dsputils.h"
#include "external/ciglet/ciglet.h"

llsm_aoptions* llsm_create_aoptions() {
  llsm_aoptions* ret = malloc(sizeof(llsm_aoptions));
  ret -> thop = 0.005;
  ret -> maxnhar = 100;
  ret -> maxnhar_e = 4;
  ret -> npsd = 64;
  ret -> nchannel = 4;
  ret -> chanfreq = calloc(3, sizeof(FP_TYPE));
  ret -> chanfreq[0] = 2000.0;
  ret -> chanfreq[1] = 4000.0;
  ret -> chanfreq[2] = 8000.0;
  ret -> noise_warp = 15000.0;
  ret -> lip_radius = 1.5;
  ret -> f0_refine = 1;
  ret -> hm_method = LLSM_AOPTION_HMCZT;
  ret -> rel_winsize = 4.0;
  return ret;
}

void llsm_delete_aoptions(llsm_aoptions* dst) {
  if(dst == NULL) return;
  free(dst -> chanfreq);
  free(dst);
}

llsm_container* llsm_aoptions_toconf(llsm_aoptions* src, FP_TYPE fnyq) {
  llsm_container* ret = llsm_create_container(10);
  llsm_container_attach(ret, LLSM_CONF_NFRM,
    llsm_create_int(0), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_THOP,
    llsm_create_fp(src -> thop), llsm_delete_fp, llsm_copy_fp);
  llsm_container_attach(ret, LLSM_CONF_MAXNHAR,
    llsm_create_int(src -> maxnhar), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_MAXNHAR_E,
    llsm_create_int(src -> maxnhar_e), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_NPSD,
    llsm_create_int(src -> npsd), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_NOSWARP,
    llsm_create_fp(src -> noise_warp), llsm_delete_fp, llsm_copy_fp);
  llsm_container_attach(ret, LLSM_CONF_FNYQ,
    llsm_create_fp(fnyq), llsm_delete_fp, llsm_copy_fp);
  llsm_container_attach(ret, LLSM_CONF_NCHANNEL,
    llsm_create_int(src -> nchannel), llsm_delete_int, llsm_copy_int);
  // llsm_container_attach(ret, LLSM_CONF_NSPEC,
  //   llsm_create_int(1024), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_LIPRADIUS,
    llsm_create_fp(src -> lip_radius), llsm_delete_fp, llsm_copy_fp);
  FP_TYPE* chanfreq = llsm_create_fparray(src -> nchannel - 1);
  memcpy(chanfreq, src -> chanfreq, sizeof(FP_TYPE) * (src -> nchannel - 1));
  llsm_container_attach(ret, LLSM_CONF_CHANFREQ,
    chanfreq, llsm_delete_fparray, llsm_copy_fparray);
  return ret;
}

llsm_soptions* llsm_create_soptions(FP_TYPE fs) {
  llsm_soptions* ret = malloc(sizeof(llsm_soptions));
  ret -> fs = fs;
  // See test/test-harmonic.c for more information.
  ret -> use_iczt = 1;
  ret -> iczt_param_a = 0.275;
  ret -> iczt_param_b = 2.26;
  return ret;
}

void llsm_delete_soptions(llsm_soptions* dst) {
  if(dst == NULL) return;
  free(dst);
}

static void llsm_analyze_harmonics(llsm_aoptions* options, FP_TYPE* x, int nx,
  FP_TYPE fs, FP_TYPE* f0, int nfrm, llsm_chunk* dst_chunk) {

  int*      tmp_nhar = calloc(nfrm, sizeof(int));
  FP_TYPE** tmp_ampl = calloc(nfrm, sizeof(FP_TYPE*));
  FP_TYPE** tmp_phse = calloc(nfrm, sizeof(FP_TYPE*));

  llsm_harmonic_analysis(x, nx, fs, f0, nfrm, options -> thop,
    options -> rel_winsize, options -> maxnhar, options -> hm_method,
    tmp_nhar, tmp_ampl, tmp_phse);

  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] == 0) continue;
    llsm_hmframe* hm = llsm_create_hmframe(tmp_nhar[i]);
    memcpy(hm -> ampl, tmp_ampl[i], tmp_nhar[i] * sizeof(FP_TYPE));
    memcpy(hm -> phse, tmp_phse[i], tmp_nhar[i] * sizeof(FP_TYPE));
    llsm_container_attach(dst_chunk -> frames[i], LLSM_FRAME_HM,
      hm, llsm_delete_hmframe, llsm_copy_hmframe);
  }

  free2d(tmp_ampl, nfrm); free2d(tmp_phse, nfrm); free(tmp_nhar);
}

static FP_TYPE* llsm_synthesize_harmonic_frame_auto(llsm_soptions* options,
  FP_TYPE* ampl, FP_TYPE* phse, int nhar, FP_TYPE f0, int nx) {
  FP_TYPE* ret = NULL;
  if(options == NULL || (! options -> use_iczt))
    ret = llsm_synthesize_harmonic_frame(ampl, phse, nhar, f0, nx);
  else {
    if(log_1(nx) * options -> iczt_param_a <
       log_1(nhar) - options -> iczt_param_b)
      ret = llsm_synthesize_harmonic_frame_iczt(ampl, phse, nhar, f0, nx);
    else
      ret = llsm_synthesize_harmonic_frame(ampl, phse, nhar, f0, nx);
  }
  return ret;
}

static FP_TYPE* llsm_synthesize_harmonics(llsm_soptions* options,
  llsm_chunk* chunk, FP_TYPE* f0, int nfrm, FP_TYPE thop, FP_TYPE fs, int ny) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  int nwin = round(thop * 2.0 * fs);
  FP_TYPE* w = hanning(nwin);
  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] == 0) continue; // skip unvoiced frames
    llsm_hmframe* hm = llsm_container_get(chunk -> frames[i], LLSM_FRAME_HM);
    FP_TYPE* yi = llsm_synthesize_harmonic_frame_auto(options,
      hm -> ampl, hm -> phse, hm -> nhar, f0[i] / fs, nwin);
    for(int j = 0; j < nwin; j ++) {
      yi[j] *= w[j];
      int idx = round((i - 1) * thop * fs + j);
      if(idx >= 0 && idx < ny)
        y[idx] += yi[j];
    }
    free(yi);
  }
  free(w);
  return y;
}

static FP_TYPE* llsm_synthesize_noise_envelope(llsm_soptions* options,
  llsm_chunk* chunk, int channel, FP_TYPE* f0, int nfrm, FP_TYPE thop,
  FP_TYPE fs, int ny) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  int nwin = round(thop * 2.0 * fs);
  FP_TYPE* w = hanning(nwin);
  llsm_hmframe* unvoiced_hm = llsm_create_hmframe(0);
  for(int i = 0; i < nfrm; i ++) {
    llsm_nmframe* nm = llsm_container_get(chunk -> frames[i], LLSM_FRAME_NM);
    llsm_hmframe* hm = f0[i] > 0 ? nm -> eenv[channel] : unvoiced_hm;
    FP_TYPE* yi = llsm_synthesize_harmonic_frame_auto(options,
      hm -> ampl, hm -> phse, hm -> nhar, f0[i] / fs, nwin);
    // Make sure the envelope is positive.
    FP_TYPE y_min = minfp(yi, nwin);
    FP_TYPE offset = nm -> edc[channel];
    if(y_min < - nm -> edc[channel]) offset = -y_min;
    for(int j = 0; j < nwin; j ++) {
      yi[j] = (yi[j] + offset + 1e-8) * w[j];
      int idx = round((i - 1) * thop * fs + j);
      if(idx >= 0 && idx < ny)
        y[idx] += yi[j];
    }
    free(yi);
  }
  llsm_delete_hmframe(unvoiced_hm);
  free(w);
  return y;
}

static void llsm_analyze_noise_psd(llsm_aoptions* options, FP_TYPE* x_res,
  int nx, FP_TYPE fs, int nfrm, llsm_chunk* dst_chunk) {
  int nwin = round(options -> thop * 4 * fs);
  int nfft = pow(2, ceil(log2(nwin)));
  int nspec = nfft / 2 + 1;
  FP_TYPE* frame_psd = calloc(nspec, sizeof(FP_TYPE));
  FP_TYPE* warp_axis = llsm_warp_frequency(0, fs / 2.0, options -> npsd,
    options -> noise_warp);

  // For each frame, compute PSD, frequency warp, and save into dst_chunk.
  for(int i = 0; i < nfrm; i ++) {
    llsm_nmframe* dst_nm = llsm_container_get(dst_chunk -> frames[i],
      LLSM_FRAME_NM);
    int center = round(i * options -> thop * fs);
    FP_TYPE* xfrm = fetch_frame(x_res, nx, center, nwin);
    llsm_estimate_psd(xfrm, nwin, nfft, frame_psd);
    FP_TYPE* wpsd = llsm_spectral_mean(frame_psd, nspec - 1, fs / 2.0,
      warp_axis, options -> npsd);
    // The PSD is squared and hence 10 * log10(.)
    for(int j = 0; j < options -> npsd; j ++)
      dst_nm -> psd[j] = 10.0 * log10(wpsd[j]);
    free(xfrm); free(wpsd);
  }
  free(warp_axis);
  free(frame_psd);
}

static void llsm_analyze_noise_envelope(llsm_aoptions* options,
  FP_TYPE* x, FP_TYPE* x_res, int nx, FP_TYPE fs, FP_TYPE* f0,
  int nfrm, llsm_chunk* dst_chunk) {
  
  int*      tmp_nhar = calloc(nfrm, sizeof(int));
  FP_TYPE** tmp_ampl = calloc(nfrm, sizeof(FP_TYPE*));
  FP_TYPE** tmp_phse = calloc(nfrm, sizeof(FP_TYPE*));

  FP_TYPE*  tmp_dc   = calloc(nfrm, sizeof(FP_TYPE));
  int*      center   = calloc(nfrm, sizeof(int));
  int*      nwin     = calloc(nfrm, sizeof(int));
  for(int i = 0; i < nfrm; i ++) {
    center[i] = round(i * options -> thop * fs);
    nwin[i] = round((f0[i] == 0 ? options -> thop * 2 : 2.0 / f0[i]) * fs);
  }

  for(int c = 0; c < options -> nchannel; c ++) {
    FP_TYPE fmin = c == 0 ? 0 : options -> chanfreq[c - 1];
    FP_TYPE fmax = c == options -> nchannel - 1 ?
                   fs / 2.0 : options -> chanfreq[c];
    // Trick: extract envelope from the original waveform in high frequencies
    //   where the residual is often smeared due to harmonic analysis errors.
    FP_TYPE* ce = llsm_subband_energy(
      fmin > 6000.0 ? x : x_res, nx, fmin / fs * 2.0, fmax / fs * 2.0);
    // Perform harmonic analysis on a squared signal and extract the lower
    //   harmonics is roughly equivalent to modeling the RMS envelope.
    llsm_harmonic_analysis(ce, nx, fs, f0, nfrm, options -> thop,
      options -> rel_winsize, options -> maxnhar_e, options -> hm_method,
      tmp_nhar, tmp_ampl, tmp_phse);
    llsm_compute_dc(ce, nx, center, nwin, nfrm, tmp_dc);
    // Store the results.
    for(int i = 0; i < nfrm; i ++) {
      llsm_nmframe* dst_nm = llsm_container_get(dst_chunk -> frames[i],
        LLSM_FRAME_NM);
      dst_nm -> edc[c] = tmp_dc[i];
      if(f0[i] == 0) continue;
      llsm_hmframe* hm = llsm_create_hmframe(tmp_nhar[i]);
      memcpy(hm -> ampl, tmp_ampl[i], tmp_nhar[i] * sizeof(FP_TYPE));
      memcpy(hm -> phse, tmp_phse[i], tmp_nhar[i] * sizeof(FP_TYPE));
      llsm_copy_hmframe_inplace(dst_nm -> eenv[c], hm);
      llsm_delete_hmframe(hm);
    }

    for(int i = 0; i < nfrm; i ++) {
      free(tmp_ampl[i]); tmp_ampl[i] = NULL;
      free(tmp_phse[i]); tmp_phse[i] = NULL;
    }
    free(ce);
  }

  free2d(tmp_ampl, nfrm); free2d(tmp_phse, nfrm); free(tmp_nhar);
  free(tmp_dc); free(center); free(nwin);
}

static void llsm_analyze_noise(llsm_aoptions* options, FP_TYPE* x,
  FP_TYPE* x_res, int nx, FP_TYPE fs, FP_TYPE* f0, int nfrm,
  llsm_chunk* dst_chunk) {
  llsm_analyze_noise_psd(options, x_res, nx, fs, nfrm, dst_chunk);
  llsm_analyze_noise_envelope(options, x, x_res, nx, fs, f0, nfrm, dst_chunk);
}

llsm_chunk* llsm_analyze(llsm_aoptions* options, FP_TYPE* x, int nx,
  FP_TYPE fs, FP_TYPE* f0, int nfrm, FP_TYPE** x_ap) {

  llsm_container* conf = llsm_aoptions_toconf(options, fs / 2.0);
  ((int*)llsm_container_get(conf, LLSM_CONF_NFRM))[0] = nfrm;
  llsm_chunk* ret = llsm_create_chunk(conf, 1);
  llsm_delete_container(conf); // conf gets copied into ret; no longer needed.
  conf = ret -> conf;

  if(options -> f0_refine)
    llsm_refine_f0(x, nx, fs, f0, nfrm, options -> thop);

  // set F0 for all the frames
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* if0 = llsm_container_get(ret -> frames[i], LLSM_FRAME_F0);
    if0[0] = f0[i];
  }

  // harmonic analysis and residual extraction
  llsm_analyze_harmonics(options, x, nx, fs, f0, nfrm, ret);
  FP_TYPE* x_sin = llsm_synthesize_harmonics(NULL, ret, f0, nfrm,
    options -> thop, fs, nx);
  FP_TYPE* x_res = calloc(nx, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) x_res[i] = x[i] - x_sin[i];
  free(x_sin);
  if(x_ap != NULL) *x_ap = x_res;

  // noise analysis
  llsm_analyze_noise(options, x, x_res, nx, fs, f0, nfrm, ret);

  if(x_ap == NULL)
    free(x_res);
  return ret;
}

int llsm_conf_checklayer0(llsm_container* src) {
  int* nfrm = llsm_container_get(src, LLSM_CONF_NFRM);
  FP_TYPE* thop = llsm_container_get(src, LLSM_CONF_THOP);
  int* npsd = llsm_container_get(src, LLSM_CONF_NPSD);
  FP_TYPE* noswarp = llsm_container_get(src, LLSM_CONF_NOSWARP);
  FP_TYPE* fnyq = llsm_container_get(src, LLSM_CONF_FNYQ);
  int* nchannel = llsm_container_get(src, LLSM_CONF_NCHANNEL);
  FP_TYPE* chanfreq = llsm_container_get(src, LLSM_CONF_CHANFREQ);
  if(nfrm == NULL || thop == NULL || npsd == NULL || noswarp == NULL ||
     fnyq == NULL || nchannel == NULL || chanfreq == NULL) return 0;
  return 1;
}

static int llsm_synthesis_check_integrity(llsm_chunk* src) {
  if(! llsm_conf_checklayer0(src -> conf)) return 0;
  int* nfrm = llsm_container_get(src -> conf, LLSM_CONF_NFRM);
  for(int i = 0; i < *nfrm; i ++)
    if(! llsm_frame_checklayer0(src -> frames[i]))
      return 0;
  return 1;
}

static FP_TYPE* llsm_synthesize_noise_excitation(llsm_soptions* options,
  llsm_chunk* src, FP_TYPE* f0, int nfrm, FP_TYPE thop, FP_TYPE fs, int ny) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  FP_TYPE* chanfreq = llsm_container_get(src -> conf, LLSM_CONF_CHANFREQ);
  int nchannel = *((int*)llsm_container_get(src -> conf, LLSM_CONF_NCHANNEL));
  for(int c = 0; c < nchannel; c ++) {
    FP_TYPE fmin = c == 0 ? 0 : chanfreq[c - 1];
    FP_TYPE fmax = c == nchannel - 1 ? fs / 2.0 : chanfreq[c];
    if(fmin >= fs / 2.0) break;
    FP_TYPE* x = llsm_generate_bandlimited_noise(ny, fmin / fs * 2.0,
      fmax / fs * 2.0);
    FP_TYPE* env = llsm_synthesize_noise_envelope(options, src, c, f0, nfrm,
      thop, fs, ny);
    for(int i = 0; i < ny; i ++) {
      x[i] *= sqrt(env[i]);
      y[i] += x[i];
    }
    free(env);
    free(x);
  }
  return y;
}

static FP_TYPE* llsm_filter_noise(llsm_chunk* src, int nfrm, FP_TYPE thop,
  FP_TYPE fs, FP_TYPE* x, int nx) {
  const int nfade = 16;
  int nwin = round(thop * fs * 2);
  FP_TYPE* w = hanning(nwin);
  FP_TYPE wsqr = 0;
  for(int i = 0; i < nwin; i ++)
    wsqr += w[i] * w[i];

  // at least 20% padding
  int nfft = pow(2, ceil(log2(nwin * 1.2 + nfade * 2)));
  int nspec = nfft / 2 + 1;
  FP_TYPE* psd = calloc(nspec, sizeof(FP_TYPE));
  FP_TYPE* fftbuffer = calloc(nfft * 4, sizeof(FP_TYPE));
  FP_TYPE* x_re = fftbuffer;
  FP_TYPE* x_im = fftbuffer + nfft;

  int npsd = *((int*)llsm_container_get(src -> conf, LLSM_CONF_NPSD));
  FP_TYPE fnyq = *((FP_TYPE*)llsm_container_get(src -> conf, LLSM_CONF_FNYQ));
  FP_TYPE noswarp = *((FP_TYPE*)llsm_container_get(src -> conf,
    LLSM_CONF_NOSWARP));
  FP_TYPE* warp_axis = llsm_warp_frequency(0, fnyq, npsd, noswarp);

  FP_TYPE* y = calloc(nx, sizeof(FP_TYPE));
  // STFT -> PSD -> warp -> diff -> filter -> ISTFT
  for(int i = 0; i < nfrm; i ++) {
    llsm_nmframe* nm = llsm_container_get(src -> frames[i], LLSM_FRAME_NM);

    // STFT
    int center = round(i * thop * fs);
    FP_TYPE* xfrm = fetch_frame(x, nx, center, nwin);
    for(int j = 0; j < nwin; j ++) xfrm[j] *= w[j];
    for(int j = 0; j < nfft; j ++) x_re[j] = 0;
    for(int j = 0; j < nwin; j ++) x_re[j - nwin / 2 + nfft / 2] = xfrm[j];
    fft(x_re, NULL, x_re, x_im, nfft, fftbuffer + nfft * 2);
    
    // PSD -> warp -> diff
    llsm_fft_to_psd(x_re, x_im, nfft, wsqr, psd);
    FP_TYPE* env = llsm_spectral_mean(psd, nspec, fs / 2.0, warp_axis, npsd);
    for(int j = 0; j < npsd; j ++)
      env[j] = pow(10.0, nm -> psd[j] / 20.0) / sqrt(env[j] + 1e-8);

    // filter
    FP_TYPE* H = llsm_spectrum_from_envelope(warp_axis, env, npsd, nspec - 1,
      fs / 2.0);
    for(int j = 0; j < nspec - 1; j ++) {
      x_re[j] *= H[j]; x_im[j] *= H[j];
    }
    x_re[nspec - 1] = x_re[nspec - 2]; complete_symm(x_re, nfft);
    x_im[nspec - 1] = x_im[nspec - 2]; complete_asymm(x_im, nfft);

    // ISTFT
    ifft(x_re, x_im, x_re, NULL, nfft, fftbuffer + nfft * 2);
    for(int j = 0; j < nfade; j ++) {
      x_re[j] *= (FP_TYPE)j / nfade;
      x_re[nfft - j - 1] *= 1.0 - (FP_TYPE)j / nfade;
    }
    for(int j = 0; j < nfft; j ++) {
      int idx = center + j - nfft / 2;
      if(idx >= 0 && idx < nx)
        y[idx] += x_re[j];
    }

    free(H); free(env);
    free(xfrm);
  }
  
  free(fftbuffer); free(psd);
  free(warp_axis);
  free(w);
  return y;
}

llsm_output* llsm_synthesize(llsm_soptions* options, llsm_chunk* src) {
  if(! llsm_synthesis_check_integrity(src)) return NULL;
  int nfrm;
  FP_TYPE thop = *((FP_TYPE*)llsm_container_get(src -> conf, LLSM_CONF_THOP));
  FP_TYPE fs = options -> fs;
  FP_TYPE* f0 = llsm_chunk_getf0(src, & nfrm);

  int ny = round((nfrm + 1) * thop * fs);
  llsm_output* ret = malloc(sizeof(llsm_output));
  ret -> ny = ny;
  ret -> fs = fs;

  FP_TYPE* y_sin = llsm_synthesize_harmonics(options, src, f0, nfrm,
    thop, fs, ny);
  ret -> y_sin = y_sin;

  FP_TYPE* y_exc = llsm_synthesize_noise_excitation(options, src, f0, nfrm,
    thop, fs, ny);
  FP_TYPE* y_nos = llsm_filter_noise(src, nfrm, thop, fs, y_exc, ny);
  ret -> y_noise = y_nos;

  ret -> y = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < ny; i ++)
    ret -> y[i] = y_sin[i] + y_nos[i];

  free(y_exc);
  free(f0);
  return ret;
}

void llsm_delete_output(llsm_output* dst) {
  if(dst == NULL) return;
  free(dst -> y);
  free(dst -> y_sin);
  free(dst -> y_noise);
  free(dst);
}

FP_TYPE* llsm_chunk_getf0(llsm_chunk* src, int* dst_nfrm) {
  int* nfrm = llsm_container_get(src -> conf, LLSM_CONF_NFRM);
  if(nfrm == NULL) return NULL;
  FP_TYPE* f0 = calloc(*nfrm, sizeof(FP_TYPE));
  *dst_nfrm = *nfrm;
  for(int i = 0; i < *nfrm; i ++) {
    FP_TYPE* if0 = llsm_container_get(src -> frames[i], LLSM_FRAME_F0);
    if(if0 != NULL)
      f0[i] = if0[0];
  }
  return f0;
}

void llsm_chunk_phasesync_rps(llsm_chunk* dst, int layer1_based) {
  int* nfrm = llsm_container_get(dst -> conf, LLSM_CONF_NFRM);
  if(nfrm == NULL) return;
  for(int i = 0; i < *nfrm; i ++)
    llsm_frame_phasesync_rps(dst -> frames[i], layer1_based);
}

void llsm_chunk_phasepropagate(llsm_chunk* dst, int sign) {
  int nfrm = 0;
  FP_TYPE* f0 = llsm_chunk_getf0(dst, & nfrm);
  FP_TYPE* thop = llsm_container_get(dst -> conf, LLSM_CONF_THOP);
  if(thop == NULL || f0 == NULL) return;
  FP_TYPE* delta_phase = cumsum(f0, nfrm);
  for(int i = 0; i < nfrm; i ++) {
    delta_phase[i] *= *thop * sign * 2.0 * M_PI;
    llsm_frame_phaseshift(dst -> frames[i], delta_phase[i]);
  }
  free(delta_phase);
  free(f0);
}

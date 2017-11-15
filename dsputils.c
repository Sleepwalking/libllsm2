#include "dsputils.h"
#include "external/ciglet/ciglet.h"

#include "filter-coef.h"

static int get_chebyshev_filter(FP_TYPE cutoff, char* type,
  FP_TYPE** dst_a, FP_TYPE** dst_b) {

  int index = max(0, round(cutoff * 2.0 / step_freq - 1));
  if(index >= filter_number) index = filter_number - 1;
  int order = coef_size;
  *dst_a = calloc(order, sizeof(FP_TYPE));
  *dst_b = calloc(order, sizeof(FP_TYPE));
  const FP_TYPE* a_line, *b_line;
  if(! strcmp(type, "lowpass")) {
    a_line = cheby_l_a + index * coef_size;
    b_line = cheby_l_b + index * coef_size;
  } else {
    a_line = cheby_h_a + index * coef_size;
    b_line = cheby_h_b + index * coef_size;
  }
  for(int i = 0; i < order; i ++) {
    (*dst_a)[i] = a_line[i];
    (*dst_b)[i] = b_line[i];
  }
  return order;
}

static FP_TYPE* chebyfilt(FP_TYPE* x, int nx, FP_TYPE c1, FP_TYPE c2) {
  c1 = max(0.0, c1);
  c2 = min(1.0, c2);
  if(c1 != 0 && c2 != 1.0) {
    FP_TYPE* x1 = chebyfilt(x , nx, c1, 1.0);
    FP_TYPE* y  = chebyfilt(x1, nx, 0.0, c2);
    free(x1);
    return y;
  }

  FP_TYPE* a, *b;
  int order = 0;
  if(c1 == 0)
    order = get_chebyshev_filter(c2, "lowpass", & a, & b);
  else
    order = get_chebyshev_filter(c1, "highpass", & a, & b);
  FP_TYPE* y = filtfilt(b, order, a, order, x, nx);
  free(a); free(b);
  return y;
}

void llsm_refine_f0(FP_TYPE* x, int nx, FP_TYPE fs, FP_TYPE* f0, int nfrm,
  FP_TYPE thop) {
  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] == 0) continue;
    FP_TYPE favg = 0;
    for(int j = 1; j <= 3; j ++) { // j-th harmonic
      ifdetector* ifd = create_ifdetector(f0[i] / fs * j, f0[i] / fs);
      FP_TYPE* xfrm = fetch_frame(x, nx, round(i * thop * fs), ifd -> nh);
      FP_TYPE f_j = ifdetector_estimate(ifd, xfrm, ifd -> nh);
      favg += f_j / j / 3.0;
      free(xfrm);
      delete_ifdetector(ifd);
    }
    f0[i] = favg * fs;
  }
}

void llsm_compute_spectrogram(FP_TYPE* x, int nx, int* center, int* winsize,
  int nfrm, int nfft, char* wintype, FP_TYPE** dst_spec, FP_TYPE** dst_phse) {
  // Call stft_forward once just to get the standard normalization factor;
  //   then scale the factor for each frame.
  int standard_winsize = 1024;
  FP_TYPE standard_normalizer = 0;
  cig_stft_forward(x, nx, center, & standard_winsize, 1, nfft, wintype, 0, 2,
    NULL, & standard_normalizer, dst_spec, dst_phse);
  standard_normalizer *= 0.5;

  // The actual STFT analysis.
  cig_stft_forward(x, nx, center, winsize, nfrm, nfft, wintype, 0, 2,
    NULL, NULL, dst_spec, dst_phse);
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE normalizer = standard_winsize / standard_normalizer / winsize[i];
    for(int j = 0; j < nfft / 2 + 1; j ++)
      dst_spec[i][j] *= normalizer;
  }
}

void llsm_compute_dc(FP_TYPE* x, int nx, int* center, int* winsize, int nfrm,
  FP_TYPE* dst_dc) {
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* xfrm = fetch_frame(x, nx, center[i], winsize[i]);
    dst_dc[i] = meanfp(xfrm, winsize[i]);
    free(xfrm);
  }
}

void llsm_harmonic_peakpicking(FP_TYPE* spectrum, FP_TYPE* phase,
  int nfft, FP_TYPE fs, int nhar, FP_TYPE f0,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse) {
  const FP_TYPE tolerance = 0.3;
  for(int i = 1; i <= nhar; i ++) {
    int l_idx = round(f0 * (i - tolerance) / fs * nfft);
    int u_idx = round(f0 * (i + tolerance) / fs * nfft);
    l_idx = max(1, l_idx);
    u_idx = min(nfft / 2 - 1, u_idx);

    int peak_bin = cig_find_peak(spectrum, l_idx, u_idx, 1);
    FP_TYPE peak_freq, peak_ampl;
    peak_ampl = qifft(spectrum, peak_bin, & peak_freq);
    dst_ampl[i - 1] = exp_3(peak_ampl);
    dst_phse[i - 1] = linterp(phase[(int)peak_freq],
      phase[(int)peak_freq + 1], fmod(peak_freq, 1.0));
  }
}

static int f0_to_nhar(FP_TYPE f0, FP_TYPE fs) {
  return floor(fs / f0 / 2);
}

void llsm_harmonic_analysis(FP_TYPE* x, int nx, FP_TYPE fs, FP_TYPE* f0,
  int nfrm, FP_TYPE thop, FP_TYPE rel_winsize, int maxnhar, int method,
  int* dst_nhar, FP_TYPE** dst_ampl, FP_TYPE** dst_phse) {
  int nfft = llsm_get_fftsize(f0, nfrm, fs, rel_winsize);
  int nspec = nfft / 2 + 1;

  int nvfrm = 0; // number of voiced frames
  for(int i = 0; i < nfrm; i ++) nvfrm += f0[i] > 0;
  int* index_vfrm = calloc(nvfrm, sizeof(int));
  int* winsize    = calloc(nvfrm, sizeof(int));
  int* center     = calloc(nvfrm, sizeof(int));
  nvfrm = 0;
  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] > 0) {
      index_vfrm[nvfrm] = i;
      winsize   [nvfrm] = round(fs / f0[i] * rel_winsize / 2) * 2;
      center    [nvfrm] = round(i * thop * fs);
      nvfrm ++;
    }
  }

  FP_TYPE** spec_magn = malloc2d(nvfrm, nspec, sizeof(FP_TYPE));
  FP_TYPE** spec_phse = malloc2d(nvfrm, nspec, sizeof(FP_TYPE));
  llsm_compute_spectrogram(x, nx, center, winsize, nvfrm, nfft, "blackman",
    spec_magn, spec_phse);
  // convert from linear to log magnitude while taking care of underflow
  for(int i = 0; i < nvfrm; i ++)
    for(int j = 0; j < nspec; j ++)
      spec_magn[i][j] = log(spec_magn[i][j] + 1e-8);

  for(int i = 0; i < nvfrm; i ++) {
    int idx = index_vfrm[i];
    dst_nhar[idx] = min(f0_to_nhar(f0[idx], fs), maxnhar);
    dst_ampl[idx] = calloc(dst_nhar[idx], sizeof(FP_TYPE));
    dst_phse[idx] = calloc(dst_nhar[idx], sizeof(FP_TYPE));
    llsm_harmonic_peakpicking(spec_magn[i], spec_phse[i], nfft, fs,
      dst_nhar[idx], f0[idx], dst_ampl[idx], dst_phse[idx]);
  }

  free(index_vfrm); free(winsize); free(center);
  free2d(spec_magn, nvfrm);
  free2d(spec_phse, nvfrm);
}

FP_TYPE* llsm_subband_energy(FP_TYPE* x, int nx, FP_TYPE fmin, FP_TYPE fmax) {
  FP_TYPE* x_filt = chebyfilt(x, nx, fmin, fmax);
  for(int i = 0; i < nx; i ++)
    x_filt[i] *= x_filt[i];
  return x_filt;
}

void llsm_fft_to_psd(FP_TYPE* X_re, FP_TYPE* X_im, int nfft, FP_TYPE wsqr,
  FP_TYPE* dst_psd) {
  // Take the power, normalize.
  for(int i = 0; i < nfft / 2 + 1; i ++) {
    dst_psd[i] = X_re[i] * X_re[i] + X_im[i] * X_im[i];
    dst_psd[i] /= wsqr;
  }
}

void llsm_estimate_psd(FP_TYPE* x, int nx, int nfft, FP_TYPE* dst_psd) {
  // Allocate window and FFT buffer.
  FP_TYPE* window = blackman(nx);
  FP_TYPE* fftbuff = calloc(nfft * 4, sizeof(FP_TYPE));
  FP_TYPE* x_re = fftbuff + 0;
  FP_TYPE* x_im = fftbuff + nfft;

  // Window the signal and compute normalization factor.
  FP_TYPE win_power = 0;
  for(int i = 0; i < nx; i ++) {
    x_re[i] = window[i] * x[i];
    win_power += window[i] * window[i];
  }

  fft(x_re, NULL, x_re, x_im, nfft, fftbuff + nfft * 2);
  llsm_fft_to_psd(x_re, x_im, nfft, win_power, dst_psd);

  free(window);
  free(fftbuff);
}

FP_TYPE* llsm_warp_frequency(FP_TYPE fmin, FP_TYPE fmax, int n,
  FP_TYPE warp_const) {
  FP_TYPE* freq = calloc(n, sizeof(FP_TYPE));
  FP_TYPE wmin = 5000.0 * log(1.0 + fmin / warp_const);
  FP_TYPE wmax = 5000.0 * log(1.0 + fmax / warp_const);
  for(int i = 0; i < n; i ++)
    freq[i] = warp_const * (
      exp(((FP_TYPE)i / n * (wmax - wmin) + wmin) / 5000.0) - 1.0);
  return freq;
}

FP_TYPE* llsm_spectral_mean(FP_TYPE* spectrum, int nspec, FP_TYPE fnyq,
  FP_TYPE* freq, int nfreq) {
  FP_TYPE* env = calloc(nfreq, sizeof(FP_TYPE));
  for(int i = 0; i < nfreq; i ++) {
    FP_TYPE fprev = i == 0 ? 0 : freq[i - 1];
    FP_TYPE fnext = i == nfreq - 1 ? freq[i] * 2 - freq[i - 1] : freq[i + 1];
    FP_TYPE fcenter = freq[i];
    int idxl = floor(fprev / fnyq * nspec);
    int idxc = round(fcenter / fnyq * nspec);
    int idxh = ceil(fnext / fnyq * nspec);
    idxl = max(0, idxl); idxl = min(nspec - 1, idxl);
    idxc = max(0, idxc); idxc = min(nspec - 1, idxc);
    idxh = max(0, idxh); idxh = min(nspec - 1, idxh);
    if(i > 0 && idxh == idxl)
      env[i] = env[i - 1];
    else {
      // Compute the mean weighted by a triangular spectral filter.
      FP_TYPE acc = 0;
      for(int j = idxl; j < idxc; j ++) {
        FP_TYPE coef = 1.0 - (j - idxl) / (idxc - idxl);
        env[i] += spectrum[j] * coef;
        acc += coef;
      }
      for(int j = idxc; j < idxh; j ++) {
        FP_TYPE coef = 1.0 - (idxh - j) / (idxh - idxc);
        env[i] += spectrum[j] * coef;
        acc += coef;
      }
      env[i] /= acc;
    }
  }
  return env;
}

FP_TYPE* llsm_spectrum_from_envelope(FP_TYPE* freq, FP_TYPE* ampl, int nfreq,
  int nspec, FP_TYPE fnyq) {
  FP_TYPE* faxis = calloc(nspec, sizeof(FP_TYPE));
  for(int i = 0; i < nspec; i ++)
    faxis[i] = (FP_TYPE)i * fnyq / nspec;
  FP_TYPE* spectrum = interp1(freq, ampl, nfreq, faxis, nspec);
  free(faxis);
  return spectrum;
}

int llsm_get_fftsize(FP_TYPE* f0, int nfrm, FP_TYPE fs, FP_TYPE rel_winsize) {
  FP_TYPE minf0 = 1000;
  for(int i = 0; i < nfrm; i ++)
    if(f0[i] > 0 && f0[i] < minf0)
      minf0 = f0[i];
  // Window size has to be odd.
  int max_winsize = round(fs / minf0 * rel_winsize / 2) * 2;
  return pow(2, ceil(log2(max_winsize)));
}

FP_TYPE* llsm_synthesize_harmonic_frame(FP_TYPE* ampl, FP_TYPE* phse, int nhar,
  FP_TYPE f0, int nx) {
  FP_TYPE* freq = malloc(nhar * sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++)
    freq[i] = f0 * (i + 1.0);
  FP_TYPE* y = gensins(freq, ampl, phse, nhar, 1.0, nx);
  free(freq);
  return y;
}

FP_TYPE* llsm_generate_white_noise(int nx) {
  FP_TYPE* ret = calloc(nx, sizeof(FP_TYPE));
  int ntemplate = min(20000, nx);
  for(int i = 0; i < ntemplate; i ++)
    ret[i] = randn(0, 1);
  for(int i = ntemplate; i < nx; i ++)
    ret[i] = ret[(i - ntemplate) % ntemplate];
  return ret;
}

static FP_TYPE* stretch_stationary_noise(FP_TYPE* x, int nx, int ny,
  int overlap) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < min(nx, ny); i ++) y[i] = x[i];
  if(ny <= nx) return y;
  int head = nx;
  while(1) {
    for(int i = 0; i < overlap; i ++) {
      FP_TYPE r = (FP_TYPE)i / overlap;
      y[head - overlap + i] *= 1.0 - r;
      y[head - overlap + i] += x[i] * r;
      y[head - overlap + i] /= sqrt(2 * r * (r - 1) + 1);
    }
    for(int i = 0; i < nx - overlap; i ++) {
      if(head + i >= ny) return y;
      y[head + i] = x[i + overlap];
    }
    head += nx - overlap;
  }
  return y;
}

FP_TYPE* llsm_generate_bandlimited_noise(int nx, FP_TYPE fmin, FP_TYPE fmax) {
  int ntemplate = min(20000, nx);
  int extension = 128;
  FP_TYPE* template_white = llsm_generate_white_noise(ntemplate + extension);
  FP_TYPE* template_colored = chebyfilt(template_white, ntemplate + extension,
    fmin, fmax);
  FP_TYPE* y = stretch_stationary_noise(template_colored, ntemplate, nx, 128);
  free(template_white); free(template_colored);
  return y;
}

void llsm_lipradiation(FP_TYPE* freq, int nfreq, FP_TYPE radius,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse) {
  FP_TYPE Rr = 128.0 / 9.0 / M_PI / M_PI;
  FP_TYPE Lr = 8.0 * radius / 100.0 / 3.0 / M_PI / 340.0;
  for(int i = 0; i < nfreq; i ++) {
    cplx iresp = c_mul(c_cplx(0, 1),
      c_div(c_cplx(freq[i] * 2.0 * M_PI * Lr * Rr, 0),
            c_cplx(Rr, freq[i] * 2.0 * M_PI * Lr)
      ));
    if(dst_ampl) dst_ampl[i] = c_abs(iresp);
    if(dst_phse) dst_phse[i] = c_arg(iresp);
  }
}

void llsm_lipfilter(FP_TYPE radius, FP_TYPE f0, int nhar,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse, int inverse) {
  FP_TYPE* tmp = calloc(nhar * 3, sizeof(FP_TYPE));
  FP_TYPE* freq = tmp;
  FP_TYPE* lip_ampl = tmp + nhar;
  FP_TYPE* lip_phse = tmp + nhar * 2;
  for(int i = 0; i < nhar; i ++) freq[i] = f0 * (1.0 + i);
  llsm_lipradiation(freq, nhar, radius, lip_ampl, lip_phse);
  if(inverse)
    for(int i = 0; i < nhar; i ++) {
      if(dst_ampl != NULL) dst_ampl[i] /= lip_ampl[i];
      if(dst_phse != NULL) dst_phse[i] -= lip_phse[i];
    }
  else
    for(int i = 0; i < nhar; i ++) {
      if(dst_ampl != NULL) dst_ampl[i] *= lip_ampl[i];
      if(dst_phse != NULL) dst_phse[i] += lip_phse[i];
    }
  free(tmp);
}

FP_TYPE* llsm_harmonic_spectrum(FP_TYPE* ampl, int nhar, FP_TYPE f0,
  int nfft) {
  FP_TYPE* X = calloc(nfft / 2 + 1, sizeof(FP_TYPE));
  int nX = nfft / 2 + 1;
  int T = 3.0 / f0;
  int width = ceil(f0 * nfft * 1.5);
  for(int i = 0; i < nhar; i ++) {
    FP_TYPE ifreq = f0 * (1.0 + i);
    int center = round(ifreq * nfft);
    for(int j = max(0, center - width);
            j < min(nX, center + width + 1); j ++) {
      FP_TYPE omega = ((FP_TYPE)j / nfft - ifreq) * 2.0 * M_PI;
      FP_TYPE resp = 0.5 * safe_aliased_sinc(T, omega) +
        0.25 * safe_aliased_sinc(T, omega - 2.0 * M_PI / T) +
        0.25 * safe_aliased_sinc(T, omega + 2.0 * M_PI / T);
      X[j] = fmax(X[j], resp * ampl[i]);
    }
  }
  // normalization
  for(int i = 0; i < nfft / 2 + 1; i ++) X[i] *= f0;
  return X;
}

static FP_TYPE compress_logspectrum(FP_TYPE x) {
  if(x > -10) return x;
  return (x + 10.0) / 2 - 10.0;
}

static FP_TYPE decompress_logspectrum(FP_TYPE x) {
  if(x > -10) return x;
  return (x + 10.0) * 2 - 10.0;
}

FP_TYPE* llsm_harmonic_envelope(FP_TYPE* ampl, int nhar, FP_TYPE f0,
  int nfft) {
  FP_TYPE* compressed_ampl = calloc(nhar, sizeof(FP_TYPE));
  FP_TYPE peak = log(maxfp(ampl, nhar));
  for(int i = 0; i < nhar; i ++)
    compressed_ampl[i] = exp(compress_logspectrum(log(ampl[i]) - peak));
  FP_TYPE* X = llsm_harmonic_spectrum(compressed_ampl, nhar, f0, nfft);
  FP_TYPE* full_spectrum = cig_spec2env(X, nfft, f0, nhar, NULL);
  free(X);
  for(int i = 0; i < nfft / 2 + 1; i ++)
    full_spectrum[i] = decompress_logspectrum(full_spectrum[i]) + peak;
  free(compressed_ampl);
  return full_spectrum;
}

FP_TYPE* llsm_harmonic_minphase(FP_TYPE* ampl, int nhar) {
  // Interpolate the harmonics to form a spectral envelope; compute the
  //   minimum-phase response; subsample the phase response at harmonic
  //   frequencies.
  int nfft = max(64, pow(2, ceil(log2(nhar) + 2)));
  FP_TYPE* har_idx  = calloc(nhar + 1, sizeof(FP_TYPE));
  FP_TYPE* har_ampl = calloc(nhar + 1, sizeof(FP_TYPE));
  FP_TYPE* fft_idx  = calloc(nfft / 2 + 1, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) {
    har_idx [i + 1] = (i + 1.0) / (nhar + 1.0) * nfft / 2.0;
    har_ampl[i + 1] = log(ampl[i] + 1e-10);
  }
  har_ampl[0] = har_ampl[1];
  for(int i = 0; i < nfft / 2 + 1; i ++) fft_idx[i] = i;
  FP_TYPE* spectrum = interp1u(0, har_idx[nhar] * 2 - har_idx[nhar - 1],
    har_ampl, nhar + 1, fft_idx, nfft / 2 + 1);
  FP_TYPE* spectrum_phase = minphase(spectrum, nfft);
  FP_TYPE* har_phse = interp1u(0, nfft / 2 + 1,
    spectrum_phase, nfft / 2 + 1, har_idx, nhar + 1);
  for(int i = 1; i < nhar; i ++)
    har_phse[i - 1] = har_phse[i];
  free(har_idx); free(har_ampl); free(fft_idx);
  free(spectrum); free(spectrum_phase);
  return har_phse;
}

typedef struct {
  FP_TYPE** power;     // squared amplitude responses
  FP_TYPE* param;      // the parameter for each response
  int nhar;            // number of harmonics
  int nresp;           // number of cached responses
} cached_glottal_model;

llsm_cached_glottal_model* llsm_create_cached_glottal_model(FP_TYPE* param,
  int nparam, int nhar) {
  cached_glottal_model* ret = malloc(sizeof(cached_glottal_model));
  ret -> nresp = nparam;
  ret -> nhar = nhar;
  ret -> power = calloc(nparam, sizeof(FP_TYPE*));
  ret -> param = calloc(nparam, sizeof(FP_TYPE));
  FP_TYPE f0 = 200.0; // the shape of LF model is f0-independent
  FP_TYPE* freq = calloc(nhar, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) freq[i] = f0 * (1.0 + i);
  for(int i = 0; i < nparam; i ++) {
    ret -> param[i] = param[i];
    lfmodel lf = lfmodel_from_rd(param[i], 1.0 / f0, 1.0);
    ret -> power[i] = lfmodel_spectrum(lf, freq, nhar, NULL);
    for(int j = 0; j < nhar; j ++) {
      ret -> power[i][j] /= j + 1.0;
      ret -> power[i][j] *= ret -> power[i][j];
    }
  }
  free(freq);
  return (llsm_cached_glottal_model*)ret;
}

void llsm_delete_cached_glottal_model(llsm_cached_glottal_model* dst_) {
  if(dst_ == NULL) return;
  cached_glottal_model* dst = (cached_glottal_model*)dst_;
  free2d(dst -> power, dst -> nresp);
  free(dst -> param);
  free(dst);
}

FP_TYPE llsm_spectral_glottal_fitting(FP_TYPE* ampl, int nhar,
  llsm_cached_glottal_model* model_) {
  cached_glottal_model* model = (cached_glottal_model*)model_;
  nhar = min(nhar, model -> nhar);
  FP_TYPE* power = calloc(nhar, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) {
    power[i] = ampl[i] * ampl[i];
  }
  FP_TYPE* distance = calloc(model -> nresp, sizeof(FP_TYPE));
  FP_TYPE* power_model = calloc(nhar, sizeof(FP_TYPE));
  for(int i = 0; i < model -> nresp; i ++) {
    FP_TYPE lgavg = 0;
    for(int j = 0; j < nhar; j ++) lgavg += log(power[j]) / nhar;
    FP_TYPE gain = power[0] / model -> power[i][0];
    for(int j = 0; j < nhar; j ++)
      power_model[j] = model -> power[i][j] * gain;
    distance[i] = exp(itakura_saito(power, power_model, nhar));
  }
  free(power_model);
  free(power);

  int valley = find_minima(distance, 0, model -> nresp - 1);
  FP_TYPE param_refined = model -> param[valley];
  if(valley > 0 && valley < model -> nresp - 1) {
    qifft(distance, valley, & param_refined);
    param_refined = linterp(model -> param[(int)param_refined],
      model -> param[(int)param_refined + 1], fmod(param_refined, 1.0));
  }
  free(distance);
  return param_refined;
}

// https://www.dsprelated.com/showarticle/1068.php
FP_TYPE* llsm_smoothing_filter(FP_TYPE* x, int nx, int order) {
  FP_TYPE* y = calloc(nx, sizeof(FP_TYPE));
  if(nx < order) {
    memcpy(y, x, nx * sizeof(FP_TYPE));
    return y;
  }
  FP_TYPE mean0 = meanfp(x, order);
  FP_TYPE mean1 = meanfp(x + nx - order, order);
  for(int i = 0; i < order / 2; i ++) {
    y[i] = mean0;
    y[nx - i - 1] = mean1;
  }
  for(int i = order / 2; i < nx - order / 2; i ++) {
    int idx_l = i - order / 2;
    int idx_u = idx_l + order;
    FP_TYPE i_mean = meanfp(x + idx_l, order);
    int npos = 0; int nneg = 0;
    FP_TYPE d_total = 0;
    for(int j = idx_l; j < idx_u; j ++) {
      npos += x[j] >= i_mean;
      nneg += x[j] <= i_mean;
      d_total += max(0, x[j] - i_mean);
    }
    y[i] = i_mean + (npos - nneg) * d_total / order / order;
  }
  return y;
}

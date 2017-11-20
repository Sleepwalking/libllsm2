/** @file */

#ifndef LLSM_DSPUTILS_H
#define LLSM_DSPUTILS_H

/** @brief A simple F0 refinment algorithm; overwrites the input. */
void llsm_refine_f0(FP_TYPE* x, int nx, FP_TYPE fs, FP_TYPE* f0, int nfrm,
  FP_TYPE thop);

/** @brief Compute the magnitude and phase spectrogram from a waveform. */
void llsm_compute_spectrogram(FP_TYPE* x, int nx, int* center, int* winsize,
  int nfrm, int nfft, char* wintype, FP_TYPE** dst_spec, FP_TYPE** dst_phse);

/** @brief Compute the short-time mean of a waveform. */
void llsm_compute_dc(FP_TYPE* x, int nx, int* center, int* winsize, int nfrm,
  FP_TYPE* dst_dc);

/** @brief Pick the spectral peak around integer multiples of F0. */
void llsm_harmonic_peakpicking(FP_TYPE* spectrum, FP_TYPE* phase,
  int nfft, FP_TYPE fs, int nhar, FP_TYPE f0,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse);

/** @brief Analyze one frame of a harmonic signal using Chirp-Z Transform. */
void llsm_harmonic_czt(FP_TYPE* x, int nx, FP_TYPE f0, FP_TYPE fs, int nhar,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse);

/** @brief Estimate the amplitude and phase of harmonics on all voiced
 *    frames; allocate memory and store the results into the designated
 *    pointer arrays. */
void llsm_harmonic_analysis(FP_TYPE* x, int nx, FP_TYPE fs, FP_TYPE* f0,
  int nfrm, FP_TYPE thop, FP_TYPE rel_winsize, int maxnhar, int method,
  int* dst_nhar, FP_TYPE** dst_ampl, FP_TYPE** dst_phse);

/** @brief Extract waveform energy from a subband. */
FP_TYPE* llsm_subband_energy(FP_TYPE* x, int nx, FP_TYPE fmin, FP_TYPE fmax);

/** @brief Convert FFT coefficients to power spectral density. */
void llsm_fft_to_psd(FP_TYPE* X_re, FP_TYPE* X_im, int nfft, FP_TYPE wsqr,
  FP_TYPE* dst_psd);

/** @brief Estimate power spectral density for one frame. */
void llsm_estimate_psd(FP_TYPE* x, int nx, int nfft, FP_TYPE* dst_psd);

/** @brief Generate a list of non-uniformly spaced frequency values. */
FP_TYPE* llsm_warp_frequency(FP_TYPE fmin, FP_TYPE fmax, int n,
  FP_TYPE warp_const);

/** @brief Compute the mean of a spectrum in the vicinity of a list of
 *    frequencies. */
FP_TYPE* llsm_spectral_mean(FP_TYPE* spectrum, int nspec, FP_TYPE fnyq,
  FP_TYPE* freq, int nfreq);

/** @brief Interpolate a list of (freq, ampl) points on a uniform axis. */
FP_TYPE* llsm_spectrum_from_envelope(FP_TYPE* freq, FP_TYPE* ampl, int nfreq,
  int nspec, FP_TYPE fnyq);

/** @brief Get the minimum FFT size (power-of-2) for harmonic analysis. */
int llsm_get_fftsize(FP_TYPE* f0, int nfrm, FP_TYPE fs, FP_TYPE rel_winsize);

/** @brief Generate a stationary harmonic signal. */
FP_TYPE* llsm_synthesize_harmonic_frame(FP_TYPE* ampl, FP_TYPE* phse, int nhar,
  FP_TYPE f0, int nx);

/** @brief Generate a Gaussian white noise (mu = 0, sigma = 1) of length nx. */
FP_TYPE* llsm_generate_white_noise(int nx);

/** @brief Generate a bandlimited Gaussian noise of length nx. */
FP_TYPE* llsm_generate_bandlimited_noise(int nx, FP_TYPE fmin, FP_TYPE fmax);

/** @brief Generate a lip radiation frequency response at designated
 *    frequencies. */
void llsm_lipradiation(FP_TYPE* freq, int nfreq, FP_TYPE radius,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse);

/** @brief Apply or remove lip radiation from a harmonic model. */
void llsm_lipfilter(FP_TYPE radius, FP_TYPE f0, int nhar,
  FP_TYPE* dst_ampl, FP_TYPE* dst_phse, int inverse);

/** @brief Generate a 3-period hanning-windowed spectrum from a zero-phase
 *    harmonic model. */
FP_TYPE* llsm_harmonic_spectrum(FP_TYPE* ampl, int nhar, FP_TYPE f0, int nfft);

/** @brief Estimate the spectral envelope of a zero-phase harmonic model. */
FP_TYPE* llsm_harmonic_envelope(FP_TYPE* ampl, int nhar, FP_TYPE f0, int nfft);

/** @brief Recover the phase of a minimum-phase harmonic model from its
 *    amplitude component. */
FP_TYPE* llsm_harmonic_minphase(FP_TYPE* ampl, int nhar);

/** @defgroup group_glottal_analysis Glottal Analysis Utilities
 *  @{ */
typedef void llsm_cached_glottal_model;

/** @brief Create a set of cached glottal model spectra. The implementation is
 *    hidden. */
llsm_cached_glottal_model* llsm_create_cached_glottal_model(FP_TYPE* param,
  int nparam, int nhar);

/** @brief Delete and free the cached glottal model. */
void llsm_delete_cached_glottal_model(llsm_cached_glottal_model* dst);

/** @brief Estimate the glottal model parameter by amplitude-only spectral
 *    fitting. */
FP_TYPE llsm_spectral_glottal_fitting(FP_TYPE* ampl, int nhar,
  llsm_cached_glottal_model* model);
/** @} */

/** @brief An improved moving average filter insensitive to impulse-like
 *    distortions. */
FP_TYPE* llsm_smoothing_filter(FP_TYPE* x, int nx, int order);

#endif

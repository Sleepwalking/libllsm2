#include "../llsm.h"
#include "../dsputils.h"
#include "../external/ciglet/ciglet.h"
#include <assert.h>

#define assert_equal(a, b) \
  assert(approx_equal(a, b))

static int approx_equal(FP_TYPE a, FP_TYPE b) {
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * 1e-6);
}

static FP_TYPE* make_spectrum(int size, int use_db, int no_modulation) {
  FP_TYPE* x = calloc(size, sizeof(FP_TYPE));
  for(int i = 0; i < size; i ++) {
    x[i] = exp(- (FP_TYPE)i / size * 10.0);
    if(! no_modulation)
      x[i] *= exp(cos(i * 0.1) * 3.0 - 3.0);
    if(use_db)
      x[i] = 20.0 * log10(x[i]);
  }
  return x;
}

// compare the difference between the original spectrum and its reconstruction
static void test_spectral_envelope() {
  FP_TYPE fnyq = 22050.0;
  FP_TYPE* spec1 = make_spectrum(1024, 0, 0);
  FP_TYPE* mod1  = make_spectrum(1024, 0, 1);
  FP_TYPE* warp_axis = llsm_warp_frequency(0, fnyq, 100, 15000.0);
  FP_TYPE* env = llsm_spectral_mean(spec1, 1024, fnyq, warp_axis, 100);
  FP_TYPE* spec2 = llsm_spectrum_from_envelope(warp_axis, env, 100,
    1024, fnyq);
  for(int i = 0; i < 1024; i ++) {
    assert(mod1[i] > (spec1[i] - spec2[i]));
    assert(mod1[i] > (spec2[i] - spec1[i]));
  }
  free(env);
  free(warp_axis);
  free(spec1); free(spec2); free(mod1);
}

static void test_harmonic_analysis(int method) {
  // Generate a linear chirp signal for testing.
  int nx = 100000;
  FP_TYPE fs = 20000;
  FP_TYPE f0_start = 100;
  FP_TYPE f0_end = 200;
  FP_TYPE ampl0_start = 0;
  FP_TYPE ampl0_end = 1.0;
  FP_TYPE* x = calloc(nx, sizeof(FP_TYPE));

  FP_TYPE thop = 0.005;
  int nfrm = floor((FP_TYPE)nx / fs / thop);
  FP_TYPE* ampl0_truth = calloc(nfrm, sizeof(FP_TYPE));

  FP_TYPE* f0 = calloc(nfrm, sizeof(FP_TYPE));

  for(int i = 0; i < nfrm; i ++) {
    int center = round(i * thop * fs);
    FP_TYPE rate = (FP_TYPE)center / nx;
    ampl0_truth[i] = ampl0_start + (ampl0_end - ampl0_start) * rate;
    f0[i] = f0_start + (f0_end - f0_start) * rate;
  }

  FP_TYPE phase = 0;
  for(int i = 0; i < nx; i ++) {
    FP_TYPE f0_inst = f0_start + (f0_end - f0_start) * i / nx;
    FP_TYPE ampl0 = ampl0_start + (ampl0_end - ampl0_start) * i / nx;
    phase = wrap(phase + f0_inst / fs * 2.0 * 3.1415927);
    x[i] = ampl0 * sin(phase) + 0.5 * sin(phase * 2.0) + 0.25 * sin(phase * 3.0);
  }
  
  // Perform harmonic analysis with known F0.
  int* nhar = calloc(nfrm, sizeof(FP_TYPE));
  FP_TYPE** ampl = calloc(nfrm, sizeof(FP_TYPE*));
  FP_TYPE** phse = calloc(nfrm, sizeof(FP_TYPE*));

  llsm_harmonic_analysis(x, nx, fs, f0, nfrm, thop, 4.0, 3, method,
    nhar, ampl, phse);

  // Check the errors.
  FP_TYPE* ampl_error = calloc(nfrm, sizeof(FP_TYPE));

  for(int i = 0; i < nfrm; i ++)
    ampl_error[i] = ampl[i][0] - ampl0_truth[i];
  FP_TYPE err_mean = meanfp(ampl_error, nfrm);
  FP_TYPE err_std  = sqrt(varfp(ampl_error, nfrm));
  assert(fabs(err_mean) < 0.01);
  assert(fabs(err_std) < 0.01);
  printf("Harmonic amplitude analysis error (h0): mean = %.2e std = %.2e\n",
    err_mean, err_std);

  for(int i = 0; i < nfrm; i ++)
    ampl_error[i] = ampl[i][1] - 0.5;
  err_mean = meanfp(ampl_error, nfrm);
  err_std  = sqrt(varfp(ampl_error, nfrm));
  assert(fabs(err_mean) < 0.01);
  assert(fabs(err_std) < 0.01);
  printf("Harmonic amplitude analysis error (h1): mean = %.2e std = %.2e\n",
    err_mean, err_std);

  for(int i = 0; i < nfrm; i ++)
    ampl_error[i] = ampl[i][2] - 0.25;
  err_mean = meanfp(ampl_error, nfrm);
  err_std  = sqrt(varfp(ampl_error, nfrm));
  assert(fabs(err_mean) < 0.01);
  assert(fabs(err_std) < 0.01);
  printf("Harmonic amplitude analysis error (h2): mean = %.2e std = %.2e\n",
    err_mean, err_std);

  free(ampl_error);
  
  FP_TYPE* phase_error = calloc(nfrm - 1, sizeof(FP_TYPE));
  for(int i = 1; i < nfrm; i ++) {
    FP_TYPE phase_inc = f0[i] * 2.0 * 3.1415927 * thop;
    phase_error[i - 1] = phase_diff(phse[i - 1][0] + phase_inc, phse[i][0]);
  }
  // Note: a better way is to take the circular mean/variance. However when
  //   error is close to zero this doesn't really matter.
  err_mean = meanfp(phase_error, nfrm - 1);
  err_std = sqrt(varfp(phase_error, nfrm - 1));
  assert(fabs(err_mean) < 0.1);
  assert(fabs(err_std) < 0.1);
  printf("Harmonic phase analysis error (h0): mean = %.2e std = %.2e\n",
    err_mean, err_std);
  free(phase_error);

  free2d(ampl, nfrm); free2d(phse, nfrm); free(nhar);
  free(ampl0_truth); free(f0);
  free(x);
}

static void test_glottal_model() {
  int nhar = 20;
  int nparam = 64;
  FP_TYPE* param_list = linspace(0.02, 3.0, nparam);
  llsm_cached_glottal_model* cgm1 = llsm_create_cached_glottal_model(
    param_list, nparam, nhar);

  FP_TYPE f0 = 200.0;
  FP_TYPE* freq = calloc(nhar, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++)
    freq[i] = f0 * (i + 1.0);

  // Generate a bunch of glottal model spectra, run the fitting algorithm on
  //   them, and compare the estimated value with the ground truth.
  for(int k = 0; k < 500; k ++) {
    FP_TYPE test_param = 0.3 + (2.5 - 0.3) / 500.0 * k;
    FP_TYPE scaling_const = rand() * 5;
    lfmodel lf = lfmodel_from_rd(test_param, 1.0 / f0, 1.0);
    FP_TYPE* ampl = lfmodel_spectrum(lf, freq, nhar, NULL);
    for(int j = 0; j < nhar; j ++) {
      ampl[j] /= j + 1.0;
      ampl[j] *= scaling_const;
    }

    FP_TYPE est_param = llsm_spectral_glottal_fitting(ampl, nhar, cgm1);
    assert(fabs(est_param - test_param) < 0.02);
    free(ampl);
  }
  free(freq);
  llsm_delete_cached_glottal_model(cgm1);
  free(param_list);
}

int main() {
  test_spectral_envelope();
  test_harmonic_analysis(LLSM_AOPTION_HMPP);
  test_glottal_model();
  return 0;
}

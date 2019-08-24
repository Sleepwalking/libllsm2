#include <ciglet/ciglet.h>
#include <assert.h>

// F. Perez-Cruz, "Kullback-Leibler divergence estimation of continuous
//   distributions," 2008 IEEE International Symposium on Information
//   Theory, 2008.
// Estimates the KL divergence between X and its approximation Y.
static inline FP_TYPE empirical_kld(FP_TYPE* x, int nx, FP_TYPE* y, int ny) {
  FP_TYPE* xsorted = sort(x, nx, NULL);
  FP_TYPE* ysorted = sort(y, ny, NULL);
  ysorted[0] = xsorted[0];
  ysorted[ny - 1] = xsorted[nx - 1];
  FP_TYPE D = 0;

  int xi = 1; int yi = 0;
  for(int i = 0; i < nx; i ++) {
    while(yi < ny - 1 && ysorted[yi] < xsorted[i]) yi ++;
    yi = max(yi, 1);
    xi = max(i, 1);
    FP_TYPE d_xi = xsorted[xi] - xsorted[xi - 1];
    FP_TYPE d_yi = ysorted[yi] - ysorted[yi - 1];
    d_xi = max(d_xi, 1e-10);
    d_yi = max(d_yi, 1e-10);
    D += log(ny * d_yi / nx / d_xi);
  }
  free(xsorted); free(ysorted);

  return D / nx - 1.0;
}

// This is a test for the helper function empirical_kld.
static inline void test_empirical_kld() {
  int nx = 100000;
  int ny = 50000;
  FP_TYPE* x_samples = calloc(nx, sizeof(FP_TYPE));
  FP_TYPE* y_samples = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++)
    x_samples[i] = randn(1.0, 1.0);
  for(int i = 0; i < ny; i ++)
    y_samples[i] = randn(1.0, 1.0);

  FP_TYPE kld1 = empirical_kld(x_samples, nx, y_samples, ny);

  for(int i = 0; i < ny; i ++)
    y_samples[i] = randn(0.0, 1.0);

  FP_TYPE kld2 = empirical_kld(x_samples, nx, y_samples, ny);
  assert(kld2 > kld1);

  for(int i = 0; i < ny; i ++)
    y_samples[i] = randn(1.0, 3.0);

  FP_TYPE kld3 = empirical_kld(x_samples, nx, y_samples, ny);
  assert(kld3 > kld1);

  for(int i = 0; i < ny; i ++)
    y_samples[i] = randu();

  FP_TYPE kld4 = empirical_kld(x_samples, nx, y_samples, ny);
  assert(kld4 > kld1);

  printf("\tKL Divergences: %.3f(~0) %.3f(~%.3f) %.3f(~%.3f) %.3f\n",
    kld1,
    kld2, log(1.0 / 1.0) + (1.0 + 1.0) / 2.0 - 0.5,
    kld3, log(sqrt(3.0) / 1.0) + (1.0 + 0.0) / 6.0 - 0.5,
    kld4);

  free(x_samples); free(y_samples);
}

static inline void dither(FP_TYPE* dst, FP_TYPE* src, int nx) {
  for(int i = 0; i < nx; i ++)
    dst[i] = src[i] + randn(0, 1.0) * 1e-4;
}

static inline void verify_data_distribution(FP_TYPE* x, int nx,
  FP_TYPE* x_approx, int ny) {
  FP_TYPE* x_dithered = calloc(nx, sizeof(FP_TYPE));
  FP_TYPE* y_dithered = calloc(ny, sizeof(FP_TYPE));
  // This is to prevent numerical overflow when taking the log of a very small
  //   number in empirical_kld.
  dither(x_dithered, x, nx);
  dither(y_dithered, x_approx, ny);

  FP_TYPE kld = empirical_kld(x_dithered, nx, y_dithered, ny);
  printf("\tKL Divergence between the original waveform distribution and its "
    "reconstruction: %f\n", kld);
  assert(kld < 0.05);

  for(int i = 1; i < nx; i ++)
    x_dithered[i] = x[i] - x[i - 1] + randn(0, 1.0) * 1e-4;
  for(int i = 1; i < ny; i ++)
    y_dithered[i] = x_approx[i] - x_approx[i - 1] + randn(0, 1.0) * 1e-4;
  kld = empirical_kld(x_dithered, nx, y_dithered, ny);
  printf("\tKL Divergence between the original waveform distribution and its "
    "reconstruction (1st derivative): %f\n", kld);
  assert(kld < 0.05);

  for(int i = 1; i < nx - 1; i ++)
    x_dithered[i] = x[i + 1] - 2.0 * x[i] + x[i - 1] + randn(0, 1.0) * 1e-4;
  for(int i = 1; i < ny - 1; i ++)
    y_dithered[i] = x_approx[i + 1] - 2.0 * x_approx[i] + x_approx[i - 1] +
      randn(0, 1.0) * 1e-4;
  kld = empirical_kld(x_dithered, nx, y_dithered, ny);
  printf("\tKL Divergence between the original waveform distribution and its "
    "reconstruction (2nd derivative): %f\n", kld);
  assert(kld < 0.05);

  free(x_dithered); free(y_dithered);
}

static inline FP_TYPE spectral_correlation(FP_TYPE** X, FP_TYPE** Y,
  int m, int n) {
  FP_TYPE* X_flat = flatten(X, m, n, sizeof(FP_TYPE));
  FP_TYPE* Y_flat = flatten(Y, m, n, sizeof(FP_TYPE));
  FP_TYPE cc = corr(X_flat, Y_flat, m * n);
  free(X_flat); free(Y_flat);
  return cc;
}

static inline void verify_spectral_distribution(FP_TYPE* x, int nx,
  FP_TYPE* x_approx, int ny) {
  int nhop = 512;
  int x_nfrm = nx / nhop;
  int y_nfrm = ny / nhop;
  int min_nfrm = min(x_nfrm, y_nfrm);
  int hop_fc = 4; int zp_fc = 1;
  int nfft = nhop * hop_fc;
  FP_TYPE** X = malloc2d(x_nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  FP_TYPE** Y = malloc2d(y_nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  stft(x       , nx, nhop, x_nfrm, hop_fc, zp_fc, NULL, NULL, X, NULL);
  stft(x_approx, ny, nhop, y_nfrm, hop_fc, zp_fc, NULL, NULL, Y, NULL);

  FP_TYPE* X_flat = flatten(X, x_nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  FP_TYPE* Y_flat = flatten(Y, y_nfrm, nfft / 2 + 1, sizeof(FP_TYPE));
  
  FP_TYPE cc = spectral_correlation(X, Y, min_nfrm, nfft / 2 + 1);
  printf("\tCorrelation coefficient between the original spectrum "
    "and its reconstruction: %f\n", cc);
  assert(cc > 0.95);

  dither(X_flat, X_flat, x_nfrm * (nfft / 2 + 1));
  dither(Y_flat, Y_flat, y_nfrm * (nfft / 2 + 1));
  FP_TYPE kld = empirical_kld(X_flat, x_nfrm * (nfft / 2 + 1), Y_flat,
    y_nfrm * (nfft / 2 + 1));
  printf("\tKL Divergence between the original spectral distribution and its "
    "reconstruction: %f\n", kld);
  free(X_flat); free(Y_flat);
  assert(kld < 0.05);

  // 1st derivative
  for(int j = 0; j < nfft / 2 + 1; j ++)
    for(int i = 0; i < x_nfrm - 1; i ++)
      X[i][j] = X[i + 1][j] - X[i][j];
  for(int j = 0; j < nfft / 2 + 1; j ++)
    for(int i = 0; i < y_nfrm - 1; i ++)
      Y[i][j] = Y[i + 1][j] - Y[i][j];

  X_flat = flatten(X, x_nfrm - 1, nfft / 2 + 1, sizeof(FP_TYPE));
  Y_flat = flatten(Y, y_nfrm - 1, nfft / 2 + 1, sizeof(FP_TYPE));
  dither(X_flat, X_flat, (x_nfrm - 1) * (nfft / 2 + 1));
  dither(Y_flat, Y_flat, (y_nfrm - 1) * (nfft / 2 + 1));
  kld = empirical_kld(X_flat, (x_nfrm - 1) * (nfft / 2 + 1), Y_flat,
    (y_nfrm - 1) * (nfft / 2 + 1));
  printf("\tKL Divergence between the original spectral distribution and its "
    "reconstruction (1st derivative): %f\n", kld);
  assert(kld < 0.05);

  free(X_flat); free(Y_flat);
  free2d(X, x_nfrm); free2d(Y, y_nfrm);
}

#include "../llsm.h"
#include "../dsputils.h"
#include "../external/ciglet/ciglet.h"
#include "verify-utils.h"

#include <sys/time.h>
static double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + (t.tv_usec / 1000000.0)) * 1000.0;
}

static FP_TYPE compare_runtime(FP_TYPE* ampl, FP_TYPE* phse,
  int nhar, int nx, int niter) {
  double t0 = get_time();
  for(int i = 0; i < niter; i ++)
    free(llsm_synthesize_harmonic_frame(ampl, phse,
      nhar, 0.001, nx));
  double t1 = get_time();
  double t_gensins = t1 - t0;

  t0 = get_time();
  for(int i = 0; i < niter; i ++)
    free(llsm_synthesize_harmonic_frame_iczt(ampl, phse,
      nhar, 0.001, nx));
  t1 = get_time();
  double t_iczt = t1 - t0;
  return t_iczt / t_gensins;
}

int main(int argc, char** argv) {
  FP_TYPE ampl[1000];
  FP_TYPE phse[1000];
  for(int i = 0; i < 1000; i ++) {
    ampl[i] = randn(0, 1.0);
    phse[i] = randn(0, 100);
  }
  
  FP_TYPE* y1 = llsm_synthesize_harmonic_frame_iczt(ampl, phse,
    100, 0.01, 1024);
  FP_TYPE* y2 = llsm_synthesize_harmonic_frame(ampl, phse,
    100, 0.01, 1024);
  for(int i = 0; i < 1024; i ++) y1[i] -= y2[i];
  FP_TYPE std = sqrt(varfp(y1, 1024));
  FP_TYPE energy = sqrt(varfp(y2, 1024));
  fprintf(stderr, "SNR = %f\n", 20.0 * log10(std / energy));
  free(y1); free(y2);

  // On Interl i7-7500U, it was found that CZT is as fast as cig_gensins when
  //   log(nhar) = log(nx) * 0.275 + 2.26.
  // That means if log(nx) * 0.275 + 2.26 < log(nhar), then CZT is faster.
  // The parameters does not seem to be strongly dependent on platform, though
  //   it hasn't been extensively tested.
  if(argc > 1 && (! strcmp(argv[1], "profile"))) {
    // First load the stuffs into the cache.
    for(int i = 0; i < 1000; i ++) {
      free(llsm_synthesize_harmonic_frame(ampl, phse,
        100, 0.001, 1024));
    }

    FP_TYPE nx_min = log(10);
    FP_TYPE nx_max = log(2048);
    FP_TYPE nhar_min = log(10);
    FP_TYPE nhar_max = log(512);
    int nstep = 20;
    for(int i = 0; i < nstep; i ++)
      printf("%f,", round(exp(linterp(nx_min, nx_max, (FP_TYPE)i / nstep))));
    printf("\n");
    for(int i = 0; i < nstep; i ++)
      printf("%f,", round(exp(linterp(nhar_min, nhar_max,
        (FP_TYPE)i / nstep))));
    printf("\n");
    for(int i = 0; i < nstep; i ++) {
      for(int j = 0; j < nstep; j ++) {
        int nx = round(exp(linterp(nx_min, nx_max, (FP_TYPE)i / nstep)));
        int nhar = round(exp(linterp(nhar_min, nhar_max, (FP_TYPE)j / nstep)));
        FP_TYPE ratio = compare_runtime(ampl, phse, nhar, nx, 200);
        printf("%f,", ratio);
      }
      printf("\n");
    }
  }

  return 0;
}

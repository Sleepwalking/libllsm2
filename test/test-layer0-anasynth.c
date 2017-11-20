#include "../llsm.h"
#include "../external/libpyin/pyin.h"
#include "verify-utils.h"

#include <sys/time.h>
static double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + (t.tv_usec / 1000000.0)) * 1000.0;
}

int main(int argc, char** argv) {
  int fs = 0;
  int nbit = 0;
  int nx = 0;
  FP_TYPE* x = wavread("test/arctic_a0001.wav", & fs, & nbit, & nx);

  int nhop = 128;
  int nfrm = 0;
  pyin_config param = pyin_init(nhop);
  param.fmin = 50.0;
  param.fmax = 500.0;
  param.trange = 24;
  param.bias = 2;
  param.nf = ceil(fs * 0.025);
  FP_TYPE* f0 = pyin_analyze(param, x, nx, fs, & nfrm);

  llsm_aoptions* opt_a = llsm_create_aoptions();
  opt_a -> thop = (FP_TYPE)nhop / fs;
  if(argc > 1 && (! strcmp(argv[1], "czt")))
    opt_a -> hm_method = LLSM_AOPTION_HMCZT;
  llsm_soptions* opt_s = llsm_create_soptions(fs);
  llsm_chunk* chunk = llsm_analyze(opt_a, x, nx, fs, f0, nfrm, NULL);

  double t0 = get_time();
  llsm_output* out = llsm_synthesize(opt_s, chunk);
  double t1 = get_time();
  wavwrite(out -> y_noise, out -> ny, opt_s -> fs, 24,
    "test/test-layer0-noise.wav");
  wavwrite(out -> y_sin  , out -> ny, opt_s -> fs, 24,
    "test/test-layer0-sin.wav");
  wavwrite(out -> y      , out -> ny, opt_s -> fs, 24,
    "test/test-layer0.wav");
  printf("Synthesize speed: %fx real-time.\n", 1000.0 / (t1 - t0));

  verify_data_distribution(x, nx, out -> y, out -> ny);
  verify_spectral_distribution(x, nx, out -> y, out -> ny);
  llsm_delete_output(out);

  llsm_chunk_phasesync_rps(chunk, 0);
  llsm_chunk_phasepropagate(chunk, 1);

  out = llsm_synthesize(opt_s, chunk);
  wavwrite(out -> y_noise, out -> ny, opt_s -> fs, 24,
    "test/test-layer0-phasesync-noise.wav");
  wavwrite(out -> y_sin  , out -> ny, opt_s -> fs, 24,
    "test/test-layer0-phasesync-sin.wav");
  wavwrite(out -> y      , out -> ny, opt_s -> fs, 24,
    "test/test-layer0-phasesync.wav");

  verify_data_distribution(x, nx, out -> y, out -> ny);
  verify_spectral_distribution(x, nx, out -> y, out -> ny);
  llsm_delete_output(out);

  llsm_delete_chunk(chunk);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0);
  free(x);
  return 0;
}

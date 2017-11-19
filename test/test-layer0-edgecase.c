#include "../llsm.h"
#include "verify-utils.h"

int main() {
  int fs = 0;
  int nbit = 0;
  int nx = 0;
  FP_TYPE* x = wavread("test/arctic_a0001.wav", & fs, & nbit, & nx);

  FP_TYPE nhop = 100.5;
  int nfrm = nx / nhop;
  // Make all frames unvoiced. Make sure it won't crash.
  FP_TYPE* f0 = calloc(nfrm, sizeof(FP_TYPE));

  llsm_aoptions* opt_a = llsm_create_aoptions();
  opt_a -> thop = nhop / fs;
  llsm_soptions* opt_s = llsm_create_soptions(fs);
  llsm_chunk* chunk = llsm_analyze(opt_a, x, nx, fs, f0, nfrm, NULL);

  llsm_output* out = llsm_synthesize(opt_s, chunk);
  wavwrite(out -> y_noise, out -> ny, opt_s -> fs, 24,
    "test/test-layer0-edgecase-noise.wav");
  wavwrite(out -> y_sin  , out -> ny, opt_s -> fs, 24,
    "test/test-layer0-edgecase-sin.wav");
  wavwrite(out -> y      , out -> ny, opt_s -> fs, 24,
    "test/test-layer0-edgecase.wav");

  for(int i = 0; i < nfrm; i ++)
    assert(out -> y_sin[i] == 0);

  llsm_delete_output(out);
  llsm_delete_chunk(chunk);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0);
  free(x);
  return 0;
}

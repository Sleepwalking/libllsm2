#include "../llsm.h"
#include "../external/libpyin/pyin.h"
#include "verify-utils.h"

int main() {
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
  llsm_soptions* opt_s = llsm_create_soptions(fs);
  llsm_chunk* chunk = llsm_analyze(opt_a, x, nx, fs, f0, nfrm, NULL);

  llsm_chunk_tolayer1(chunk, 2048);

  llsm_chunk* reconstructed = llsm_create_chunk(chunk -> conf, 0);
  llsm_coder* coder = llsm_create_coder(chunk -> conf, 64, 5);
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* enc = llsm_coder_encode(coder, chunk -> frames[i]);
    reconstructed -> frames[i] = llsm_coder_decode_layer1(coder, enc);
    free(enc);
  }
  llsm_delete_coder(coder);

  llsm_chunk_tolayer0(reconstructed);
  llsm_chunk_phasepropagate(reconstructed, 1);

  llsm_output* out = llsm_synthesize(opt_s, reconstructed);
  wavwrite(out -> y, out -> ny, opt_s -> fs, 24, "test/test-coder.wav");
  printf("Checking the reconstruction against the input...\n");
  verify_data_distribution(x, nx, out -> y, out -> ny);
  llsm_delete_output(out);

  llsm_delete_chunk(reconstructed);

  llsm_delete_chunk(chunk);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0);
  free(x);
  return 0;
}

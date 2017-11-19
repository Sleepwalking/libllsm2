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
  llsm_chunk_tolayer0(chunk);

  llsm_output* out = llsm_synthesize(opt_s, chunk);
  wavwrite(out -> y_noise, out -> ny, opt_s -> fs, 24,
    "test/test-layer1-noise.wav");
  wavwrite(out -> y_sin  , out -> ny, opt_s -> fs, 24,
    "test/test-layer1-sin.wav");
  wavwrite(out -> y      , out -> ny, opt_s -> fs, 24,
    "test/test-layer1.wav");

  verify_data_distribution(x, nx, out -> y, out -> ny);
  verify_spectral_distribution(x, nx, out -> y, out -> ny);
  llsm_delete_output(out);

  // Shift pitch by 1.5x.
  // Alternatively, use llsm_chunk_phasesync_rps.
  llsm_chunk_phasepropagate(chunk, -1);
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* f0_i = llsm_container_get(chunk -> frames[i], LLSM_FRAME_F0);
    f0_i[0] *= 1.5;
    // Compensate for the amplitude gain.
    FP_TYPE* vt_magn = llsm_container_get(chunk -> frames[i],
      LLSM_FRAME_VTMAGN);
    if(vt_magn != NULL) {
      int nspec = llsm_fparray_length(vt_magn);
      for(int j = 0; j < nspec; j ++)
        vt_magn[j] -= log(1.5);
    }
  }
  llsm_chunk_tolayer0(chunk);
  llsm_chunk_phasepropagate(chunk, 1);

  out = llsm_synthesize(opt_s, chunk);
  wavwrite(out -> y_noise, out -> ny, opt_s -> fs, 24,
    "test/test-layer1-pitchshift-noise.wav");
  wavwrite(out -> y_sin  , out -> ny, opt_s -> fs, 24,
    "test/test-layer1-pitchshift-sin.wav");
  wavwrite(out -> y      , out -> ny, opt_s -> fs, 24,
    "test/test-layer1-pitchshift.wav");
  llsm_delete_output(out);

  llsm_delete_chunk(chunk);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0);
  free(x);
  return 0;
}

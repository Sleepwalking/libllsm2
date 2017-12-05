#include "../llsmrt.h"
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
  opt_a -> npsd = 128;
  opt_a -> maxnhar = 400;
  opt_a -> maxnhar_e = 5;
  llsm_soptions* opt_s = llsm_create_soptions(fs);
  llsm_chunk* chunk = llsm_analyze(opt_a, x, nx, fs, f0, nfrm, NULL);

  llsm_rtsynth_buffer* rtbuffer = llsm_create_rtsynth_buffer(opt_s,
    chunk -> conf, 10, 4096);
  int latency = llsm_rtsynth_buffer_getlatency(rtbuffer);
  FP_TYPE* y = calloc(nx, sizeof(FP_TYPE));

  double t0 = get_time();
  int count = 0;
  for(int i = 0; i < nfrm; i ++) {
    llsm_rtsynth_buffer_feed(rtbuffer, chunk -> frames[i]);
    while(count < nx) {
      FP_TYPE tmp = 0;
      int status = llsm_rtsynth_buffer_fetch(rtbuffer, & tmp);
      if(status == 0) break;
      if(count >= latency) {
        y[count - latency] = tmp;
      }
      count ++;
    }
  }
  llsm_delete_rtsynth_buffer(rtbuffer);

  double t1 = get_time();
  wavwrite(y, nx, opt_s -> fs, 24, "test/test-llsmrt.wav");
  printf("Synthesis speed (llsmrt): %f ms, %fx real-time.\n", t1 - t0,
    1000.0 / (t1 - t0) * ((FP_TYPE)nx / opt_s -> fs));

  t0 = get_time();
  llsm_output* out = llsm_synthesize(opt_s, chunk);
  t1 = get_time();
  printf("Synthesis speed (llsm): %f ms, %fx real-time.\n", t1 - t0,
    1000.0 / (t1 - t0) * ((FP_TYPE)nx / opt_s -> fs));

  verify_data_distribution(out -> y, out -> ny, y, nx - latency);
  verify_spectral_distribution(out -> y, out -> ny, y, nx - latency);

  llsm_delete_output(out);
  free(y);

  llsm_delete_chunk(chunk);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0);
  free(x);
  return 0;
}

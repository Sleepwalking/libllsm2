#define _DEFAULT_SOURCE

#include "../llsmrt.h"
#include "../external/libpyin/pyin.h"
#include "verify-utils.h"

#ifdef USE_PTHREAD
#include <pthread.h>
#include <unistd.h>
#endif

#include <sys/time.h>
static double get_time() {
  struct timeval t;
  gettimeofday(&t, NULL);
  return (t.tv_sec + (t.tv_usec / 1000000.0)) * 1000.0;
}

llsm_rtsynth_buffer* rtbuffer = NULL;
int nfrm = 0;
int nx = 0;
int latency = 0;
llsm_chunk* chunk = NULL;
FP_TYPE* y = NULL;
int synth_finish = 0;

#ifdef USE_PTHREAD
static void* synthesis_thread(void* ptr) {
  printf("Synthesis thread: starting.\n");
  for(int i = 0; i < nfrm; i ++) {
    // llsm_rtsynth_buffer_feed blocks when the buffer is full.
    llsm_rtsynth_buffer_feed(rtbuffer, chunk -> frames[i]);
  }
  printf("Synthesis thread: task completed.\n");
  synth_finish = 1;
  return NULL;
}

static void* reading_thread(void* ptr) {
  printf("Reading thread: starting.\n");
  int count = 0;
  while(count < nx) {
    FP_TYPE tmp = 0;
    int status = llsm_rtsynth_buffer_fetch(rtbuffer, & tmp);
    if(synth_finish) break;
    if(status == 0) {
      // llsm_rtsynth_buffer_fetch does not block.
      // So give it a bit of time waiting for more samples to come.
      usleep(10);
      continue;
    }
    if(count >= latency) {
      y[count - latency] = tmp;
    }
    count ++;
  }
  printf("Reading thread: task completed.\n");
  return NULL;
}
#endif

int main(int argc, char** argv) {
  int fs = 0;
  int nbit = 0;
  FP_TYPE* x = wavread("test/arctic_a0001.wav", & fs, & nbit, & nx);

  int nhop = 128;
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
  chunk = llsm_analyze(opt_a, x, nx, fs, f0, nfrm, NULL);
  llsm_chunk_tolayer1(chunk, 2048);
  llsm_chunk_phasesync_rps(chunk, 1);
  opt_s -> use_l1 = 1;
  
  // Keep switching between pulse-by-pulse and harmonic synthesis modes.
  for(int i = 0; i < nfrm; i ++) {
    llsm_container_attach(chunk -> frames[i], LLSM_FRAME_HM, NULL, NULL, NULL);
    /*
    // Enable this to test real-time switching between different synthesis
    //   modes on different pitches.
    FP_TYPE* f0_i = llsm_container_get(chunk -> frames[i], LLSM_FRAME_F0);
    f0_i[0] *= 1.5;
    // Compensate for the amplitude gain.
    FP_TYPE* vt_magn = llsm_container_get(chunk -> frames[i],
      LLSM_FRAME_VTMAGN);
    if(vt_magn != NULL) {
      int nspec = llsm_fparray_length(vt_magn);
      for(int j = 0; j < nspec; j ++)
        vt_magn[j] -= 20.0 * log10(1.5);
    }*/
    if(i % 100 > 50) {
      llsm_container_attach(chunk -> frames[i], LLSM_FRAME_PBPSYN,
        llsm_create_int(1), llsm_delete_int, llsm_copy_int);
    }
  }
  llsm_chunk_phasepropagate(chunk, 1);

  rtbuffer = llsm_create_rtsynth_buffer(opt_s, chunk -> conf, 4096);
  latency = llsm_rtsynth_buffer_getlatency(rtbuffer);
  y = calloc(nx, sizeof(FP_TYPE));

  double t0 = get_time();
# ifdef USE_PTHREAD
  printf("Main thread: launching threads.\n");
  pthread_t threads[2];
  pthread_create(& threads[0], NULL, & synthesis_thread, NULL);
  pthread_create(& threads[1], NULL, & reading_thread, NULL);
  
  pthread_join(threads[0], NULL);
  pthread_join(threads[1], NULL);
# else
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
# endif
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

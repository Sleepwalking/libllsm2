#include "../llsm.h"
#include "../buffer.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define to_fp(x) ((FP_TYPE*)x)

#define assert_equal(a, b) \
  assert(approx_equal(a, b))

static int approx_equal(FP_TYPE a, FP_TYPE b) {
  return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * 1e-6);
}

void test_container() {
  // test: creation and attach
  llsm_container* c1 = llsm_create_container(10);
  llsm_container_attach(c1, 0 , llsm_create_fp(5.0) , free, llsm_copy_fp);
  llsm_container_attach(c1, 1 , llsm_create_fp(10.0), NULL, llsm_copy_fp);
  assert_equal(to_fp(c1 -> members[0])[0], 5.0);
  assert_equal(to_fp(c1 -> members[1])[0], 10.0);

  llsm_container_attach(c1, 15, llsm_create_fp(50.0), free, NULL);
  assert(c1 -> nmember >= 16);
  assert_equal(to_fp(c1 -> members[15])[0], 50.0);

  // test: copy-constructing
  llsm_container* c2 = llsm_copy_container(c1);
  assert_equal(to_fp(c2 -> members[0])[0], 5.0);
  assert_equal(to_fp(c2 -> members[1])[0], 10.0);
  assert_equal(to_fp(c2 -> members[15])[0], 50.0);
  to_fp(c2 -> members[15])[0] = 45.0;
  assert_equal(to_fp(c1 -> members[15])[0], 45.0);

  llsm_container* c3 = llsm_create_container(5);
  llsm_container_attach(c3, 0, llsm_create_fp(-5.0), free, llsm_copy_fp);
  llsm_copy_container_inplace(c3, c2);
  assert_equal(to_fp(c3 -> members[0])[0], 5.0);
  assert_equal(to_fp(c3 -> members[1])[0], 10.0);
  assert_equal(to_fp(c3 -> members[15])[0], 45.0);
  to_fp(c3 -> members[15])[0] = 50.0;
  assert_equal(to_fp(c1 -> members[15])[0], 50.0);

  // test: removal
  llsm_container_remove(c1, 0);
  assert(c1 -> members[0] == NULL);
  assert_equal(to_fp(c2 -> members[0])[0], 5.0);

  free(c1 -> members[1]);
  free(c2 -> members[1]);
  free(c3 -> members[1]);
  llsm_delete_container(c1);
  llsm_delete_container(c2);
  llsm_delete_container(c3);
}

// Use valgrind to check if there's any memory leak.
void test_container_fparray() {
  llsm_container* c1 = llsm_create_container(1);
  llsm_container_attach(c1, 0, llsm_create_fparray(0),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c1, 0, llsm_create_fparray(1),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c1, 1, llsm_create_fparray(1),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c1, 2, llsm_create_fparray(0),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c1, 3, llsm_create_fparray(1),
    llsm_delete_fparray, llsm_copy_fparray);
  
  llsm_container* c2 = llsm_create_container(1);
  llsm_copy_container_inplace(c2, c1);
  llsm_container_attach(c2, 0, llsm_create_fparray(1),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c2, 1, llsm_create_fparray(0),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c2, 2, llsm_create_fparray(1),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(c2, 3, llsm_create_fparray(0),
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_delete_container(c1);
  llsm_delete_container(c2);
}

void test_hmframe() {
  // test: creation
  llsm_hmframe* h1 = llsm_create_hmframe(3);
  h1 -> ampl[0] = 1.0; h1 -> phse[0] = 1.0;
  h1 -> ampl[1] = 0.5; h1 -> phse[1] = -0.5;
  h1 -> ampl[2] = 0.2; h1 -> phse[2] = 2.5;

  // test: copy-constructing
  llsm_hmframe* h2 = llsm_copy_hmframe(h1);
  assert(h2 -> nhar == 3);
  assert_equal(h2 -> ampl[0], 1.0);
  assert_equal(h2 -> phse[0], 1.0);
  assert_equal(h2 -> ampl[2], 0.2);
  assert_equal(h2 -> phse[2], 2.5);

  // test: rotate and rotate back
  llsm_hmframe_phaseshift(h2, 3.14);
  llsm_hmframe_phaseshift(h2, 3.14);
  llsm_hmframe_phaseshift(h2, -6.28);
  assert_equal(h2 -> phse[0], 1.0);
  assert_equal(h2 -> phse[1], -0.5);
  assert_equal(h2 -> phse[2], 2.5);

  llsm_delete_hmframe(h1);
  llsm_delete_hmframe(h2);
}

void test_nmframe() {
  // test: creation
  llsm_nmframe* n1 = llsm_create_nmframe(3, 2, 20);
  for(int i = 0; i < 20; i ++)
    n1 -> psd[i] = i - 10.0;
  for(int i = 0; i < 3; i ++) {
    n1 -> edc[i] = i * 0.1;
    n1 -> eenv[i] -> ampl[0] = 1.0;
    n1 -> eenv[i] -> ampl[1] = 0.5;
  }

  // test: copy-constructing
  llsm_nmframe* n2 = llsm_copy_nmframe(n1);
  assert(n2 -> npsd == 20);
  assert(n2 -> nchannel == 3);
  for(int i = 0; i < 20; i ++)
    assert_equal(n2 -> psd[i], i - 10.0);
  for(int i = 0; i < 3; i ++) {
    assert_equal(n2 -> edc[i], i * 0.1);
    assert_equal(n2 -> eenv[i] -> ampl[0], 1.0);
    assert_equal(n2 -> eenv[i] -> ampl[1], 0.5);
  }

  llsm_delete_nmframe(n1);
  llsm_delete_nmframe(n2);
}

void test_chunk() {
  llsm_aoptions* opt = llsm_create_aoptions();
  llsm_container* conf = llsm_aoptions_toconf(opt, 22050);
  int npsd = ((int*)llsm_container_get(conf, LLSM_CONF_NPSD))[0];
  ((int*)llsm_container_get(conf, LLSM_CONF_NFRM))[0] = 100;
  llsm_chunk* chunk1 = llsm_create_chunk(conf, 1);
  for(int i = 0; i < 100; i ++) {
    llsm_container* currframe = chunk1 -> frames[i];
    llsm_nmframe* nm = llsm_container_get(currframe, LLSM_FRAME_NM);
    for(int j = 0; j < npsd; j ++)
      nm -> psd[j] = sin(j * 0.1);
  }

  llsm_chunk* chunk2 = llsm_copy_chunk(chunk1);
  for(int i = 0; i < 100; i ++) {
    llsm_container* currframe = chunk2 -> frames[i];
    llsm_nmframe* nm = llsm_container_get(currframe, LLSM_FRAME_NM);
    for(int j = 0; j < npsd; j ++)
      assert_equal(nm -> psd[j], sin(j * 0.1));
  }

  llsm_delete_chunk(chunk1);
  llsm_delete_chunk(chunk2);
  llsm_delete_container(conf);
  llsm_delete_aoptions(opt);
}

void test_buffer() {
  llsm_ringbuffer* rb1 = llsm_create_ringbuffer(4096);
  FP_TYPE x[100];
  FP_TYPE y[200];
  for(int i = 0; i < 100; i ++)
    x[i] = (FP_TYPE)rand() / RAND_MAX - 0.5;
  for(int i = 0; i < 100; i ++) {
    // The following two realizations of appending 100 samples should
    //   give equivalent results.
    if(i % 2 == 0) {
      llsm_ringbuffer_appendchunk(rb1, 100, x);
    } else {
      llsm_ringbuffer_forward(rb1, 100);
      llsm_ringbuffer_writechunk(rb1, -100, 100, x);
    }
    // The following two realizations of appending 100 zeros should
    //   give equivalent results.
    if(i % 3 == 0) {
      llsm_ringbuffer_appendblank(rb1, 100);
    } else {
      for(int j = 0; j < 100; j ++)
        llsm_ringbuffer_append(rb1, 0);
    }
    // The following two realizations of reading 200 samples should
    //   give equivalent results.
    if(i % 4 == 0) {
      for(int j = 0; j < 200; j ++)
        y[j] = llsm_ringbuffer_read(rb1, j - 200);
    } else {
      llsm_ringbuffer_readchunk(rb1, -200, 200, y);
    }

    for(int j = 0; j < 100; j ++)
      assert(y[j] == x[j]);
    for(int j = 100; j < 200; j ++)
      assert(y[j] == 0);

    if(i <= 5) continue;

    llsm_ringbuffer_readchunk(rb1, -1300, 200, y);
    for(int j = 0; j < 100; j ++)
      assert(y[j] == 0);
    for(int j = 100; j < 200; j ++)
      assert(y[j] == x[j - 100]);
  }
  llsm_delete_ringbuffer(rb1);
}

int main() {
  test_container();
  test_container_fparray();
  test_hmframe();
  test_nmframe();
  test_chunk();
  test_buffer();
  return 0;
}

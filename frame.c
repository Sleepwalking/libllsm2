#include "llsm.h"
#include "external/ciglet/ciglet.h"

llsm_hmframe* llsm_create_hmframe(int nhar) {
  llsm_hmframe* ret = malloc(sizeof(llsm_hmframe));
  ret -> ampl = calloc(nhar, sizeof(FP_TYPE));
  ret -> phse = calloc(nhar, sizeof(FP_TYPE));
  ret -> nhar = nhar;
  return ret;
}

llsm_hmframe* llsm_copy_hmframe(llsm_hmframe* src) {
  llsm_hmframe* ret = llsm_create_hmframe(src -> nhar);
  llsm_copy_hmframe_inplace(ret, src);
  return ret;
}

void llsm_copy_hmframe_inplace(llsm_hmframe* dst, llsm_hmframe* src) {
  int memsize = sizeof(FP_TYPE) * src -> nhar;
  if(dst -> nhar < src -> nhar) {
    dst -> ampl = realloc(dst -> ampl, memsize);
    dst -> phse = realloc(dst -> phse, memsize);
  }
  memcpy(dst -> ampl, src -> ampl, memsize);
  memcpy(dst -> phse, src -> phse, memsize);
  dst -> nhar = src -> nhar;
}

void llsm_delete_hmframe(llsm_hmframe* dst) {
  if(dst == NULL) return;
  free(dst -> ampl);
  free(dst -> phse);
  free(dst);
}

void llsm_hmframe_phaseshift(llsm_hmframe* dst, FP_TYPE theta) {
  for(int i = 0; i < dst -> nhar; i ++)
    dst -> phse[i] = wrap(dst -> phse[i] + theta * (i + 1.0));
}

llsm_nmframe* llsm_create_nmframe(int nchannel, int nhar_e, int npsd) {
  llsm_nmframe* ret = malloc(sizeof(llsm_nmframe));
  ret -> eenv = calloc(nchannel, sizeof(llsm_hmframe*));
  ret -> edc = calloc(nchannel, sizeof(FP_TYPE));
  ret -> psd = calloc(npsd, sizeof(FP_TYPE));
  ret -> npsd = npsd;
  ret -> nchannel = nchannel;

  // set the default PSD to -120 dB.
  for(int i = 0; i < npsd; i ++)
    ret -> psd[i] = -120.0;

  for(int i = 0; i < nchannel; i ++) {
    ret -> eenv[i] = llsm_create_hmframe(nhar_e);

    // set the default edc to 1e-5.
    ret -> edc[i] = 1e-5;
  }
  return ret;
}

void llsm_copy_nmframe_inplace(llsm_nmframe* dst, llsm_nmframe* src) {
  if(dst -> npsd < src -> npsd)
    dst -> psd = realloc(dst -> psd, sizeof(FP_TYPE) * src -> npsd);
  memcpy(dst -> psd, src -> psd, sizeof(FP_TYPE) * src -> npsd);
  dst -> npsd = src -> npsd;
  
  if(dst -> nchannel < src -> nchannel) {
    dst -> edc = realloc(dst -> edc, sizeof(FP_TYPE) * src -> nchannel);
    dst -> eenv = realloc(dst -> eenv,
      sizeof(llsm_hmframe*) * src -> nchannel);
    for(int i = dst -> nchannel; i < src -> nchannel; i ++)
      dst -> eenv[i] = llsm_create_hmframe(0);
  } else if(dst -> nchannel > src -> nchannel) {
    for(int i = src -> nchannel; i < dst -> nchannel; i ++)
      llsm_delete_hmframe(dst -> eenv[i]);
  }
  for(int i = 0; i < src -> nchannel; i ++) {
    dst -> edc[i] = src -> edc[i];
    llsm_copy_hmframe_inplace(dst -> eenv[i], src -> eenv[i]);
  }
  dst -> nchannel = src -> nchannel;
}

llsm_nmframe* llsm_copy_nmframe(llsm_nmframe* src) {
  llsm_nmframe* ret = llsm_create_nmframe(src -> nchannel, 0, src -> npsd);
  llsm_copy_nmframe_inplace(ret, src);
  return ret;
}

void llsm_delete_nmframe(llsm_nmframe* dst) {
  if(dst == NULL) return;
  for(int i = 0; i < dst -> nchannel; i ++)
    llsm_delete_hmframe(dst -> eenv[i]);
  free(dst -> eenv);
  free(dst -> edc);
  free(dst -> psd);
  free(dst);
}

static FP_TYPE* copy_fp(FP_TYPE* src) {
  FP_TYPE* ret = malloc(sizeof(FP_TYPE));
  ret[0] = src[0];
  return ret;
}

llsm_container* llsm_create_frame(int nhar, int nchannel, int nhar_e,
  int npsd) {
  llsm_container* ret = llsm_create_container(3);
  llsm_hmframe* hm = llsm_create_hmframe(nhar);
  llsm_nmframe* nm = llsm_create_nmframe(nchannel, nhar_e, npsd);
  FP_TYPE* f0 = malloc(sizeof(FP_TYPE));
  f0[0] = 0;
  llsm_container_attach(ret, LLSM_FRAME_F0, f0, free, copy_fp);
  llsm_container_attach(ret, LLSM_FRAME_HM, hm, llsm_delete_hmframe,
    llsm_copy_hmframe);
  llsm_container_attach(ret, LLSM_FRAME_NM, nm, llsm_delete_nmframe,
    llsm_copy_nmframe);
  return ret;
}

void llsm_frame_phaseshift(llsm_container* dst, FP_TYPE theta) {
  llsm_hmframe* hm = llsm_container_get(dst, LLSM_FRAME_HM);
  if(hm != NULL)
    llsm_hmframe_phaseshift(hm, theta);
  llsm_nmframe* nm = llsm_container_get(dst, LLSM_FRAME_NM);
  if(nm != NULL)
    for(int i = 0; i < nm -> nchannel; i ++)
      llsm_hmframe_phaseshift(nm -> eenv[i], theta);
  FP_TYPE* vs_phse = llsm_container_get(dst, LLSM_FRAME_VSPHSE);
  if(vs_phse != NULL) {
    int nhar = llsm_fparray_length(vs_phse);
    for(int i = 0; i < nhar; i ++)
      vs_phse[i] = wrap(vs_phse[i] + theta * (i + 1.0));
  }
}

void llsm_frame_phasesync_rps(llsm_container* dst, int layer1_based) {
  llsm_hmframe* hm = llsm_container_get(dst, LLSM_FRAME_HM);
  FP_TYPE* vs_phse = llsm_container_get(dst, LLSM_FRAME_VSPHSE);
  FP_TYPE phase_ref = 0;
  if(layer1_based && vs_phse != NULL && llsm_fparray_length(vs_phse) > 0) {
    phase_ref = vs_phse[0];
  } else if(hm != NULL && hm -> nhar > 0) {
    phase_ref = hm -> phse[0];
  }
  llsm_frame_phaseshift(dst, -phase_ref);
}

int llsm_frame_checklayer0(llsm_container* src) {
  FP_TYPE* f0 = llsm_container_get(src, LLSM_FRAME_F0);
  llsm_nmframe* hm = llsm_container_get(src, LLSM_FRAME_HM);
  llsm_nmframe* nm = llsm_container_get(src, LLSM_FRAME_NM);
  if(f0 == NULL || nm == NULL) return 0;
  if(*f0 != 0 && hm == NULL) return 0;
  return 1;
}

int llsm_frame_checklayer1(llsm_container* src) {
  FP_TYPE* f0 = llsm_container_get(src, LLSM_FRAME_F0);
  FP_TYPE* rd = llsm_container_get(src, LLSM_FRAME_RD);
  FP_TYPE* spec_env = llsm_container_get(src, LLSM_FRAME_VTMAGN);
  FP_TYPE* vs_phse = llsm_container_get(src, LLSM_FRAME_VSPHSE);
  llsm_nmframe* nm = llsm_container_get(src, LLSM_FRAME_NM);
  if(f0 == NULL || rd == NULL || spec_env == NULL || vs_phse == NULL ||
     nm == NULL) return 0;
  return 1;
}

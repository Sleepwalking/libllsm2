/*
  libllsm2 - Low Level Speech Model (version 2)
  ===
  Copyright (c) 2017 Kanru Hua.

  libllsm2 is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  libllsm2 is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with libllsm. If not, see <http://www.gnu.org/licenses/>.
*/

#include "llsm.h"
#include "dsputils.h"
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

FP_TYPE* llsm_hmframe_harpsd(llsm_hmframe* src, int db_scale) {
  FP_TYPE* psd = calloc(src -> nhar, sizeof(FP_TYPE));
  for(int i = 0; i < src -> nhar; i ++) {
    psd[i] = src -> ampl[i] * src -> ampl[i] * 0.5;
    if(db_scale) psd[i] = 10.0 * log10(psd[i]);
  }
  return psd;
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

FP_TYPE* llsm_frame_compute_snr(llsm_container* src, llsm_container* conf,
  int as_aperiodicity) {
  FP_TYPE* f0 = llsm_container_get(src, LLSM_FRAME_F0);
  llsm_hmframe* hm = llsm_container_get(src, LLSM_FRAME_HM);
  llsm_nmframe* nm = llsm_container_get(src, LLSM_FRAME_NM);
  FP_TYPE* fnyq = llsm_container_get(conf, LLSM_CONF_FNYQ);
  FP_TYPE* noswarp = llsm_container_get(conf, LLSM_CONF_NOSWARP);
  if(f0 == NULL || hm == NULL || nm == NULL) return NULL;
  if(fnyq == NULL || noswarp == NULL) return NULL;
  int nfft = max(64, pow(2, ceil(log2(hm -> nhar) + 2)));
  FP_TYPE* spec_env = llsm_harmonic_envelope(hm -> ampl, hm -> nhar,
    *f0 / *fnyq / 2.0, nfft);
  for(int i = 0; i < nfft / 2 + 1; i ++) {
    spec_env[i] = pow(10.0, spec_env[i] / 20.0); // dB to magnitude
    spec_env[i] *= spec_env[i] * 0.5; // magnitude to variance
  }
  FP_TYPE* warp_axis = llsm_warp_frequency(0, *fnyq, nm -> npsd, *noswarp);
  FP_TYPE* spec_warp = llsm_spectral_mean(spec_env, nfft / 2 + 1, *fnyq,
    warp_axis, nm -> npsd);
  for(int i = 0; i < nm -> npsd; i ++) {
    if(as_aperiodicity) {
      FP_TYPE snr = spec_warp[i] / pow(10.0, nm -> psd[i] / 10.0);
      spec_warp[i] = 1.0 / (1.0 + snr);
    } else
      spec_warp[i] = 10.0 * log10(spec_warp[i]) - nm -> psd[i];
  }
  free(warp_axis);
  free(spec_env);
  return spec_warp;
}

int llsm_frame_checklayer0(llsm_container* src) {
  FP_TYPE* f0 = llsm_container_get(src, LLSM_FRAME_F0);
  llsm_hmframe* hm = llsm_container_get(src, LLSM_FRAME_HM);
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

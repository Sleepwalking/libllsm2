/*
  libllsm2 - Low Level Speech Model (version 2)
  ===
  Copyright (c) 2017-2019 Kanru Hua.

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

// Ooura-fft
void ddct(int, int, FP_TYPE *);

typedef struct {
  int order_spec;
  int order_bap;
  int nfullspec;
  int nchannel;
  int nhar_e;
  int npsd;
  FP_TYPE fnyq;
  FP_TYPE liprad;
  FP_TYPE* psdaxis;
  FP_TYPE* melaxis;
  FP_TYPE* faxis;
  FP_TYPE* apaxis;
} llsm_coder_;

llsm_coder* llsm_create_coder(llsm_container* conf, int order_spec,
  int order_bap) {
  llsm_coder_* ret = malloc(sizeof(llsm_coder_));
  FP_TYPE* fnyq = llsm_container_get(conf, LLSM_CONF_FNYQ);
  int* nchannel = llsm_container_get(conf, LLSM_CONF_NCHANNEL);
  int* nhar_e = llsm_container_get(conf, LLSM_CONF_MAXNHAR_E);
  int* npsd = llsm_container_get(conf, LLSM_CONF_NPSD);
  int* nspec = llsm_container_get(conf, LLSM_CONF_NSPEC);
  FP_TYPE* liprad = llsm_container_get(conf, LLSM_CONF_LIPRADIUS);
  ret -> order_spec = order_spec;
  ret -> order_bap = order_bap;
  ret -> nfullspec = (nspec[0] - 1) * 2;
  ret -> nchannel = nchannel[0];
  ret -> nhar_e = nhar_e[0];
  ret -> npsd = npsd[0];
  ret -> fnyq = fnyq[0];
  ret -> liprad = liprad[0];
  ret -> faxis = calloc(ret -> nfullspec, sizeof(FP_TYPE));
  for(int i = 0; i < ret -> nfullspec; i ++)
    ret -> faxis[i] = fnyq[0] * 2 * i / ret -> nfullspec;
  ret -> psdaxis = linspace(0, ret -> fnyq, ret -> npsd);
  FP_TYPE mel_ceil = freq2mel(ret -> fnyq);
  FP_TYPE mel_floor = freq2mel(50);
  ret -> melaxis = calloc(nspec[0], sizeof(FP_TYPE));
  for(int i = 0; i < nspec[0]; i ++)
    ret -> melaxis[i] = mel2freq(
      mel_floor + (mel_ceil - mel_floor) * i / nspec[0]);
  ret -> apaxis = linspace(0, fnyq[0], order_bap + 1);
  return ret;
}

void llsm_delete_coder(llsm_coder* dst_) {
  llsm_coder_* dst = (llsm_coder_*)dst_;
  if(dst == NULL) return;
  free(dst -> psdaxis);
  free(dst -> melaxis);
  free(dst -> faxis);
  free(dst -> apaxis);
  free(dst);
}

FP_TYPE* llsm_coder_encode(llsm_coder* c_, llsm_container* src) {
  llsm_coder_* c = (llsm_coder_*)c_;
  int ns = c -> nfullspec / 2 + 1;
  FP_TYPE* enc = calloc(c -> order_spec + c -> order_bap + 3, sizeof(FP_TYPE));

  FP_TYPE* f0 = llsm_container_get(src, LLSM_FRAME_F0);
  llsm_nmframe* nm = llsm_container_get(src, LLSM_FRAME_NM);
  enc[0] = f0[0] > 0; // voicing
  enc[1] = f0[0];     // f0

  // from scaled frequency axis to full frequency axis
  FP_TYPE* spec_psd = interp1(c -> psdaxis, nm -> psd, nm -> npsd,
    c -> faxis, ns);
  // from intensity to power
  for(int j = 0; j < ns; j ++)
    spec_psd[j] = exp_2(spec_psd[j] * 2.30258 / 10.0);

  if(f0[0] > 0) {
    FP_TYPE* rd = llsm_container_get(src, LLSM_FRAME_RD);
    FP_TYPE* vtmagn = llsm_container_get(src, LLSM_FRAME_VTMAGN);
    enc[2] = rd[0];
    // spectral synthesis
    lfmodel gfm = lfmodel_from_rd(rd[0], 1.0 / f0[0], 1.0);
    FP_TYPE* lfmagnresp = lfmodel_spectrum(gfm, c -> faxis, ns, NULL);
    FP_TYPE* lfmagnf0 = lfmodel_spectrum(gfm, f0, 1, NULL);
    FP_TYPE* spec_env = calloc(ns, sizeof(FP_TYPE));
    for(int j = 1; j < ns; j ++) {
      spec_env[j] = exp_2(vtmagn[j] * 2.30258 / 20.0)
        * lfmagnresp[j] / lfmagnf0[0] * f0[0] / c -> faxis[j];
    }
    spec_env[0] = spec_env[1];
    llsm_lipfilter(c -> liprad, c -> fnyq / ns, ns, spec_env, NULL, 0);
    // magnitude to PSD (power distributed over the spacing between harmonics)
    for(int j = 1; j < ns; j ++)
      spec_env[j] *= spec_env[j] * 44100 / 4 / f0[0];
    // total power (harmonics + noise)
    for(int j = 0; j < ns; j ++)
      spec_psd[j] += spec_env[j];
    // convert PSD into BAP
    for(int j = 0; j < c -> order_bap; j ++) {
      int n0 = j * (ns - 1) / c -> order_bap;
      int n1 = (j + 1) * (ns - 1) / c -> order_bap;
      FP_TYPE apsum = 0;
      for(int k = n0; k < n1; k ++)
        apsum += 1 - spec_env[k] / spec_psd[k];
      enc[3 + c -> order_spec + j] = apsum / (n1 - n0);
    }
    free(lfmagnf0);
    free(lfmagnresp);
    free(spec_env);
  } else {
    // AP = 1.0
    for(int j = 0; j < c -> order_bap; j ++)
      enc[3 + c -> order_spec + j] = 1.0;
  }

  // from power to log intensity
  for(int j = 0; j < ns; j ++)
    spec_psd[j] = log_2(spec_psd[j]) * 0.5;

  // from full linear frequency axis to full mel-frequency axis
  FP_TYPE* mel_psd = interp1(c -> faxis, spec_psd, ns, c -> melaxis, ns);

  // DCT
  ddct(ns - 1, -1, mel_psd);
  // IDCT to low-order spectrum
  mel_psd[0] *= 0.5;
  ddct(c -> order_spec, 1, mel_psd);
  // normalize
  for(int j = 0; j < c -> order_spec; j ++) {
    mel_psd[j] *= 2.0 / (ns - 1);
    enc[3 + j] = mel_psd[j];
  }

  free(spec_psd);
  free(mel_psd);

  return enc;
}

static llsm_container* llsm_coder_decode(llsm_coder* c_, FP_TYPE* src,
  int use_layer1) {
  llsm_coder_* c = (llsm_coder_*)c_;
  int ns = c -> nfullspec / 2 + 1;
  int voicing = src[0] > 0.5;
  FP_TYPE f0 = max(20.0, src[1]);
  FP_TYPE rd = min(3.0, max(0.02, src[2]));
  int nhar = voicing ? c -> fnyq / f0 : 0;
  llsm_container* ret = llsm_create_frame(
    nhar, c -> nchannel, c -> nhar_e, c -> npsd);
  llsm_nmframe* nm = llsm_container_get(ret, LLSM_FRAME_NM);
  llsm_container_attach(ret, LLSM_FRAME_RD,
    llsm_create_fp(rd), llsm_delete_fp, llsm_copy_fp);
  llsm_container_attach(ret, LLSM_FRAME_F0,
    llsm_create_fp(f0 * voicing), llsm_delete_fp, llsm_copy_fp);

  FP_TYPE* src_spec = src + 3;
  FP_TYPE* src_bap = src + 3 + c -> order_spec;
  FP_TYPE* mel_psd = calloc(ns, sizeof(FP_TYPE));
  FP_TYPE* bap_pad = calloc(c -> order_bap + 1, sizeof(FP_TYPE));
  // undo the IDCT
  for(int j = 0; j < c -> order_spec; j ++)
    mel_psd[j] = src_spec[j] * 0.5 * (ns - 1) * 2.0 / c -> order_spec;
  for(int j = c -> order_spec; j < ns; j ++) mel_psd[j] = 0;
  ddct(c -> order_spec, -1, mel_psd);
  // IDCT to full-order spectrum
  mel_psd[0] *= 0.5;
  ddct(ns - 1, 1, mel_psd);
  for(int j = 0; j < ns - 1; j ++)
    mel_psd[j] *= 2.0 / (ns - 1);
  mel_psd[ns - 1] = mel_psd[ns - 2];
  // band aperiodicity to full aperiodicity
  for(int j = 0; j < c -> order_bap; j ++)
    bap_pad[j + 1] = src_bap[j];
  bap_pad[0] = voicing ? 0 : 1;
  FP_TYPE* full_psd = interp1(c -> melaxis, mel_psd, ns, c -> faxis, ns);
  FP_TYPE* full_ap = interp1(c -> apaxis, bap_pad, c -> order_bap + 1,
    c -> faxis, ns);
  for(int j = 0; j < ns; j ++) {
    if(voicing) { // post-processing trick to reduce low-frequency noise
      FP_TYPE fj = j * c -> fnyq / ns;
      if(fj < 500)
        full_ap[j] = 1e-3;
      else if(fj < 2000)
        full_ap[j] = 1e-3 + (full_ap[j] - 1e-3) * (fj - 500) / 1500;
    }
    // decompose PSD into harmonic magnitude and noise power
    full_psd[j] = exp_2(2.0 * full_psd[j]);
    FP_TYPE sum_psd = full_psd[j];
    FP_TYPE per_psd = sum_psd * (1.0 - full_ap[j]);
    full_psd[j] = sqrt(per_psd * f0 * 4 / 44100);
    full_ap[j]  = sum_psd * full_ap[j];
  }
  FP_TYPE* full_spec  = full_psd; full_psd = NULL;
  FP_TYPE* full_noise = full_ap;  full_ap = NULL;

  // power to log intensity
  free(nm -> psd);
  nm -> psd = interp1(c -> faxis, full_noise, ns, c -> psdaxis, nm -> npsd);
  for(int j = 0; j < nm -> npsd; j ++)
    nm -> psd[j] = 10.0 / 2.30258 * log_2(nm -> psd[j]);

  if(nhar > 0 && use_layer1) {
    llsm_container_remove(ret, LLSM_FRAME_HM);
    lfmodel gfm = lfmodel_from_rd(rd, 1.0 / f0, 1.0);
    FP_TYPE* lfmagnresp = lfmodel_spectrum(gfm, c -> faxis, ns, NULL);
    FP_TYPE* lfmagnf0 = lfmodel_spectrum(gfm, & f0, 1, NULL);
    llsm_lipfilter(c -> liprad, c -> fnyq / ns, ns, full_spec, NULL, 1);
    // magnitude to log
    for(int j = 1; j < ns; j ++)
      full_spec[j] = 20.0 / 2.30258 * log_2(full_spec[j]
        * c -> faxis[j] / f0 * lfmagnf0[0] / lfmagnresp[j]);
    full_spec[0] = full_spec[1];
    FP_TYPE* vtmagn = llsm_create_fparray(ns);
    FP_TYPE* vsphse = llsm_create_fparray(nhar);
    llsm_container_attach(ret, LLSM_FRAME_VTMAGN, vtmagn,
      llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(ret, LLSM_FRAME_VSPHSE, vsphse,
      llsm_delete_fparray, llsm_copy_fparray);
    for(int j = 0; j < ns; j ++) vtmagn[j] = full_spec[j];
    FP_TYPE* harfreq = linspace(0, nhar * f0, nhar + 1);
    free(lfmodel_spectrum(gfm, harfreq + 1, nhar, vsphse));
    free(harfreq);
    free(lfmagnresp);
    free(lfmagnf0);
  }
  if(nhar > 0 && ! use_layer1) {
    llsm_hmframe* hm = llsm_create_hmframe(nhar);
    llsm_container_attach(ret, LLSM_FRAME_HM, hm,
      llsm_delete_hmframe, llsm_copy_hmframe);
    FP_TYPE* harfreq = linspace(0, nhar * f0, nhar + 1);
    FP_TYPE* ampl = interp1(c -> faxis, full_spec, ns, harfreq + 1, nhar);
    for(int i = 0; i < nhar; i ++)
      hm -> ampl[i] = ampl[i];
    llsm_lipfilter(c -> liprad, f0, nhar, ampl, NULL, 1);
    lfmodel gfm = lfmodel_from_rd(rd, 1.0 / f0, 1.0);
    FP_TYPE* vsphse = calloc(nhar, sizeof(FP_TYPE));
    // recover vocal tract magnitude response
    FP_TYPE* lfmagnresp = lfmodel_spectrum(gfm, harfreq + 1, nhar, vsphse);
    for(int i = 0; i < nhar; i ++) {
      FP_TYPE vs_ampl = lfmagnresp[i] / (i + 1.0) / lfmagnresp[0];
      ampl[i] /= vs_ampl;
    }
    // compute radiated phase
    FP_TYPE* vtphse = llsm_harmonic_minphase(ampl, nhar);
    llsm_lipfilter(c -> liprad, f0, nhar, NULL, vtphse, 0);
    for(int i = 0; i < nhar; i ++)
      hm -> phse[i] = vtphse[i] + vsphse[i];
    free(harfreq);
    free(ampl);
    free(lfmagnresp);
    free(vtphse);
  }

  free(mel_psd);
  free(bap_pad);
  free(full_spec);
  free(full_noise);
  return ret;
}

llsm_container* llsm_coder_decode_layer1(llsm_coder* c, FP_TYPE* src) {
  return llsm_coder_decode(c, src, 1);
}

llsm_container* llsm_coder_decode_layer0(llsm_coder* c, FP_TYPE* src) {
  return llsm_coder_decode(c, src, 0);
}

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
  int npsd;
  FP_TYPE fnyq;
  FP_TYPE liprad;
  FP_TYPE* freqwarp;
  FP_TYPE* freq2mel;
  FP_TYPE* faxis;
} llsm_coder_;

llsm_coder* llsm_create_coder(llsm_container* conf, int order_spec,
  int order_bap) {
  llsm_coder_* ret = malloc(sizeof(llsm_coder_));
  FP_TYPE* fnyq = llsm_container_get(conf, LLSM_CONF_FNYQ);
  int* npsd = llsm_container_get(conf, LLSM_CONF_NPSD);
  int* nspec = llsm_container_get(conf, LLSM_CONF_NSPEC);
  FP_TYPE* noswarp = llsm_container_get(conf, LLSM_CONF_NOSWARP);
  FP_TYPE* liprad = llsm_container_get(conf, LLSM_CONF_LIPRADIUS);
  ret -> order_spec = order_spec;
  ret -> order_bap = order_bap;
  ret -> nfullspec = (nspec[0] - 1) * 2;
  ret -> npsd = npsd[0];
  ret -> fnyq = fnyq[0];
  ret -> liprad = liprad[0];
  ret -> faxis = calloc(ret -> nfullspec, sizeof(FP_TYPE));
  for(int i = 0; i < ret -> nfullspec; i ++)
    ret -> faxis[i] = fnyq[0] * 2 * i / ret -> nfullspec;
  ret -> freqwarp = llsm_warp_frequency(0, ret -> fnyq, ret -> npsd, *noswarp);
  FP_TYPE mel_ceil = freq2mel(ret -> fnyq);
  FP_TYPE mel_floor = freq2mel(50);
  ret -> freq2mel = calloc(nspec[0], sizeof(FP_TYPE));
  for(int i = 0; i < nspec[0]; i ++)
    ret -> freq2mel[i] = mel2freq(
      mel_floor + (mel_ceil - mel_floor) * i / nspec[0]);
  return ret;
}

void llsm_delete_coder(llsm_coder* dst_) {
  llsm_coder_* dst = (llsm_coder_*)dst_;
  if(dst == NULL) return;
  free(dst -> freqwarp);
  free(dst -> freq2mel);
  free(dst -> faxis);
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

  // from warpped frequency axis to full linear frequency axis
  FP_TYPE* spec_psd = interp1(c -> freqwarp, nm -> psd, nm -> npsd,
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
  FP_TYPE* mel_psd = interp1(c -> faxis, spec_psd, ns, c -> freq2mel, ns);

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

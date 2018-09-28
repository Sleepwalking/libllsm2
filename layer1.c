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

static int llsm_layer0to1_check_integrity(llsm_chunk* src) {
  int* nfrm = llsm_container_get(src -> conf, LLSM_CONF_NFRM);
  FP_TYPE* thop = llsm_container_get(src -> conf, LLSM_CONF_THOP);
  FP_TYPE* fnyq = llsm_container_get(src -> conf, LLSM_CONF_FNYQ);
  FP_TYPE* liprad = llsm_container_get(src -> conf, LLSM_CONF_LIPRADIUS);
  if(nfrm == NULL || thop == NULL || fnyq == NULL || liprad == NULL) 
    return 0;
  for(int i = 0; i < *nfrm; i ++)
    if(! llsm_frame_checklayer0(src -> frames[i]))
      return 0;
  return 1;
}

static int llsm_layer1to0_check_integrity(llsm_container* conf) {
  int* nfrm = llsm_container_get(conf, LLSM_CONF_NFRM);
  FP_TYPE* fnyq = llsm_container_get(conf, LLSM_CONF_FNYQ);
  FP_TYPE* liprad = llsm_container_get(conf, LLSM_CONF_LIPRADIUS);
  int* nspec = llsm_container_get(conf, LLSM_CONF_NSPEC);
  if(nfrm == NULL || fnyq == NULL || liprad == NULL || nspec == NULL)
    return 0;
  return 1;
}

static FP_TYPE* llsm_analyze_rd(llsm_chunk* src) {
  int nfrm = *((int*)llsm_container_get(src -> conf, LLSM_CONF_NFRM));
  FP_TYPE thop = *((FP_TYPE*)llsm_container_get(src -> conf, LLSM_CONF_THOP));
  FP_TYPE lip_radius = *((FP_TYPE*)llsm_container_get(src -> conf,
    LLSM_CONF_LIPRADIUS));
  
  int ncandidate = 64;
  FP_TYPE* rd_list = linspace(0.02, 3.0, ncandidate);
  llsm_cached_glottal_model* cgm = llsm_create_cached_glottal_model(
    rd_list, ncandidate, 80);
  free(rd_list);
  FP_TYPE* rd = calloc(nfrm, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE f0 = *((FP_TYPE*)llsm_container_get(src -> frames[i],
      LLSM_FRAME_F0));
    if(f0 == 0) continue;
    llsm_hmframe* hm = llsm_container_get(src -> frames[i], LLSM_FRAME_HM);
    int nhar = min(hm -> nhar, round(8000.0 / f0));
    FP_TYPE* ampl = calloc(nhar, sizeof(FP_TYPE));
    memcpy(ampl, hm -> ampl, nhar * sizeof(FP_TYPE));
    llsm_lipfilter(lip_radius, f0, nhar, ampl, NULL, 1);
    rd[i] = llsm_spectral_glottal_fitting(ampl, nhar, cgm);
    free(ampl);
  }
  llsm_delete_cached_glottal_model(cgm);

  FP_TYPE* rd_cont = interp_in_blank(rd, nfrm, 0);
  FP_TYPE* rd_smooth = llsm_smoothing_filter(rd_cont, nfrm,
    round(0.02 / thop));

  free(rd); free(rd_cont);
  return rd_smooth;
}

// Note: the reason we don't expose this function in llsm.h just like its
//   counterpart llsm_frame_tolayer0, is that layer0-to-layer1 conversion
//   cannot be done in a frame-by-frame manner, unless the Rd parameter is
//   known in advance.
static void llsm_frame_tolayer1(llsm_container* dst, FP_TYPE lip_radius,
  FP_TYPE fnyq, int nfft) {
  int nspec = nfft / 2 + 1;
  FP_TYPE rd = *((FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_RD));
  FP_TYPE f0 = *((FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_F0));
  llsm_hmframe* hm = llsm_container_get(dst, LLSM_FRAME_HM);
  FP_TYPE* ampl = calloc(hm -> nhar, sizeof(FP_TYPE));
  FP_TYPE* phse = calloc(hm -> nhar, sizeof(FP_TYPE));
  FP_TYPE* freq = calloc(hm -> nhar, sizeof(FP_TYPE));
  memcpy(ampl, hm -> ampl, hm -> nhar * sizeof(FP_TYPE));
  memcpy(phse, hm -> phse, hm -> nhar * sizeof(FP_TYPE));
  for(int i = 0; i < hm -> nhar; i ++) freq[i] = f0 * (i + 1.0);

  // Generate a LF pulse, and normalize.
  lfmodel glott = lfmodel_from_rd(rd, 1.0 / f0, 1.0);
  FP_TYPE* vs_ampl = lfmodel_spectrum(glott, freq, hm -> nhar, NULL);
  for(int i = 1; i < hm -> nhar; i ++) vs_ampl[i] /= (1.0 + i) * vs_ampl[0];
  vs_ampl[0] = 1.0;

  llsm_lipfilter(lip_radius, f0, hm -> nhar, ampl, phse, 1);
  for(int i = 0; i < hm -> nhar; i ++) ampl[i] /= vs_ampl[i];

  FP_TYPE* vt_phse = llsm_harmonic_minphase(ampl, hm -> nhar);
  FP_TYPE* vs_phse = calloc(hm -> nhar, sizeof(FP_TYPE));
  for(int i = 0; i < hm -> nhar; i ++) vs_phse[i] = phse[i] - vt_phse[i];

  // The spectral envelope after removing lip and glottal responses.
  FP_TYPE* spec_env = llsm_harmonic_envelope(ampl, hm -> nhar,
    f0 / fnyq / 2.0, nfft);
  
  FP_TYPE* arr_vs_phse = llsm_create_fparray(hm -> nhar);
  FP_TYPE* arr_spec_env = llsm_create_fparray(nspec);
  memcpy(arr_spec_env, spec_env, nspec * sizeof(FP_TYPE));
  memcpy(arr_vs_phse, vs_phse, hm -> nhar * sizeof(FP_TYPE));
  llsm_container_attach(dst, LLSM_FRAME_VTMAGN, arr_spec_env,
    llsm_delete_fparray, llsm_copy_fparray);
  llsm_container_attach(dst, LLSM_FRAME_VSPHSE, arr_vs_phse,
    llsm_delete_fparray, llsm_copy_fparray);
  
  free(ampl); free(phse); free(freq);
  free(vt_phse); free(vs_ampl); free(vs_phse); free(spec_env);
}

void llsm_chunk_tolayer1(llsm_chunk* dst, int nfft) {
  if(! llsm_layer0to1_check_integrity(dst)) return;
  int nfrm = *((int*)llsm_container_get(dst -> conf, LLSM_CONF_NFRM));
  FP_TYPE fnyq = *((FP_TYPE*)llsm_container_get(dst -> conf, LLSM_CONF_FNYQ));
  FP_TYPE lip_radius = *((FP_TYPE*)llsm_container_get(dst -> conf,
    LLSM_CONF_LIPRADIUS));

  llsm_container_attach(dst -> conf, LLSM_CONF_NSPEC,
    llsm_create_int(nfft / 2 + 1), llsm_delete_int, llsm_copy_int);

  FP_TYPE* rd = llsm_analyze_rd(dst);
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE f0 = *((FP_TYPE*)llsm_container_get(dst -> frames[i],
      LLSM_FRAME_F0));
    llsm_container_attach(dst -> frames[i], LLSM_FRAME_RD,
      llsm_create_fp(rd[i]), llsm_delete_fp, llsm_copy_fp);
    if(f0 == 0) continue;
    llsm_frame_tolayer1(dst -> frames[i], lip_radius, fnyq, nfft);
  }
  free(rd);
}

void llsm_frame_tolayer0(llsm_container* dst, llsm_container* conf) {
  if(! llsm_layer1to0_check_integrity(conf)) return;
  if(! llsm_frame_checklayer1(dst)) return;
  FP_TYPE fnyq = *((FP_TYPE*)llsm_container_get(conf, LLSM_CONF_FNYQ));
  FP_TYPE lip_radius = *((FP_TYPE*)llsm_container_get(conf,
    LLSM_CONF_LIPRADIUS));
  int nspec = *((int*)llsm_container_get(conf, LLSM_CONF_NSPEC));
  FP_TYPE* f0 = llsm_container_get(dst, LLSM_FRAME_F0);
  FP_TYPE* rd = llsm_container_get(dst, LLSM_FRAME_RD);
  FP_TYPE* spec_env = llsm_container_get(dst, LLSM_FRAME_VTMAGN);
  FP_TYPE* vs_phse = llsm_container_get(dst, LLSM_FRAME_VSPHSE);
  if(*f0 == 0) return;

  int nhar = llsm_fparray_length(vs_phse);
  int* maxnhar = llsm_container_get(conf, LLSM_CONF_MAXNHAR);
  if(maxnhar != NULL) nhar = min(nhar, *maxnhar);

  FP_TYPE* freq = calloc(nhar, sizeof(FP_TYPE));
  for(int i = 0; i < nhar; i ++) freq[i] = *f0 * (i + 1.0);

  lfmodel glott = lfmodel_from_rd(*rd, 1.0 / *f0, 1.0);
  FP_TYPE* vs_ampl = lfmodel_spectrum(glott, freq, nhar, NULL);
  for(int i = 1; i < nhar; i ++) vs_ampl[i] /= (1.0 + i) * vs_ampl[0];
  vs_ampl[0] = 1.0;

  FP_TYPE* faxis = linspace(0, fnyq, nspec);
  FP_TYPE* vt_ampl = interp1(faxis, spec_env, nspec, freq, nhar);
  for(int i = 0; i < nhar; i ++)
    vt_ampl[i] = exp(vt_ampl[i] / 20.0 * 2.3025851); // db2mag
  FP_TYPE* vt_phse = llsm_harmonic_minphase(vt_ampl, nhar);
  
  llsm_hmframe* hm = llsm_create_hmframe(nhar);
  for(int i = 0; i < nhar; i ++) {
    hm -> ampl[i] = vt_ampl[i] * vs_ampl[i];
    hm -> phse[i] = vt_phse[i] + vs_phse[i];
  }
  llsm_lipfilter(lip_radius, *f0, nhar, hm -> ampl, hm -> phse, 0);
  llsm_container_attach(dst, LLSM_FRAME_HM, hm, llsm_delete_hmframe,
    llsm_copy_hmframe);

  free(freq);
  free(faxis);
  free(vt_ampl); free(vt_phse); free(vs_ampl);
}

void llsm_chunk_tolayer0(llsm_chunk* dst) {
  int nfrm = *((int*)llsm_container_get(dst -> conf, LLSM_CONF_NFRM));
  for(int i = 0; i < nfrm; i ++)
    llsm_frame_tolayer0(dst -> frames[i], dst -> conf);
}

FP_TYPE* llsm_make_filtered_pulse(llsm_container* src, lfmodel source,
  FP_TYPE phase_correction, int pre_rotate, int size, FP_TYPE fnyq,
  FP_TYPE lip_radius, FP_TYPE fs) {
  FP_TYPE* buffer = calloc(size * 6, sizeof(FP_TYPE));
  FP_TYPE* lfphseresp = buffer + size * 2;
  FP_TYPE* real_resp  = buffer + size * 3;
  FP_TYPE* imag_resp  = buffer + size * 4;
  FP_TYPE* freq_axis  = buffer + size * 5;
  int halfsize = size / 2 + 1;
  for(int i = 0; i < halfsize; i ++)
    freq_axis[i] = i * fs / size;
  
  FP_TYPE* f0 = llsm_container_get(src, LLSM_FRAME_F0);
  FP_TYPE* vsphse = llsm_container_get(src, LLSM_FRAME_VSPHSE);
  int nhar = llsm_fparray_length(vsphse);
  FP_TYPE* freq_har = calloc(nhar + 1, sizeof(FP_TYPE));
  FP_TYPE* phse_har = calloc(nhar + 1, sizeof(FP_TYPE));
  for(int i = 0; i <= nhar; i ++) freq_har[i] = i * f0[0];
  // First, we compute the difference between LF-model phase and the actual
  //   phase for each harmonic. This difference will be freq-interpolated and
  //   added back to the phase spectrum so that the pulse-by-pulse synthesized
  //   speech matches the result from harmonic models.
  free(lfmodel_spectrum(source, freq_har + 1, nhar, phse_har + 1));
  FP_TYPE vsshift = vsphse[0] - (phse_har[1] - 0.5 * M_PI);
  for(int i = 1; i <= nhar; i ++) {
    phse_har[i] -= 0.5 * M_PI;
    phse_har[i] = wrap(vsphse[i - 1] - phse_har[i] - vsshift * i);
  }
  
  // Again, to keep PbP synthesis consistent with HM, we need to compute
  //   vocal tract phase directly from its harmonic representation, albeit
  //   at a cost of slightly breaking minimum phase property (w.r.t the
  //   full-sized spectrum).
  FP_TYPE* vtmagn = llsm_container_get(src, LLSM_FRAME_VTMAGN);
  int nspec = llsm_fparray_length(vtmagn);
  FP_TYPE* vtaxis = linspace(0, fnyq, nspec);
  FP_TYPE* vtamplhar = interp1(vtaxis, vtmagn, nspec, freq_har + 1, nhar);
  for(int i = 0; i < nhar; i ++)
    vtamplhar[i] = exp(vtamplhar[i] / 20.0 * 2.3025851); // db2mag
  FP_TYPE* vt_phse = llsm_harmonic_minphase(vtamplhar, nhar);
  // accumulate VT phase to the phase delta vector
  for(int i = 0; i < nhar; i ++)
    phse_har[i + 1] += vt_phse[i];
  free(vtamplhar); free(vt_phse);
  
  // Now the harmonic phase delta will be expanded into a full-sized phase
  //   envelope. Phase interpolation is error-prone but in this context the
  //   errors won't matter (too much) thanks to the error-cancelling effect
  //   of PSOLA.
  FP_TYPE* phse_har_unwrap = unwrap(phse_har, nhar + 1);
  FP_TYPE* phse_delta = interp1(freq_har, phse_har_unwrap, nhar + 1,
    freq_axis, halfsize);
  free(phse_har);
  free(phse_har_unwrap);
  
  // From this point we will move from harmonic reprensentations to full-sized
  //   spectra. A spectrum generated from LF model is first integrated (to
  //   become glottal flow velocity); then various phase corrections are
  //   applied. The amplitude is normalized by the first glottal harmonic and
  //   appropriately scaled for IFFT.
  FP_TYPE* lfmagnf0 = lfmodel_spectrum(source, f0, 1, NULL);
  FP_TYPE* lfmagnresp = lfmodel_spectrum(
    source, freq_axis, halfsize, lfphseresp);
  lfmagnresp[0] = 0;
  lfphseresp[0] = 0;
  for(int i = 1; i < halfsize; i ++) {
    lfmagnresp[i] *= (fnyq / freq_axis[i]) / lfmagnf0[0];
    lfphseresp[i] += (phase_correction - pre_rotate) * i * 2 * M_PI / size;
    lfphseresp[i] += phse_delta[i] - 0.5 * M_PI;
  }
  free(lfmagnf0);
  free(phse_delta);
  
  // Apply the lip radiation filter.
  llsm_lipfilter(lip_radius, fs / size, halfsize, lfmagnresp, lfphseresp, 0);
  
  // Apply the vocal tract magnitude filter (whose phase part has already been
  //   addressed). Convert the result into complex form.
  FP_TYPE* vtmagn_scaled = interp1(vtaxis, vtmagn, nspec, freq_axis, halfsize);
  for(int i = 0; i < halfsize; i ++)
    vtmagn_scaled[i] *= 2.3025851 / 20.0; // db2log
  for(int i = 0; i < halfsize; i ++) {
    lfmagnresp[i] *= exp(vtmagn_scaled[i]);
    real_resp[i] = lfmagnresp[i] * cos_2(lfphseresp[i]);
    imag_resp[i] = lfmagnresp[i] * sin_2(lfphseresp[i]);
  }
  free(vtaxis); free(vtmagn_scaled);
  free(freq_har);
  
  // Recover the negative part using real-dft symmetry and inverse transform.
  complete_symm (real_resp, size);
  complete_asymm(imag_resp, size);
  ifft(real_resp, imag_resp, real_resp, NULL, size, buffer);
  
  // Some tricks to reduce glitches at boundaries.
  int fadein = min(256, pre_rotate);
  int fadeout = min(256, size);
  FP_TYPE* y = calloc(size, sizeof(FP_TYPE));
  for(int i = 0; i < size; i ++)
    y[i] = real_resp[i];
  for(int i = 0; i < fadein; i ++)
    y[i] *= (FP_TYPE)i / fadein;
  for(int i = size - fadeout; i < size; i ++)
    y[i] *= (FP_TYPE)(size - i) / fadeout;
  free(lfmagnresp);
  free(buffer);

  return y;
}

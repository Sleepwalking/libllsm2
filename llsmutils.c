/*
  libllsm2 - Low Level Speech Model (version 2)
  ===
  Copyright (c) 2017-2018 Kanru Hua.

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
#include "llsmutils.h"

FP_TYPE* llsm_synthesize_harmonic_frame_auto(llsm_soptions* options,
  FP_TYPE* ampl, FP_TYPE* phse, int nhar, FP_TYPE f0, int nx) {
  FP_TYPE* ret = NULL;
  if(options == NULL || (! options -> use_iczt))
    ret = llsm_synthesize_harmonic_frame(ampl, phse, nhar, f0, nx);
  else {
    if(log_1(nx) * options -> iczt_param_a <
       log_1(nhar) - options -> iczt_param_b)
      ret = llsm_synthesize_harmonic_frame_iczt(ampl, phse, nhar, f0, nx);
    else
      ret = llsm_synthesize_harmonic_frame(ampl, phse, nhar, f0, nx);
  }
  return ret;
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

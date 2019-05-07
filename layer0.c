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
#include "llsmutils.h"
#include "constants.h"
#include "external/ciglet/ciglet.h"

llsm_aoptions* llsm_create_aoptions() {
  llsm_aoptions* ret = malloc(sizeof(llsm_aoptions));
  ret -> thop = 0.005;
  ret -> maxnhar = 100;
  ret -> maxnhar_e = 4;
  ret -> npsd = 256;
  ret -> nchannel = 4;
  ret -> chanfreq = calloc(3, sizeof(FP_TYPE));
  ret -> chanfreq[0] = 2000.0;
  ret -> chanfreq[1] = 4000.0;
  ret -> chanfreq[2] = 8000.0;
  ret -> lip_radius = 1.5;
  ret -> f0_refine = 1;
  ret -> hm_method = LLSM_AOPTION_HMCZT;
  ret -> rel_winsize = 4.0;
  return ret;
}

void llsm_delete_aoptions(llsm_aoptions* dst) {
  if(dst == NULL) return;
  free(dst -> chanfreq);
  free(dst);
}

llsm_container* llsm_aoptions_toconf(llsm_aoptions* src, FP_TYPE fnyq) {
  llsm_container* ret = llsm_create_container(10);
  llsm_container_attach(ret, LLSM_CONF_NFRM,
    llsm_create_int(0), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_THOP,
    llsm_create_fp(src -> thop), llsm_delete_fp, llsm_copy_fp);
  llsm_container_attach(ret, LLSM_CONF_MAXNHAR,
    llsm_create_int(src -> maxnhar), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_MAXNHAR_E,
    llsm_create_int(src -> maxnhar_e), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_NPSD,
    llsm_create_int(src -> npsd), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_FNYQ,
    llsm_create_fp(fnyq), llsm_delete_fp, llsm_copy_fp);
  llsm_container_attach(ret, LLSM_CONF_NCHANNEL,
    llsm_create_int(src -> nchannel), llsm_delete_int, llsm_copy_int);
  // llsm_container_attach(ret, LLSM_CONF_NSPEC,
  //   llsm_create_int(1024), llsm_delete_int, llsm_copy_int);
  llsm_container_attach(ret, LLSM_CONF_LIPRADIUS,
    llsm_create_fp(src -> lip_radius), llsm_delete_fp, llsm_copy_fp);
  FP_TYPE* chanfreq = llsm_create_fparray(src -> nchannel - 1);
  memcpy(chanfreq, src -> chanfreq, sizeof(FP_TYPE) * (src -> nchannel - 1));
  llsm_container_attach(ret, LLSM_CONF_CHANFREQ,
    chanfreq, llsm_delete_fparray, llsm_copy_fparray);
  return ret;
}

llsm_soptions* llsm_create_soptions(FP_TYPE fs) {
  llsm_soptions* ret = malloc(sizeof(llsm_soptions));
  ret -> fs = fs;
  // See test/test-harmonic.c for more information.
  ret -> use_iczt = 1;
  ret -> use_l1 = 0;
  ret -> iczt_param_a = 0.275;
  ret -> iczt_param_b = 2.26;
  return ret;
}

void llsm_delete_soptions(llsm_soptions* dst) {
  if(dst == NULL) return;
  free(dst);
}

static void llsm_analyze_harmonics(llsm_aoptions* options, FP_TYPE* x, int nx,
  FP_TYPE fs, FP_TYPE* f0, int nfrm, llsm_chunk* dst_chunk) {

  int*      tmp_nhar = calloc(nfrm, sizeof(int));
  FP_TYPE** tmp_ampl = calloc(nfrm, sizeof(FP_TYPE*));
  FP_TYPE** tmp_phse = calloc(nfrm, sizeof(FP_TYPE*));

  llsm_harmonic_analysis(x, nx, fs, f0, nfrm, options -> thop,
    options -> rel_winsize, options -> maxnhar, options -> hm_method,
    tmp_nhar, tmp_ampl, tmp_phse);

  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] == 0) continue;
    llsm_hmframe* hm = llsm_create_hmframe(tmp_nhar[i]);
    memcpy(hm -> ampl, tmp_ampl[i], tmp_nhar[i] * sizeof(FP_TYPE));
    memcpy(hm -> phse, tmp_phse[i], tmp_nhar[i] * sizeof(FP_TYPE));
    llsm_container_attach(dst_chunk -> frames[i], LLSM_FRAME_HM,
      hm, llsm_delete_hmframe, llsm_copy_hmframe);
  }

  free2d(tmp_ampl, nfrm); free2d(tmp_phse, nfrm); free(tmp_nhar);
}

static FP_TYPE* llsm_synthesize_harmonics_l0(llsm_soptions* options,
  llsm_chunk* chunk, FP_TYPE* f0, int nfrm, FP_TYPE thop, FP_TYPE fs, int ny) {
  const int maxnhar = 2048;
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  int nwin = round(thop * fs) * 2;
  FP_TYPE* w = hanning(nwin);
  FP_TYPE* phase = calloc(maxnhar, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] == 0) continue; // skip unvoiced frames
    llsm_hmframe* hm = llsm_container_get(chunk -> frames[i], LLSM_FRAME_HM);
    FP_TYPE rawidx = i * thop * fs;
    int baseidx = round(rawidx);
    FP_TYPE phase_correction = (rawidx - baseidx) * 2 * M_PI / fs * f0[i];
    int nhar = min(maxnhar, hm -> nhar);
    for(int k = 0; k < nhar; k ++)
      phase[k] = hm -> phse[k] - phase_correction * (k + 1.0);
    FP_TYPE* yi = llsm_synthesize_harmonic_frame_auto(options,
      hm -> ampl, phase, nhar, f0[i] / fs, nwin);
    for(int j = 0; j < nwin; j ++) {
      yi[j] *= w[j];
      int idx = baseidx + j - nwin / 2;
      if(idx >= 0 && idx < ny)
        y[idx] += yi[j];
    }
    free(yi);
  }
  free(phase);
  free(w);
  return y;
}

static FP_TYPE* llsm_synthesize_harmonics(llsm_soptions* options,
  llsm_chunk* chunk, FP_TYPE* f0, int nfrm, FP_TYPE thop, FP_TYPE fs, int ny) {
  if(! options -> use_l1) {
    return llsm_synthesize_harmonics_l0(
      options, chunk, f0, nfrm, thop, fs, ny);
  }

  const int maxnhar = 2048;
  FP_TYPE* y_hm  = calloc(ny, sizeof(FP_TYPE)); // harmonic model
  FP_TYPE* y_pbp = calloc(ny, sizeof(FP_TYPE)); // pulse-by-pulse synthesis
  FP_TYPE* y_mix = calloc(ny, sizeof(FP_TYPE));
  int nwin = round(thop * fs) * 2;
  FP_TYPE* w = hanning(nwin);
  FP_TYPE* fnyq = llsm_container_get(chunk -> conf, LLSM_CONF_FNYQ);
  FP_TYPE* liprad = llsm_container_get(chunk -> conf, LLSM_CONF_LIPRADIUS);

  FP_TYPE pulse_previous = 0; // aligned to zero-time phase of glottal flow
  int pbp_periods = 0;
  int pbp_periods_thrd = 3;
  FP_TYPE pbp_switch_rate = 0;
  FP_TYPE pbp_switch_state = 0;
  int baseidx_prev = 0;

  for(int i = 0; i < nfrm; i ++) {
    if(f0[i] == 0) continue; // skip unvoiced frames
    int baseidx = i * thop * fs;
    llsm_container* src_frame = chunk -> frames[i];
    FP_TYPE* vsphse = llsm_container_get(src_frame, LLSM_FRAME_VSPHSE);
    FP_TYPE* vtmagn = llsm_container_get(src_frame, LLSM_FRAME_VTMAGN);
    FP_TYPE* rd = llsm_container_get(src_frame, LLSM_FRAME_RD);
    int* pbpsyn = llsm_container_get(src_frame, LLSM_FRAME_PBPSYN);
    llsm_pbpeffect* pbpeff = llsm_container_get(src_frame, LLSM_FRAME_PBPEFF);
    if(vsphse == NULL || vtmagn == NULL || rd == NULL) continue;
    int pbp_on = pbpsyn != NULL && pbpsyn[0] == 1;
    int nspec = llsm_fparray_length(vtmagn);
    // update locations of pulses locked onto the first source harmonic
    FP_TYPE len_period = fs / f0[i];
    FP_TYPE t_period = 1.0 / f0[i];
    lfmodel source_model = lfmodel_from_rd(*rd, t_period, 1.0);
    FP_TYPE source_p0 = 0;
    free(lfmodel_spectrum(source_model, & f0[i], 1, & source_p0));
    source_p0 -= 0.5 * M_PI; // integrate (flow derivative to flow velocity)
    FP_TYPE p0 = wrap(vsphse[0]);
    FP_TYPE p0_dist = phase_diff(source_p0, p0);
    if(p0_dist < 0) p0_dist += 2.0 * M_PI;
    // the next position where a glottal flow cycle begins
    FP_TYPE pulse_projected = baseidx + p0_dist / 2 / M_PI * len_period;
    // reset the pulse tracker after unvoiced part
    int len_reset = max(len_period, thop * fs) * 2;
    if(pulse_projected - pulse_previous > len_reset)
      pulse_previous = pulse_projected - len_reset;
    int num_periods = round((pulse_projected - pulse_previous) / len_period);
    len_period = (pulse_projected - pulse_previous) / num_periods;
    int pulse_size = pow(2, ceil(log2(max(len_period * 2, nspec))));

    // pulse-by-pulse synthesis
    if(pbp_on || pbp_periods > 0) {
      if(num_periods > 0) {
        FP_TYPE* offsets = calloc(num_periods, sizeof(FP_TYPE));
        lfmodel* sources = calloc(num_periods, sizeof(lfmodel));
        for(int j = 0; j < num_periods; j ++) {
          FP_TYPE delta_t = 0;
          if(pbpeff != NULL) {
            llsm_gfm g = llsm_lfmodel_to_gfm(source_model);
            pbpeff -> modifier(& g, & delta_t, pbpeff -> info, src_frame);
            sources[j] = llsm_gfm_to_lfmodel(g);
          } else
            sources[j] = source_model;
          offsets[j] = pulse_previous + j * len_period + delta_t * fs;
        }
        int pulse_base = offsets[0];
        for(int j = 0; j < num_periods; j ++) offsets[j] -= pulse_base;
        FP_TYPE* y = llsm_make_filtered_pulse(src_frame, sources, offsets,
          num_periods, len_period, pulse_size, *fnyq, *liprad, fs);
        for(int k = 0; k < pulse_size; k ++) {
          int idx = pulse_base + k - len_period;
          if(idx >= 0 && idx < ny) y_pbp[idx] += y[k];
        }
        free(y);
        free(offsets);
        free(sources);
        pbp_periods += pbp_on ? num_periods : -num_periods;
        pbp_periods = min(pbp_periods, pbp_periods_thrd);
        pbp_periods = max(pbp_periods, 0);
      }
    }
    pulse_previous = pulse_projected;

    pbp_switch_rate = 1.0 / min(len_period, thop * fs);
    int require_hm = 0;
    if(pbp_on && pbp_periods == pbp_periods_thrd) {
      for(int j = baseidx_prev; j < baseidx; j ++) {
        if(pbp_switch_state < 1.0) {
          pbp_switch_state += pbp_switch_rate;
          require_hm = 1;
        }
        y_mix[j] = pbp_switch_state;
      }
    } else
    if(! pbp_on && pbp_periods == 0) {
      for(int j = baseidx_prev; j < baseidx; j ++) {
        if(pbp_switch_state > 0) {
          pbp_switch_state -= pbp_switch_rate;
          require_hm = 1;
        }
        y_mix[j] = pbp_switch_state;
      }
    } else {
      for(int j = baseidx_prev; j < baseidx; j ++)
        y_mix[j] = pbp_switch_state;
    }
    baseidx_prev = baseidx;

    if(pbp_on && pbp_periods == pbp_periods_thrd && (! require_hm)) continue;

    // harmonic model synthesis
    if(llsm_container_get(src_frame, LLSM_FRAME_HM) == NULL)
      llsm_frame_tolayer0(src_frame, chunk -> conf);
    llsm_hmframe* hm = llsm_container_get(src_frame, LLSM_FRAME_HM);
    if(hm == NULL) continue;
    int nhar = min(maxnhar, hm -> nhar);
    FP_TYPE* yi = llsm_synthesize_harmonic_frame_auto(options,
      hm -> ampl, hm -> phse, nhar, f0[i] / fs, nwin);
    for(int j = 0; j < nwin; j ++) {
      yi[j] *= w[j];
      int idx = baseidx + j - nwin / 2;
      if(idx >= 0 && idx < ny)
        y_hm[idx] += yi[j];
    }
    free(yi);
  }
  free(w);

  for(int i = 0; i < ny; i ++)
    y_mix[i] = y_hm[i] * (1.0 - y_mix[i]) + y_pbp[i] * y_mix[i];

  free(y_hm);
  free(y_pbp);
  return y_mix;
}

static FP_TYPE* llsm_synthesize_noise_envelope(llsm_soptions* options,
  llsm_chunk* chunk, int channel, FP_TYPE* f0, int nfrm, FP_TYPE thop,
  FP_TYPE fs, int ny) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  int nwin = round(thop * 2.0 * fs);
  FP_TYPE* w = hanning(nwin);
  llsm_hmframe* unvoiced_hm = llsm_create_hmframe(0);
  for(int i = 0; i < nfrm; i ++) {
    llsm_nmframe* nm = llsm_container_get(chunk -> frames[i], LLSM_FRAME_NM);
    llsm_hmframe* hm = f0[i] > 0 ? nm -> eenv[channel] : unvoiced_hm;
    FP_TYPE* yi = llsm_synthesize_harmonic_frame_auto(options,
      hm -> ampl, hm -> phse, hm -> nhar, f0[i] / fs, nwin);
    // Make sure the envelope is positive.
    FP_TYPE offset = nm -> edc[channel];
    for(int j = 0; j < nwin; j ++)
      yi[j] = max(yi[j] + offset, 1e-8);
    for(int j = 0; j < nwin; j ++) {
      yi[j] *= w[j];
      int idx = round((i - 1) * thop * fs + j);
      if(idx >= 0 && idx < ny)
        y[idx] += yi[j];
    }
    free(yi);
  }
  llsm_delete_hmframe(unvoiced_hm);
  free(w);
  return y;
}

static void llsm_analyze_noise_psd(llsm_aoptions* options, FP_TYPE* x,
  FP_TYPE* x_res, int nx, FP_TYPE fs, int nfrm, llsm_chunk* dst_chunk) {
  int nwin = round(options -> thop * 4 * fs);
  int nfft = pow(2, ceil(log2(nwin)));
  int nspec = nfft / 2 + 1;

  // compute spectral envelope (harmonic + noise power)
  int nfft_spgm = pow(2, ceil(log2(0.03 * fs)));
  FP_TYPE** spgm = malloc2d(nfrm, nfft_spgm / 2 + 1, sizeof(FP_TYPE));
  int* center = calloc(nfrm, sizeof(int));
  int* winsize_spgm = calloc(nfrm, sizeof(int));
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* f0 = llsm_container_get(dst_chunk -> frames[i], LLSM_FRAME_F0);
    winsize_spgm[i] = f0 == NULL || f0[0] == 0 ? nwin : fs / f0[0] * 3;
    center[i] = round(i * options -> thop * fs);
  }
  llsm_compute_spectrogram(x, nx, center, winsize_spgm, nfrm, nfft_spgm,
    "hanning", spgm, NULL);
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* f0 = llsm_container_get(dst_chunk -> frames[i], LLSM_FRAME_F0);
    FP_TYPE f0_scaled = (f0 == NULL || f0[0] == 0 ? 200 : f0[0]) / fs;
    FP_TYPE* env = spec2env(spgm[i], nfft_spgm, f0_scaled, NULL);
    for(int j = 0; j < nspec; j ++) {
      int idx = j * nfft_spgm / nfft;
      spgm[i][j] = env[idx] * 2; // magnitude to power (log)
    }
    free(env);
  }
  free(winsize_spgm);

  // transposed PSD spectrogram
  FP_TYPE** spgm_psd = malloc2d(nspec, nfrm, sizeof(FP_TYPE));
  // PSD residual vectors
  FP_TYPE** spgm_res = malloc2d(nfrm, nspec, sizeof(FP_TYPE));
  // compute noise PSD
  FP_TYPE* psdvec = calloc(nspec, sizeof(FP_TYPE));
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* xfrm = fetch_frame(x_res, nx, center[i], nwin);
    llsm_estimate_psd(xfrm, nwin, nfft, psdvec);
    for(int j = 0; j < nspec; j ++)
      spgm_psd[j][i] = log(max(1e-10, psdvec[j]));
    free(xfrm);
  }
  FP_TYPE* Q = calloc(nfrm, sizeof(FP_TYPE)); // process variance
  FP_TYPE* R = calloc(nfrm, sizeof(FP_TYPE)); // observation variance
  FP_TYPE* P = calloc(nfrm, sizeof(FP_TYPE)); // forward posterior
  for(int i = 0; i < nfrm; i ++) R[i] = LOGCHI2VAR;
  for(int j = 0; j < nspec; j ++) {
    // moving statistics -> process variance
    for(int i = 0; i < nfrm; i ++) {
      FP_TYPE m1 = 0;
      FP_TYPE m2 = 0;
      for(int k = -1; k <= 1; k ++) {
        int idx = min(nfrm - 1, max(0, i + k));
        m1 += spgm[idx][j];
        m2 += spgm[idx][j] * spgm[idx][j];
      }
      Q[i] = max(1e-8, m2 / 3 - m1 * m1 / 9);
    }
    // smoothen PSD along time and extract the residual
    FP_TYPE* y = kalmanf1d(spgm_psd[j], Q, R, nfrm, P, NULL);
    FP_TYPE* s = kalmans1d(y, P, Q, nfrm);
    for(int i = 0; i < nfrm; i ++) {
      spgm_res[i][j] = spgm_psd[j][i] - s[i];
      spgm_psd[j][i] = s[i] + EULERGAMMA; // bias removal
    }
    free(y); free(s);
  }
  free(P); free(Q); free(R);

  FP_TYPE* dst_axis = linspace(0, fs / 2.0, options -> npsd);
  for(int i = 0; i < nfrm; i ++) {
    llsm_nmframe* dst_nm = llsm_container_get(
      dst_chunk -> frames[i], LLSM_FRAME_NM);
    for(int j = 0; j < nspec; j ++) psdvec[j] = spgm_psd[j][i];
    FP_TYPE* dst_psd = interp1u(
      0, fs / 2.0, psdvec, nspec, dst_axis, options -> npsd);
    FP_TYPE* dst_res = interp1u(
      0, fs / 2.0, spgm_res[i], nspec, dst_axis, options -> npsd);
    FP_TYPE* resvec = llsm_create_fparray(options -> npsd);
    // The PSD is squared and hence 10 * log10(.)
    // -120 dB noise floor for underflow protection.
    for(int j = 0; j < options -> npsd; j ++) {
      resvec[j] = LOG2IN(dst_res[j]);
      dst_psd[j] = exp(dst_psd[j]);
      dst_nm -> psd[j] = 10.0 * log10(dst_psd[j] * 44100 / fs + 1e-12);
    }
    llsm_container_attach(dst_chunk -> frames[i], LLSM_FRAME_PSDRES,
      resvec, llsm_delete_fparray, llsm_copy_fparray);
    free(dst_psd); free(dst_res);
  }
  free(dst_axis);
  free2d(spgm_psd, nspec);
  free2d(spgm_res, nfrm);
  free2d(spgm, nfrm);
  free(center);
  free(psdvec);
}

static void llsm_analyze_noise_envelope(llsm_aoptions* options,
  FP_TYPE* x, FP_TYPE* x_res, int nx, FP_TYPE fs, FP_TYPE* f0,
  int nfrm, llsm_chunk* dst_chunk) {

  int*      tmp_nhar = calloc(nfrm, sizeof(int));
  FP_TYPE** tmp_ampl = calloc(nfrm, sizeof(FP_TYPE*));
  FP_TYPE** tmp_phse = calloc(nfrm, sizeof(FP_TYPE*));

  FP_TYPE*  tmp_dc   = calloc(nfrm, sizeof(FP_TYPE));
  int*      center   = calloc(nfrm, sizeof(int));
  int*      nwin     = calloc(nfrm, sizeof(int));
  for(int i = 0; i < nfrm; i ++) {
    center[i] = round(i * options -> thop * fs);
    nwin[i] = round((f0[i] == 0 ? options -> thop * 2 : 2.0 / f0[i]) * fs);
  }

  for(int c = 0; c < options -> nchannel; c ++) {
    FP_TYPE fmin = c == 0 ? 0 : options -> chanfreq[c - 1];
    FP_TYPE fmax = c == options -> nchannel - 1 ?
                   fs / 2.0 : options -> chanfreq[c];
    // Trick: extract envelope from the original waveform in high frequencies
    //   where the residual is often smeared due to harmonic analysis errors.
    FP_TYPE* ce = llsm_subband_energy(
      fmin > 6000.0 ? x : x_res, nx, fmin / fs, fmax / fs);
    // Perform harmonic analysis on a squared signal and extract the lower
    //   harmonics is roughly equivalent to modeling the RMS envelope.
    llsm_harmonic_analysis(ce, nx, fs, f0, nfrm, options -> thop,
      options -> rel_winsize, options -> maxnhar_e, options -> hm_method,
      tmp_nhar, tmp_ampl, tmp_phse);
    llsm_compute_dc(ce, nx, center, nwin, nfrm, tmp_dc);
    // Store the results.
    for(int i = 0; i < nfrm; i ++) {
      llsm_nmframe* dst_nm = llsm_container_get(dst_chunk -> frames[i],
        LLSM_FRAME_NM);
      dst_nm -> edc[c] = tmp_dc[i];
      if(f0[i] == 0) continue;
      llsm_hmframe* hm = llsm_create_hmframe(tmp_nhar[i]);
      memcpy(hm -> ampl, tmp_ampl[i], tmp_nhar[i] * sizeof(FP_TYPE));
      memcpy(hm -> phse, tmp_phse[i], tmp_nhar[i] * sizeof(FP_TYPE));
      llsm_copy_hmframe_inplace(dst_nm -> eenv[c], hm);
      llsm_delete_hmframe(hm);
    }

    for(int i = 0; i < nfrm; i ++) {
      free(tmp_ampl[i]); tmp_ampl[i] = NULL;
      free(tmp_phse[i]); tmp_phse[i] = NULL;
    }
    free(ce);
  }

  free2d(tmp_ampl, nfrm); free2d(tmp_phse, nfrm); free(tmp_nhar);
  free(tmp_dc); free(center); free(nwin);
}

static void llsm_analyze_noise(llsm_aoptions* options, FP_TYPE* x,
  FP_TYPE* x_res, int nx, FP_TYPE fs, FP_TYPE* f0, int nfrm,
  llsm_chunk* dst_chunk) {
  llsm_analyze_noise_psd(options, x, x_res, nx, fs, nfrm, dst_chunk);
  llsm_analyze_noise_envelope(options, x, x_res, nx, fs, f0, nfrm, dst_chunk);
}

llsm_chunk* llsm_analyze(llsm_aoptions* options, FP_TYPE* x, int nx,
  FP_TYPE fs, FP_TYPE* f0, int nfrm, FP_TYPE** x_ap) {

  llsm_container* conf = llsm_aoptions_toconf(options, fs / 2.0);
  ((int*)llsm_container_get(conf, LLSM_CONF_NFRM))[0] = nfrm;
  llsm_chunk* ret = llsm_create_chunk(conf, 1);
  llsm_delete_container(conf); // conf gets copied into ret; no longer needed.
  conf = ret -> conf;

  if(options -> f0_refine)
    llsm_refine_f0(x, nx, fs, f0, nfrm, options -> thop);

  // set F0 for all the frames
  for(int i = 0; i < nfrm; i ++) {
    FP_TYPE* if0 = llsm_container_get(ret -> frames[i], LLSM_FRAME_F0);
    if0[0] = f0[i];
  }

  // harmonic analysis and residual extraction
  llsm_analyze_harmonics(options, x, nx, fs, f0, nfrm, ret);
  FP_TYPE* x_sin = llsm_synthesize_harmonics_l0(NULL, ret, f0, nfrm,
    options -> thop, fs, nx);
  FP_TYPE* x_res = calloc(nx, sizeof(FP_TYPE));
  for(int i = 0; i < nx; i ++) x_res[i] = x[i] - x_sin[i];
  free(x_sin);
  if(x_ap != NULL) *x_ap = x_res;

  // noise analysis
  llsm_analyze_noise(options, x, x_res, nx, fs, f0, nfrm, ret);

  if(x_ap == NULL)
    free(x_res);
  return ret;
}

int llsm_conf_checklayer0(llsm_container* src) {
  int* nfrm = llsm_container_get(src, LLSM_CONF_NFRM);
  FP_TYPE* thop = llsm_container_get(src, LLSM_CONF_THOP);
  int* npsd = llsm_container_get(src, LLSM_CONF_NPSD);
  FP_TYPE* fnyq = llsm_container_get(src, LLSM_CONF_FNYQ);
  int* nchannel = llsm_container_get(src, LLSM_CONF_NCHANNEL);
  FP_TYPE* chanfreq = llsm_container_get(src, LLSM_CONF_CHANFREQ);
  if(nfrm == NULL || thop == NULL || npsd == NULL ||
     fnyq == NULL || nchannel == NULL || chanfreq == NULL) return 0;
  return 1;
}

static int llsm_synthesis_check_integrity(llsm_chunk* src) {
  if(! llsm_conf_checklayer0(src -> conf)) return 0;
  int* nfrm = llsm_container_get(src -> conf, LLSM_CONF_NFRM);
  for(int i = 0; i < *nfrm; i ++)
    if(! llsm_frame_checklayer0(src -> frames[i]) &&
       ! llsm_frame_checklayer1(src -> frames[i]))
      return 0;
  return 1;
}

static FP_TYPE* llsm_synthesize_noise_excitation(llsm_soptions* options,
  llsm_chunk* src, FP_TYPE* f0, int nfrm, FP_TYPE thop, FP_TYPE fs, int ny) {
  FP_TYPE* y = calloc(ny, sizeof(FP_TYPE));
  FP_TYPE* chanfreq = llsm_container_get(src -> conf, LLSM_CONF_CHANFREQ);
  int nchannel = *((int*)llsm_container_get(src -> conf, LLSM_CONF_NCHANNEL));
  for(int c = 0; c < nchannel; c ++) {
    FP_TYPE fmin = c == 0 ? 0 : chanfreq[c - 1];
    FP_TYPE fmax = c == nchannel - 1 ? fs / 2.0 : chanfreq[c];
    if(fmin >= fs / 2.0) break;
    FP_TYPE* x = llsm_generate_bandlimited_noise(ny, fmin / fs, fmax / fs);
    FP_TYPE* env = llsm_synthesize_noise_envelope(options, src, c, f0, nfrm,
      thop, fs, ny);
    for(int i = 0; i < ny; i ++) {
      x[i] *= sqrt(env[i]);
      y[i] += x[i];
    }
    free(env);
    free(x);
  }
  return y;
}

static FP_TYPE* llsm_filter_noise(llsm_chunk* src, int nfrm, FP_TYPE thop,
  FP_TYPE fs, FP_TYPE* x, int nx) {
  const int nfade = 16;
  int nwin = round(thop * fs * 2);
  FP_TYPE* w = hanning(nwin);
  FP_TYPE wsqr = 0;
  for(int i = 0; i < nwin; i ++)
    wsqr += w[i] * w[i];

  // at least 20% padding
  int nfft = pow(2, ceil(log2(nwin * 1.2 + nfade * 2)));
  int nspec = nfft / 2 + 1;
  FP_TYPE* psd = calloc(nspec, sizeof(FP_TYPE));
  FP_TYPE* fftbuffer = calloc(nfft * 4, sizeof(FP_TYPE));
  FP_TYPE* x_re = fftbuffer;
  FP_TYPE* x_im = fftbuffer + nfft;

  int npsd = *((int*)llsm_container_get(src -> conf, LLSM_CONF_NPSD));
  FP_TYPE fnyq = *((FP_TYPE*)llsm_container_get(src -> conf, LLSM_CONF_FNYQ));

  FP_TYPE* y = calloc(nx, sizeof(FP_TYPE));
  FP_TYPE* src_axis = linspace(0, fnyq, npsd);
  FP_TYPE* src_psd = calloc(npsd, sizeof(FP_TYPE));
  // STFT -> PSD -> diff -> filter -> ISTFT
  for(int i = 0; i < nfrm; i ++) {
    llsm_nmframe* nm = llsm_container_get(src -> frames[i], LLSM_FRAME_NM);
    FP_TYPE* resvec = llsm_container_get(src -> frames[i], LLSM_FRAME_PSDRES);

    // STFT
    int center = round(i * thop * fs);
    FP_TYPE* xfrm = fetch_frame(x, nx, center, nwin);
    for(int j = 0; j < nwin; j ++) xfrm[j] *= w[j];
    for(int j = 0; j < nfft; j ++) x_re[j] = 0;
    for(int j = 0; j < nwin; j ++) x_re[j - nwin / 2 + nfft / 2] = xfrm[j];
    fft(x_re, NULL, x_re, x_im, nfft, fftbuffer + nfft * 2);

    // PSD -> diff
    llsm_fft_to_psd(x_re, x_im, nfft, wsqr, psd);
    FP_TYPE* env = moving_avg(psd, nspec, 3);
    for(int j = 0; j < npsd; j ++) src_psd[j] = nm -> psd[j];
    if(resvec != NULL)
    for(int j = 0; j < npsd; j ++)
      src_psd[j] += resvec[j] - LOG2IN(LOGRESBIAS);
    FP_TYPE* H = llsm_spectrum_from_envelope(
      src_axis, src_psd, npsd, nspec - 1, fs / 2.0);
    for(int j = 0; j < nspec; j ++)
      H[j] = exp(DB2LOG(H[j])) / sqrt(env[j] * 44100 / fs + 1e-8);

    // filter
    for(int j = 0; j < nspec - 1; j ++) {
      x_re[j] *= H[j]; x_im[j] *= H[j];
    }
    x_re[nspec - 1] = x_re[nspec - 2]; complete_symm(x_re, nfft);
    x_im[nspec - 1] = x_im[nspec - 2]; complete_asymm(x_im, nfft);

    // ISTFT
    ifft(x_re, x_im, x_re, NULL, nfft, fftbuffer + nfft * 2);
    for(int j = 0; j < nfade; j ++) {
      x_re[j] *= (FP_TYPE)j / nfade;
      x_re[nfft - j - 1] *= 1.0 - (FP_TYPE)j / nfade;
    }
    for(int j = 0; j < nfft; j ++) {
      int idx = center + j - nfft / 2;
      if(idx >= 0 && idx < nx)
        y[idx] += x_re[j];
    }

    free(H); free(env);
    free(xfrm);
  }

  free(fftbuffer); free(psd);
  free(src_axis); free(src_psd);
  free(w);
  return y;
}

llsm_output* llsm_synthesize(llsm_soptions* options, llsm_chunk* src) {
  if(! llsm_synthesis_check_integrity(src)) return NULL;
  int nfrm;
  FP_TYPE thop = *((FP_TYPE*)llsm_container_get(src -> conf, LLSM_CONF_THOP));
  FP_TYPE fs = options -> fs;
  FP_TYPE* f0 = llsm_chunk_getf0(src, & nfrm);

  int ny = round((nfrm + 1) * thop * fs);
  llsm_output* ret = malloc(sizeof(llsm_output));
  ret -> ny = ny;
  ret -> fs = fs;

  FP_TYPE* y_sin = llsm_synthesize_harmonics(options, src, f0, nfrm,
    thop, fs, ny);
  ret -> y_sin = y_sin;

  FP_TYPE* y_exc = llsm_synthesize_noise_excitation(options, src, f0, nfrm,
    thop, fs, ny);
  FP_TYPE* y_nos = llsm_filter_noise(src, nfrm, thop, fs, y_exc, ny);
  ret -> y_noise = y_nos;

  ret -> y = calloc(ny, sizeof(FP_TYPE));
  for(int i = 0; i < ny; i ++)
    ret -> y[i] = y_sin[i] + y_nos[i];

  free(y_exc);
  free(f0);
  return ret;
}

void llsm_delete_output(llsm_output* dst) {
  if(dst == NULL) return;
  free(dst -> y);
  free(dst -> y_sin);
  free(dst -> y_noise);
  free(dst);
}

FP_TYPE* llsm_chunk_getf0(llsm_chunk* src, int* dst_nfrm) {
  int* nfrm = llsm_container_get(src -> conf, LLSM_CONF_NFRM);
  if(nfrm == NULL) return NULL;
  FP_TYPE* f0 = calloc(*nfrm, sizeof(FP_TYPE));
  *dst_nfrm = *nfrm;
  for(int i = 0; i < *nfrm; i ++) {
    FP_TYPE* if0 = llsm_container_get(src -> frames[i], LLSM_FRAME_F0);
    if(if0 != NULL)
      f0[i] = if0[0];
  }
  return f0;
}

void llsm_chunk_phasesync_rps(llsm_chunk* dst, int layer1_based) {
  int* nfrm = llsm_container_get(dst -> conf, LLSM_CONF_NFRM);
  if(nfrm == NULL) return;
  for(int i = 0; i < *nfrm; i ++)
    llsm_frame_phasesync_rps(dst -> frames[i], layer1_based);
}

void llsm_chunk_phasepropagate(llsm_chunk* dst, int sign) {
  int nfrm = 0;
  FP_TYPE* f0 = llsm_chunk_getf0(dst, & nfrm);
  FP_TYPE* thop = llsm_container_get(dst -> conf, LLSM_CONF_THOP);
  if(thop == NULL || f0 == NULL) return;
  FP_TYPE* delta_phase = cumsum(f0, nfrm);
  for(int i = 0; i < nfrm; i ++) {
    delta_phase[i] *= *thop * sign * 2.0 * M_PI;
    llsm_frame_phaseshift(dst -> frames[i], delta_phase[i]);
  }
  free(delta_phase);
  free(f0);
}

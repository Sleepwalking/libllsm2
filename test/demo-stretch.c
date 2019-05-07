#include "../llsm.h"
#include "../external/libpyin/pyin.h"
#include "../external/ciglet/ciglet.h"

// circular interpolation of two radian values
static FP_TYPE linterpc(FP_TYPE a, FP_TYPE b, FP_TYPE ratio) {
  FP_TYPE ax = cos_2(a);
  FP_TYPE ay = sin_2(a);
  FP_TYPE bx = cos_2(b);
  FP_TYPE by = sin_2(b);
  FP_TYPE cx = linterp(ax, bx, ratio);
  FP_TYPE cy = linterp(ay, by, ratio);
  return atan2(cy, cx);
}

static void interp_nmframe(llsm_nmframe* dst, llsm_nmframe* src,
  FP_TYPE ratio, int dst_voiced, int src_voiced) {
  for(int i = 0; i < dst -> npsd; i ++)
    dst -> psd[i] = linterp(dst -> psd[i], src -> psd[i], ratio);

  for(int b = 0; b < dst -> nchannel; b ++) {
    llsm_hmframe* srceenv = src -> eenv[b];
    llsm_hmframe* dsteenv = dst -> eenv[b];
    dst -> edc[b] = linterp(dst -> edc[b], src -> edc[b], ratio);
    int b_minnhar = min(srceenv -> nhar, dsteenv -> nhar);
    int b_maxnhar = max(srceenv -> nhar, dsteenv -> nhar);
    if(dsteenv -> nhar < b_maxnhar) {
      dsteenv -> ampl = realloc(dsteenv -> ampl, sizeof(FP_TYPE) * b_maxnhar);
      dsteenv -> phse = realloc(dsteenv -> phse, sizeof(FP_TYPE) * b_maxnhar);
    }
    for(int i = 0; i < b_minnhar; i ++) {
      dsteenv -> ampl[i] =
        linterp(dsteenv -> ampl[i], srceenv -> ampl[i], ratio);
      dsteenv -> phse[i] =
        linterpc(dsteenv -> phse[i], srceenv -> phse[i], ratio);
    }
    if(b_maxnhar == srceenv -> nhar) {
      for(int i = b_minnhar; i < b_maxnhar; i ++) {
        dsteenv -> ampl[i] = srceenv -> ampl[i];
        dsteenv -> phse[i] = srceenv -> phse[i];
      }
    }
    dsteenv -> nhar = b_maxnhar;
  }
}

#define LOG2DB (20.0 / 2.3025851)
#define mag2db(x) (log_2(x) * LOG2DB)

// dst <- (dst &> src)
static void interp_llsm_frame(llsm_container* dst, llsm_container* src,
  FP_TYPE ratio) {
# define EPS 1e-8
  FP_TYPE dst_f0 = *((FP_TYPE*)llsm_container_get(dst, LLSM_FRAME_F0));
  FP_TYPE src_f0 = *((FP_TYPE*)llsm_container_get(src, LLSM_FRAME_F0));
  llsm_nmframe* dst_nm = llsm_container_get(dst, LLSM_FRAME_NM);
  llsm_nmframe* src_nm = llsm_container_get(src, LLSM_FRAME_NM);
  FP_TYPE* src_rd = llsm_container_get(src, LLSM_FRAME_RD);
  FP_TYPE* dst_rd = llsm_container_get(dst, LLSM_FRAME_RD);
  FP_TYPE* dst_vsphse = llsm_container_get(dst, LLSM_FRAME_VSPHSE);
  FP_TYPE* src_vsphse = llsm_container_get(src, LLSM_FRAME_VSPHSE);
  FP_TYPE* dst_vtmagn = llsm_container_get(dst, LLSM_FRAME_VTMAGN);
  FP_TYPE* src_vtmagn = llsm_container_get(src, LLSM_FRAME_VTMAGN);

  // always take the frequency of the voiced frame
  llsm_container* voiced = dst_f0 <= 0 && src_f0 <= 0 ? NULL :
    (src_f0 > 0 ? src : dst);
  int bothvoiced = dst_f0 > 0 && src_f0 > 0;

  int dstnhar = dst_vsphse == NULL ? 0 : llsm_fparray_length(dst_vsphse);
  int srcnhar = src_vsphse == NULL ? 0 : llsm_fparray_length(src_vsphse);
  int maxnhar = max(dstnhar, srcnhar);
  int minnhar = min(dstnhar, srcnhar);

  if(! bothvoiced && voiced == src) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(src_f0),
      llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(*src_rd),
      llsm_delete_fp, llsm_copy_fp);
  } else
  if(voiced == NULL) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(0),
      llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(1.0),
      llsm_delete_fp, llsm_copy_fp);
  }
  int nspec = dst_vtmagn != NULL ? llsm_fparray_length(dst_vtmagn) :
    (src_vtmagn != NULL ? llsm_fparray_length(src_vtmagn) : 0);

  if(bothvoiced) {
    llsm_container_attach(dst, LLSM_FRAME_F0, llsm_create_fp(
      linterp(dst_f0, src_f0, ratio)), llsm_delete_fp, llsm_copy_fp);
    llsm_container_attach(dst, LLSM_FRAME_RD, llsm_create_fp(
      linterp(*dst_rd, *src_rd, ratio)), llsm_delete_fp, llsm_copy_fp);

    FP_TYPE* vsphse = llsm_create_fparray(maxnhar);
    FP_TYPE* vtmagn = llsm_create_fparray(nspec);
    for(int i = 0; i < minnhar; i ++)
      vsphse[i] = linterpc(dst_vsphse[i], src_vsphse[i], ratio);
    for(int i = 0; i < nspec; i ++)
      vtmagn[i] = linterp(dst_vtmagn[i], src_vtmagn[i], ratio);
    if(dstnhar < srcnhar)
      for(int i = minnhar; i < maxnhar; i ++)
        vsphse[i] = src_vsphse[i];

    dst_vsphse = vsphse;
    dst_vtmagn = vtmagn;
    llsm_container_attach(dst, LLSM_FRAME_VSPHSE, dst_vsphse,
      llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(dst, LLSM_FRAME_VTMAGN, dst_vtmagn,
      llsm_delete_fparray, llsm_copy_fparray);
  } else if(voiced == src) {
    dst_vsphse = llsm_copy_fparray(src_vsphse);
    dst_vtmagn = llsm_copy_fparray(src_vtmagn);
    llsm_container_attach(dst, LLSM_FRAME_VSPHSE, dst_vsphse,
      llsm_delete_fparray, llsm_copy_fparray);
    llsm_container_attach(dst, LLSM_FRAME_VTMAGN, dst_vtmagn,
      llsm_delete_fparray, llsm_copy_fparray);
    FP_TYPE fade = mag2db(max(EPS, ratio));
    for(int i = 0; i < nspec; i ++) dst_vtmagn[i] += fade;
  } else {
    FP_TYPE fade = mag2db(max(EPS, 1.0 - ratio));
    for(int i = 0; i < nspec; i ++) dst_vtmagn[i] += fade;
  }
  for(int i = 0; i < nspec; i ++) dst_vtmagn[i] = max(-80, dst_vtmagn[i]);

  interp_nmframe(dst_nm, src_nm, ratio, dst_f0 > 0, src_f0 > 0);
# undef EPS
}

int main() {
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
  opt_a -> hm_method = LLSM_AOPTION_HMCZT;
  llsm_soptions* opt_s = llsm_create_soptions(fs);
  llsm_chunk* chunk = llsm_analyze(opt_a, x, nx, fs, f0, nfrm, NULL);

  llsm_output* out0 = llsm_synthesize(opt_s, chunk);
  wavwrite(out0 -> y, out0 -> ny, opt_s -> fs, 24, "test/demo-stretch-orig.wav");
  wavwrite(out0 -> y_sin, out0 -> ny, opt_s -> fs, 24, "test/demo-stretch-orig-sin.wav");
  wavwrite(out0 -> y_noise, out0 -> ny, opt_s -> fs, 24, "test/demo-stretch-orig-noise.wav");
  llsm_delete_output(out0);

  llsm_chunk_tolayer1(chunk, 2048);
  llsm_chunk_phasepropagate(chunk, -1);

  int nfrm_new = nfrm * 2;
  llsm_container* conf_new = llsm_copy_container(chunk -> conf);
  llsm_container_attach(conf_new, LLSM_CONF_NFRM,
    llsm_create_int(nfrm_new), llsm_delete_int, llsm_copy_int);
  llsm_chunk* chunk_new = llsm_create_chunk(conf_new, 0);
  llsm_delete_container(conf_new);

  for(int i = 0; i < nfrm_new; i ++) {
    FP_TYPE mapped = (FP_TYPE)i * nfrm / nfrm_new;
    int base = mapped;
    FP_TYPE ratio = mapped - base;
    int residx = base + rand() % 5 - 2;
    residx = max(0, min(nfrm - 1, residx));
    base = min(base, nfrm - 2);
    chunk_new -> frames[i] = llsm_copy_container(chunk -> frames[base]);
    interp_llsm_frame(
      chunk_new -> frames[i], chunk -> frames[base + 1], ratio);
    FP_TYPE* resvec = llsm_container_get(chunk -> frames[residx],
      LLSM_FRAME_PSDRES);
    llsm_container_attach(chunk_new -> frames[i], LLSM_FRAME_PSDRES,
      llsm_copy_fparray(resvec), llsm_delete_fparray, llsm_copy_fparray);
  }
  llsm_chunk_tolayer0(chunk_new);
  llsm_chunk_phasepropagate(chunk_new, 1);

  llsm_output* out = llsm_synthesize(opt_s, chunk_new);
  wavwrite(out -> y, out -> ny, opt_s -> fs, 24, "test/demo-stretch.wav");
  wavwrite(out -> y_noise, out -> ny, opt_s -> fs, 24, "test/demo-stretch-noise.wav");
  llsm_delete_output(out);

  llsm_delete_chunk(chunk);
  llsm_delete_chunk(chunk_new);
  llsm_delete_aoptions(opt_a);
  llsm_delete_soptions(opt_s);
  free(f0);
  free(x);
  return 0;
}

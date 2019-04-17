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

/** @file */

#ifndef LLSM_H
#define LLSM_H

/** @brief Function pointer to destructors (e.g. llsm_delete_container). */
typedef void (*llsm_fdestructor)(void*);

/** @brief Function pointer to copy constructors (e.g. llsm_copy_container). */
typedef void* (*llsm_fcopy)(void*);

/** @defgroup group_utils Container-related Utilities
 *  @{ */
FP_TYPE* llsm_create_fp(FP_TYPE x);
int*     llsm_create_int(int x);
FP_TYPE* llsm_create_fparray(int size);
FP_TYPE* llsm_copy_fp(FP_TYPE* src);
int*     llsm_copy_int(int* src);
FP_TYPE* llsm_copy_fparray(FP_TYPE* src);
void     llsm_delete_fp(FP_TYPE* dst);
void     llsm_delete_int(int* dst);
void     llsm_delete_fparray(FP_TYPE* dst);
int      llsm_fparray_length(FP_TYPE* src);
/** @} */

/** @defgroup group_container llsm_container
 *  @{ */
/** @brief A generic container that can store multiple structures of
 *    different types. */
typedef struct {
  void** members;
  llsm_fdestructor* destructors;
  llsm_fcopy* copyctors;
  int nmember;
} llsm_container;

/** @brief Create an empty container object. */
llsm_container* llsm_create_container(int nmember);
/** @brief Copy-construct a container object from an existing one. */
llsm_container* llsm_copy_container(llsm_container* src);
/** @brief In-place version of llsm_copy_container. */
void llsm_copy_container_inplace(llsm_container* dst, llsm_container* src);
/** @brief Delete and free a container.
 *
 *  For each member, if a destructor is specified, call the destructor before
 *    deleting the container. */
void llsm_delete_container(llsm_container* dst);

/** @brief Get the member at index from a container. */
void* llsm_container_get(llsm_container* src, int index);
/** @brief Attach (shallow-copy) an object to a container.
 *
 *  If destructor is NULL, the added member will not be deleted when
 *    llsm_delete_container or llsm_container_remove is called.
 *  If copyctor is NULL, the added member will be shallow-copied when
 *    llsm_copy_container is called.
 *  If index is greater than the size of the existing list, the container
 *    will be automatically expanded. */
#define llsm_container_attach(dst, index, ptr, dtor, copyctor) \
        llsm_container_attach_(dst, index, ptr, (llsm_fdestructor)dtor, \
          (llsm_fcopy)copyctor)
void llsm_container_attach_(llsm_container* dst, int index, void* ptr,
  llsm_fdestructor dtor, llsm_fcopy copyctor);
/** @brief Remove a member from a container.
 *
 *  If a destructor is specified, call the destructor to delete and free
 *    the member and then set the member to NULL. */
void llsm_container_remove(llsm_container* dst, int index);
/** @} */

/** @defgroup group_frame_index Indexing Macros for LLSM Frame
 *  @brief List of macros indicating indices of parameters in a LLSM frame.
 *  @{ */
#define LLSM_FRAME_F0        0  /**< fundamental frequency (FP_TYPE) */
#define LLSM_FRAME_HM        1  /**< harmonic model (llsm_hmframe) */
#define LLSM_FRAME_NM        2  /**< noise model (llsm_nmframe) */
#define LLSM_FRAME_RD       10  /**< Rd parameter (FP_TYPE) */
#define LLSM_FRAME_VTMAGN   11  /**< vocal tract magnitude response
                                     (FP_TYPE*, dB) */
#define LLSM_FRAME_VSPHSE   12  /**< vocal source harmonic phase (FP_TYPE*) */
/** @} */

/** @defgroup group_config_index Indexing Macros for LLSM Configuration
 *  @brief List of macros indicating indices of attributes in the model
 *     configuration.
 *  @{ */
#define LLSM_CONF_NFRM       0  /**< number of frames (int) */
#define LLSM_CONF_THOP       1  /**< time interval (FP_TYPE, sec) */
#define LLSM_CONF_MAXNHAR    2  /**< maximum number of harmonics (int) */
#define LLSM_CONF_MAXNHAR_E  3  /**< maximum number of harmonics for noise
                                     envelope (int) */
#define LLSM_CONF_NPSD       4  /**< size of noise PSD vector (int) */
#define LLSM_CONF_NOSWARP    5  /**< noise PSD warping constant (FP_TYPE) */
#define LLSM_CONF_FNYQ       6  /**< Nyquist frequency (FP_TYPE, Hz) */
#define LLSM_CONF_NCHANNEL   7  /**< number of noise channels (int) */
#define LLSM_CONF_CHANFREQ   8  /**< frequencies of noise channels
                                     (FP_TYPE*, Hz) */
#define LLSM_CONF_NSPEC     10  /**< size of magnitude response (int) */
#define LLSM_CONF_LIPRADIUS 11  /**< assumed radius of lip opening
                                     (FP_TYPE, cm)*/
/** @} */

/** @defgroup group_hmframe llsm_hmframe
 *  @{ */
/** @brief Harmonic model parameters (in one frame). */
typedef struct {
  FP_TYPE* ampl;       /**< harmonic amplitude (linear) vector */
  FP_TYPE* phse;       /**< harmonic phase (radian) vector */
  int      nhar;       /**< number of harmonics */
} llsm_hmframe;

/** @brief Create an empty harmonic model frame with nhar harmonics. */
llsm_hmframe* llsm_create_hmframe(int nhar);
/** @brief Copy-construct a harmonic model frame from an existing one. */
llsm_hmframe* llsm_copy_hmframe(llsm_hmframe* src);
/** @brief In-place version of llsm_copy_hmframe. */
void llsm_copy_hmframe_inplace(llsm_hmframe* dst, llsm_hmframe* src);
/** @brief Delete and free a harmonic model frame. */
void llsm_delete_hmframe(llsm_hmframe* dst);
/** @brief Rotate the phases by (theta * 1-based index of the harmonic). */
void llsm_hmframe_phaseshift(llsm_hmframe* dst, FP_TYPE theta);
/** @brief Compute the noise-equivalent harmonic power spectral density. */
FP_TYPE* llsm_hmframe_harpsd(llsm_hmframe* src, int db_scale);
/** @} */

/** @defgroup group_nmframe llsm_nmframe
 *  @{ */
/** @brief Noise model parameters (in one frame). */
typedef struct {
  llsm_hmframe** eenv; /**< the harmonic model describing the noise envelope
                            in each channel */
  FP_TYPE* edc;        /**< the short-time mean of the noise envelope in each
                            channel */
  FP_TYPE* psd;        /**< power spectral density (dB) vector */
  int npsd;            /**< size of psd */
  int nchannel;        /**< number of channels */
} llsm_nmframe;

/** @brief Create an empty noise model frame. */
llsm_nmframe* llsm_create_nmframe(int nchannel, int nhar_e, int npsd);
/** @brief Copy-construct a harmonic model frame from an existing one. */
llsm_nmframe* llsm_copy_nmframe(llsm_nmframe* src);
/** @brief In-place version of llsm_copy_npframe. */
void llsm_copy_nmframe_inplace(llsm_nmframe* dst, llsm_nmframe* src);
/** @brief Delete and free a noise model frame. */
void llsm_delete_nmframe(llsm_nmframe* dst);
/** @} */

/** @defgroup group_frame LLSM Frame
 *  @brief A LLSM frame is essentially a container with a harmonic model and
 *    a noise model inside.
 *  @{ */
/** @brief Create an empty LLSM frame. */
llsm_container* llsm_create_frame(int nhar, int nchannel, int nhar_e,
  int npsd);
/** @brief Build the layer 0 harmonic model representation from an existing
 *    layer 1 representation. */
void llsm_frame_tolayer0(llsm_container* dst, llsm_container* conf);
/** @brief An extension of llsm_hmframe_phaseshift to LLSM frames. */
void llsm_frame_phaseshift(llsm_container* dst, FP_TYPE theta);
/** @brief Convert from absolute phase to relative phase shift (RPS). */
void llsm_frame_phasesync_rps(llsm_container* dst, int layer1_based);
/** @brief Compute the Signal-to-Noise Ratio from the layer 0 representation;
      return SNR (dB) or Aperiodicity (linear) on a warped frequency axis. */
FP_TYPE* llsm_frame_compute_snr(llsm_container* src, llsm_container* conf,
  int as_aperiodicity);
/** @brief Verify if the frame contains the information necessary for layer 0
 *    representation. */
int llsm_frame_checklayer0(llsm_container* src);
/** @brief Verify if the frame contains the information necessary for layer 1
 *    representation. */
int llsm_frame_checklayer1(llsm_container* src);
/** @} */

/** @brief Verify if the configuration contains the information necessary for
 *    layer 0 representation. */
int llsm_conf_checklayer0(llsm_container* src);
/** @brief Verify if the configuration contains the information necessary for
 *    layer 1 representation. */
int llsm_conf_checklayer1(llsm_container* src);

/** @brief Synthesis results. */
typedef struct {
  int ny;              /**< size of the output waveform */
  FP_TYPE fs;          /**< sampling rate (Hz) */
  FP_TYPE* y;          /**< output waveform */
  FP_TYPE* y_sin;      /**< sinusoidal component of the output waveform */
  FP_TYPE* y_noise;    /**< noise component of the output waveform */
} llsm_output;

/** @brief Delete and free the synthesis results. */
void llsm_delete_output(llsm_output* dst);

/** @defgroup group_aoptions llsm_aoptions
 *  @{ */
/** @brief Options for the analysis routine. */
typedef struct {
  FP_TYPE thop;        /**< hop time (seconds) */
  int maxnhar;         /**< maximum number of harmonics */
  int maxnhar_e;       /**< maximum number of harmonics for noise envelopes */
  int npsd;            /**< size of the PSD vector */
  int nchannel;        /**< number of channels for noise modeling */
  FP_TYPE* chanfreq;   /**< channel frequencies for noise modeling */
  FP_TYPE noise_warp;  /**< spectral warping factor for noise modeling */
  FP_TYPE lip_radius;  /**< default lip radius (cm) */

  int f0_refine;       /**< flag for enabling F0 refminement */
  int hm_method;       /**< method for harmonic analysis */
  FP_TYPE rel_winsize; /**< the ratio of window size to period length */
} llsm_aoptions;

/** @brief Create default analysis options. */
llsm_aoptions* llsm_create_aoptions();
/** @brief Delete and free analysis options. */
void llsm_delete_aoptions(llsm_aoptions* dst);
/** @brief Create a model configuration from analysis options. */
llsm_container* llsm_aoptions_toconf(llsm_aoptions* src, FP_TYPE fnyq);

#define LLSM_AOPTION_HMPP  0 /**< Peaking-Picking method for harmonic
                                analysis. */
#define LLSM_AOPTION_HMCZT 1 /**< Chirp-Z Transform for harmonic analysis. */

/** @} */

/** @defgroup group_soptions llsm_soptions
 *  @{ */
/** @brief Options for the synthesis routine. */
typedef struct {
  FP_TYPE fs;           /**< output sampling rate (Hz) */
  int use_iczt;         /**< automatically switch to ICZT based harmonic
                             signal generation if it's predicted to be faster
                             than the recurrent method */
  FP_TYPE iczt_param_a; /**< the slope parameter for switching on/off ICZT */
  FP_TYPE iczt_param_b; /**< the offset parameter for switching on/off ICZT */
} llsm_soptions;

/** @brief Create default synthesis options. */
llsm_soptions* llsm_create_soptions(FP_TYPE fs);
/** @brief Delete and free synthesis options. */
void llsm_delete_soptions(llsm_soptions* dst);
/** @} */

/** @defgroup group_chunk llsm_chunk
 *  @{ */
/** @brief A LLSM parameter chunk consisting of an array of LLSM frames. */
typedef struct {
  llsm_container* conf;
  llsm_container** frames;
} llsm_chunk;

/** @brief Create an empty parameter chunk from model configurations. */
llsm_chunk* llsm_create_chunk(llsm_container* conf, int init_frames);
/** @brief Copy-construct a parameter chunk from an existing one. */
llsm_chunk* llsm_copy_chunk(llsm_chunk* src);
/** @brief Delete and free a parameter chunk. */
void llsm_delete_chunk(llsm_chunk* dst);

/** @brief Build the layer 1 representation from an existing layer 0
 *    representation. */
void llsm_chunk_tolayer1(llsm_chunk* dst, int nfft);
/** @brief Build the layer 0 harmonic model representation from an existing
 *    layer 1 representation. */
void llsm_chunk_tolayer0(llsm_chunk* dst);
/** @brief An extension of llsm_frame_phasesync_rps to LLSM chunks. */
void llsm_chunk_phasesync_rps(llsm_chunk* dst, int layer1_based);
/** @brief Add or subtract the integration of F0 to/from the phase vectors. */
void llsm_chunk_phasepropagate(llsm_chunk* dst, int sign);
/** @brief Get F0 and number of frames from a parameter chunk. */
FP_TYPE* llsm_chunk_getf0(llsm_chunk* src, int* dst_nfrm);

/** @brief Perform layer 0 analysis on a speech signal. */
llsm_chunk* llsm_analyze(llsm_aoptions* options, FP_TYPE* x, int nx,
  FP_TYPE fs, FP_TYPE* f0, int nfrm, FP_TYPE** x_ap);
/** @brief Generate speech from a LLSM parameter chunk. */
llsm_output* llsm_synthesize(llsm_soptions* options, llsm_chunk* src);
/** @} */

/** @defgroup group_coder LLSM Coder
 *  @{ */
/** @brief Temporary data for (lossy) encoding and decoding of LLSM frames.
      The implementation is not visible to users. */
typedef void llsm_coder;

/** @brief Create a coder for a certain model configuration. The coder can be
      used for both encoding and decoding. The dimensionality of encoded
      frames is order_spec + order_bap + 3. */
llsm_coder* llsm_create_coder(llsm_container* conf, int order_spec,
  int order_bap);
/** @brief Delete and free an LLSM coder. */
void llsm_delete_coder(llsm_coder* dst);
/** @brief Convert an LLSM frame into a fixed-dimensional vector. */
FP_TYPE* llsm_coder_encode(llsm_coder* c, llsm_container* src);
/** @brief Reconstruct an LLSM frame (layer 1 representation) from a
      fixed-dimensional vector. */
llsm_container* llsm_coder_decode_layer1(llsm_coder* c, FP_TYPE* src);
/** @brief Reconstruct an LLSM frame (layer 0 representation) from a
      fixed-dimensional vector, without going through layer 1. */
llsm_container* llsm_coder_decode_layer0(llsm_coder* c, FP_TYPE* src);

/** @} */

#endif

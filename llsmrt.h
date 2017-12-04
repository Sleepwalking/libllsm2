/** @file */

#ifndef LLSM_LLSMRT_H
#define LLSM_LLSMRT_H

#include "llsm.h"

/** @defgroup group_llsmrt LLSM Real-time Synthesis
 *  @{ */
/** @brief A type hiding the implementation details of a LLSM real-time
 *    synthesis buffer. */
typedef void llsm_rtsynth_buffer;

/** @brief Create a real-time synthesis buffer based on LLSM configuraiton. */
llsm_rtsynth_buffer* llsm_create_rtsynth_buffer(llsm_soptions* options,
  llsm_container* conf, int capacity_frames, int capacity_samples);
/** @brief Delete and free a real-time synthesis buffer and all the frames
 *    stored in it. */
void llsm_delete_rtsynth_buffer(llsm_rtsynth_buffer* dst);
/** @brief Append a LLSM frame to the real-time synthesis buffer. */
void llsm_rtsynth_buffer_feed(llsm_rtsynth_buffer* dst, llsm_container* frame);
/** @brief Get one sample from the real-time synthesis buffer; returns 1 on
 *    success. */
int llsm_rtsynth_buffer_fetch(llsm_rtsynth_buffer* src, FP_TYPE* dst);

/** #} */

#endif

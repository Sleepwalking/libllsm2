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
  llsm_container* conf, int capacity_samples);
/** @brief Delete and free a real-time synthesis buffer and all the frames
 *    stored in it. */
void llsm_delete_rtsynth_buffer(llsm_rtsynth_buffer* dst);

/** @brief Get the real-time synthesis latency (in samples). */
int llsm_rtsynth_buffer_getlatency(llsm_rtsynth_buffer* src);
/** @brief Get the number of output samples available. */
int llsm_rtsynth_buffer_numoutput(llsm_rtsynth_buffer* src);
/** @brief Append a LLSM frame to the real-time synthesis buffer. */
void llsm_rtsynth_buffer_feed(llsm_rtsynth_buffer* dst, llsm_container* frame);
/** @brief Get one sample from the real-time synthesis buffer; returns 1 on
 *    success. */
int llsm_rtsynth_buffer_fetch(llsm_rtsynth_buffer* src, FP_TYPE* dst);

/** #} */

#endif

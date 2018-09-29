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

/** @file Some LLSM-specific utilities shared across the library. Users
 *    should not include this file. */

#ifndef LLSM_UTILS_H
#define LLSM_UTILS_H

#include "llsm.h"
#include "external/ciglet/ciglet.h"

/** @brief Automatically pick the most efficient method (from sinusoid bank
 *    and CZT) for harmonic synthesis . */
FP_TYPE* llsm_synthesize_harmonic_frame_auto(llsm_soptions* options,
  FP_TYPE* ampl, FP_TYPE* phse, int nhar, FP_TYPE f0, int nx);

/** @brief Generate a pulse period from a LF model and filter it by layer 1
 *    parameters, while keeping the phases coherent with the harmonic model. */
FP_TYPE* llsm_make_filtered_pulse(llsm_container* src, lfmodel source,
  FP_TYPE phase_correction, int pre_rotate, int size, FP_TYPE fnyq,
  FP_TYPE lip_radius, FP_TYPE fs);

#endif

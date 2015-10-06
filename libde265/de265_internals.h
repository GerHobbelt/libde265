/*
* H.265 video codec.
* Copyright (c) 2013-2014 struktur AG, Dirk Farin <farin@struktur.de>
*
* This file is part of libde265.
*
* libde265 is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as
* published by the Free Software Foundation, either version 3 of
* the License, or (at your option) any later version.
*
* libde265 is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU Lesser General Public License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with libde265.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DE265_DEBUG_H
#define DE265_DEBUG_H

#ifdef __cplusplus
extern "C" {
#endif

#include <libde265/de265-version.h>
#include "de265.h"

#if defined(_MSC_VER) && !defined(LIBDE265_STATIC_BUILD)
#ifdef LIBDE265_EXPORTS
#define LIBDE265_API __declspec(dllexport)
#else
#define LIBDE265_API __declspec(dllimport)
#endif
#elif HAVE_VISIBILITY
#ifdef LIBDE265_EXPORTS
#define LIBDE265_API __attribute__((__visibility__("default")))
#else
#define LIBDE265_API
#endif
#else
#define LIBDE265_API
#endif

#if __GNUC__
#define LIBDE265_DEPRECATED __attribute__((deprecated))
#elif defined(_MSC_VER)
#define LIBDE265_DEPRECATED __declspec(deprecated)
#else
#define LIBDE265_DEPRECATED
#endif

#if defined(_MSC_VER)
#define LIBDE265_INLINE __inline
#else
#define LIBDE265_INLINE inline
#endif

/// Get the number of CTBs in the image
LIBDE265_API int de265_internals_get_numCTB(const struct de265_image*);

/* Get the slice index of each CTB in the image. Take cake that the array idxArray is at least
 * as long as there are CTBs (see de265_internals_get_numCTB()).
 */
LIBDE265_API void de265_internals_get_CTB_sliceIdx(const struct de265_image*, uint16_t *idxArray);

#ifdef __cplusplus
}
#endif

#endif

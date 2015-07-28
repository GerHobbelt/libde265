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

#ifndef FALLBACK_UPSAMPLING_H
#define FALLBACK_UPSAMPLING_H

#include <stddef.h>
#include <stdint.h>

// Scalability extensions - Inter layer upsampling process for a whole image
void resampling_process_of_luma_sample_values_fallback (uint8_t *src, ptrdiff_t srcstride, int src_size[2],
                                                        uint8_t *dst, ptrdiff_t dststride, int dst_size[2],
                                                        int position_params[10]);
void resampling_process_of_chroma_sample_values_fallback (uint8_t *src, ptrdiff_t srcstride, int src_size[2],
                                                          uint8_t *dst, ptrdiff_t dststride, int dst_size[2],
                                                          int position_params[10]);

// Inter layer upsampling for a block. The src pointer must already be padded if necessary.
// (So src-4 and src-4*src_stride must be accessable)
void resampling_process_of_luma_block_fallback_8bit (const uint8_t *src,  ptrdiff_t src_stride, int16_t src_height,
                                                     int16_t *dst, ptrdiff_t dst_stride, int dst_width, int dst_height,
                                                     int x, int y, const int *position_params);
void resampling_process_of_luma_block_fallback_16bit (const uint16_t *src,  ptrdiff_t src_stride, int16_t src_height,
                                                     int16_t *dst, ptrdiff_t dst_stride, int dst_width, int dst_height,
                                                     int x, int y, const int *position_params);
void resampling_process_of_chroma_block_fallback_8bit (const uint8_t *src,  ptrdiff_t src_stride, int16_t src_height,
                                                       int16_t *dst, ptrdiff_t dst_stride, int dst_width, int dst_height,
                                                       int x, int y, const int *position_params);
void resampling_process_of_chroma_block_fallback_16bit (const uint16_t *src,  ptrdiff_t src_stride, int16_t src_height,
                                                       int16_t *dst, ptrdiff_t dst_stride, int dst_width, int dst_height,
                                                       int x, int y, const int *position_params);


#endif
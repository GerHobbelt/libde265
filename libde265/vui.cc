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

#include "vui.h"
#include "vps.h"

const char* get_video_format_name(enum VideoFormat format)
{
  switch (format) {
  case VideoFormat_Component: return "component";
  case VideoFormat_PAL:       return "PAL";
  case VideoFormat_NTSC:      return "NTSC";
  case VideoFormat_SECAM:     return "SECAM";
  case VideoFormat_MAC:       return "MAC";
  default:                    return "unspecified";
  }
}


const char* get_SAR_Indicator_name(enum SAR_Inidcator format)
{
  switch (format)
  {
  case UNSPECIFIED: return "UNSPECIFIED ";
  case SAR_1_1: return "SAR_1_1";
  case SAR_12_11: return "SAR_12_11";
  case SAR_10_11: return "SAR_10_11";
  case SAR_16_11: return "SAR_16_11";
  case SAR_40_33: return "SAR_40_33";
  case SAR_24_11: return "SAR_24_11";
  case SAR_20_11: return "SAR_20_11";
  case SAR_32_11: return "SAR_32_11";
  case SAR_80_33: return "SAR_80_33";
  case SAR_18_11: return "SAR_18_11";
  case SAR_15_11: return "SAR_15_11";
  case SAR_64_33: return "SAR_64_33";
  case SAR_160_99: return "SAR_160_99";
  case SAR_4_3: return "SAR_4_3";
  case SAR_3_2: return "SAR_3_2";
  case SAR_2_1: return "SAR_2_1";
  case SAR_RESERVED: return "SAR_RESERVED";
  case SAR_EXTENDED: return "SAR_EXTENDED";
  default: return "";
  }
}

video_usability_information::video_usability_information()
{
  aspect_ratio_info_present_flag = false;
  sar_width = 0;
  sar_height = 0;


  // --- overscan ---

  overscan_info_present_flag = false;
  overscan_appropriate_flag = false;


  // --- video signal type ---

  video_signal_type_present_flag = false;
  video_format = VideoFormat_Unspecified;
  video_full_range_flag = false;
  colour_description_present_flag = false;
  colour_primaries = 2;
  transfer_characteristics = 2;
  matrix_coeffs = 2;

  // --- chroma / interlaced ---

  chroma_loc_info_present_flag = false;
  chroma_sample_loc_type_top_field = 0;
  chroma_sample_loc_type_bottom_field = 0;

  neutral_chroma_indication_flag = false;
  field_seq_flag = false;
  frame_field_info_present_flag = false;

  // --- default display window ---

  default_display_window_flag = false;
  def_disp_win_left_offset = 0;
  def_disp_win_right_offset = 0;
  def_disp_win_top_offset = 0;
  def_disp_win_bottom_offset = 0;


  // --- timing ---

  vui_timing_info_present_flag = false;
  vui_num_units_in_tick = 0;
  vui_time_scale = 0;

  vui_poc_proportional_to_timing_flag = false;
  vui_num_ticks_poc_diff_one_minus1 = 0;


  // --- hrd parameters ---

  vui_hrd_parameters_present_flag = false;
  //hrd_parameters vui_hrd_parameters;


  // --- bitstream restriction ---

  bitstream_restriction_flag = false;
  tiles_fixed_structure_flag = false;
  motion_vectors_over_pic_boundaries_flag = true;
  restricted_ref_pic_lists_flag = false;
  min_spatial_segmentation_idc = 0;
  max_bytes_per_pic_denom = 2;
  max_bits_per_min_cu_denom = 1;
  log2_max_mv_length_horizontal = 15;
  log2_max_mv_length_vertical = 15;
}

de265_error video_usability_information::read(bitreader* br,
                                              int sps_max_sub_layers_minus1)
{
  aspect_ratio_info_present_flag = get_bits(br,1);
  if (aspect_ratio_info_present_flag) {
    int code = get_bits(br,8);
    if (code <= 16) aspect_ratio_idc = (SAR_Inidcator) code;
    else if (code < 255 ) aspect_ratio_idc = SAR_RESERVED;
    else aspect_ratio_idc = SAR_EXTENDED;
    if (aspect_ratio_idc == SAR_EXTENDED) {
      sar_width  = get_bits(br,16);
      sar_height = get_bits(br,16);
    }
  }

  overscan_info_present_flag = get_bits(br,1);
  if (overscan_info_present_flag) {
    overscan_appropriate_flag = get_bits(br,1);
  }

  video_signal_type_present_flag = get_bits(br,1);
  if (video_signal_type_present_flag) {
    int video_format_idc = get_bits(br, 3);
    if (video_format_idc > 5) {
      video_format_idc = VideoFormat_Unspecified;
    }
    video_format = (VideoFormat)video_format_idc;
    video_full_range_flag = get_bits(br,1);
    colour_description_present_flag = get_bits(br,1);
    if (colour_description_present_flag) {
      colour_primaries = get_bits(br,8);
      transfer_characteristics = get_bits(br,8);
      matrix_coeffs = get_bits(br,8);
    }
  }

  chroma_loc_info_present_flag = get_bits(br,1);
  if( chroma_loc_info_present_flag ) {
    chroma_sample_loc_type_top_field = get_uvlc(br);
    chroma_sample_loc_type_bottom_field = get_uvlc(br);
  }

  neutral_chroma_indication_flag = get_bits(br,1);
  field_seq_flag = get_bits(br,1);
  frame_field_info_present_flag = get_bits(br,1);
  default_display_window_flag = get_bits(br,1);
  if( default_display_window_flag ) {
    def_disp_win_left_offset = get_uvlc(br);
    def_disp_win_right_offset = get_uvlc(br);
    def_disp_win_top_offset = get_uvlc(br);
    def_disp_win_bottom_offset = get_uvlc(br);
  }

  vui_timing_info_present_flag = get_bits(br,1);
  if( vui_timing_info_present_flag ) {
    vui_num_units_in_tick = get_bits(br,32);
    vui_time_scale = get_bits(br,32);
    vui_poc_proportional_to_timing_flag = get_bits(br,1);
    if (vui_poc_proportional_to_timing_flag) {
      vui_num_ticks_poc_diff_one_minus1 = get_uvlc(br);
    }
    vui_hrd_parameters_present_flag = get_bits(br,1);
    if (vui_hrd_parameters_present_flag) {
      de265_error err = hrd.read(br, true, sps_max_sub_layers_minus1 );
      if (err != DE265_OK) return err;
    }
  }

  bitstream_restriction_flag = get_bits(br,1);
  if( bitstream_restriction_flag ) {
    tiles_fixed_structure_flag = get_bits(br,1);
    motion_vectors_over_pic_boundaries_flag = get_bits(br,1);
    restricted_ref_pic_lists_flag = get_bits(br,1);
    min_spatial_segmentation_idc = get_uvlc(br);
    max_bytes_per_pic_denom = get_uvlc(br);
    max_bits_per_min_cu_denom = get_uvlc(br);
    log2_max_mv_length_horizontal = get_uvlc(br);
    log2_max_mv_length_vertical = get_uvlc(br);
  }

  return DE265_OK;
}

de265_error vps_vui_video_signal_info::read_video_signal_info(bitreader* br)
{
  video_vps_format = get_bits(br,3);
  video_full_range_vps_flag = get_bits(br,1);
  colour_primaries_vps = get_bits(br,8);
  transfer_characteristics_vps = get_bits(br,8);
  matrix_coeffs_vps = get_bits(br,8);

  return DE265_OK;
}

de265_error vps_vui_bsp_hrd_params::read_vps_vui_bsp_hrd_params(bitreader* br,
                                                                video_parameter_set* vps)
{
  video_parameter_set_extension *vps_ext = &vps->vps_extension;

  vps_num_add_hrd_params = get_uvlc(br);
  for( int i = vps->vps_num_hrd_parameters;
        i < vps->vps_num_hrd_parameters + vps_num_add_hrd_params; i++ ) {
    if (i > 0) {
      cprms_add_present_flag[ i ] = get_bits(br,1);
    }
    num_sub_layer_hrd_minus1[i] = get_uvlc(br);

    hrd_params[i].read(br, cprms_add_present_flag[i], num_sub_layer_hrd_minus1[i]);
  }

  if( vps->vps_num_hrd_parameters + vps_num_add_hrd_params > 0 ) {
    for( int h = 1; h < vps_ext->NumOutputLayerSets; h++ ) {
      num_signalled_partitioning_schemes[ h ] = get_uvlc(br);
      for( int j = 1; j < num_signalled_partitioning_schemes[ h ] + 1; j++ ) {
        num_partitions_in_scheme_minus1[ h ][ j ] = get_uvlc(br);
        for( int k = 0; k <= num_partitions_in_scheme_minus1[ h ][ j ]; k++ ) {
          for( int r = 0; r < vps_ext->NumLayersInIdList[ vps_ext->OlsIdxToLsIdx[ h ] ]; r++ ) {
            layer_included_in_partition_flag[ h ][ j ][ k ][ r ] = get_bits(br,1);
          }
        }
      }

      for( int i = 0; i < num_signalled_partitioning_schemes[ h ] + 1; i++ ) {
        for( int t = 0; t <= vps_ext->MaxSubLayersInLayerSetMinus1[ vps_ext->OlsIdxToLsIdx[ h ] ]; t++ ) {
          num_bsp_schedules_minus1[ h ][ i ][ t ] = get_uvlc(br);
          for (int j = 0; j <= num_bsp_schedules_minus1[h][i][t]; j++) {
            for( int k = 0; k <= num_partitions_in_scheme_minus1[ h ][ i ]; k++ ) {
              if( vps->vps_num_hrd_parameters + vps_num_add_hrd_params > 1 ) {
                int nr_bits = ceil_log2( vps->vps_num_hrd_parameters + vps_num_add_hrd_params );
                bsp_hrd_idx[ h ][ i ][ t ][ j ][ k ] = get_bits(br,nr_bits);
              }
              bsp_sched_idx[ h ][ i ][ t ][ j ][ k ] = get_uvlc(br);
            }
          }
        }
      }
    }
  }

  return DE265_OK;
}

de265_error vps_vui::read_vps_vui( bitreader* br,
                                   video_parameter_set* vps)
{
  video_parameter_set_extension * vps_ext = &vps->vps_extension;

  cross_layer_pic_type_aligned_flag = get_bits(br,1);
  cross_layer_irap_aligned_flag = vps_ext->vps_vui_present_flag;
  if (!cross_layer_pic_type_aligned_flag) {
    cross_layer_irap_aligned_flag = get_bits(br,1);
  }
  if (cross_layer_irap_aligned_flag) {
all_layers_idr_aligned_flag = get_bits(br, 1);
  }
  bit_rate_present_vps_flag = get_bits(br, 1);
  pic_rate_present_vps_flag = get_bits(br, 1);
  if (bit_rate_present_vps_flag || pic_rate_present_vps_flag) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i < vps_ext->NumLayerSets; i++) {
      for (int j = 0; j <= vps_ext->MaxSubLayersInLayerSetMinus1[i]; j++) {
        if (bit_rate_present_vps_flag) {
          bit_rate_present_flag[i][j] = get_bits(br, 1);
        }
        if (pic_rate_present_vps_flag) {
          pic_rate_present_flag[i][j] = get_bits(br, 1);
        }
        if (bit_rate_present_flag[i][j]) {
          avg_bit_rate[i][j] = get_bits(br, 16);
          max_bit_rate[i][j] = get_bits(br, 16);
        }
        if (pic_rate_present_flag[i][j]) {
          constant_pic_rate_idc[i][j] = get_bits(br, 2);
          avg_pic_rate[i][j] = get_bits(br, 16);
        }
      }
    }
  }

  video_signal_info_idx_present_flag = get_bits(br, 1);
  if (video_signal_info_idx_present_flag) {
    vps_num_video_signal_info_minus1 = get_bits(br, 4);
  }
  for (int i = 0; i <= vps_num_video_signal_info_minus1; i++) {
    video_signal_info[i].read_video_signal_info(br);
  }

  if (video_signal_info_idx_present_flag && vps_num_video_signal_info_minus1 > 0) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= vps->MaxLayersMinus1; i++) {
      vps_video_signal_info_idx[i] = get_bits(br, 4);
    }
  }

  tiles_not_in_use_flag = get_bits(br, 1);

  if (!tiles_not_in_use_flag) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= vps->MaxLayersMinus1; i++) {
      tiles_in_use_flag[i] = get_bits(br, 1);
      if (tiles_in_use_flag[i]) {
        loop_filter_not_across_tiles_flag[i] = get_bits(br, 1);
      }
    }
    for (int i = vps->vps_base_layer_internal_flag ? 1 : 2; i <= vps->MaxLayersMinus1; i++) {
      for (int j = 0; j < vps_ext->NumDirectRefLayers[vps_ext->layer_id_in_nuh[i]]; j++) {
        int layerIdx = vps_ext->LayerIdxInVps[vps_ext->IdDirectRefLayer[vps_ext->layer_id_in_nuh[i]][j]];
        if (tiles_in_use_flag[i] && tiles_in_use_flag[layerIdx]) {
          tile_boundaries_aligned_flag[i][j] = get_bits(br, 1);
        }
      }
    }
  }

  wpp_not_in_use_flag = get_bits(br, 1);
  if (!wpp_not_in_use_flag) {
    for (int i = vps->vps_base_layer_internal_flag ? 0 : 1; i <= vps->MaxLayersMinus1; i++) {
      wpp_in_use_flag[i] = get_bits(br, 1);
    }
  }
  single_layer_for_non_irap_flag = get_bits(br, 1);
  higher_layer_irap_skip_flag = get_bits(br, 1);
  ilp_restricted_ref_layers_flag = get_bits(br, 1);

  if (ilp_restricted_ref_layers_flag) {
    for (int i = 1; i <= vps->MaxLayersMinus1; i++) {
      for (int j = 0; j < vps_ext->NumDirectRefLayers[vps_ext->layer_id_in_nuh[i]]; j++) {
        if (vps->vps_base_layer_internal_flag ||
          vps_ext->IdDirectRefLayer[vps_ext->layer_id_in_nuh[i]][j] > 0) {
          min_spatial_segment_offset_plus1[i][j] = get_uvlc(br);
          if (min_spatial_segment_offset_plus1[i][j] > 0) {
            ctu_based_offset_enabled_flag[i][j] = get_bits(br, 1);
            if (ctu_based_offset_enabled_flag[i][j]) {
              min_horizontal_ctu_offset_plus1[i][j] = get_uvlc(br);
            }
          }
        }
      }
    }
  }

  vps_vui_bsp_hrd_present_flag = get_bits(br, 1);
  if (vps_vui_bsp_hrd_present_flag) {
    bsp_hrd_params.read_vps_vui_bsp_hrd_params(br, vps);
  }

  for (int i = 1; i <= vps->MaxLayersMinus1; i++) {
    if (vps_ext->NumDirectRefLayers[vps_ext->layer_id_in_nuh[i]] == 0) {
      base_layer_parameter_set_compatibility_flag[i] = get_bits(br, 1);
    }
  }

  return DE265_OK;
}

de265_error hrd_parameters::read( bitreader* br,
                                  bool commonInfPresentFlag,
                                  int maxNumSubLayersMinus1)
{
  if (commonInfPresentFlag) {
    nal_hrd_parameters_present_flag = get_bits(br,1);
    vcl_hrd_parameters_present_flag = get_bits(br,1);
    if (nal_hrd_parameters_present_flag || vcl_hrd_parameters_present_flag) {
      sub_pic_hrd_params_present_flag = get_bits(br,1);
      if (sub_pic_hrd_params_present_flag) {
          tick_divisor_minus2 = get_bits(br,8);
          du_cpb_removal_delay_increment_length_minus1 = get_bits(br,5);
          sub_pic_cpb_params_in_pic_timing_sei_flag = get_bits(br,1);
          dpb_output_delay_du_length_minus1 = get_bits(br,5);
      }
      bit_rate_scale = get_bits(br,4);
      cpb_size_scale = get_bits(br,4);
      if (sub_pic_hrd_params_present_flag) {
        cpb_size_du_scale = get_bits(br,4);
      }
      initial_cpb_removal_delay_length_minus1 = get_bits(br,5);
      au_cpb_removal_delay_length_minus1 = get_bits(br,5);
      dpb_output_delay_length_minus1 = get_bits(br,5);
    }
  }

  for( int i = 0; i <= maxNumSubLayersMinus1; i++ ) {
    fixed_pic_rate_general_flag[ i ] = get_bits(br,1);
    if (!fixed_pic_rate_general_flag[i]) {
      fixed_pic_rate_within_cvs_flag[ i ] = get_bits(br,1);
    }
    if (fixed_pic_rate_within_cvs_flag[i]) {
      elemental_duration_in_tc_minus1[i] = get_uvlc(br);
    }
    else {
      low_delay_hrd_flag[ i ] = get_bits(br,1);
    }
    if (!low_delay_hrd_flag[i]) {
      cpb_cnt_minus1[i] = get_uvlc(br);
    }
    if (nal_hrd_parameters_present_flag) {
      sub_layer_hrd[i].read(br, this, i);
    }
    if (vcl_hrd_parameters_present_flag) {
      sub_layer_hrd[i].read(br, this, i);
    }
  }

  return DE265_OK;
}

de265_error sub_layer_hrd_parameters::read( bitreader* br,
                                            hrd_parameters *hrd,
                                            int subLayerId)
{
  int CpbCnt = hrd->cpb_cnt_minus1[subLayerId];
  for( int i = 0; i <= CpbCnt; i++ ) {
    bit_rate_value_minus1[ i ] = get_uvlc(br);
    cpb_size_value_minus1[ i ] = get_uvlc(br);
    if( hrd->sub_pic_hrd_params_present_flag ) {
      cpb_size_du_value_minus1[ i ] = get_uvlc(br);
      bit_rate_du_value_minus1[ i ] = get_uvlc(br);
    }
    cbr_flag[ i ] = get_bits(br,1);
  }

  return DE265_OK;
}

void video_usability_information::dump(int fd) const
{
  //#if (_MSC_VER >= 1500)
  //#define LOG0(t) loginfo(LogHeaders, t)
  //#define LOG1(t,d) loginfo(LogHeaders, t,d)
  //#define LOG2(t,d1,d2) loginfo(LogHeaders, t,d1,d2)
  //#define LOG3(t,d1,d2,d3) loginfo(LogHeaders, t,d1,d2,d3)

  FILE* fh;
  if (fd == 1) fh = stdout;
  else if (fd == 2) fh = stderr;
  else { return; }

#define LOG0(t) log2fh(fh, t)
#define LOG1(t,d) log2fh(fh, t,d)
#define LOG2(t,d1,d2) log2fh(fh, t,d1,d2)
#define LOG3(t,d1,d2,d3) log2fh(fh, t,d1,d2,d3)

  LOG0("----------------- VUI -----------------\n");
  LOG1("aspect_ratio_info_present_flag : %d\n", aspect_ratio_info_present_flag);
  if (aspect_ratio_info_present_flag) {
    LOG1("aspect_ratio_idc           : %d\n", get_SAR_Indicator_name(aspect_ratio_idc));
    LOG2("sample aspect ratio        : %d:%d\n", sar_width, sar_height);
  }
  
  LOG1("overscan_info_present_flag : %d\n", overscan_info_present_flag);
  LOG1("overscan_appropriate_flag  : %d\n", overscan_appropriate_flag);

  LOG1("video_signal_type_present_flag: %d\n", video_signal_type_present_flag);
  if (video_signal_type_present_flag) {
    LOG1("  video_format                : %s\n", get_video_format_name(video_format));
    LOG1("  video_full_range_flag       : %d\n", video_full_range_flag);
    LOG1("  colour_description_present_flag : %d\n", colour_description_present_flag);
    LOG1("  colour_primaries            : %d\n", colour_primaries);
    LOG1("  transfer_characteristics    : %d\n", transfer_characteristics);
    LOG1("  matrix_coeffs               : %d\n", matrix_coeffs);
  }

  LOG1("chroma_loc_info_present_flag: %d\n", chroma_loc_info_present_flag);
  if (chroma_loc_info_present_flag) {
    LOG1("  chroma_sample_loc_type_top_field   : %d\n", chroma_sample_loc_type_top_field);
    LOG1("  chroma_sample_loc_type_bottom_field: %d\n", chroma_sample_loc_type_bottom_field);
  }

  LOG1("neutral_chroma_indication_flag: %d\n", neutral_chroma_indication_flag);
  LOG1("field_seq_flag                : %d\n", field_seq_flag);
  LOG1("frame_field_info_present_flag : %d\n", frame_field_info_present_flag);

  LOG1("default_display_window_flag   : %d\n", default_display_window_flag);
  LOG1("  def_disp_win_left_offset    : %d\n", def_disp_win_left_offset);
  LOG1("  def_disp_win_right_offset   : %d\n", def_disp_win_right_offset);
  LOG1("  def_disp_win_top_offset     : %d\n", def_disp_win_top_offset);
  LOG1("  def_disp_win_bottom_offset  : %d\n", def_disp_win_bottom_offset);

  LOG1("vui_timing_info_present_flag  : %d\n", vui_timing_info_present_flag);
  if (vui_timing_info_present_flag) {
    LOG1("  vui_num_units_in_tick       : %d\n", vui_num_units_in_tick);
    LOG1("  vui_time_scale              : %d\n", vui_time_scale);
  }

  LOG1("vui_poc_proportional_to_timing_flag : %d\n", vui_poc_proportional_to_timing_flag);
  LOG1("vui_num_ticks_poc_diff_one          : %d\n", vui_num_ticks_poc_diff_one_minus1 + 1);

  LOG1("vui_hrd_parameters_present_flag : %d\n", vui_hrd_parameters_present_flag);
  if (vui_hrd_parameters_present_flag) {
    //hrd_parameters vui_hrd_parameters;
  }


  // --- bitstream restriction ---

  LOG1("bitstream_restriction_flag         : %d\n", bitstream_restriction_flag);
  if (bitstream_restriction_flag) {
    LOG1("  tiles_fixed_structure_flag       : %d\n", tiles_fixed_structure_flag);
    LOG1("  motion_vectors_over_pic_boundaries_flag : %d\n", motion_vectors_over_pic_boundaries_flag);
    LOG1("  restricted_ref_pic_lists_flag    : %d\n", restricted_ref_pic_lists_flag);
    LOG1("  min_spatial_segmentation_idc     : %d\n", min_spatial_segmentation_idc);
    LOG1("  max_bytes_per_pic_denom          : %d\n", max_bytes_per_pic_denom);
    LOG1("  max_bits_per_min_cu_denom        : %d\n", max_bits_per_min_cu_denom);
    LOG1("  log2_max_mv_length_horizontal    : %d\n", log2_max_mv_length_horizontal);
    LOG1("  log2_max_mv_length_vertical      : %d\n", log2_max_mv_length_vertical);
  }

#undef LOG0
#undef LOG1
#undef LOG2
#undef LOG3
  //#endif
}
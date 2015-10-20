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

#include "de265_internals.h"
#include "image.h"

LIBDE265_API void de265_internals_get_CTB_Info_Layout(const struct de265_image *img, int *widthInUnits, int *heightInUnits, int *log2UnitSize)
{
	img->internals_get_CTB_Info_Layout(widthInUnits, heightInUnits, log2UnitSize);
}

LIBDE265_API void de265_internals_get_CTB_sliceIdx(const struct de265_image *img, uint16_t *idxArray)
{
	img->internals_get_sliceIdx(idxArray);
}

LIBDE265_API void de265_internals_get_CB_Info_Layout(const struct de265_image *img, int *widthInUnits, int *heightInUnits, int *log2UnitSize)
{
	img->internals_get_CB_Info_Layout(widthInUnits, heightInUnits, log2UnitSize);
}

LIBDE265_API void de265_internals_get_CB_info(const struct de265_image *img, uint16_t *idxArray)
{
	img->internals_get_CB_info(idxArray);
}

LIBDE265_API void de265_internals_get_PB_Info_layout(const struct de265_image *img, int *widthInUnits, int *heightInUnits, int *log2UnitSize)
{
	img->internals_get_PB_Info_layout(widthInUnits, heightInUnits, log2UnitSize);
}

LIBDE265_API void de265_internals_get_PB_info(const struct de265_image *img, int16_t *refPOC0, int16_t *refPOC1, int16_t *x0, int16_t *y0, int16_t *x1, int16_t *y1)
{
	img->internals_get_PB_info(refPOC0, refPOC1, x0, y0, x1, y1);
}

LIBDE265_API void de265_internals_get_IntraDir_Info_layout(const struct de265_image *img, int *widthInUnits, int *heightInUnits, int *log2UnitSize)
{
	img->internals_get_IntraDir_Info_Layer(widthInUnits, heightInUnits, log2UnitSize);
}

LIBDE265_API void de265_internals_get_intraDir_info(const struct de265_image *img, uint8_t *intraDir, uint8_t *intraDirChroma)
{
	img->internals_get_intraDir_info(intraDir, intraDirChroma);
}

LIBDE265_API void de265_internals_get_TUInfo_Info_layout(const struct de265_image *img, int *widthInUnits, int *heightInUnits, int *log2UnitSize)
{
	img->internals_get_TUInfo_Info_layout(widthInUnits, heightInUnits, log2UnitSize);
}

LIBDE265_API void de265_internals_get_TUInfo_info(const struct de265_image *img, uint8_t *tuInfo)
{
	img->internals_get_TUInfo_info(tuInfo);
}
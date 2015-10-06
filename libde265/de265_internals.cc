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

LIBDE265_API int de265_internals_get_numCTB(const struct de265_image* img)
{
	return img->number_of_ctbs();
}

LIBDE265_API void de265_internals_get_CTB_sliceIdx(const struct de265_image *img, uint16_t *idxArray)
{
	img->internals_get_sliceIdx(idxArray);
}
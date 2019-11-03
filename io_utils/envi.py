#!/usr/bin/python3

#   Copyright 2016-2019 Adamo Ferro
#
#   This file is part of SOPA.
#
#   SOPA is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   SOPA is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with SOPA. If not, see <http://www.gnu.org/licenses/>.
#
#   The use of SOPA or part of it for the creation of any sub-product
#   (e.g., scientific papers, posters, images, other softwares)
#   must be acknowledged.

from osgeo import gdal, gdalconst
import numpy as np

NP2GDAL_CONVERSION = {      # source: https://borealperspectives.org/2014/01/16/data-type-mapping-when-using-pythongdal-to-write-numpy-arrays-to-geotiff/
  "uint8": 1,
  "int8": 1,
  "uint16": 2,
  "int16": 3,
  "uint32": 4,
  "int32": 5,
  "float32": 6,
  "float64": 7,
  "complex64": 10,
  "complex128": 11,
}

def read(filename):
    # gdal.Driver.Register()

    img = gdal.Open(filename, gdalconst.GA_ReadOnly)

    if img is None:
        # error occurred
        print("ERROR")
        return None
    else:
        cols = img.RasterXSize
        rows = img.RasterYSize

        band = img.GetRasterBand(1)
        np_array = band.ReadAsArray(0, 0, cols, rows)

        return np_array


def write(filename, np_array, complex_output=False):
    if (np_array.dtype.name == "complex64" or np_array.dtype.name == "complex128") and not complex_output:
        np_array_tmp = np.abs(np_array).astype("float64")**2.
    else:
        np_array_tmp = np_array
    try:
        driver = gdal.GetDriverByName('ENVI')
        ds = driver.Create(filename, np_array_tmp.shape[1], np_array_tmp.shape[0], 1, NP2GDAL_CONVERSION[np_array_tmp.dtype.name])
        ds.GetRasterBand(1).WriteArray(np_array_tmp)
        ds.FlushCache()             # write to disk.
    except:
        print("ERROR: problem writing image to disk.")
        return False
    return True

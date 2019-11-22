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

# TODO add limit checks

from dem.dem_base import DEMBase
import numpy as np
from io_utils import envi


class MOLA128(DEMBase):
    """MOLA MEGDR 1/128 degrees DEM"""

    def __init__(self, dem_name="MOLA128", filename_base="", dummy_value=0, max_lat=88, min_long=0):
        super().__init__(dem_name, filename_base, dummy_value)
        self._MOLA_MAX_LAT = max_lat        # degrees
        self._MOLA_MIN_LONG = min_long      # degrees
        self._MOLA_SAMPLING = 1./128.       # degrees per pixel

        self._mask = None

        self.deg_sampling = self._MOLA_SAMPLING
        self.offset_to_be_added = 3396000       # as per MOLA labels

    def read(self):
        if self.filename_base != "":
            try:
                self.dem = envi.read(self.filename_base)
                self.dem_shape = self.dem.shape
                return True
            except IOError:
                print("ERROR: cannot read MOLA128 DEM file.")
        else:
            print("ERROR: header filename not set.")
        return False

    def reset(self):
        del self.dem

    def get_dem_radius_from_lat_lon(self, lat, lon):
        if self.dem is not None:
            mola_x, mola_y = self.to_XY_from_LL(lat, lon)
            return self._XY_vectors_to_radius(mola_x, mola_y)
        else:
            print("ERROR: MOLA DEM not yet loaded.")
            return None

    def get_dem_radius_and_latlon_grid_from_latlon_area(self, lat_n, lat_s, lon_e, lon_w):
        if self.dem is not None:
            lat_vector = np.arange(lat_s, lat_n+self.deg_sampling, self.deg_sampling)
            lon_vector = np.arange(lon_w, lon_e+self.deg_sampling, self.deg_sampling)
            lon_vector[lon_vector < 0] = lon_vector[lon_vector < 0] + 360
            latlon_mesh_lat, latlon_mesh_lon = np.meshgrid(lat_vector, lon_vector)
            mola_x, mola_y = self.to_XY_from_LL(latlon_mesh_lat, latlon_mesh_lon)
            return latlon_mesh_lat, latlon_mesh_lon, self._XY_vectors_to_radius(mola_x, mola_y)
        else:
            print("ERROR: MOLA DEM not yet loaded.")
            return None

    def to_XY_from_LL(self, lat, long, return_as_float=False):
        mola_x = (long - self._MOLA_MIN_LONG) / self.deg_sampling
        mola_y = (self._MOLA_MAX_LAT - lat) / self.deg_sampling

        if not return_as_float:
            mola_x = np.round(mola_x).astype("int")
            mola_y = np.round(mola_y).astype("int")

        if len(mola_x.shape) == 0:
            mola_x = np.array([mola_x])
            mola_y = np.array([mola_y])

        mola_x[mola_x == self.dem_shape[1]] = 0

        self._mask = (mola_x < 0) | (mola_y < 0) | (mola_x >= self.dem_shape[1]) | (mola_y >= self.dem_shape[0])
        mola_x[mola_x < 0] = 0
        mola_y[mola_y < 0] = 0
        mola_x[mola_x >= self.dem_shape[1]] = 0
        mola_y[mola_y >= self.dem_shape[0]] = 0

        return (mola_x, mola_y)

    def to_LL_from_XY(self, x, y):
        lat = -y * self.deg_sampling + self._MOLA_MAX_LAT
        lon = x * self.deg_sampling + self._MOLA_MIN_LONG
        return (lat, lon)

    def XY_ranges_to_radius(self, x_min, x_max, y_min, y_max):
        output_dem = self.dem[y_min:y_max, x_min:x_max].copy() + self.offset_to_be_added
        return output_dem

    def _XY_vectors_to_radius(self, x_vect, y_vect):
        # get vector containing MOLA radius corresponding to each frame to be processed (in the same order)
        output_dem = self.dem[y_vect, x_vect].copy() + self.offset_to_be_added
        output_dem[self._mask] = self.dummy_value
        return output_dem

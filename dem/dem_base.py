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

class DEMBase():
    """docstring for DEMBase."""

    def __init__(self, dem_name="no name", filename_base="", dummy_value=0):
        self.dem_name = dem_name
        self.filename_base = filename_base
        self.dummy_value = dummy_value
        self.dem = None
        self.offset_to_be_added = 0

    def __str__(self):
        return self.dem_name + ", base filename: " + self.filename_base

    def read(self):
        raise NotImplementedError

    def get_dem_radius_from_lat_lon(self, lat, lon):
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

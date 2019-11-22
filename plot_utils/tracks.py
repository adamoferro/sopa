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

from matplotlib import pyplot as pp
import numpy as np


class TrackPlotOnDEM(object):
    """docstring"""

    def __init__(self, dem_obj, track_lats, track_lons, title=None):
        self._dem_obj = dem_obj
        self._track_lats = track_lats
        self._track_lons = track_lons
        self._n_labels = 10
        if title is not None:
            self._title = title
        else:
            self._title = "Spacecraft ground track on DEM [km]"

    def _calculate_xy_limits(self, lats, lons):
        sc_track_x, sc_track_y = self._dem_obj.to_XY_from_LL(lats, lons, return_as_float=True)
        min_x = int(np.floor(np.min(sc_track_x)))
        max_x = int(np.ceil(np.max(sc_track_x)))
        min_y = int(np.floor(np.min(sc_track_y)))
        max_y = int(np.ceil(np.max(sc_track_y)))
        return min_x, max_x, min_y, max_y, sc_track_x, sc_track_y

    def _create_labels(self, max_lat, min_lat, max_long, min_long):
        lat_range = np.round(np.arange(max_lat, min_lat, -self._dem_obj.deg_sampling)*100)/100
        lat_range_len = len(lat_range)
        long_range = np.round(np.arange(min_long, max_long, self._dem_obj.deg_sampling)*100)/100
        long_range_len = len(long_range)
        step_x = int(long_range_len / (self._n_labels - 1))
        step_y = int(lat_range_len / (self._n_labels - 1))
        labels_x_pos = np.arange(0, long_range_len, step_x)
        labels_y_pos = np.arange(0, lat_range_len, step_y)
        labels_x = long_range[::step_x]
        labels_y = lat_range[::step_y]
        return labels_x, labels_y, labels_x_pos, labels_y_pos

    def plot(self):
        min_x, max_x, min_y, max_y, sc_track_x, sc_track_y = self._calculate_xy_limits(self._track_lats, self._track_lons)

        min_lat, min_long = self._dem_obj.to_LL_from_XY(min_x, max_y)
        max_lat, max_long = self._dem_obj.to_LL_from_XY(max_x, min_y)

        labels_x, labels_y, labels_x_pos, labels_y_pos = self._create_labels(max_lat, min_lat, max_long, min_long)

        pp.figure(figsize=(12, 6))
        pp.imshow(self._dem_obj.XY_ranges_to_radius(min_x, max_x+1, min_y, max_y+1)/1000, interpolation="none")
        pp.plot(sc_track_x-min_x, sc_track_y-min_y, "r")
        pp.axis("tight")
        pp.title(self._title, fontweight="bold")
        pp.xlabel("Long", fontweight="bold")
        pp.xticks(labels_x_pos, labels_x)
        pp.ylabel("Lat", fontweight="bold")
        pp.yticks(labels_y_pos, labels_y)
        pp.colorbar()
        pp.show()


class TrackAndFirstReturnsPlotOnDEM(TrackPlotOnDEM):
    """docstring"""

    def __init__(self, dem_obj, track_lats, track_lons, fr_lats, fr_lons, title=None):
        super().__init__(dem_obj, track_lats, track_lons, title=None)
        self._fr_lats = fr_lats
        self._fr_lons = fr_lons
        if title is not None:
            self._title = title
        else:
            self._title = "Spacecraft ground track (yellow) vs.\nsimulated first return positions (magenta) on DEM [km]"

    def plot(self):
        min_x, max_x, min_y, max_y, fr_x, fr_y = self._calculate_xy_limits(self._fr_lats, self._fr_lons)
        _, _, _, _, sc_track_x, sc_track_y = self._calculate_xy_limits(self._track_lats, self._track_lons)

        sc_track_x_ok_ids = (sc_track_x >= min_x) & (sc_track_x <= max_x)
        sc_track_y_ok_ids = (sc_track_y >= min_y) & (sc_track_y <= max_y)
        sc_track_ok_ids = sc_track_x_ok_ids & sc_track_y_ok_ids
        sc_track_x = sc_track_x[sc_track_ok_ids]
        sc_track_y = sc_track_y[sc_track_ok_ids]

        min_lat, min_long = self._dem_obj.to_LL_from_XY(min_x, max_y)
        max_lat, max_long = self._dem_obj.to_LL_from_XY(max_x, min_y)

        labels_x, labels_y, labels_x_pos, labels_y_pos = self._create_labels(max_lat, min_lat, max_long, min_long)

        pp.figure(figsize=(12, 6))
        pp.imshow(self._dem_obj.XY_ranges_to_radius(min_x, max_x+1, min_y, max_y+1)/1000, interpolation="none")
        pp.plot(sc_track_x-min_x, sc_track_y-min_y, "y")
        pp.plot(fr_x-min_x, fr_y-min_y, "+", color="magenta", markersize=3)
        pp.axis("tight")
        pp.title(self._title, fontweight="bold")
        pp.xlabel("Long", fontweight="bold")
        pp.xticks(labels_x_pos, labels_x)
        pp.ylabel("Lat", fontweight="bold")
        pp.yticks(labels_y_pos, labels_y)
        pp.colorbar()
        pp.show()

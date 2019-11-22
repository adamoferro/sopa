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

import matplotlib
from matplotlib import pyplot as pp
import numpy as np


class _BaseRadargramPlot(object):
    def __init__(self, image, db_output=False, min_value=-1, max_value=-1, method=None, top_point=0, bottom_point=-1, title=""):
        self._db_output = db_output
        if self._db_output:
            self._image = 10.*np.log10(image)
        else:
            self._image = image

        if method == "hist" and (min_value == -1 or max_value == -1):
            hist_min, hist_max = self._hist_limits()

        if min_value == -1:
            if method != "hist":
                self._min_value = np.min(self._image)
            else:
                self._min_value = hist_min
        else:
            self._min_value = min_value
        if max_value == -1:
            if method != "hist":
                self._max_value = np.max(self._image)
            else:
                self._max_value = hist_max
        else:
            self._max_value = max_value
        self._top_point = top_point
        if bottom_point == -1:
            self._bottom_point = self._image.shape[0]
        else:
            self._bottom_point = bottom_point
        self._title = title

    def _hist_limits(self, ll=0.5, ul=0.995):
        lower_perc_limit = ll
        upper_perc_limit = ul
        h, be = np.histogram(self._image, bins=1000, density=True)
        bw = be[1] - be[0]
        h *= bw
        hs = np.cumsum(h)
        h_ok_ids = np.where((hs >= lower_perc_limit) & (hs <= upper_perc_limit))[0]
        return be[h_ok_ids[0]], be[h_ok_ids[-1]]

    def plot(self, cmap="gray"):
        pp.figure(figsize=(12, 6))
        pp.imshow(self._image[self._top_point:self._bottom_point, :], cmap=cmap, vmin=self._min_value, vmax=self._max_value, interpolation="none")
        pp.title(self._title, fontweight="bold")
        pp.ylabel("Samples", fontweight="bold")
        pp.xlabel("Frames", fontweight="bold")
        pp.colorbar()
        pp.show()


class RadargramPlot(_BaseRadargramPlot):
    """docstring"""

    def __init__(self, image, db_output=True, min_value=-1, max_value=-1, method="hist", top_point=0, bottom_point=-1, title=""):
        super().__init__(image, db_output, min_value, max_value, method, top_point, bottom_point, title)


class SimulationPlot(_BaseRadargramPlot):
    """docstring"""

    def __init__(self, image, min_value=-1, max_value=-1, method="hist", top_point=0, bottom_point=-1, title=""):
        super().__init__(image, False, min_value, max_value, method, top_point, bottom_point, title)


class GroundDistancePlot(_BaseRadargramPlot):
    """docstring"""

    def __init__(self, image, min_value=0, max_value=20, top_point=0, bottom_point=-1, title=""):
        super().__init__(image, False, min_value, max_value, None, top_point, bottom_point, title)

    def plot(self):
        nan_cmap = matplotlib.cm.get_cmap("gist_rainbow")
        nan_cmap.set_bad(color='black')
        super().plot(cmap=nan_cmap)

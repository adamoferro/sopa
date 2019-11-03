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

# TODO:
# - modify print with logger

from log.log import fake_logger
import gzip


class RSData():
    """docstring for RSData."""

    def __init__(self, dataset_name="no name", pri=0, presumming=0, filename_base="", use_gzip_and_iso8859_1=False, logger=None):
        self.dataset_name = dataset_name
        self.n_frames = 0
        self.pri = pri                      # in s
        self.presumming = presumming
        self.filename_base = filename_base
        self.use_gzip_and_iso8859_1 = use_gzip_and_iso8859_1

        self._orbit_data = None
        self._data = None
        self._data_type = ""

        # select the "open" function according to the selected encoding
        self.encoding = "UTF-8"
        self.open_function = open
        if use_gzip_and_iso8859_1:
            self.encoding = "ISO-8859-1"
            self.open_function = gzip.open

        if logger is None:
            self.logger = fake_logger()
        else:
            self.logger = logger

    def _correct_lon(self):
        if self._orbit_data is not None:
            self._orbit_data["Lon"][self._orbit_data["Lon"] < 0] = 360 + self._orbit_data["Lon"][self._orbit_data["Lon"] < 0]

    @property
    def data(self):
        if self._data is None:
            print("ERROR: data not yet read.")
            return None
        else:
            return self._data

    @data.setter
    def data(self, value):
        print("ERROR: cannot set data manually.")

    @property
    def orbit_data(self):
        if self._orbit_data is None:
            print("ERROR: orbit data not yet generated.")
            return None
        else:
            return self._orbit_data

    @orbit_data.setter
    def orbit_data(self, value):
        print("ERROR: cannot set orbit data manually.")

    def __str__(self):
        output_str = self.dataset_name + ", n_frames=" + str(self.n_frames) + ", PRI=" + "{:.1f}".format(self.pri*1e6) + " us, presum=" + str(self.presumming)
        if self.filename_base != "":
            output_str += ", base filename: " + self.filename_base
        return output_str

    def load(self, mode="full", force_reload=False):            # mode can be: "all", "ancillary", "data". Not all implementation may be able to provide different behaviors
        if self.filename_base != "":
            return self._load_from_file(mode, force_reload)
        else:
            print("ERROR: dataset base filename not set.")
            return False

    def set_orbit_data_from_ext(self, ext_orbit_data):
        self._orbit_data = ext_orbit_data

    def _load_from_file(self, mode, force_reload):      # "virtual" method
        raise NotImplementedError

    def generate_orbit_data(self, skip):                # generates orbit data which match 1:1 with the data
        raise NotImplementedError

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

from input.data_base import RSData
from . import pds_us_log as log_messages
from geometry.coordinate_converter import CoordinateConverter
import os
import numpy as np

# TODO use more the logger


class RDR(RSData):
    def __init__(self, dataset_name="no name", pri=0, presumming=0, filename_base="", use_gzip_and_iso8859_1=False, logger=None):
        super().__init__(dataset_name, pri, presumming, filename_base, use_gzip_and_iso8859_1, logger)

        # input file suffixes
        self._DATA_SUFFIX = "_rgram.img"
        self._ORBIT_SUFFIX = "_geom.tab"

        self._DT = 3 / 80e6

        self._ancillary_already_read = False
        self._data_already_read = False

        self._geom_data = None     # numpy array, dtype={
                                    #                "names": ("column_id", "timestamp", "Lat", "Lon", "MarsRadius", "Radius", "Vradial", "Vtang", "SZA", "PhaseDist"),
                                    #                "formats": ("i", "datetime64[ms]", "d", "d", "d", "d", "d", "d", "d", "d")
                                    #              }

        self._data_type = "focused"

    def _load_from_file(self, mode="full", force_reload=False):      # filename_base is set, checked by super
        if self._check_file_presence(mode):
            load_ok1 = True
            load_ok2 = True
            if (mode == "full" or mode == "ancillary") and ((not self._ancillary_already_read) or force_reload):
                load_ok1 = self._load_ancillary_data()
            if (mode == "full" or mode == "data") and ((not self._data_already_read) or force_reload):
                load_ok2 = self._load_rdr_data()
            return load_ok1 and load_ok2
        return False

    def _check_file_presence(self, mode="full"):
        files_ok1 = True
        files_ok2 = True
        if mode == "full" or mode == "ancillary":
            self._ORBIT_FILENAME = self.filename_base + self._ORBIT_SUFFIX
            files_ok1 = os.path.isfile(self._ORBIT_FILENAME)
        if (mode == "full" and files_ok1) or mode == "data":
            self._DATA_FILENAME = self.filename_base + self._DATA_SUFFIX
            files_ok2 = os.path.isfile(self._DATA_FILENAME)
        return (files_ok1 and files_ok2)

    def _load_ancillary_data(self):
        # orbit data parsing
        geom_data = self._read_rdr_data_geom_file()
        if geom_data is None:
            self.logger.log_messages(log_messages.ERR_MESSAGES["orbit_file_generic"])
            return False
        else:
            self._geom_data = geom_data
        self._ancillary_already_read = True
        return True

    def _read_rdr_data_geom_file(self):
        def parsetime(v):
            return np.datetime64(v)

        try:
            with self.open_function(self._ORBIT_FILENAME, mode="rt", encoding=self.encoding) as fp:
                geom_data = np.loadtxt(
                    fp,
                    dtype={
                        "names": ("column_id", "timestamp", "Lat", "Lon", "MarsRadius", "Radius", "Vradial", "Vtang", "SZA", "PhaseDist"),
                        "formats": ("i", "datetime64[ms]", "d", "d", "d", "d", "d", "d", "d", "d")
                    },
                    delimiter=',',
                    skiprows=0,
                    converters={1: parsetime},
                )
                self.n_frames = geom_data.shape[0]
                return geom_data
        except IOError:
            self.logger.log_messages(log_messages.ERR_MESSAGES["orbit_file_read"])
        except UnicodeDecodeError:
            self.logger.log_messages(log_messages.ERR_MESSAGES["orbit_file_decode"])
        return None

    def _load_rdr_data(self, start_frame=0, end_frame=-1):
        '''
        Read PDS US RDR data (filenames in the form s_xxxxxxxx_rgram.img)
        - start_frame and end_frame can be in the range [0, OBS_NUMBER_OF_FRAMES-1]
            (Python array indexing convention)
        - end_frame == -1 means "all the frames starting from start_frame"
        - if end_frame is greater than the available number of frames, only the
            available data are returned without any error or warning
        - if an error occurs, the function returns False
        '''
        N_SAMPLES = 3600
        frame_count = -1
        if start_frame < 0 or (end_frame != -1 and start_frame > end_frame):
            self.logger.log_messages(log_messages.ERR_MESSAGES["data_frame-limits"])
            return None
        elif end_frame != -1:
            frame_count = end_frame-start_frame
        try:
            with open(self._DATA_FILENAME, 'rb') as fp:
                fp.seek(start_frame*N_SAMPLES)
                rdr_data = np.fromfile(fp, dtype=np.float32, count=frame_count*N_SAMPLES).reshape((3600, -1)).astype("float64")
                self._data = rdr_data
                return True
        except IOError:
            self.logger.log_messages(log_messages.ERR_MESSAGES["data_read"])

        return False

    def generate_orbit_data(self, skip=1):
        '''
        Saves in self._orbit_data the following fields as a dictionary:
        "time_offset_s", "rxwin_time_s", "Lat", "Lon", "Radius", "Vtang", "Vradial", "X", "Y", "Z", "SZA"
        There is no need for interpolation, as geometric data are already matched 1:1 with RDR data.
        The field rxwin_time_s is arbitrarily set to 0, as it is not given in the data.
        Values of X, Y and Z are calculated, as they are not present in the data either.
        '''
        print("WARNING: orbit data obtained from RDR PDS US data may be not sufficiently precise for simulation and/or focusing.")
        geom_fields_to_include = ["Lat", "Lon", "Radius", "Vtang", "Vradial", "SZA"]
        if self._ancillary_already_read:
            n_output_frames = int(self.n_frames/skip)
            orbit_data = dict()
            orbit_data["rxwin_time_s"] = np.zeros(n_output_frames)

            for field in geom_fields_to_include:
                orbit_data[field] = self._geom_data[field][::skip].copy()
            orbit_data["Radius"] *= 1000.
            orbit_data["X"], orbit_data["Y"], orbit_data["Z"] = CoordinateConverter.to_xyz_from_latlon_and_radius(orbit_data["Lat"], orbit_data["Lon"], orbit_data["Radius"])
            orbit_data["time_offset_s"] = (self._geom_data["timestamp"][::skip].astype("float64")-self._geom_data["timestamp"][0].astype("float64"))/1000.
            orbit_data["dt"] = self._DT
            self._orbit_data = orbit_data
            self._correct_lon()

        else:
            print("ERROR: ancillary data not yet loaded.")
            return None

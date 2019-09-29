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
#

from input.data_base import RSData
from . import pds_log as log_messages
from geometry.coordinate_converter import CoordinateConverter
import os
import numpy as np


class RDR(RSData):
    def __init__(self, dataset_name="no name", pri=0, presumming=0, filename_base="", use_gzip_and_iso8859_1=False, logger=None, read_hk_data=False):
        super().__init__(dataset_name, pri, presumming, filename_base, use_gzip_and_iso8859_1, logger)

        # input file suffixes
        self._HEADER_SUFFIX = ".lbl"
        self._DATA_SUFFIX = ".dat"

        # required fields and values in the ancillary data to consider the data valid
        # note: it is mandatory to set these listis and they must have the same length
        self._FIELDS_TO_SEARCH = ["START_TIME", "FILE_RECORDS", "INSTRUMENT_MODE_ID", "MRO:PULSE_REPETITION_INTERVAL", "DATA_QUALITY_ID"]
        self._FIELDS_TO_SEARCH_DTYPE = [str, int, str, str, str, str, str]
        self._FIELDS_TO_SEARCH_CONSTRAIN = ["", "", "", ["1428 <MICROSECOND>", "1428 <MICROSECONDS>"], "\"0\""]
        self._N_FIELDS_TO_SEARCH = len(self._FIELDS_TO_SEARCH)

        self._HEADER_FIELD_N_FRAMES = "FILE_RECORDS"
        self._HEADER_FIELD_PRI = "MRO:PULSE_REPETITION_INTERVAL"
        self._HEADER_FIELD_IM = "INSTRUMENT_MODE_ID"

        # acquisition modes and related parameters
        self._IMS = [3, 6, 9, 12, 15, 18, 21, 2, 5, 8, 11, 14, 17, 20, 1, 4, 7, 10, 13, 16, 19]
        self._IM_PRESUMMING = [16, 2, 28, 4, 32, 8, 1, 28, 4, 32, 8, 1, 16, 2, 32, 8, 1, 16, 2, 28, 4]

        self._N_STORED_SAMPLES = 667
        self._DT_ADC = 3 / 80e6
        self._DT = self._DT_ADC*2
        self._DCG_DELAY = 11.98e-6      # digital chirp generator delay, used in the calculation of the rxwin opening time

        self._header_already_read = False
        self._data_already_read = False

        # output
        self._header_data = None    # dictionary, general values (e.g., PRI, presumming)
        self._whole_data = None

    def _load_from_file(self, mode=None, force_reload=False):      # filename_base is set, checked by super
        # no "mode", only "data" because everything is contained in one file
        if self._check_file_presence():
            if ((not self._data_already_read) or force_reload):
                return self._load_rdr_data()
        return False

    def _check_file_presence(self):
        self._HEADER_FILENAME = self.filename_base + self._HEADER_SUFFIX
        if os.path.isfile(self._HEADER_FILENAME):
            return self._check_file_presence_data()
        else:
            self.logger.log_messages(log_messages.ERR_MESSAGES["other_file_not-found"])
        return False

    def _check_file_presence_data(self):
        self._DATA_FILENAME = self.filename_base + self._DATA_SUFFIX
        if os.path.isfile(self._DATA_FILENAME):
            return True
        else:
            self.logger.log_messages(log_messages.ERR_MESSAGES["other_file_not-found"])
        return False

    def _decode_im(self, im):
        '''
        service function that returns a tuple containing the
        main parameters of the acquisition mode
        '''
        id = self._IMS.index(im)
        return self._IM_PRESUMMING[id]

    def _read_rdr_data_header(self):
        try:
            with self.open_function(self._HEADER_FILENAME, mode="rt", encoding=self.encoding) as fp:
                header_data = dict()
                error_field = ""
                i_field_to_search = 0
                line = fp.readline()
                while line and i_field_to_search < self._N_FIELDS_TO_SEARCH:
                    if line[0] != "/" and len(line) > 1:
                        unpacked_line = line.split('=')
                        if len(unpacked_line) == 2:
                            field, value = unpacked_line[0:2]
                            value = value.strip()
                            if field.strip() == self._FIELDS_TO_SEARCH[i_field_to_search]:
                                value = self._FIELDS_TO_SEARCH_DTYPE[i_field_to_search](value)
                                header_data[self._FIELDS_TO_SEARCH[i_field_to_search]] = value
                                if self._FIELDS_TO_SEARCH_CONSTRAIN[i_field_to_search] != "":
                                    if isinstance(self._FIELDS_TO_SEARCH_CONSTRAIN[i_field_to_search], list):
                                        if value not in self._FIELDS_TO_SEARCH_CONSTRAIN[i_field_to_search]:
                                            error_field = self._FIELDS_TO_SEARCH[i_field_to_search]
                                            break
                                    elif value != self._FIELDS_TO_SEARCH_CONSTRAIN[i_field_to_search]:
                                        error_field = self._FIELDS_TO_SEARCH[i_field_to_search]
                                        break

                                i_field_to_search += 1
                    line = fp.readline()

                if len(header_data) != self._N_FIELDS_TO_SEARCH and error_field == "":
                    self.logger.log_messages(log_messages.ERR_MESSAGES["header_file_generic"])
                else:
                    return header_data, error_field
        except IOError:
            self.logger.log_messages(log_messages.ERR_MESSAGES["header_file_read"])
        except UnicodeDecodeError:
            self.logger.log_messages(log_messages.ERR_MESSAGES["header_file_decode"])
        return None

    @staticmethod
    def _correct_range_shift(data, range_shifts):
        ref_range_shift = range_shifts.max()
        frame_shifts = np.round((ref_range_shift - range_shifts) / 2).astype("int")      # range_shifts are expressed in ADC sampling time, that it one half of RDR radargram sample time
        zero_padding_length = frame_shifts.max()
        zero_padded_data = np.zeros([data.shape[0]+zero_padding_length, data.shape[1]], dtype=data.dtype)
        for i_f in np.arange(data.shape[1], dtype="int"):
            zero_padded_data[frame_shifts[i_f]:frame_shifts[i_f]+data.shape[0], i_f] = data[:, i_f]
        return zero_padded_data

    def _load_rdr_data(self, start_frame=0, end_frame=-1):
        if not self._header_already_read:
            # header parsing
            header_data_grouped = self._read_rdr_data_header()

            if header_data_grouped is not None:
                header_data, header_error_field = header_data_grouped

                if header_error_field != "":
                    err_msg = log_messages.ERR_MESSAGES["discarded"] + " " + header_error_field + " = " + str(header_data[header_error_field])
                    self.logger.log_messages(err_msg)
                    return False
                else:
                    self._header_data = header_data
                    self.n_frames = self._header_data[self._HEADER_FIELD_N_FRAMES]
                    self.pri = int(self._header_data[self._HEADER_FIELD_PRI].split(" ")[0])/1e6

                    im = self._header_data[self._HEADER_FIELD_IM]
                    if im[:2] == "SS":
                        im = int(im[2:])
                        self.presumming = self._decode_im(im)
                    else:
                        self.logger.log_messages(log_messages.ERR_MESSAGES["header_file_no_ss"])
                        return False
            else:
                self.logger.log_messages(log_messages.ERR_MESSAGES["header_file_generic"])
                return False
            self._header_already_read = True

        dt = np.dtype([("SCET_BLOCK_WHOLE", "<u4"),         #    0 (starts at)
                ("SCET_BLOCK_FRAC", "<u2"),                 #    4
                ("SKIP1", np.void, 180),                    #    6
                ("RECEIVE_WINDOW_OPENING_TIME", "<f4"),     #  186
                ("SKIP2", np.void, 4),                      #  190
                ("ECHO_SAMPLES_REAL", "<f4", 667),          #  194
                ("ECHO_SAMPLES_IMAGINARY", "<f4", 667),     # 2862
                ("SKIP3", np.void, 30),                     # 5530
                ("RANGE_SHIFT", "<i2"),                     # 5560
                ("SKIP4", np.void, 31),                     # 5562
                ("SOLAR_LONGITUDE", "<f8"),                 # 5593
                ("ORBIT_NUMBER", "<i4"),
                ("X_MARS_SC_POSITION_VECTOR", "<f8"),       # 5606
                ("Y_MARS_SC_POSITION_VECTOR", "<f8"),
                ("Z_MARS_SC_POSITION_VECTOR", "<f8"),
                ("SPACECRAFT_ALTITUDE", "<f8"),             # 5630
                ("SUB_SC_EAST_LONGITUDE", "<f8"),
                ("SUB_SC_PLANETOCENTRIC_LATITUDE", "<f8"),
                ("SUB_SC_PLANETOGRAPHIC_LATITUDE", "<f8"),
                ("X_MARS_SC_VELOCITY_VECTOR", "<f8"),
                ("Y_MARS_SC_VELOCITY_VECTOR", "<f8"),
                ("Z_MARS_SC_VELOCITY_VECTOR", "<f8"),
                ("MARS_SC_RADIAL_VELOCITY", "<f8"),
                ("MARS_SC_TANGENTIAL_VELOCITY", "<f8"),
                ("LOCAL_TRUE_SOLAR_TIME", "<f8"),
                ("SOLAR_ZENITH_ANGLE", "<f8"),
                ("SC_PITCH_ANGLE", "<f8"),
                ("SC_YAW_ANGLE", "<f8"),
                ("SC_ROLL_ANGLE", "<f8"),
                ("MRO_SAMX_INNER_GIMBAL_ANGLE", "<f8"),
                ("MRO_SAMX_OUTER_GIMBAL_ANGLE", "<f8"),
                ("MRO_SAPX_INNER_GIMBAL_ANGLE", "<f8"),
                ("MRO_SAPX_OUTER_GIMBAL_ANGLE", "<f8"),
                ("MRO_HGA_INNER_GIMBAL_ANGLE", "<f8"),
                ("MRO_HGA_OUTER_GIMBAL_ANGLE", "<f8"),
                ("DES_TEMP", "<f4"),
                ("DES_5V", "<f4"),
                ("DES_12V", "<f4"),
                ("DES_2V5", "<f4"),
                ("RX_TEMP", "<f4"),
                ("TX_TEMP", "<f4"),
                ("TX_LEV", "<f4"),
                ("TX_CURR", "<f4"),
                ("CORRUPTED_DATA_FLAG", "<i1")])

        frame_count = -1
        if start_frame < 0 or (end_frame != -1 and start_frame > end_frame):
            self.logger.log_messages(log_messages.ERR_MESSAGES["data_frame-limits"])
            return None
        elif end_frame != -1:
            frame_count = end_frame-start_frame
        try:
            with open(self._DATA_FILENAME, 'rb') as fp:
                fp.seek(start_frame * dt.itemsize)
                raw_file = fp.read()
                raw_data_np = np.fromstring(raw_file, dtype=dt, count=frame_count)
                raw_data = (raw_data_np["ECHO_SAMPLES_REAL"] + 1j*raw_data_np["ECHO_SAMPLES_IMAGINARY"]).T

                self._whole_data = raw_data_np
                self._data = self._correct_range_shift(raw_data, raw_data_np["RANGE_SHIFT"])
                self._data_already_read = True
                return True
        except IOError:
            self.logger.log_messages(log_messages.ERR_MESSAGES["data_read"])
        return False

    def generate_orbit_data(self, skip=1):
        '''
        Saves in self._orbit_data the following fields as a dictionary,
        interpolated to match 1:skip the frames of the raw data:
        "time_offset_s", "rxwin_time_s", "Lat", "Lon", "Radius", "Vtang", "Vradial", "X", "Y", "Z", "Roll", "Pitch", "Yaw", "SZA"
        '''
        geom_fields_to_include = "SUB_SC_PLANETOCENTRIC_LATITUDE", "SUB_SC_EAST_LONGITUDE", "MARS_SC_TANGENTIAL_VELOCITY", "MARS_SC_RADIAL_VELOCITY", "X_MARS_SC_POSITION_VECTOR", "Y_MARS_SC_POSITION_VECTOR", "Z_MARS_SC_POSITION_VECTOR", "SC_ROLL_ANGLE", "SC_PITCH_ANGLE", "SC_YAW_ANGLE", "SOLAR_ZENITH_ANGLE"
        geom_fields_to_include_new_names = ["Lat", "Lon", "Vtang", "Vradial", "X", "Y", "Z", "Roll", "Pitch", "Yaw", "SZA"]
        geom_fields_to_convert_to_m = ["Radius", "X", "Y", "Z"]
        if self._data_already_read:
            orbit_data = dict()

            orbit_data["time_offset_s"] = (self._whole_data["SCET_BLOCK_WHOLE"] + self._whole_data["SCET_BLOCK_FRAC"]*(2**-16))[::skip]
            orbit_data["time_offset_s"] -= orbit_data["time_offset_s"][0]

            for i_f, field in enumerate(geom_fields_to_include):
                orbit_data[geom_fields_to_include_new_names[i_f]] = self._whole_data[field][::skip].astype("float64")
            orbit_data["Radius"] = CoordinateConverter.to_radius_from_xyz(orbit_data["X"], orbit_data["Y"], orbit_data["Z"])

            for field in geom_fields_to_convert_to_m:
                orbit_data[field] *= 1000.

            orbit_data["rxwin_time_s"] = self._whole_data["RECEIVE_WINDOW_OPENING_TIME"][::skip]*self._DT_ADC - self._DCG_DELAY + (0 if self.pri > 1500e-6 else self.pri)

            self._orbit_data = orbit_data

        else:
            print("ERROR: data not yet loaded.")
            return None

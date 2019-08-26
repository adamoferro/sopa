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

ERR_MESSAGES = dict()
ERR_MESSAGES["argparse"] = "ERROR: problem while parsing input parameters."
ERR_MESSAGES["frame_skip"] = "ERROR: frame skip must be greater than zero."
ERR_MESSAGES["fft_length"] = "ERROR: minimum allowed range FFT length is 3600."
ERR_MESSAGES["median_window_perc_length"] = "ERROR: the median window percentual length must be between 0 and 100."
ERR_MESSAGES["emi_thr_db"] = "ERROR: the dB threshold must be greater than zero."
ERR_MESSAGES["avg_win_length"] = "ERROR: the average power window length must be greater than zero."
ERR_MESSAGES["fft_max_filter_length"] = "ERROR: the FFT max filter length must be greater than zero."
ERR_MESSAGES["summary_half_n_bins"] = "ERROR: the summary half number of FFT bins must be greater than zero and smaller than (fft_max_filter_length-1)/2."
ERR_MESSAGES["header_file_generic"] = "ERROR: problem while parsing header file."
ERR_MESSAGES["header_file_read"] = "ERROR: cannot read or cannot parse header file."
ERR_MESSAGES["header_file_decode"] = "ERROR: cannot decode header file. Try --use_gzip option."
ERR_MESSAGES["header_file_no-fields-required"] = "ERROR: at least one header field must be parsed, or simply avoid calling this function."
ERR_MESSAGES["header_file_options"] = "ERROR: options not valid."
ERR_MESSAGES["aux_file_read"] = "ERROR: cannot read aux file."
ERR_MESSAGES["aux_file_generic"] = "ERROR: problem while reading aux file."
ERR_MESSAGES["orbit_file_not-found"] = "ERROR: orbit file not present."
ERR_MESSAGES["orbit_file_read"] = "ERROR: cannot read orbit file."
ERR_MESSAGES["orbit_file_decode"] = "ERROR: cannot decode orbit file. Try --use_gzip option."
ERR_MESSAGES["orbit_file_generic"] = "ERROR: problem while reading orbit file."
ERR_MESSAGES["other_file_not-found"] = "ERROR: header and/or aux and/or data file not present."
ERR_MESSAGES["discarded"] = "--> discarded, header field not compliant:"
ERR_MESSAGES["no-consistency"] = "ERROR: no consistency between header and aux data (different number of frames is reported)."
ERR_MESSAGES["raw-data_frame-limits"] = "ERROR: start_frame and/or end_frame not valid."
ERR_MESSAGES["raw-data_read"] = "ERROR: could not read raw data."
ERR_MESSAGES["hk_file_not-found"] = "ERROR: HK file not found."
ERR_MESSAGES["hk-data_read"] = "ERROR: could not read HK data."
ERR_MESSAGES["hk_file_generic"] = "ERROR: problem while reading HK file."

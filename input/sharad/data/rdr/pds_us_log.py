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
ERR_MESSAGES["orbit_file_not-found"] = "ERROR: orbit file not present."
ERR_MESSAGES["orbit_file_read"] = "ERROR: cannot read orbit file."
ERR_MESSAGES["orbit_file_decode"] = "ERROR: cannot decode orbit file. Try --use_gzip option."
ERR_MESSAGES["orbit_file_generic"] = "ERROR: problem while reading orbit file."
ERR_MESSAGES["other_file_not-found"] = "ERROR: header and/or aux and/or data file not present."
ERR_MESSAGES["data_frame-limits"] = "ERROR: start_frame and/or end_frame not valid."
ERR_MESSAGES["data_read"] = "ERROR: could not read data."

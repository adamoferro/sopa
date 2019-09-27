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
ERR_MESSAGES["header_file_generic"] = "ERROR: problem while parsing header file."
ERR_MESSAGES["header_file_read"] = "ERROR: cannot read or cannot parse header file."
ERR_MESSAGES["header_file_decode"] = "ERROR: cannot decode header file. Try --use_gzip option."
ERR_MESSAGES["header_file_no_ss"] = "ERROR: the acquisition mode is not SS."
ERR_MESSAGES["other_file_not-found"] = "ERROR: label and/or aux and/or data file not present."
ERR_MESSAGES["discarded"] = "--> discarded, header field not compliant:"

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

class _logger_base():
    def __init__(self, verbose=False):
        raise NotImplementedError

    def log_messages(self, new_message):
        pass

    def log_write(self, filename):
        pass


class logger(_logger_base):
    def __init__(self, verbose=False):
        self.LOG_MESSAGES = list()
        self.verbose = verbose

    def log_messages(self, new_message):
        self.LOG_MESSAGES.append(new_message)
        if self.verbose:
            print(new_message)

    def log_write(self, filename):
        try:
            with open(filename, 'a') as fp:
                for lm in self.LOG_MESSAGES:
                    fp.write(lm + '\n')
                self.LOG_MESSAGES = list()           # empty for next write
        except IOError:
            print("ERROR: cannot write log file or append new messages to existing log file.")


class fake_logger(_logger_base):
    def __init__(self, verbose=False):
        pass

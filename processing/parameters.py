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
# - express half aperture in seconds or meters (when expressed in frames it depends on presumming)
# - check correctness of parameters
# - add Doppler max search and adjacent Doppler windows (as in SOFA "original")


class _BaseProcessingParameters():
    def __init__(self, frame_start=0, frame_end=-1, frame_skip=1, meters_skip=None, frame_half_aperture=768, squint_angle=0, range_weighting_function="hanning", range_weighting_kaiser_beta=2):
        # frame_half_aperture               samples
        # squint_angle                      degrees
        # range_weighting_function          ["hanning","hamming","kaiser","none"]
        # range_weighting_kaiser_beta       float > 0
        # meters_skip has priority on frame_skip

        self.frame_start = max([frame_start, frame_half_aperture])
        self.frame_end = frame_end
        if meters_skip is not None:
            self.meters_skip = meters_skip
            self.frame_skip = 0
        else:
            self.frame_skip = frame_skip
            self.meters_skip = 0

        self.frame_half_aperture = frame_half_aperture
        self.squint_angle = squint_angle
        self.range_weighting_function = range_weighting_function
        self.range_weighting_kaiser_beta = range_weighting_kaiser_beta


class SimulatorParameters(_BaseProcessingParameters):
    def __init__(self, frame_start=0, frame_end=-1, frame_skip=1, meters_skip=None, frame_half_aperture=768, squint_angle=0, range_weighting_function="hanning", range_weighting_kaiser_beta=2):
        super().__init__(frame_start, frame_end, frame_skip, meters_skip, frame_half_aperture, squint_angle, range_weighting_function, range_weighting_kaiser_beta)


class FocuserParameters(_BaseProcessingParameters):
    def __init__(self, frame_start=0, frame_end=-1, frame_skip=1, meters_skip=None, frame_half_aperture=768, squint_angle=0, multilooking_doppler_bandwidth=0.8, correct_range_migration=True, compensate_doppler_phase=True, range_weighting_function="hanning", range_weighting_kaiser_beta=2, azimuth_weighting_function="hanning", azimuth_weighting_kaiser_beta=2, gauss_rec_n_filters=0, gauss_rec_alpha=0.66):
        # multilooking_doppler_bandwidth    Hz
        # correct_range_migration           boolean
        # compensate_doppler_phase          boolean
        # azimuth_weighting_function        ["hanning","hamming","kaiser","none"]
        # azimuth_weighting_kaiser_beta     float > 0
        # gauss_rec_n_filters               int >= 0
        # gauss_rec_alpha                   float > 0

        super().__init__(frame_start, frame_end, frame_skip, meters_skip, frame_half_aperture, squint_angle, range_weighting_function, range_weighting_kaiser_beta)
        self.multilooking_doppler_bandwidth = multilooking_doppler_bandwidth
        self.correct_range_migration = correct_range_migration
        self.compensate_doppler_phase = compensate_doppler_phase
        self.azimuth_weighting_function = azimuth_weighting_function
        self.azimuth_weighting_kaiser_beta = azimuth_weighting_kaiser_beta
        self.gauss_rec_n_filters = gauss_rec_n_filters
        self.gauss_rec_alpha = gauss_rec_alpha

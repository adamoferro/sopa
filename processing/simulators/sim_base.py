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

import numpy as np


class _SIMBase(object):
    """docstring for SIMBase."""
    _C = 299792458.

    def __init__(self, geom_obj, dem_obj, param_obj, n_processes=1, use_data_dt=True):
        self.DEFAULT_DT = 37.5e-9
        self.USE_DATA_DT = use_data_dt

        self.dt = self.DEFAULT_DT
        self._N_OUTPUT_SAMPLES = 667*2
        self.sim_top_offset_samples = 200

        self.geom_obj = geom_obj
        self.dem_obj = dem_obj

        self._frames_to_be_simulated_start = param_obj.frame_start
        self._frames_to_be_simulated_end = param_obj.frame_end
        if param_obj.frame_skip == 0:
            self._skip_meters = param_obj.meters_skip
            self._skip_frames = 0
        else:
            self._skip_frames = param_obj.frame_skip
            self._skip_meters = None
        self._n_half_aperture_frames = param_obj.frame_half_aperture

        if n_processes > 1:
            try:
                import ray
            except ImportError:
                n_processes = 1
                print("WARNING: Multiprocessing not possible on this machine (failed to import module ""ray""). Working with 1 process.")
        elif n_processes < 1:
            n_processes = 1

        self.n_processes = n_processes

    def _pre_simulate(self, orbit_data):
        n_input_frames = orbit_data["X"].shape[0]

        if self._frames_to_be_simulated_end >= n_input_frames:
            print("ERROR: last frame to be simulated after the end of the input data.")
            return None

        if self.USE_DATA_DT:
            if "dt" in orbit_data:
                self.dt = orbit_data["dt"]
            else:
                print("WARNING: orbit data do not contain the dt value. Using the default dt = " + str(self.DEFAULT_DT) + " s.")

        # calculate actual start and end frames according to aperture
        if self._frames_to_be_simulated_start == -1:
            self._frames_to_be_simulated_start = np.max([self._n_half_aperture_frames, 0])

        if self._frames_to_be_simulated_end == -1:
            self._frames_to_be_simulated_end = np.min([n_input_frames-self._n_half_aperture_frames-1, n_input_frames-1])

        if self._frames_to_be_simulated_end < self._frames_to_be_simulated_start:
            print("ERROR: last frame to be simulated before first frame to be simulated.")
            return None

        # calculation of range of frames to be simulated
        if self._skip_frames != 0:
            frames_to_be_simulated = np.arange(self._frames_to_be_simulated_start, self._frames_to_be_simulated_end, self._skip_frames)
        else:
            frames_to_be_simulated, _, _ = self.geom_obj.get_frame_ids_from_skip_distance(orbit_data["Lat"][self._frames_to_be_simulated_start:self._frames_to_be_simulated_end], orbit_data["Lon"][self._frames_to_be_simulated_start:self._frames_to_be_simulated_end], self._skip_meters)
            frames_to_be_simulated += self._frames_to_be_simulated_start

        # correction of the last frame to be simulated
        self._frames_to_be_simulated_end = frames_to_be_simulated[-1]
        n_frames_to_be_simulated = len(frames_to_be_simulated)

        # select only data related to frames to be processed
        input_rxwin_data_t = orbit_data["rxwin_time_s"][frames_to_be_simulated]
        input_scx_data = orbit_data["X"][frames_to_be_simulated]
        input_scy_data = orbit_data["Y"][frames_to_be_simulated]
        input_scz_data = orbit_data["Z"][frames_to_be_simulated]
        input_lat_data = orbit_data["Lat"][frames_to_be_simulated]
        input_long_data = orbit_data["Lon"][frames_to_be_simulated]
        input_scr_data = orbit_data["Radius"][frames_to_be_simulated]

        # shifts for final frame alignment
        input_scr_data_alignment = input_scr_data - self.geom_obj.get_ellipsoid_radius_from_latlon(input_lat_data, input_long_data)
        input_scr_data_shifts_alignment = max(input_scr_data_alignment)-input_scr_data_alignment
        input_scr_data_shifts_alignment_t = (input_scr_data_shifts_alignment*2. / self._C)

        # # calculate reference rxwin and sc radius to be used as reference for frame time alignment
        # sc_radius_ref = np.min(input_scr_data)
        # ellipsoid_radius = self.geom_obj.get_ellipsoid_radius_from_latlon(input_lat_data, input_long_data)
        # dem_radius = self.dem_obj.get_dem_radius_from_lat_lon(input_lat_data, input_long_data)
        # dem_radius_ref_id = np.argmax(dem_radius)
        # ellipsoid_radius_ref = ellipsoid_radius[dem_radius_ref_id]
        # input_scr_data_shifts_alignment_t = ((sc_radius_ref - input_scr_data) + (ellipsoid_radius - ellipsoid_radius_ref)) * 2. / self._C

        if np.max(input_rxwin_data_t) == 0:     # if rxwin is not set in the orbit data
            print("NOTICE: field rxwin_time_s not set. Using the default value.")

            # rxwin is calculated as the time delay that sets the "highest point" of the simulation at sim_top_offset_samples from the top border
            dem_radius = self.dem_obj.get_dem_radius_from_lat_lon(input_lat_data, input_long_data)
            t_surface_vect_ref = 2 * (input_scr_data - dem_radius) / self._C
            rxwin_ref = np.min(t_surface_vect_ref) - self.sim_top_offset_samples * self.dt
            self._N_OUTPUT_SAMPLES += int(np.max(np.round((input_scr_data_shifts_alignment_t) / self.dt)))
        else:
            rxwin_ref = np.min(input_rxwin_data_t)
            self._N_OUTPUT_SAMPLES += int(np.max(np.round((input_scr_data_shifts_alignment_t + (input_rxwin_data_t-rxwin_ref)) / self.dt)))

        self.sim_image = np.zeros([self._N_OUTPUT_SAMPLES, n_frames_to_be_simulated])   #, dtype="int16")   # see comment below
        self.uncert_image = np.zeros(self.sim_image.shape)                          #, dtype="uint8")
        self.ground_distance_min_image = np.zeros(self.sim_image.shape)
        self.ground_distance_max_image = np.zeros(self.sim_image.shape)
        self.first_return_lats = np.zeros(n_frames_to_be_simulated)
        self.first_return_lons = np.zeros(n_frames_to_be_simulated)

        return (n_frames_to_be_simulated, input_rxwin_data_t, input_scx_data, input_scy_data, input_scz_data, input_lat_data, input_long_data, input_scr_data, rxwin_ref, input_scr_data_shifts_alignment_t)

    def _define_frame_ranges(self, n_frames_to_be_simulated):
        n_frames_per_subprocess = int(np.ceil(n_frames_to_be_simulated / self.n_processes))
        sub_frame_range_starts = np.arange(0, n_frames_to_be_simulated, n_frames_per_subprocess)
        sub_frame_range_stops = sub_frame_range_starts + n_frames_per_subprocess
        sub_frame_range_stops[-1] = np.min([sub_frame_range_stops[-1], n_frames_to_be_simulated])
        return sub_frame_range_starts, sub_frame_range_stops

    def simulate(self, orbit_data, sim_callback, further_pars=[]):

        pre_simulate_output = self._pre_simulate(orbit_data)
        if pre_simulate_output is not None:
            n_frames_to_be_simulated, input_rxwin_data_t, input_scx_data, input_scy_data, input_scz_data, input_lat_data, input_long_data, input_scr_data, rxwin_ref, input_scr_data_shifts_alignment_t = pre_simulate_output
        else:
            return None

        UTMNorthing_g, UTMEasting_g, UTMZone_g, UTMZoneNumber = self.geom_obj.LL_to_UTM(input_lat_data, input_long_data)

        # frame before, with same UTMZone of central frame
        UTMNorthing_before, UTMEasting_before = self.geom_obj.LL_to_UTM_with_given_zone_number(input_lat_data[:-self.n_frames_offset_for_dir_estimation], input_long_data[:-self.n_frames_offset_for_dir_estimation], UTMZoneNumber[self.n_frames_offset_for_dir_estimation:])
        UTMNorthing_before = np.concatenate((np.ones(self.n_frames_offset_for_dir_estimation) * UTMNorthing_before[0], UTMNorthing_before), 0)
        UTMEasting_before = np.concatenate((np.ones(self.n_frames_offset_for_dir_estimation) * UTMEasting_before[0], UTMEasting_before), 0)

        # frame after, with same UTMZone of central frame
        UTMNorthing_after, UTMEasting_after = self.geom_obj.LL_to_UTM_with_given_zone_number(input_lat_data[self.n_frames_offset_for_dir_estimation:], input_long_data[self.n_frames_offset_for_dir_estimation:], UTMZoneNumber[:-self.n_frames_offset_for_dir_estimation])
        UTMNorthing_after = np.concatenate((UTMNorthing_after, np.ones(self.n_frames_offset_for_dir_estimation) * UTMNorthing_after[-1]), 0)
        UTMEasting_after = np.concatenate((UTMEasting_after, np.ones(self.n_frames_offset_for_dir_estimation) * UTMEasting_after[-1]), 0)

        local_direction_g = np.arctan2(UTMNorthing_after - UTMNorthing_before, UTMEasting_after - UTMEasting_before)

        orbit_xyz = (input_scx_data, input_scy_data, input_scz_data)

        if self.n_processes > 1:
            sub_frame_range_starts, sub_frame_range_stops = self._define_frame_ranges(n_frames_to_be_simulated)

            n_bytes_to_reserve = int((self.dem_obj.dem.nbytes + self.alt_sim_matrix.nbytes + self.act_sim_matrix.nbytes)*1.30)  # calculate the amount of RAM to reserve for shared objects as 10% more of what summed here
            import ray
            ray.init(num_cpus=self.n_processes, object_store_memory=n_bytes_to_reserve)

            @ray.remote
            def _simulate_frame_range_callback(id, n_processes, orbit_xyz, frame_range, sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, input_scr_data_shifts_alignment_t, dt, dem_obj, geom_obj, alt_sim_matrix, act_sim_matrix, gaussian_filter_sigma, further_pars):
                return sim_callback(id, n_processes, orbit_xyz, frame_range, sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, input_scr_data_shifts_alignment_t, dt, dem_obj, geom_obj, alt_sim_matrix, act_sim_matrix, gaussian_filter_sigma, further_pars)

            dem_obj_id = ray.put(self.dem_obj)
            geom_obj_id = ray.put(self.geom_obj)
            alt_sim_matrix_id = ray.put(self.alt_sim_matrix)
            act_sim_matrix_id = ray.put(self.act_sim_matrix)

            result_ids = list()
            sub_frame_ranges = list()
            for id in np.arange(self.n_processes):
                sub_frame_range = np.arange(sub_frame_range_starts[id], sub_frame_range_stops[id])
                sub_frame_ranges.append(sub_frame_range)
                result_ids.append(_simulate_frame_range_callback.remote(id, self.n_processes, orbit_xyz, sub_frame_range, self._N_OUTPUT_SAMPLES, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, input_scr_data_shifts_alignment_t, self.dt, dem_obj_id, geom_obj_id, alt_sim_matrix_id, act_sim_matrix_id, self.gaussian_filter_sigma, further_pars))

            results = ray.get(result_ids)

            for id in np.arange(self.n_processes):
                self.sim_image[:, sub_frame_ranges[id]] = results[id][0]
                self.uncert_image[:, sub_frame_ranges[id]] = results[id][1]
                self.first_return_lats[sub_frame_ranges[id]] = results[id][2]
                self.first_return_lons[sub_frame_ranges[id]] = results[id][3]
                self.ground_distance_min_image[:, sub_frame_ranges[id]] = results[id][4]
                self.ground_distance_max_image[:, sub_frame_ranges[id]] = results[id][5]

            ray.shutdown()
        else:
            self.sim_image, self.uncert_image, self.first_return_lats, self.first_return_lons, self.ground_distance_min_image, self.ground_distance_max_image = sim_callback(0, 1, orbit_xyz, np.arange(n_frames_to_be_simulated), self._N_OUTPUT_SAMPLES, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, input_scr_data_shifts_alignment_t, self.dt, self.dem_obj, self.geom_obj, self.alt_sim_matrix, self.act_sim_matrix, self.gaussian_filter_sigma, further_pars)

        # self.sim_image = self.sim_image.astype("int16")
        self.uncert_image = self.uncert_image.astype("uint8")
        return self.sim_image, self.uncert_image, self.first_return_lats, self.first_return_lons, self.ground_distance_min_image, self.ground_distance_max_image

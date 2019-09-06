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


class DRSSim(object):
    """docstring for DRSSim."""
    _C = 299792458.

    def __init__(self, geom_obj, dem_obj, n_processes=1):
        self.n_frames_offset_for_dir_estimation = 1
        self.alt_sim_spacing = 50
        self.act_sim_spacing = 50
        self.alt_sim_forward_max_distance = 500
        self.alt_sim_backward_max_distance = 500
        self.act_sim_max_distance = 30000
        self.dt = 37.5e-9
        self.sim_n_samples = 1800
        self.sim_top_offset_samples = 200

        self.geom_obj = geom_obj
        self.dem_obj = dem_obj

        alt_sim_range = np.arange(-self.alt_sim_backward_max_distance, self.alt_sim_forward_max_distance + self.alt_sim_spacing, self.alt_sim_spacing)
        act_sim_range = np.arange(-self.act_sim_max_distance, self.act_sim_max_distance + self.act_sim_spacing, self.act_sim_spacing)
        self.alt_sim_matrix = np.array([alt_sim_range]).T * np.ones([1, len(act_sim_range)])
        self.act_sim_matrix = np.array([act_sim_range]) * np.ones([len(alt_sim_range), 1])

        multiprocessing_enabled = False
        if n_processes > 1:
            try:
                import ray
                multiprocessing_enabled = True
            except:
                n_processes = 1
                print("WARNING: Multiprocessing not possible on this machine. Working with 1 process.")
        elif n_processes < 1:
            n_processes = 1

        self.n_processes = n_processes

    def simulate(self, orbit_data):
        n_frames_to_be_simulated = len(orbit_data["Lat"])

        UTMNorthing_g, UTMEasting_g, UTMZone_g, UTMZoneNumber = self.geom_obj.LL_to_UTM(orbit_data["Lat"], orbit_data["Lon"])

        # frame before, with same UTMZone of central frame
        UTMNorthing_before, UTMEasting_before = self.geom_obj.LL_to_UTM_with_given_zone_number(orbit_data["Lat"][:-self.n_frames_offset_for_dir_estimation], orbit_data["Lon"][:-self.n_frames_offset_for_dir_estimation], UTMZoneNumber[self.n_frames_offset_for_dir_estimation:])
        UTMNorthing_before = np.concatenate((np.ones(self.n_frames_offset_for_dir_estimation) * UTMNorthing_before[0], UTMNorthing_before), 0)
        UTMEasting_before = np.concatenate((np.ones(self.n_frames_offset_for_dir_estimation) * UTMEasting_before[0], UTMEasting_before), 0)

        # frame after, with same UTMZone of central frame
        UTMNorthing_after, UTMEasting_after = self.geom_obj.LL_to_UTM_with_given_zone_number(orbit_data["Lat"][self.n_frames_offset_for_dir_estimation:], orbit_data["Lon"][self.n_frames_offset_for_dir_estimation:], UTMZoneNumber[:-self.n_frames_offset_for_dir_estimation])
        UTMNorthing_after = np.concatenate((UTMNorthing_after, np.ones(self.n_frames_offset_for_dir_estimation) * UTMNorthing_after[-1]), 0)
        UTMEasting_after = np.concatenate((UTMEasting_after, np.ones(self.n_frames_offset_for_dir_estimation) * UTMEasting_after[-1]), 0)

        local_direction_g = np.arctan2(UTMNorthing_after - UTMNorthing_before, UTMEasting_after - UTMEasting_before)

        # calculate reference rxwin and sc radius to be used as reference for frame time alignment
        sc_radius_ref = np.min(orbit_data["Radius"])
        ellipsoid_radius = self.geom_obj.get_ellipsoid_radius_from_latlon(orbit_data["Lat"], orbit_data["Lon"])
        dem_radius = self.dem_obj.get_dem_radius_from_lat_lon(orbit_data["Lat"], orbit_data["Lon"])
        dem_radius_ref_id = np.argmax(dem_radius)
        ellipsoid_radius_ref = ellipsoid_radius[dem_radius_ref_id]
        sc_radius_and_ellipsoid_time_correction = ((sc_radius_ref - orbit_data["Radius"]) + (ellipsoid_radius - ellipsoid_radius_ref)) * 2. / self._C

        # rxwin is calculated as the time delay that sets the "highest point" of the simulation at sim_top_offset_samples from the top border
        t_surface_vect_ref = 2 * (orbit_data["Radius"] - dem_radius) / self._C
        rxwin_ref = np.min(t_surface_vect_ref + sc_radius_and_ellipsoid_time_correction) - self.sim_top_offset_samples * self.dt

        self.sim_image = np.zeros([self.sim_n_samples, n_frames_to_be_simulated])   #, dtype="int16")   # see comment below
        self.uncert_image = np.zeros(self.sim_image.shape)                          #, dtype="uint8")

        if self.n_processes > 1:
            n_bytes_to_reserve = int((self.dem_obj.dem.nbytes + self.alt_sim_matrix.nbytes + self.act_sim_matrix.nbytes)*1.10)  # calculate the amount of RAM to reserve for shared objects as 5% more of what summed here
            import ray
            ray.init(num_cpus=self.n_processes, object_store_memory=n_bytes_to_reserve)

            @ray.remote
            def _simulate_frame_range_callback(id, n_processes, orbit_data, frame_range, sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, sc_radius_and_ellipsoid_time_correction, dt, dem_obj, geom_obj, alt_sim_matrix, act_sim_matrix):
                return DRSSim._simulate_frame_range(id, n_processes, orbit_data, frame_range, sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, sc_radius_and_ellipsoid_time_correction, dt, dem_obj, geom_obj, alt_sim_matrix, act_sim_matrix)

            n_frames_per_subprocess = int(np.round(n_frames_to_be_simulated / self.n_processes))
            sub_frame_range_starts = np.arange(0, n_frames_to_be_simulated, n_frames_per_subprocess)
            sub_frame_range_stops = sub_frame_range_starts + n_frames_per_subprocess
            sub_frame_range_stops[-1] = np.min([sub_frame_range_stops[-1], n_frames_to_be_simulated])

            dem_obj_id = ray.put(self.dem_obj)
            geom_obj_id = ray.put(self.geom_obj)
            alt_sim_matrix_id = ray.put(self.alt_sim_matrix)
            act_sim_matrix_id = ray.put(self.act_sim_matrix)

            result_ids = list()
            sub_frame_ranges = list()
            for id in np.arange(self.n_processes):
                sub_frame_range = np.arange(sub_frame_range_starts[id], sub_frame_range_stops[id])
                sub_frame_ranges.append(sub_frame_range)
                result_ids.append(_simulate_frame_range_callback.remote(id, self.n_processes, orbit_data, sub_frame_range, self.sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, sc_radius_and_ellipsoid_time_correction, self.dt, dem_obj_id, geom_obj_id, alt_sim_matrix_id, act_sim_matrix_id))

            results = ray.get(result_ids)

            for id in np.arange(self.n_processes):
                self.sim_image[:, sub_frame_ranges[id]] = results[id][0]
                self.uncert_image[:, sub_frame_ranges[id]] = results[id][1]

            ray.shutdown()
        else:
            self.sim_image, self.uncert_image = DRSSim._simulate_frame_range(0, 1, orbit_data, np.arange(n_frames_to_be_simulated), self.sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, sc_radius_and_ellipsoid_time_correction, self.dt, self.dem_obj, self.geom_obj, self.alt_sim_matrix, self.act_sim_matrix)

        self.sim_image = self.sim_image.astype("int16")
        self.uncert_image = self.uncert_image.astype("uint8")
        return self.sim_image, self.uncert_image

    @staticmethod       # defined static to allow calling from ray avoiding the "pickling" of the object
    def _simulate_frame_range(id, n_processes, orbit_data, frame_range, sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, sc_radius_and_ellipsoid_time_correction, dt, dem_obj, geom_obj, alt_sim_matrix, act_sim_matrix):
        _C = 299792458.

        n_frames_to_be_simulated = len(frame_range)
        sim_image = np.zeros([sim_n_samples, n_frames_to_be_simulated])         #, dtype="int16")   # int16 and uint8 saves memory but considerably slows down the processing as numpy cannot handle ints in a fast way
        uncert_image = np.zeros([sim_n_samples, n_frames_to_be_simulated])      #, dtype="uint8")

        status_str_pre = '\t'+"|\t\t"*id
        status_str_post = "\t\t|"*(n_processes-id-1)
        last_status = -1
        for count, i_frame_tbs in enumerate(frame_range):
            if count == n_frames_to_be_simulated-1:
                new_status = 100
            else:
                new_status = int(count / n_frames_to_be_simulated * 10) * 10
            if new_status != last_status:
                print(status_str_pre + '{:>2}'.format(new_status) + status_str_post)
                last_status = new_status

            local_direction = local_direction_g[i_frame_tbs]
            UTMNorthing = UTMNorthing_g[i_frame_tbs]
            UTMEasting = UTMEasting_g[i_frame_tbs]
            UTMZone = UTMZone_g[i_frame_tbs]
            sc_x_cart_local = orbit_data["X"][i_frame_tbs]
            sc_y_cart_local = orbit_data["Y"][i_frame_tbs]
            sc_z_cart_local = orbit_data["Z"][i_frame_tbs]

            UTMNorthing_surface_matrix = UTMNorthing + (alt_sim_matrix * np.sin(local_direction) - act_sim_matrix * np.cos(local_direction))
            UTMEasting_surface_matrix = UTMEasting + (alt_sim_matrix * np.cos(local_direction) + act_sim_matrix * np.sin(local_direction))
            Lat_surface_matrix, Long_surface_matrix = geom_obj.UTM_to_LL(UTMNorthing_surface_matrix, UTMEasting_surface_matrix, UTMZone)
            dem_surface_radius_matrix = dem_obj.get_dem_radius_from_lat_lon(Lat_surface_matrix, Long_surface_matrix)

            dummy_check_matrix = (dem_surface_radius_matrix == dem_obj.dummy_value)
            if dummy_check_matrix.any():
                uncert_image[:, count] = 1
                dem_surface_radius_matrix[dummy_check_matrix] = 0             # radius too far for being recorded in the simulation

            x_cart_surface_matrix, y_cart_surface_matrix, z_cart_surface_matrix = geom_obj.to_xyz_from_latlon_and_radius(Lat_surface_matrix, Long_surface_matrix, dem_surface_radius_matrix)

            d_surface = geom_obj.euclidean_distance_cart(sc_x_cart_local, sc_y_cart_local, sc_z_cart_local, x_cart_surface_matrix, y_cart_surface_matrix, z_cart_surface_matrix)

            t_surface = 2. * d_surface / _C
            sample_pos_surface = np.round(((t_surface - rxwin_ref + sc_radius_and_ellipsoid_time_correction[i_frame_tbs]) / dt)).astype("int")        # use rxwin_ref instead of input_rxwin_data_t_local

            d_surface = d_surface.flatten()
            sample_pos_surface = sample_pos_surface.flatten()
            sample_pos_surface = sample_pos_surface[(sample_pos_surface >= 0) & (sample_pos_surface < sim_n_samples)]
            for i_sample_tbs in np.arange(len(sample_pos_surface)):
                sim_image[sample_pos_surface[i_sample_tbs], count] += 1

        return (sim_image, uncert_image)

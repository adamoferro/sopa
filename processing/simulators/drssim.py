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

from processing.simulators.sim_base import _SIMBase
import numpy as np
from scipy.ndimage import gaussian_filter


class dRSsim(_SIMBase):
    """docstring for dRSsim."""

    def __init__(self, geom_obj, dem_obj, param_obj, n_processes=1, use_data_dt=True):
        super().__init__(geom_obj, dem_obj, param_obj, n_processes, use_data_dt)
        self.n_frames_offset_for_dir_estimation = 1
        self.alt_sim_spacing = 50
        self.act_sim_spacing = 50
        self.gaussian_filter_sigma = 4
        self.alt_sim_forward_max_distance = 500
        self.alt_sim_backward_max_distance = 500
        self.act_sim_max_distance = 50000

        alt_sim_range = np.arange(-self.alt_sim_backward_max_distance, self.alt_sim_forward_max_distance + self.alt_sim_spacing, self.alt_sim_spacing)
        act_sim_range = np.arange(-self.act_sim_max_distance, self.act_sim_max_distance + self.act_sim_spacing, self.act_sim_spacing)
        self.alt_sim_matrix = np.array([alt_sim_range]).T * np.ones([1, len(act_sim_range)])
        self.act_sim_matrix = np.array([act_sim_range]) * np.ones([len(alt_sim_range), 1])

    def simulate(self, orbit_data):
        return super().simulate(orbit_data, dRSsim._simulate_frame_range)

    @staticmethod       # defined static to allow calling from ray avoiding the "pickling" of the object
    def _simulate_frame_range(id, n_processes, orbit_xyz, frame_range, sim_n_samples, local_direction_g, UTMNorthing_g, UTMEasting_g, UTMZone_g, rxwin_ref, input_scr_data_shifts_alignment_t, dt, dem_obj, geom_obj, alt_sim_matrix, act_sim_matrix, gaussian_filter_sigma, further_pars=[]):
        _C = 299792458.

        n_frames_to_be_simulated = len(frame_range)
        sim_image = np.zeros([sim_n_samples, n_frames_to_be_simulated])         #, dtype="int16")   # int16 and uint8 saves memory but considerably slows down the processing as numpy cannot handle ints in a fast way
        uncert_image = np.zeros(sim_image.shape, dtype="uint8")
        ground_distance_min_image = np.nan * np.ones(sim_image.shape)
        ground_distance_max_image = np.nan * np.ones(sim_image.shape)
        first_return_lats = np.zeros(n_frames_to_be_simulated)
        first_return_lons = np.zeros(n_frames_to_be_simulated)

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
            sc_x_cart_local = orbit_xyz[0][i_frame_tbs]
            sc_y_cart_local = orbit_xyz[1][i_frame_tbs]
            sc_z_cart_local = orbit_xyz[2][i_frame_tbs]

            UTMNorthing_surface_matrix = UTMNorthing + (alt_sim_matrix * np.sin(local_direction) - act_sim_matrix * np.cos(local_direction))
            UTMEasting_surface_matrix = UTMEasting + (alt_sim_matrix * np.cos(local_direction) + act_sim_matrix * np.sin(local_direction))
            Lat_surface_matrix, Long_surface_matrix = geom_obj.UTM_to_LL(UTMNorthing_surface_matrix, UTMEasting_surface_matrix, UTMZone)
            dem_surface_radius_matrix = dem_obj.get_dem_radius_from_lat_lon(Lat_surface_matrix, Long_surface_matrix)

            ground_distance_from_nadir = np.sqrt((UTMNorthing_surface_matrix-UTMNorthing)**2 + (UTMEasting_surface_matrix-UTMEasting)**2)

            dummy_check_matrix = (dem_surface_radius_matrix == dem_obj.dummy_value)

            # smooth the DEM to avoid "stripes" in the simulation
            dem_surface_radius_matrix = gaussian_filter(dem_surface_radius_matrix, gaussian_filter_sigma)

            if dummy_check_matrix.any():
                uncert_image[:, count] = 1
                dem_surface_radius_matrix[dummy_check_matrix] = 0             # radius too far for being recorded in the simulation

            x_cart_surface_matrix, y_cart_surface_matrix, z_cart_surface_matrix = geom_obj.to_xyz_from_latlon_and_radius(Lat_surface_matrix, Long_surface_matrix, dem_surface_radius_matrix)

            d_surface = geom_obj.euclidean_distance_cart(sc_x_cart_local, sc_y_cart_local, sc_z_cart_local, x_cart_surface_matrix, y_cart_surface_matrix, z_cart_surface_matrix)

            t_surface = 2. * d_surface / _C
            sample_pos_surface = np.round(((t_surface - rxwin_ref + input_scr_data_shifts_alignment_t[i_frame_tbs]) / dt)).astype("int")        # use rxwin_ref instead of input_rxwin_data_t_local

            # set out-of-range sample positions to the farthest position (useful later for detecting the first-return position)
            sample_pos_surface[(sample_pos_surface < 0) | (sample_pos_surface >= sim_n_samples)] = sim_n_samples

            # find first return coordinates
            first_return_ids = np.unravel_index(np.argmin(sample_pos_surface), sample_pos_surface.shape)
            first_return_lats[count] = Lat_surface_matrix[first_return_ids]
            first_return_lons[count] = Long_surface_matrix[first_return_ids]

            ground_distance_from_nadir = ground_distance_from_nadir.flatten()
            sample_pos_surface = sample_pos_surface.flatten()

            # find and keep only the sample_pos_surface items within the simulation range
            sample_pos_surface_ok_ids = np.where(sample_pos_surface < sim_n_samples)[0]
            sample_pos_surface = sample_pos_surface[sample_pos_surface_ok_ids]

            # increment sim_image values by 1 each time they correspond to a valid sample_pos_surface
            np.add.at(sim_image[:, count], sample_pos_surface, 1)

            # keep only ground_distance_from_nadir items corresponding to valid sample_pos_surface items (selected before)
            ground_distance_from_nadir = ground_distance_from_nadir[sample_pos_surface_ok_ids]

            # sort distances and sample positions both by increasing distance
            sorting_indexes = np.argsort(ground_distance_from_nadir)
            ground_distance_from_nadir = ground_distance_from_nadir[sorting_indexes]
            sample_pos_surface = sample_pos_surface[sorting_indexes]

            # for each simulation sample detect the corresponding ground distance
            # note: as the arrays were first ordered by ground distance, and considering that
            # "unique" selects the first instance of the array items, the "shortest" distance
            # associated to a certain sample is kept
            unique_sample_pos_surface, unique_sample_pos_surface_ids = np.unique(sample_pos_surface, return_index=True)
            unique_ground_distance_from_nadir = ground_distance_from_nadir[unique_sample_pos_surface_ids]
            ground_distance_min_image[unique_sample_pos_surface, count] = unique_ground_distance_from_nadir

            # reverse the ordered arrays and use the sample principle to save the "longest" distances for every sample
            sample_pos_surface = np.flip(sample_pos_surface)
            ground_distance_from_nadir = np.flip(ground_distance_from_nadir)
            unique_sample_pos_surface, unique_sample_pos_surface_ids = np.unique(sample_pos_surface, return_index=True)
            unique_ground_distance_from_nadir = ground_distance_from_nadir[unique_sample_pos_surface_ids]
            ground_distance_max_image[unique_sample_pos_surface, count] = unique_ground_distance_from_nadir

        return (sim_image, uncert_image, first_return_lats, first_return_lons, ground_distance_min_image, ground_distance_max_image)

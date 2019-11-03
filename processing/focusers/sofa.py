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

# TODO
# - estimate platform velocity from data instead of using reference velocity
# - calculate _SIZE_RANGE_FFT automatically? (next pow2?)
# - estimate of average number of frames corresponding to the defined squint offset using only the frames required for processing and not the whole track

import numpy as np
from scipy.signal import chirp
from scipy.ndimage.filters import gaussian_filter1d

class SOFA(object):
    """docstring for SOFA."""
    _C = 299792458.

    def __init__(self, instr_obj, geom_obj, dem_radius_profile, data_obj, param_obj, n_processes=1, read_block_size_max=5000, debug_mode=False):
        self.geom_obj = geom_obj
        self.data_obj = data_obj

        if self.data_obj._data_type != "raw":
            raise ValueError("ERROR: input data are not raw.")

        self._dem_radius_profile = dem_radius_profile

        self._frequency = instr_obj.frequency
        self._bandwidth = instr_obj.bandwidth
        self._fast_time_sampling_frequency = instr_obj.fast_time_sampling_frequency
        self._wc = 2. * np.pi * instr_obj.frequency
        self._lmbda_min = self._C / (instr_obj.frequency + instr_obj.bandwidth/2.)    # wavelength at highest frequency
        self._dt = 1. / instr_obj.fast_time_sampling_frequency
        self._chirp_duration = instr_obj.chirp_duration
        self._n_input_frames = data_obj.data.shape[1]
        self._n_input_samples = data_obj.data.shape[0]

        self._pri = data_obj.pri
        self._presumming = data_obj.presumming

        self._frames_to_be_focused_start = param_obj.frame_start
        self._frames_to_be_focused_end = param_obj.frame_end
        if param_obj.frame_skip == 0:
            self._skip_meters = param_obj.meters_skip
            self._skip_frames = 0
        else:
            self._skip_frames = param_obj.frame_skip
            self._skip_meters = None
        self._n_half_aperture_frames = param_obj.frame_half_aperture

        self._SIZE_RANGE_FFT = 4096
        self._N_OUTPUT_SAMPLES = 667*2

        self._CORRECT_RANGE_MIGRATION = param_obj.correct_range_migration
        self._COMPENSATE_DOPPLER_PHASE = param_obj.compensate_doppler_phase

        self._n_gaussian_filters = param_obj.gauss_rec_n_filters
        self._gaussian_alpha = param_obj.gauss_rec_alpha

        self._squint_angle_deg = param_obj.squint_angle                 # =0, look nadir; >0, look ahead; <0 look behind

        self._read_block_size_max = read_block_size_max
        self._debug_mode = debug_mode

        # ----------- CALCULATED PARAMETERS -----------
        self._pri_data = self._pri * self._presumming
        self._n_aperture_frames = 2*self._n_half_aperture_frames+1        # center sample is included and in the middle of aperture
        self._coherence_time = (self._n_aperture_frames-1) * self._pri_data
        self._du_global = self._pri_data * instr_obj.ref_v_platform

        # TODO
        self._focusing_ML_doppler_half_bandwidth = param_obj.multilooking_doppler_bandwidth / 2.
        self._doppler_half_window_width = np.round(self._focusing_ML_doppler_half_bandwidth * self._coherence_time).astype("int")-1
        self._doppler_half_window_width = np.max([0, self._doppler_half_window_width])      # avoid negative window widths

        if self._debug_mode:
            print("n_input_frames =", self._n_input_frames)
            print("coherence time =", self._coherence_time)
            print("Doppler half window width =", self._doppler_half_window_width)

        self._squint_angle_rad = self._squint_angle_deg * np.pi / 180.

        if n_processes > 1:
            try:
                import ray
            except ImportError:
                n_processes = 1
                print("WARNING: Multiprocessing not possible on this machine (failed to import module ""ray""). Working with 1 process.")
        elif n_processes < 1:
            n_processes = 1

        self.n_processes = n_processes

        # gaussian filter parameters
        #   power coefficients
        self._gf_power_coefficients = list()
        if self._n_gaussian_filters > 0:
            self._gf_power_coefficients.append(self._gaussian_alpha)
            for i_gf in np.arange(self._n_gaussian_filters-1):
                self._gf_power_coefficients.append((1.-np.sum(self._gf_power_coefficients))*self._gaussian_alpha)
            self._gf_power_coefficients.append(1-np.sum(self._gf_power_coefficients))

        #   std devs
        self._gf_std_devs = list()
        gf_std_dev_max = (np.round((self._focusing_ML_doppler_half_bandwidth * self._coherence_time) * 10.)) / 10.
        self._gf_std_devs = [(i_gf+1)*gf_std_dev_max/self._n_gaussian_filters for i_gf in np.arange(self._n_gaussian_filters)]

        if debug_mode:
            print("Gaussian filters number:", self._n_gaussian_filters)
            print("   power coefficients =", self._gf_power_coefficients, "(including no filter coeff.)")
            print("   std devs =", self._gf_std_devs)
        else:
            if debug_mode:
                print("NO Gaussian filtering.")

        self._range_weighting_window_original = 1.
        rwfs = param_obj.range_weighting_function
        if rwfs != "none":
            if rwfs == "kaiser":
                self._range_weighting_window_original = np.kaiser(self._SIZE_RANGE_FFT/2., param_obj.range_weighting_kaiser_beta)
            else:
                self._range_weighting_window_original = eval("np." + rwfs + "(self._SIZE_RANGE_FFT/2.)")

        self._azimuth_weighting_windows = 1.
        awfs = param_obj.azimuth_weighting_function
        if awfs != "none":
            if awfs == "kaiser":
                self._azimuth_weighting_window_original = np.kaiser(self._n_aperture_frames, param_obj.azimuth_weighting_kaiser_beta)
            else:
                self._azimuth_weighting_window_original = eval("np." + awfs + "(self._n_aperture_frames)")

    def _create_synthetic_chirp_conj(self):
        t_chirp = np.arange(0, self._chirp_duration, 1. / self._fast_time_sampling_frequency)
        range_weighting_window = np.concatenate((self._range_weighting_window_original, np.zeros(int(self._SIZE_RANGE_FFT/2))), 0)

        f0 = self._frequency + self._bandwidth/2.
        t1 = self._chirp_duration
        f1 = self._frequency - self._bandwidth/2.
        chirp2_real = chirp(t_chirp, f0, t1, f1)

        chirp_fft = (np.fft.fft(chirp2_real, self._SIZE_RANGE_FFT)) / np.sqrt(self._SIZE_RANGE_FFT)
        chirp_fft = chirp_fft*range_weighting_window
        chirp_fft_conj = np.conj(chirp_fft).T

        return chirp_fft_conj

    def focus(self):
        # Estimation of average number of frames corresponding to the defined squint offset
        input_sc_alt = self.data_obj.orbit_data["Radius"] - self._dem_radius_profile
        mean_altitude = np.average(input_sc_alt)
        squint_offset_approx_distance = -(mean_altitude*np.tan(self._squint_angle_rad))*1.01    # TODO *** 1.01 because (abs) result is tipically underestimated
        squint_offset_approx_frames = int(np.round(squint_offset_approx_distance / self._du_global))

        # calculate actual start and end frames according to aperture and squint angle
        if self._frames_to_be_focused_start == -1:
            self._frames_to_be_focused_start = np.max([self._n_half_aperture_frames-squint_offset_approx_frames+1, 0])

        if self._frames_to_be_focused_end == -1:
            self._frames_to_be_focused_end = np.min([self._n_input_frames-self._n_half_aperture_frames-1-squint_offset_approx_frames, self._n_input_frames-1])

        if self._frames_to_be_focused_end < self._frames_to_be_focused_start:
            print("ERROR: last frame to be focused before first frame to be focused.")
            return None

        if self._debug_mode:
            print("SC mean altitude =", mean_altitude)
            print("squint_offset_approx_distance =", squint_offset_approx_distance)
            print("squint_offset_approx_frames =", squint_offset_approx_frames)

        # calculation of range of frames to be focused
        if self._skip_frames != 0:      # when the skip is expressed in frames
            frames_to_be_focused = np.arange(self._frames_to_be_focused_start, self._frames_to_be_focused_end, self._skip_frames)
            min_skip_frames = self._skip_frames
            max_skip_frames = self._skip_frames
        else:                           # when the skip is expressed in meters
            frames_to_be_focused, min_skip_frames, max_skip_frames = self.geom_obj.get_frame_ids_from_skip_distance(self.data_obj.orbit_data["Lat"][self._frames_to_be_focused_start:self._frames_to_be_focused_end], self.data_obj.orbit_data["Lon"][self._frames_to_be_focused_start:self._frames_to_be_focused_end], self._skip_meters)
            frames_to_be_focused += self._frames_to_be_focused_start

        # correction of the last frame to be focused
        self._frames_to_be_focused_end = frames_to_be_focused[-1]
        n_frames_to_be_focused = len(frames_to_be_focused)

        self._block_overlapping_size = self._n_aperture_frames - max_skip_frames
        if self._block_overlapping_size < 0:
            raise ValueError("ERROR: (max) frame skip greater than aperture size is not handled.")

        if self._debug_mode:
            print("frames_to_be_focused_start =", self._frames_to_be_focused_start)
            print("frames_to_be_focused_end =", self._frames_to_be_focused_end)
            print("n_frames_to_be_focused =", n_frames_to_be_focused)
            print("last frame to be focused =", frames_to_be_focused[-1])

        # calculation of range of frames to be processed to focus the selected frames
        frames_to_be_processed_start = self._frames_to_be_focused_start - self._n_half_aperture_frames + squint_offset_approx_frames
        if frames_to_be_processed_start > self._frames_to_be_focused_start:
            frames_to_be_processed_start = self._frames_to_be_focused_start
            if self._debug_mode:
                print("\nNOTICE: first frame to be processed = first frame to be focused, as aperture is very delayed.")

        frames_to_be_processed_end = self._frames_to_be_focused_end + self._n_half_aperture_frames + squint_offset_approx_frames + 1
        if frames_to_be_processed_end < self._frames_to_be_focused_end:
            frames_to_be_processed_end = self._frames_to_be_focused_end
            if self._debug_mode:
                print("\nNOTICE: last frame to be processed = last frame to be focused, as aperture is very anticipated.\n")

        frames_to_be_processed = np.arange(frames_to_be_processed_start, frames_to_be_processed_end+1, dtype="int")
        n_frames_to_be_processed = len(frames_to_be_processed)

        if self._debug_mode:
            print("n_frames_to_be_processed =", n_frames_to_be_processed)
            print("frames_to_be_processed_start =", frames_to_be_processed_start)
            print("frames_to_be_processed_end =", frames_to_be_processed_end)

        if frames_to_be_processed_start < 0:
            print("ERROR: first frame to be processed is before the beginning of input raw data.")
            return None

        if frames_to_be_processed_end >= self._n_input_frames:
            print("ERROR: the last frame to be processed is after the end of input raw data.")
            return None

        # calculation of ids of frames to be focused wrt frames_to_be_processed
        frames_to_be_focused_within_frames_to_be_processed = frames_to_be_focused - frames_to_be_processed_start

        # select only data related to frames to be processed
        input_rxwin_data_t = self.data_obj.orbit_data["rxwin_time_s"][frames_to_be_processed]
        input_scx_data = self.data_obj.orbit_data["X"][frames_to_be_processed]
        input_scy_data = self.data_obj.orbit_data["Y"][frames_to_be_processed]
        input_scz_data = self.data_obj.orbit_data["Z"][frames_to_be_processed]
        input_lat_data = self.data_obj.orbit_data["Lat"][frames_to_be_processed]
        input_long_data = self.data_obj.orbit_data["Lon"][frames_to_be_processed]
        input_scr_data = self.data_obj.orbit_data["Radius"][frames_to_be_processed]
        dem_radius_profile = self._dem_radius_profile[frames_to_be_processed]

        input_theta = np.arctan2(input_scy_data, input_scx_data)
        input_phi = np.arccos(input_scz_data/input_scr_data)

        # shifts for final frame alignment
        input_scr_data_alignment = input_scr_data - self.geom_obj.get_ellipsoid_radius_from_latlon(input_lat_data, input_long_data)
        input_scr_data_shifts_alignment = max(input_scr_data_alignment[frames_to_be_focused_within_frames_to_be_processed])-input_scr_data_alignment
        input_scr_data_shifts_alignment_t = (input_scr_data_shifts_alignment*2. / self._C)

        input_rxwin_data_t_shifts = input_rxwin_data_t-min(input_rxwin_data_t[frames_to_be_focused_within_frames_to_be_processed])

        # frame sample offsets for final frame alignment
        sample_start_shifts = np.round((input_rxwin_data_t_shifts + input_scr_data_shifts_alignment_t) / self._dt).astype("int")
        sample_start_shifts_max = np.max(sample_start_shifts[frames_to_be_focused_within_frames_to_be_processed])

        self._N_OUTPUT_SAMPLES += sample_start_shifts_max       # update N_OUTPUT_SAMPLES
        OUTPUT_SAMPLES = np.arange(self._N_OUTPUT_SAMPLES)

        # create here azimuth weighting windows because _N_OUTPUT_SAMPLES has been updated
        self._azimuth_weighting_windows = self._azimuth_weighting_window_original * np.ones([self._N_OUTPUT_SAMPLES, 1])

        delay_start_shifts = input_rxwin_data_t_shifts + input_scr_data_shifts_alignment_t

        # create complex phasor
        data_samples = np.arange(self._n_input_samples)
        phasor_samples = np.arange(self._SIZE_RANGE_FFT)

        # baseband conversion
        phasor_const = 2. * np.pi * (-self._frequency)
        phasor = np.exp(1j * (phasor_const * (phasor_samples / self._fast_time_sampling_frequency)))

        chirp_fft_conj = self._create_synthetic_chirp_conj()

        t_frame_base = self._dt * OUTPUT_SAMPLES
        fftfreqs = np.fft.fftfreq(self._SIZE_RANGE_FFT, self._dt)

        # calculate block size for parallel processing
        read_block_size_optimal_approx = np.ceil(n_frames_to_be_processed / self.n_processes) + self._block_overlapping_size
        if read_block_size_optimal_approx > self._read_block_size_max:
            read_block_size_optimal_approx = self._read_block_size_max
        n_frames_to_be_focused_within_block_max = int(np.ceil((read_block_size_optimal_approx - self._n_half_aperture_frames*2) / min_skip_frames) + 1)
        n_blocks = int(np.ceil(n_frames_to_be_focused / n_frames_to_be_focused_within_block_max))
        n_frames_to_be_focused_in_last_block = n_frames_to_be_focused - (n_blocks-1)*n_frames_to_be_focused_within_block_max
        n_frames_to_be_focused_per_block = np.concatenate((n_frames_to_be_focused_within_block_max*np.ones(n_blocks-1), [n_frames_to_be_focused_in_last_block])).astype("int")
        n_frames_to_be_focused_per_block_cum = np.cumsum(n_frames_to_be_focused_per_block)

        if self._debug_mode:
            msg = ""
            if read_block_size_optimal_approx == self._read_block_size_max:
                msg = "(MAXIMUM ALLOWED SIZE)"
            print("read_block_size_optimal_approx =", read_block_size_optimal_approx, msg)
            print("n_frames_to_be_focused_within_block_max =", n_frames_to_be_focused_within_block_max)
            print("N BLOCKS =", n_blocks)

        # initialize output images
        output_image = np.zeros([self._N_OUTPUT_SAMPLES, n_frames_to_be_focused])

        # create processing input block parameters
        input_parameters = list()
        for block_id in np.arange(n_blocks):
            frames_to_be_focused_in_this_block = frames_to_be_focused[int(block_id * n_frames_to_be_focused_within_block_max):int(block_id * n_frames_to_be_focused_within_block_max + n_frames_to_be_focused_per_block[block_id])]
            block_frame_to_be_focused_start = frames_to_be_focused_in_this_block[0]
            block_frame_to_be_focused_end = frames_to_be_focused_in_this_block[-1]
            input_data_frame_start = block_frame_to_be_focused_start - self._n_half_aperture_frames + squint_offset_approx_frames
            input_data_frame_end = np.min([block_frame_to_be_focused_end + self._n_half_aperture_frames + 1 + squint_offset_approx_frames, frames_to_be_processed_end])
            input_parameters.append((input_data_frame_start, input_data_frame_end, frames_to_be_focused_in_this_block, frames_to_be_processed_start, self._SIZE_RANGE_FFT, data_samples, phasor, self._C, self._wc, input_rxwin_data_t, chirp_fft_conj, self._N_OUTPUT_SAMPLES, OUTPUT_SAMPLES, self._dt, input_scr_data, input_theta, input_phi, self._n_half_aperture_frames, self._n_aperture_frames, self._CORRECT_RANGE_MIGRATION, input_scx_data, input_scy_data, input_scz_data, self._COMPENSATE_DOPPLER_PHASE, self._lmbda_min, self._azimuth_weighting_windows, self._n_gaussian_filters, self._gf_std_devs, self._gf_power_coefficients, self._doppler_half_window_width, squint_offset_approx_frames, t_frame_base, dem_radius_profile, fftfreqs, delay_start_shifts))

        if self.n_processes > 1:
            n_bytes_to_reserve = int((self.data_obj.data.nbytes)*1.25)  # calculate the amount of RAM to reserve for shared objects as 75% more of raw data
            import ray
            ray.init(num_cpus=self.n_processes, redis_max_memory=int(n_bytes_to_reserve*.85), object_store_memory=n_bytes_to_reserve)

            @ray.remote
            def _block_focus_callback(id, input_parameters, do, go):
                return SOFA._block_focus(id, *input_parameters, do, go)

            geom_obj_id = ray.put(self.geom_obj)
            data_obj_id = ray.put(self.data_obj)

            result_ids = list()
            for block_id in np.arange(n_blocks):
                result_ids.append(_block_focus_callback.remote(block_id, input_parameters[block_id], data_obj_id, geom_obj_id))

            results = ray.get(result_ids)

            output_image[:, 0:n_frames_to_be_focused_per_block[0]] = results[0][0]
            for block_id in np.arange(1, n_blocks):
                output_image[:, n_frames_to_be_focused_per_block_cum[block_id-1]:n_frames_to_be_focused_per_block_cum[block_id-1] + n_frames_to_be_focused_per_block[block_id]] = results[block_id][0]

            ray.shutdown()
        else:
            output_image[:, 0:n_frames_to_be_focused_per_block[0]] = SOFA._block_focus(0, *(input_parameters[0]), self.data_obj, self.geom_obj)[0]
            for block_id in np.arange(1, n_blocks):
                output_image[:, n_frames_to_be_focused_per_block_cum[block_id-1]:n_frames_to_be_focused_per_block_cum[block_id-1] + n_frames_to_be_focused_per_block[block_id]] = SOFA._block_focus(block_id, *(input_parameters[block_id]), self.data_obj, self.geom_obj)[0]

        return output_image

    @staticmethod
    def _block_focus(block_id, input_data_frame_start, input_data_frame_end, frames_to_be_focused_in_this_block, frames_to_be_processed_start, sizeFFT1, data_samples, phasor, c, wc, input_rxwin_data_t, chirp_fft_conj, N_OUTPUT_SAMPLES, OUTPUT_SAMPLES, dt_compr, input_scr_data, input_theta, input_phi, n_half_aperture_frames, n_aperture_frames, CORRECT_RANGE_MIGRATION, input_scx_data, input_scy_data, input_scz_data, USE_DOPPLER_CORRECTION, lmbda_min, azimuth_weighting_windows, N_GAUSSIAN_FILTERS, std_devs, gf_power_coefficients, doppler_half_window_width, squint_offset_approx_frames, t_frame_base, dem_radius_profile, fftfreqs, delay_start_shifts, data_obj, geom_obj):

        input_data = data_obj.data[:, input_data_frame_start:input_data_frame_end]
        read_block_size_real = input_data.shape[1]

        # NOTE: frames to be focused could be OUTSIDE the block because of squint!!!
        frames_to_be_focused_within_block = frames_to_be_focused_in_this_block
        n_frames_to_be_focused_within_block = len(frames_to_be_focused_within_block)

        print("PROCESSING BLOCK", block_id)
        print("  FOCUSING FROM ABS ", frames_to_be_focused_within_block[0], "TO", frames_to_be_focused_within_block[-1])
        print("  n_frames_to_be_focused_within_block =", n_frames_to_be_focused_within_block)

        output_image = np.zeros([N_OUTPUT_SAMPLES, n_frames_to_be_focused_within_block])
        signals = np.zeros([sizeFFT1, read_block_size_real], dtype="complex")

        # TODO *** check if possible to convert in matrix calculations without killing RAM (as already happened...)
        for i_frame in np.arange(read_block_size_real):
            signals[data_samples, i_frame] = input_data[:, i_frame]     # DELETED * phasor[data_samples]

        print("...matched filtering...")
        signals_fft = np.fft.fft(signals, axis=0)      # spectrum: 0-positive ... negative-0

        # FOR VERSION (slower and less RAM)   TODO *** try to convert in matrix format
        signals_fft_matched = np.zeros([len(chirp_fft_conj), read_block_size_real], dtype="complex")
        for i_frame in np.arange(read_block_size_real):
            signals_fft_matched[:, i_frame] = signals_fft[:, i_frame] * chirp_fft_conj

        signals = None
        signals_fft = None

        i_frame_focusing = 0
        angle_between_aperture_start_and_ref_sample = 0
        angle_between_aperture_center_and_ref_sample = 0
        angle_between_aperture_end_and_ref_sample = 0

        # aperture has final length = 2*n_half_aperture_frames+1
        aperture_range_offset = np.arange(squint_offset_approx_frames-n_half_aperture_frames, squint_offset_approx_frames+n_half_aperture_frames+1, dtype="int")
        for ref_frame_abs in frames_to_be_focused_within_block:
            ref_frame = ref_frame_abs - frames_to_be_processed_start        # referred to the frames_to_be_processed interval
            aperture_frames = ref_frame + aperture_range_offset

            if CORRECT_RANGE_MIGRATION:
                # range migration is corrected wrt DEM radius: not necessary to correct for every single sample, error is negligible TODO *** check (already done, but investigate more)
                r_sample_ref = dem_radius_profile[ref_frame]
                x_sample_ref = r_sample_ref*(np.cos(input_theta[ref_frame]))*(np.sin(input_phi[ref_frame]))
                y_sample_ref = r_sample_ref*(np.sin(input_theta[ref_frame]))*(np.sin(input_phi[ref_frame]))
                z_sample_ref = r_sample_ref*np.cos(input_phi[ref_frame])

                sc_distance_to_ref_sample_within_aperture = np.sqrt((input_scx_data[aperture_frames]-x_sample_ref)**2+(input_scy_data[aperture_frames]-y_sample_ref)**2+(input_scz_data[aperture_frames]-z_sample_ref)**2)
                sc_distance_to_ref_sample_within_aperture_ref = np.sqrt((input_scx_data[ref_frame]-x_sample_ref)**2+(input_scy_data[ref_frame]-y_sample_ref)**2+(input_scz_data[ref_frame]-z_sample_ref)**2)

                # t of beginning of frames wrt tx -rxwin
                t_within_frames_within_aperture = sc_distance_to_ref_sample_within_aperture*2./c-input_rxwin_data_t[aperture_frames]
                t_within_frames_within_aperture_ref = sc_distance_to_ref_sample_within_aperture_ref*2./c-input_rxwin_data_t[ref_frame]

                sc_position_wrt_ref_sample = (input_scx_data[ref_frame]-x_sample_ref, input_scy_data[ref_frame]-y_sample_ref, input_scz_data[ref_frame]-z_sample_ref)
                sc_aperture_start_position_wrt_ref_sample = (input_scx_data[ref_frame+squint_offset_approx_frames-n_half_aperture_frames]-x_sample_ref, input_scy_data[ref_frame+squint_offset_approx_frames-n_half_aperture_frames]-y_sample_ref, input_scz_data[ref_frame+squint_offset_approx_frames-n_half_aperture_frames]-z_sample_ref)
                sc_aperture_center_position_wrt_ref_sample = (input_scx_data[ref_frame+squint_offset_approx_frames]-x_sample_ref, input_scy_data[ref_frame+squint_offset_approx_frames]-y_sample_ref, input_scz_data[ref_frame+squint_offset_approx_frames]-z_sample_ref)
                sc_aperture_end_position_wrt_ref_sample = (input_scx_data[ref_frame+squint_offset_approx_frames+n_half_aperture_frames]-x_sample_ref, input_scy_data[ref_frame+squint_offset_approx_frames+n_half_aperture_frames]-y_sample_ref, input_scz_data[ref_frame+squint_offset_approx_frames+n_half_aperture_frames]-z_sample_ref)

                angle_between_aperture_start_and_ref_sample += geom_obj.angle_between(sc_position_wrt_ref_sample, sc_aperture_start_position_wrt_ref_sample)
                angle_between_aperture_center_and_ref_sample += geom_obj.angle_between(sc_position_wrt_ref_sample, sc_aperture_center_position_wrt_ref_sample)
                angle_between_aperture_end_and_ref_sample += geom_obj.angle_between(sc_position_wrt_ref_sample, sc_aperture_end_position_wrt_ref_sample)

                t_sample_compr_ref = (t_frame_base.T+input_rxwin_data_t[ref_frame])
                r_sample_ref = np.array([input_scr_data[ref_frame]-t_sample_compr_ref*c/2.]).T
                x_sample_ref = r_sample_ref*(np.cos(input_theta[ref_frame]))*(np.sin(input_phi[ref_frame]))
                y_sample_ref = r_sample_ref*(np.sin(input_theta[ref_frame]))*(np.sin(input_phi[ref_frame]))
                z_sample_ref = r_sample_ref*np.cos(input_phi[ref_frame])

                sc_distance_to_ref_sample_within_aperture_matrix = np.sqrt((input_scx_data[aperture_frames]*np.ones([N_OUTPUT_SAMPLES, 1])-x_sample_ref*np.ones([1, n_aperture_frames]))**2+(input_scy_data[aperture_frames]*np.ones([N_OUTPUT_SAMPLES, 1])-y_sample_ref*np.ones([1, n_aperture_frames]))**2+(input_scz_data[aperture_frames]*np.ones([N_OUTPUT_SAMPLES, 1])-z_sample_ref*np.ones([1, n_aperture_frames]))**2)

                delays_within_aperture = -(t_within_frames_within_aperture-t_within_frames_within_aperture_ref)

                if USE_DOPPLER_CORRECTION:
                    # phase correction "US style", only referred to a reference sample
                    doppler_phase_corr_within_aperture = -4.*np.pi*sc_distance_to_ref_sample_within_aperture_matrix/lmbda_min
                else:
                    doppler_phase_corr_within_aperture = np.zeros([N_OUTPUT_SAMPLES, n_aperture_frames])
            else:
                delays_within_aperture = np.zeros(n_aperture_frames)
                doppler_phase_corr_within_aperture = np.zeros([N_OUTPUT_SAMPLES, n_aperture_frames])

            signals_matched_within_aperture = np.zeros([N_OUTPUT_SAMPLES, n_aperture_frames], dtype="complex")

            for i_aperture_frame in aperture_frames-aperture_frames[0]:
                delay_vector_fft = np.exp(-1j*2.*np.pi*fftfreqs*(delays_within_aperture[i_aperture_frame]+delay_start_shifts[ref_frame]))
                signals_matched_within_aperture[:, i_aperture_frame] = np.fft.ifft(signals_fft_matched[:, i_aperture_frame+aperture_frames[0]+frames_to_be_processed_start-input_data_frame_start]*delay_vector_fft)[OUTPUT_SAMPLES]

            # Doppler correction on whole aperture
            signals_matched_within_aperture *= np.exp(1j*doppler_phase_corr_within_aperture)
            signals_matched_within_aperture *= azimuth_weighting_windows

            signals_matched_within_aperture_azimuth_fft = np.fft.fftshift(np.fft.fft(signals_matched_within_aperture, axis=1), axes=(1,))
            signals_matched_within_aperture_azimuth_fft = np.abs(signals_matched_within_aperture_azimuth_fft)**2.

            if N_GAUSSIAN_FILTERS > 0:
                signals_matched_within_aperture_upsampled_azimuth_fft_g = list()

                for std_dev_tmp in std_devs:
                    signals_matched_within_aperture_upsampled_azimuth_fft_g.append(gaussian_filter1d(signals_matched_within_aperture_azimuth_fft, std_dev_tmp, axis=1))

                signals_matched_within_aperture_azimuth_fft *= gf_power_coefficients[0]
                for i_std_dev in np.arange(N_GAUSSIAN_FILTERS):
                    signals_matched_within_aperture_azimuth_fft += gf_power_coefficients[N_GAUSSIAN_FILTERS-1-i_std_dev]*signals_matched_within_aperture_upsampled_azimuth_fft_g[i_std_dev]

            zero_doppler_search_half_window_width = 0
            frame_power_max_index = (n_half_aperture_frames-zero_doppler_search_half_window_width)+np.argmax(signals_matched_within_aperture_azimuth_fft[:, n_half_aperture_frames-zero_doppler_search_half_window_width:n_half_aperture_frames+zero_doppler_search_half_window_width+1], axis=1)
            frame_multilooking_start = frame_power_max_index-doppler_half_window_width
            output_frame_power = np.zeros(N_OUTPUT_SAMPLES)
            for i_output_sample in OUTPUT_SAMPLES:
                output_frame_power[i_output_sample] = np.sum(signals_matched_within_aperture_azimuth_fft[i_output_sample, frame_multilooking_start[i_output_sample]:frame_multilooking_start[i_output_sample]+doppler_half_window_width*2+1])/(doppler_half_window_width*2+1)

            output_image[OUTPUT_SAMPLES, i_frame_focusing] = output_frame_power[:]

            i_frame_focusing += 1

        angle_between_aperture_start_and_ref_sample /= n_frames_to_be_focused_within_block
        angle_between_aperture_center_and_ref_sample /= n_frames_to_be_focused_within_block
        angle_between_aperture_end_and_ref_sample /= n_frames_to_be_focused_within_block

        if abs(squint_offset_approx_frames) < n_half_aperture_frames:
            angle_between_aperture_start_and_ref_sample *= -1
        aperture_angle = np.abs(angle_between_aperture_end_and_ref_sample-angle_between_aperture_start_and_ref_sample)

        return output_image, angle_between_aperture_center_and_ref_sample, aperture_angle

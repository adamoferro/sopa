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

import sys
import numpy as np
from log.log import logger
from input.sharad.instrument import SHARAD
from input.sharad.data.edr.pds import EDR as EDR_PDS
from geometry.coordinate_converter import MARS_IAU2000
from dem.mola.mola128 import MOLA128
from processing.parameters import FocuserParameters
from processing.simulators.sfrssim import sfRSsim
# from processing.simulators.drssim import dRSsim         # uncomment if needed
from processing.focusers.sofa import SOFA
from io_utils import envi
from plot_utils.tracks import TrackPlotOnDEM, TrackAndFirstReturnsPlotOnDEM
from plot_utils.radargrams import SimulationPlot, GroundDistancePlot, RadargramPlot


def main(argv=None):
    """
    This example:
    1) creates an instance of a SHARAD object
    2) defines a coordinate converter based on Mars ellipsoid characteristics
    3) loads MOLA128 radius data
    4) creates an object containing the parameters used in the focuser (that
       will be used also by the simulator)
    5) defines a square facets simulator using the pre-defined coordinate converter, the
       MOLA128 DEM and the parameter object
    6) loads a PDS Italian EDR product
    7) simulates the radargram using the orbit data and saves the result on disk
    8) creates some plots: simulation, first return positions, minimum ground distance
       of clutter contributions for each simulation sample
    9) focuses the data, saves the result on disk and plots the focused radargram

    Examples of how to read other data formats are commented at the end of the code.
    """

    MOLA128_FILENAME_BASE = "/yourpath/mola-radius-data"

    EDR_PDS_INPUT_FILENAME_BASE = "/yourpath/e_1260201_001_ss19_700_a"   # OBS folder path+base file name (without suffixes)
    EDR_PDS_INPUT_DATASET_NAME = "edr1260201"
    EDR_PDS_OUTPUT_SIM_FILENAME_BASE = "/yourpath/PDS-EDR-1260201_SIM"
    EDR_PDS_OUTPUT_FOC_FILENAME_BASE = "/yourpath/PDS-EDR-1260201_FOC"

    FRAME_START = 80000
    FRAME_END = 100000
    METERS_SKIP = 450
    FRAME_HALF_APERTURE = 128       # use 758 for "standard" focusing with presumming = 4
    SQUINT_ANGLE = 0

    lg = logger(verbose=True)

    # create SHARAD instrument instance
    i = SHARAD()

    # Mars ellipsoid characteristics
    cc = MARS_IAU2000()

    # the MOLA128 radius data are contained in a single ENVI file
    m = MOLA128(filename_base=MOLA128_FILENAME_BASE)        #mola_mars_data.hdr
    m.read()

    # definition of focuser and simulator parameters
    pars = FocuserParameters(frame_start=FRAME_START, frame_end=FRAME_END, meters_skip=METERS_SKIP, frame_half_aperture=FRAME_HALF_APERTURE, squint_angle=SQUINT_ANGLE)

    # definition of a square facets RS simulator using the Mars coordinate converter and the MOLA128 DEM
    s = sfRSsim(geom_obj=cc, dem_obj=m, param_obj=pars, n_processes=6)       # simulation parameters are contained in the DRSSim class

    # ...or dRSsim:
    # s = dRSsim(geom_obj=cc, dem_obj=m, param_obj=pars, n_processes=6)       # simulation parameters are contained in the DRSSim class

    # EDR PDS
    fn_base = EDR_PDS_INPUT_FILENAME_BASE
    d = EDR_PDS(dataset_name=EDR_PDS_INPUT_DATASET_NAME, filename_base=fn_base, logger=lg)
    d.load(mode="full")
    d.generate_orbit_data()
    if d.orbit_data is not None:

        # plot the track path on the DEM
        tpod = TrackPlotOnDEM(m, d.orbit_data["Lat"], d.orbit_data["Lon"], title="Spacecraft ground track on DEM (full track)")
        tpod.plot()

        # simulate surface clutter radargram
        sim_image, uncert_image, first_return_lats, first_return_lons, ground_distance_min_image, ground_distance_max_image = s.simulate(d.orbit_data)

        # save results on disk
        envi.write(EDR_PDS_OUTPUT_SIM_FILENAME_BASE, sim_image)
        envi.write(EDR_PDS_OUTPUT_SIM_FILENAME_BASE + "_uncert", uncert_image)

        # plot the simulated radargram
        # only the top part of the simulation is shown (the top point is calculated dinamically)
        top_point = np.min(np.argmin(np.isnan(ground_distance_min_image).astype(int), axis=0)) - 50
        bottom_point = top_point + 1000
        sp = SimulationPlot(sim_image, top_point=top_point, bottom_point=bottom_point, title="Radargram simulation with sfRSsim\n(square facets and cos^exp function)")
        sp.plot()

        # plot the track path and the location of the first simulated returns on the DEM
        tfrpod = TrackAndFirstReturnsPlotOnDEM(m, d.orbit_data["Lat"], d.orbit_data["Lon"], first_return_lats, first_return_lons)
        tfrpod.plot()

        # plot the image showing the minimum ground distance from nadir of each simulated range sample
        gdp = GroundDistancePlot(ground_distance_min_image/1000., top_point=top_point, bottom_point=bottom_point, title="Distance [km] of the closest facet contributing\nto each radargram sample simulation")
        gdp.plot()

        # get the DEM radius profile related to the spacecraft track, to be used later by the focuser
        dem_radius_profile = m.get_dem_radius_from_lat_lon(d.orbit_data["Lat"], d.orbit_data["Lon"])

        # free the memory occupied by the DEM (not mandatory, but it may be useful if system RAM is limited)
        del m

        if d.data is not None:
            # create a SOFA instance and focus the input radargram
            f = SOFA(i, cc, dem_radius_profile, d, pars, n_processes=4, debug_mode=True)
            focused_rdr = f.focus()
            if focused_rdr is not None:
                envi.write(EDR_PDS_OUTPUT_FOC_FILENAME_BASE + "_ha-" + str(pars.frame_half_aperture) + "_s-" + str(pars.squint_angle), focused_rdr)

                # plot the focused radargram
                rp = RadargramPlot(focused_rdr, top_point=top_point, bottom_point=bottom_point, title="Radargram")
                rp.plot()

    # ------ OTHER INPUT EXAMPLES ------
    # NOTES:
    # 1) focusing is possible only on EDR-type input formats
    # 2) the code for simulating the orbits is the same as in the EDR PDS example
    #
    # # EDR CO-SHARPS
    # from input.sharad.data.edr.cosharps import EDR as EDR_COSHARPS
    # fn_base = "/yourpath/OBS_1260201000_1"
    # d = EDR_COSHARPS(dataset_name="OBS_1260201000_1", filename_base=fn_base, logger=lg)
    # d.load(mode="ancillary")            # use "ancillary" if you only want to check the orbit or to simulate, reading is faster
    # d.generate_orbit_data()
    #
    # # RDR PDS
    # from input.sharad.data.rdr.pds import RDR as RDR_PDS
    # fn_base = "/yourpath/r_1260201_001_ss19_700_a"
    # d = RDR_PDS(dataset_name="rdr1260201", filename_base=fn_base, logger=lg)
    # d.load()                        # note 1: d.data is np.complex64 for PDS RDR data
    #                                 # note 2: no "mode" is provided, as all data are contained in one unique file
    # envi.write("/yourpath/test_pds_rdr_data", d.data)      # data is automatically converted into power. If complex output is needed, use the complex_output=True flag
    # d.generate_orbit_data()
    #
    # # RDR PDS US
    # # Note: this data are already undersampled in the along-track direction
    # from input.sharad.data.rdr.pds_us import RDR as RDR_PDS_US
    # fn_base = "/yourpath/s_01260201"
    # d = RDR_PDS_US(dataset_name="s_01260201", filename_base=fn_base, logger=lg)
    # d.load(mode="full")
    # d.generate_orbit_data()
    # envi.write("/yourpath/test_pds-us-rdr_data", d.data)


if __name__ == "__main__":
    sys.exit(main())

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
from log.log import logger
from input.sharad.data.edr.pds import EDR as EDR_PDS
from input.sharad.data.edr.cosharps import EDR as EDR_COSHARPS
from input.sharad.data.rdr.pds import RDR as RDR_PDS
from input.sharad.data.rdr.pds_us import RDR as RDR_PDS_US
from geometry.coordinate_converter import CoordinateConverter
from dem.mola.mola128 import MOLA128
from simulators.drssim import dRSsim
from io_utils import envi


def main(argv=None):
    """
    This example:
    1) defines a coordinate converter based on Mars ellipsoid characteristics
    2) loads MOLA128 radius data
    3) defines a dummy simulator using the pre-defined coordinate converter and
       the MOLA128 DEM
    4)
        a) loads the ancillary part of a PDS Italian EDR product
        b) simulates the radargram using the orbit data (one frame every 26)
        c) saves the simulation and the image representing uncertain frames on disk
    5) loads the ancillary part of a CO-SHARPS raw product
        a) interpolates the geometric data in order to get the orbit information
           of one frame every 26
        b) simulates the radargram using the orbit data (one frame every 26)
        c) saves the simulation and the image representing uncertain frames on disk
    6) repeats the same steps using RDR PDS Italian data, also saving the RDR focused
       data on disk
    7) repeats the same steps using RDR PDS US data, also saving the RDR focused
       data on disk
    """

    lg = logger(verbose=True)

    # Mars ellipsoid characteristics
    cc = CoordinateConverter("Mars IAU 2000", 3396190., 3376200.)

    # the MOLA128 radius data are contained in a single ENVI file
    m = MOLA128(filename_base="/yourpath/r128")        # MOLA radius data
    m.read()

    # definition of a dummy RS simulator using the Mars coordinate converter and the MOLA128 DEM
    s = dRSsim(geom_obj=cc, dem_obj=m, n_processes=6)       # simulation parameters are contained in the dRSsim class

    # change fn_base with your OBS folder path+base file name (without suffixes)
    # EDR PDS
    fn_base = "/yourpath/e_1260201_001_ss19_700_a"
    d = EDR_PDS(dataset_name="edr1260201", filename_base=fn_base, logger=lg)
    d.load(mode="ancillary")
    d.generate_orbit_data(skip=26)
    od_pds_edr = d.orbit_data
    if od_pds_edr is not None:
        sim_image, uncert_image = s.simulate(od_pds_edr)
        envi.write("/yourpath/test_pds-edr_sim", sim_image)
        envi.write("/yourpath/test_pds-edr_sim_uncert", uncert_image)

    # EDR CO-SHARPS
    fn_base = "/yourpath/OBS_1260201000_1"
    d = EDR_COSHARPS(dataset_name="OBS_1260201000_1", filename_base=fn_base, logger=lg)
    d.load(mode="ancillary")
    d.generate_orbit_data(skip=26)
    od_cosharps_edr = d.orbit_data
    if od_cosharps_edr is not None:
        sim_image, uncert_image = s.simulate(od_cosharps_edr)
        envi.write("/yourpath/test_cosharps-edr_sim", sim_image)
        envi.write("/yourpath/test_cosharps-edr_sim_uncert", uncert_image)

    # RDR PDS
    fn_base = "/yourpath/r_1260201_001_ss19_700_a"
    d = RDR_PDS(dataset_name="rdr1260201", filename_base=fn_base, logger=lg)
    d.load()                        # note 1: d.data is np.complex64 for PDS RDR data
                                    # note 2: no "mode" is provided, as all data are contained in one unique file
    envi.write("/yourpath/test_pds_rdr_data", d.data)      # data is automatically converted into power. If complex output is needed, use the complex_output=True flag
    d.generate_orbit_data(skip=26)
    od_pds_rdr = d.orbit_data
    if od_pds_rdr is not None:
        sim_image, uncert_image = s.simulate(od_pds_rdr)
        envi.write("/yourpath/test_pds-rdr_sim", sim_image)
        envi.write("/yourpath/test_pds-rdr_sim-uncert", uncert_image)

    # RDR PDS US
    fn_base = "/yourpath/s_01260201"
    d = RDR_PDS_US(dataset_name="s_01260201", filename_base=fn_base, logger=lg)
    d.load(mode="full")
    envi.write("/yourpath/test_pds-us-rdr_data", d.data)

    d.generate_orbit_data()
    od_pds_us_rdr = d.orbit_data
    if od_pds_us_rdr is not None:
        sim_image, uncert_image = s.simulate(od_pds_us_rdr)
        envi.write("/yourpath/test_pds-us-rdr_sim", sim_image)
        envi.write("/yourpath/test_pds-us-rdr_sim_uncert", uncert_image)


if __name__ == "__main__":
    sys.exit(main())

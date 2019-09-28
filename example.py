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
from input.sharad.data.rdr.pds_us import RDR as RDR_PDS_US
from geometry.coordinate_converter import CoordinateConverter
from dem.mola.mola128 import MOLA128
from simulators.drssim import dRSsim
from io_utils import envi


def main(argv=None):
    """
    This example:
    - defines a coordinate converter based on Mars ellipsoid characteristics
    - loads MOLA128 radius data
    - defines a dummy simulator using the pre-defined coordinate converter and
      the MOLA128 DEM
    - loads the ancillary part of a CO-SHARPS EDR product
    - interpolates the geometric data in order to get the orbit information
      of one frame every 26
    - simulates the orbit
    - saves the simulation and the image representing uncertain frames on disk
    - repeats the same steps using RDR PDS US data, also saving the RDR focused
      data on disk
    """

    lg = logger(verbose=True)

    # Mars ellipsoid characteristics
    cc = CoordinateConverter("Mars IAU 2000", 3396190., 3376200.)

    # the MOLA128 radius data are contained in a single ENVI file
    m = MOLA128(filename_base="/yourpath/r128.hdr")        #mola_mars_data.hdr
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

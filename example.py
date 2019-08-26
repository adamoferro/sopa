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
from input.sharad.data.edr.cosharps import EDR
from geometry.coordinate_converter import CoordinateConverter
from dem.mola.mola128 import MOLA128
from simulators.drssim import DRSSim
import spectral.io.envi as envi                 # used for read/write of ENVI files


def main(argv=None):
    """
    This example:
    - loads the ancillary part of a COSHARPS EDR product;
    - interpolates the geometric data in order to get the orbit information
      of one frame every 26
    - loads MOLA128 data
    - simulates the orbit
    - saves the simulation and the image representing uncertain frames on disk
    """

    lg = logger(verbose=True)

    # change fn_base with your OBS folder path+base file name (without suffixes)
    fn_base = "/yourpath/OBS_1260201000_1/OBS_1260201000_1"
    d = EDR(dataset_name="OBS_1260201000_1", filename_base=fn_base, logger=lg)
    d.load(mode="ancillary")
    od = d.generate_orbit_data(skip=26)

    if od is not None:

        # # check interpolation
        # from matplotlib import pyplot as pp
        # pp.plot(od["time_offset_s"]+d._geom_data["timestamp"][0].astype("float64")/1000, od["Radius"])
        # pp.plot(d._geom_data["timestamp"].astype("float64")/1000, d._geom_data["Radius"]*1000, "r+", markersize=10)
        # pp.show()
        # exit()

        # Mars ellipsoid characteristics
        cc = CoordinateConverter("Mars IAU 2000", 3396190., 3376200.)

        # the MOLA128 data are contained in a single ENVI file
        d = MOLA128(filename_base="/yourpath/mola_mars_data.hdr")
        d.read()

        s = DRSSim(geom_obj=cc, dem_obj=d, n_processes=6)       # simulation parameters are contained in the DRSSim class
        sim_image, uncert_image = s.simulate(od)

        envi.save_image("/yourpath/test.hdr", sim_image, force=True)
        envi.save_image("/yourpath/test_uncert.hdr", uncert_image, force=True)


if __name__ == "__main__":
    sys.exit(main())

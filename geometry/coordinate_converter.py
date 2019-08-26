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
import math


class CoordinateConverter():
    def __init__(self, ellipsoid_name="", eer=0, epr=0):
        self._ELLIPSOID_EQUATOR_RADIUS = eer
        self._ELLIPSOID_POLAR_RADIUS = epr
        self._FLATTENING = (self._ELLIPSOID_EQUATOR_RADIUS - self._ELLIPSOID_POLAR_RADIUS) / self._ELLIPSOID_EQUATOR_RADIUS
        self._ECCENTRICITY_SQUARED = (eer**2 - epr**2) / eer**2
        self._deg2rad = math.pi / 180
        self._rad2deg = 180.0 / math.pi

    def euclidean_distance_cart(self, x1, y1, z1, x2, y2, z2):
        return np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

    def to_radius_from_xyz(self, xcart, ycart, zcart):
        return np.sqrt(xcart * xcart + ycart * ycart + zcart * zcart)

    def to_xyz_from_latlon_and_radius(self, latitude, longitude, radius):
        latitude = latitude * math.pi / 180.
        longitude_w = (longitude - 180.) * math.pi / 180.
        xcart = radius / np.sqrt((1 + np.tan(longitude_w)**2) * (1 + np.tan(latitude)**2))
        xcart[np.abs(longitude) <= 270] *= -1
        xcart[np.abs(longitude) <= 90] *= -1
        ycart = xcart * np.tan(longitude_w)
        zcart = np.sqrt(xcart**2 + ycart**2) * np.tan(latitude)
        return xcart, ycart, zcart

    def get_ellipsoid_radius_from_latlon(self, latitude, longitude):
        latitude = latitude * math.pi / 180.
        ellipsoid_local_radius = np.sqrt(self._ELLIPSOID_EQUATOR_RADIUS**2 / (1 + (1. / ((1 - self._FLATTENING)**2) - 1)*(np.sin(latitude)**2)))
        return ellipsoid_local_radius

    def get_ellipsoid_xyz_from_latlon(self, latitude, longitude):
        radius = self.get_ellipsoid_radius_from_latlong(latitude, longitude)
        latitude = latitude * math.pi / 180.
        longitude_w = (longitude - 180.) * math.pi / 180.
        xcart, ycart, zcart = self.to_xyz_from_latlon_and_radius(latitude, longitude_w, radius)
        return xcart, ycart, zcart

    # by David Wolever @ http://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python
    # extended for matrix use by Adamo Ferro
    def _unit_vector(self, vector):
        """ Returns the unit vector of the vector.  """
        return vector / np.linalg.norm(vector, axis=0)

    def _angle_between(self, v1, v2):
        """ Returns the angle in radians between vectors 'v1' and 'v2'::

                >>> angle_between((1, 0, 0), (0, 1, 0))
                1.5707963267948966
                >>> angle_between((1, 0, 0), (1, 0, 0))
                0.0
                >>> angle_between((1, 0, 0), (-1, 0, 0))
                3.141592653589793

        """
        v1_u = self._unit_vector(v1)
        v2_u = self._unit_vector(v2)

        if len(np.array(v2).shape) <= 1:
            angle = np.arccos(np.dot(v1_u, v2_u))
            if np.isnan(angle):
                if (v1_u == v2_u).all():
                    return 0.0
                else:
                    return np.pi
        else:
            angle = np.zeros([v2.shape[1], v2.shape[2]])
            for i_y in np.arange(angle.shape[0]):
                for i_x in np.arange(angle.shape[1]):
                    angle[i_y, i_x] = np.dot(v1_u[:, i_y, i_x], v2_u[:, i_y, i_x])
            angle = np.arccos(angle)
            nans = np.isnan(angle)
            angle[nans] = (((v1_u[:, nans] == v2_u[:,nans]).all()) * np.pi)
        return angle

    def LL_to_UTM(self, LatInput, LongInput):
        """
        converts lat/long to UTM coords.  Equations from USGS Bulletin 1532
        East Longitudes are positive, West longitudes are negative.
        North latitudes are positive, South latitudes are negative
        Lat and Long are in decimal degrees
        Written by Chuck Gantz- chuck.gantz@globalstar.com
        Ported to python3+numpy by Adamo Ferro
        """

        if len(LatInput.shape) == 0:
            Lat = np.array([LatInput])
            Long = np.array([LongInput])
        else:
            Lat = LatInput
            Long = LongInput

        a = self._ELLIPSOID_EQUATOR_RADIUS
        eccSquared = self._ECCENTRICITY_SQUARED
        k0 = 0.9996

        LongOrigin=0
        eccPrimeSquared=0
        N=0
        T=0
        C=0
        A=0
        M=0

    # Make sure the longitude is between -180.00 .. 179.9
        LongTemp = (Long+180)-((Long+180)/360).astype("int")*360-180 #-180.00 .. 179.9;

        LatRad = Lat * self._deg2rad
        LongRad = LongTemp * self._deg2rad
        LongOriginRad = 0
        ZoneNumber = 0

        ZoneNumber = ((LongTemp + 180)/6).astype("int") + 1

        #
        # if( Lat >= 56.0 && Lat < 64.0 && LongTemp >= 3.0 && LongTemp < 12.0 )
        # 		ZoneNumber = 32;
        #
        #  # Special zones for Svalbard
        # 	if( Lat >= 72.0 && Lat < 84.0 )
        # 	{
        # 	  if(      LongTemp >= 0.0  && LongTemp <  9.0 ) ZoneNumber = 31;
        # 	  else if( LongTemp >= 9.0  && LongTemp < 21.0 ) ZoneNumber = 33;
        # 	  else if( LongTemp >= 21.0 && LongTemp < 33.0 ) ZoneNumber = 35;
        # 	  else if( LongTemp >= 33.0 && LongTemp < 42.0 ) ZoneNumber = 37;
        # 	 }

        LongOrigin = (ZoneNumber - 1)*6 - 180 + 3  #+3 puts origin in middle of zone
        LongOriginRad = LongOrigin * self._deg2rad;

        #compute the UTM Zone from the latitude and longitude
        UTMLetters = self.UTMLetterDesignator(Lat)

        UTMZone = [str(ZoneNumber[i_z])+UTMLetters[i_z] for i_z in np.arange(len(ZoneNumber))]

        eccPrimeSquared = (eccSquared)/(1-eccSquared)

        N = a/np.sqrt(1-eccSquared*np.sin(LatRad)*np.sin(LatRad))
        T = np.tan(LatRad)*np.tan(LatRad)
        C = eccPrimeSquared*np.cos(LatRad)*np.cos(LatRad)
        A = np.cos(LatRad)*(LongRad-LongOriginRad)

        M = a*((1	- eccSquared/4		- 3*eccSquared*eccSquared/64	- 5*eccSquared*eccSquared*eccSquared/256)*LatRad
    				- (3*eccSquared/8	+ 3*eccSquared*eccSquared/32	+ 45*eccSquared*eccSquared*eccSquared/1024)*np.sin(2*LatRad)
    									+ (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024)*np.sin(4*LatRad)
    									- (35*eccSquared*eccSquared*eccSquared/3072)*np.sin(6*LatRad))

        UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
    					+ (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
    					+ 500000.0)

        UTMNorthing = (k0*(M+N*np.tan(LatRad)*(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
    				 + (61-58*T+T*T+600*C-330*eccPrimeSquared)*A*A*A*A*A*A/720)));

        UTMNorthing[Lat < 0] += 10000000.0 #10000000 meter offset for southern hemisphere
    #    if Lat < 0:
    #        UTMNorthing += 10000000.0; #10000000 meter offset for southern hemisphere

        return (UTMNorthing, UTMEasting, UTMZone, ZoneNumber)

    def LL_to_UTM_with_given_zone_number(self, LatInput, LongInput, ZoneNumber):

        # same as LLtoUTM, but the zone is given by the user
        if len(LatInput.shape)==0:
            Lat=np.array([LatInput])
            Long=np.array([LongInput])
        else:
            Lat=LatInput
            Long=LongInput

        a = self._ELLIPSOID_EQUATOR_RADIUS
        eccSquared = self._ECCENTRICITY_SQUARED
        k0 = 0.9996

        LongOrigin=0
        eccPrimeSquared=0
        N=0
        T=0
        C=0
        A=0
        M=0

    #Make sure the longitude is between -180.00 .. 179.9
        LongTemp = (Long+180)-((Long+180)/360).astype("int")*360-180; # -180.00 .. 179.9;
        LongOriginRad=0
        LongOrigin = (ZoneNumber - 1)*6 - 180 + 3  #+3 puts origin in middle of zone
        LongOriginRad = LongOrigin * self._deg2rad

        LatRad = Lat*self._deg2rad


        LongTemp[(LongTemp-LongOrigin)>180]-=360
        LongTemp[(LongTemp-LongOrigin)<-180]+=360

        LongRad = LongTemp*self._deg2rad;

        eccPrimeSquared = (eccSquared)/(1-eccSquared);

        N = a/np.sqrt(1-eccSquared*np.sin(LatRad)*np.sin(LatRad));
        T = np.tan(LatRad)*np.tan(LatRad);
        C = eccPrimeSquared*np.cos(LatRad)*np.cos(LatRad);
        A = np.cos(LatRad)*(LongRad-LongOriginRad);

        M = a*((1	- eccSquared/4		- 3*eccSquared*eccSquared/64	- 5*eccSquared*eccSquared*eccSquared/256)*LatRad
                                - (3*eccSquared/8	+ 3*eccSquared*eccSquared/32	+ 45*eccSquared*eccSquared*eccSquared/1024)*np.sin(2*LatRad)
                                                                        + (15*eccSquared*eccSquared/256 + 45*eccSquared*eccSquared*eccSquared/1024)*np.sin(4*LatRad)
                                                                        - (35*eccSquared*eccSquared*eccSquared/3072)*np.sin(6*LatRad));

        UTMEasting = (k0*N*(A+(1-T+C)*A*A*A/6
                                        + (5-18*T+T*T+72*C-58*eccPrimeSquared)*A*A*A*A*A/120)
                                        + 500000.0);

        UTMNorthing = (k0*(M+N*np.tan(LatRad)*(A*A/2+(5-T+9*C+4*C*C)*A*A*A*A/24
                                 + (61-58*T+T*T+600*C-330*eccPrimeSquared)*A*A*A*A*A*A/720)));

        UTMNorthing[Lat<0]+= 10000000.0 #10000000 meter offset for southern hemisphere

        return (UTMNorthing, UTMEasting)

    def UTMLetterDesignator(self, LatDataInput):
        """
        This routine determines the correct UTM letter designator for the given latitude
        returns 'Z' if latitude is outside the UTM limits of 84N to 80S
        Written by Chuck Gantz- chuck.gantz@globalstar.com
        Ported to python3+numpy by Adamo Ferro
        """

        Letter_designator_list=list()

        if len(LatDataInput.shape)==0:
            LatData=np.array([LatDataInput])
        else:
            LatData=LatDataInput


        for Lat in LatData:
            LetterDesignator=""
    #	if((84 >= Lat) && (Lat >= 72)) LetterDesignator = 'X';
            if Lat >= 72:
                LetterDesignator = 'X'			#*** modified ***
            elif 72 > Lat and Lat >= 64:
                LetterDesignator = 'W'
            elif 64 > Lat and Lat >= 56:
                LetterDesignator = 'V'
            elif 56 > Lat and Lat >= 48:
                LetterDesignator = 'U'
            elif 48 > Lat and Lat >= 40:
                LetterDesignator = 'T'
            elif 40 > Lat and Lat >= 32:
                LetterDesignator = 'S'
            elif 32 > Lat and Lat >= 24:
                LetterDesignator = 'R'
            elif 24 > Lat and Lat >= 16:
                LetterDesignator = 'Q'
            elif 16 > Lat and Lat >= 8:
                LetterDesignator = 'P'
            elif 8 > Lat and Lat >= 0:
                LetterDesignator = 'N'
            elif 0 > Lat and Lat >= -8:
                LetterDesignator = 'M'
            elif -8> Lat and Lat >= -16:
                LetterDesignator = 'L'
            elif -16 > Lat and Lat >= -24:
                LetterDesignator = 'K'
            elif -24 > Lat and Lat >= -32:
                LetterDesignator = 'J'
            elif -32 > Lat and Lat >= -40:
                LetterDesignator = 'H'
            elif -40 > Lat and Lat >= -48:
                LetterDesignator = 'G'
            elif -48 > Lat and Lat >= -56:
                LetterDesignator = 'F'
            elif -56 > Lat and Lat >= -64:
                LetterDesignator = 'E'
            elif -64 > Lat and Lat >= -72:
                LetterDesignator = 'D'
        #	else if((-72 > Lat) && (Lat >= -80)) LetterDesignator = 'C';	#*** modified ***
            elif -72 > Lat:
                LetterDesignator = 'C'
        #	else LetterDesignator = 'Z'; //This is here as an error flag to show that the Latitude is outside the UTM limits
            Letter_designator_list.append(LetterDesignator)
        return Letter_designator_list

    def UTM_to_LL(self, UTMNorthing, UTMEasting, fixed_UTMZone):
        """
        converts UTM coords to lat/long.  Equations from USGS Bulletin 1532
        East Longitudes are positive, West longitudes are negative.
        North latitudes are positive, South latitudes are negative
        Lat and Long are in decimal degrees.
        Written by Chuck Gantz- chuck.gantz@globalstar.com

        Ported to python3 by Adamo Ferro
        Modified in order to return only positive Longitudes between 0-360 degs
        """

        k0 = 0.9996
        a = self._ELLIPSOID_EQUATOR_RADIUS
        eccSquared = self._ECCENTRICITY_SQUARED
        eccPrimeSquared=0
        e1 = (1-np.sqrt(1-eccSquared))/(1+np.sqrt(1-eccSquared))
        N1=0
        T1=0
        C1=0
        R1=0
        D=0
        M=0
        LongOrigin=0
        mu=0
        phi1=0
        phi1Rad=0
        x=0
        y=0
        ZoneNumber=0
        ZoneLetter=""
        NorthernHemisphere=0 #1 for northern hemispher, 0 for southern

        x = UTMEasting - 500000.0 #remove 500,000 meter offset for longitude
        y = UTMNorthing

        ZoneNumber = int(fixed_UTMZone[:-1])
        ZoneLetter=fixed_UTMZone[-1]
        if ord(ZoneLetter) - ord('N') >= 0:
            NorthernHemisphere = 1 #point is in northern hemisphere
        else:
            NorthernHemisphere = 0 #point is in southern hemisphere
            y -= 10000000.0 #remove 10,000,000 meter offset used for southern hemisphere

        LongOrigin = (ZoneNumber - 1)*6 - 180 + 3   #+3 puts origin in middle of zone

        eccPrimeSquared = (eccSquared)/(1-eccSquared)

        M = y / k0
        mu = M/(a*(1-eccSquared/4-3*eccSquared*eccSquared/64-5*eccSquared*eccSquared*eccSquared/256))

        phi1Rad = mu+ (3*e1/2-27*e1*e1*e1/32)*np.sin(2*mu) + (21*e1*e1/16-55*e1*e1*e1*e1/32)*np.sin(4*mu)+(151*e1*e1*e1/96)*np.sin(6*mu)
        phi1 = phi1Rad*self._rad2deg

        N1 = a/np.sqrt(1-eccSquared*np.sin(phi1Rad)*np.sin(phi1Rad))
        T1 = np.tan(phi1Rad)*np.tan(phi1Rad)
        C1 = eccPrimeSquared*np.cos(phi1Rad)*np.cos(phi1Rad)
        R1 = a*(1-eccSquared)/(1-eccSquared*np.sin(phi1Rad)*np.sin(phi1Rad))**(1.5)
        D = x/(N1*k0)

        Lat = phi1Rad - (N1*np.tan(phi1Rad)/R1)*(D*D/2-(5+3*T1+10*C1-4*C1*C1-9*eccPrimeSquared)*D*D*D*D/24
                                        +(61+90*T1+298*C1+45*T1*T1-252*eccPrimeSquared-3*C1*C1)*D*D*D*D*D*D/720)
        Lat = Lat * self._rad2deg

        Long = (D-(1+2*T1+C1)*D*D*D/6+(5-2*C1+28*T1-3*C1*C1+8*eccPrimeSquared+24*T1*T1)
                                        *D*D*D*D*D/120)/np.cos(phi1Rad);
        Long = LongOrigin + Long * self._rad2deg;

        if len(np.array(Long).shape)>0:
            Long[Long<0]=360+Long[Long<0]
        elif Long<0:
            Long=360+Long

        return (Lat, Long)

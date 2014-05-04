from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from datetime import *
from sat_math_utils import *

def main():
    # A random TLE I found
    line1 = '1 25544U 98067A   13042.10021568  .00009549  00000-0  16043-3 0  1955'
    line2 = '2 25544 051.6484 014.6335 0012183 260.9055 192.8707 15.52430028815125'

    # Initialize the satellite propagation model
    satellite = twoline2rv(line1, line2, wgs72)

    now = datetime.today()

    # get current position and velocity in ECI
    position_eci, velocity_eci = satellite.propagate(now.year, now.month, now.day, now.hour, now.minute, now.second)

    # Using GMST, transform ECI to various other frames
    gmst =  current_GMST(now)
    position_sp = spherical (position_eci)
    velocity_sp = spherical (velocity_eci)

    path_points_eci = []
    path_points_eci.append( [ position_eci[0],
                                position_eci[1],
                                position_eci[2] ])
    path_points_sp = []
    path_points_sp.append( [ position_sp[0],
                             position_sp[1],
                             position_sp[2] ] )
    print path_points_sp
    print path_points_eci

if __name__ == '__main__':
    main()

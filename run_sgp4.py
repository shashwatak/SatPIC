from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4
from datetime import *
from sat_math_utils import *
import csv
from jdcal import *

NUM_SECONDS_IN_DAY = 86400

def main():
    tle_list = []
    with open('data/zarya_tle', 'rb') as tlefile:
        for _ in range(50):
            line1 = tlefile.readline()
            line2 = tlefile.readline()
            tle_list.append((line1, line2))

    # line1 = '1 25544U 98067A   13042.10021568  .00009549  00000-0  16043-3 0  1955'
    # line2 = '2 25544 051.6484 014.6335 0012183 260.9055 192.8707 15.52430028815125'

    phi_file = open('phi.csv', 'wb')
    theta_file = open('theta.csv', 'wb')
    phi_writer = csv.writer(phi_file)
    theta_writer = csv.writer(theta_file)

    for tle in tle_list:
        line1, line2 = tle

        # Initialize the satellite propagation model
        satellite = twoline2rv(line1, line2, wgs84)

        # get current position and velocity in ECI
        position_eci, velocity_eci = sgp4(satellite, 0)

        # transform rectangular coordinates to spherical coordinates
        position_sp = spherical(position_eci)
        velocity_sp = spherical(velocity_eci)

        year, month, day, time_frac = jd2gcal(MJD_0, satellite.jdsatepoch-MJD_0)

        seconds_of_day = time_frac * NUM_SECONDS_IN_DAY
        minutes_of_day = seconds_of_day / 60
        hours = minutes_of_day / 60
        minutes = minutes_of_day % 60
        seconds = seconds_of_day % 60

        d = datetime(year, month, day, int(hours), int(minutes), int(seconds))
        time_string = '{:%Y-%m-%d %H:%M:%S}'.format(d)

        phi_writer.writerow([time_string, position_sp[0], '0'])
        theta_writer.writerow([time_string, position_sp[1], '0'])

if __name__ == '__main__':
    main()

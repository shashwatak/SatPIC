from sgp4.earth_gravity import wgs84
from sgp4.io import twoline2rv
from sgp4.propagation import sgp4
from datetime import *
from sat_math_utils import *
import csv
from jdcal import *

NUM_SECONDS_IN_DAY = 86400
TIME_TEMPLATE = '{:%Y-%m-%d %H:%M:%S}'

phi_file = open('phi.csv', 'wb')
theta_file = open('theta.csv', 'wb')
phi_csv_writer = csv.writer(phi_file)
theta_csv_writer = csv.writer(theta_file)

def main():
    tle_list = []
    with open('data/zarya_tle', 'rb') as tlefile:
        for _ in range(50):
            line1 = tlefile.readline()
            line2 = tlefile.readline()
            tle_list.append((line1, line2))

    phi_file.write('dttm,value,label\n')
    theta_file.write('dttm,value,label\n')

    prev_sat = None
    prev_tle = None
    prev_date = None

    for tle in tle_list:
        line1, line2 = tle

        # Initialize the satellite propagation model
        satellite = twoline2rv(line1, line2, wgs84)
        date = date_from_jds(satellite.jdsatepoch)

        if prev_tle is None:
            prev_tle = tle
            prev_sat = satellite
            prev_date = date
            continue

        time_delta_minutes = (date - prev_date).total_seconds()/60
        time_delta_half = time_delta_minutes / 2

        # First we catch the previous satellite up to the midway point
        propogate_over_range(prev_sat, prev_date, 0, int(time_delta_half))

        # Then we run the current satellite backwards to th halfway point
        propogate_over_range(satellite, date, -int(time_delta_half), 0)

        prev_tle = tle
        prev_sat = satellite
        prev_date = date

    phi_file.close()
    theta_file.close()

def propogate_over_range(satellite, date, min_time, max_time):
    # First we catch the previous satellite up to the midway point
    for tsince in range(min_time, max_time):
        # get current position and velocity in ECI
        position_eci, velocity_eci = sgp4(satellite, tsince)

        interval_date = increment_date_by_minutes(date, tsince)
        time_string = TIME_TEMPLATE.format(interval_date)

        # transform rectangular coordinates to spherical coordinates
        position_sp = spherical(position_eci)
        velocity_sp = spherical(velocity_eci)

        phi_csv_writer.writerow([time_string, position_sp[0], '0'])
        theta_csv_writer.writerow([time_string, position_sp[1], '0'])

def date_from_jds(jds):
    year, month, day, time_frac = jd2gcal(MJD_0, jds-MJD_0)
    seconds_of_day = time_frac * NUM_SECONDS_IN_DAY
    minutes_of_day = seconds_of_day / 60
    hours = minutes_of_day / 60
    minutes = minutes_of_day % 60
    seconds = seconds_of_day % 60
    date = datetime(year, month, day, int(hours), int(minutes), int(seconds))
    return date

def increment_date_by_minutes(date, minutes):
    time_difference = timedelta(minutes=minutes)
    return date + time_difference



if __name__ == '__main__':
    main()

import numpy as np
from sgp4.api import Satrec, jday
from sgp4.earth_gravity import wgs72

def get_satellite_position_velocity(satellite, date_time, val):
    jd, fr = jday(date_time.year, date_time.month, date_time.day, 
                  date_time.hour, date_time.minute, date_time.second)
    e, r, v = satellite.sgp4(jd, fr)
    if e != 0:
        raise RuntimeError(f'Error propagating satellite: {e}')
    return np.array(r), np.array(v)

def check_collision(pos1, pos2, threshold=1.0):  # threshold in kilometers
    distance = np.linalg.norm(pos1 - pos2)
    return distance < threshold

def check_specific_satellites_collision(sat_tle, all_sats, start_time, end_time, time_step):
    positions = []
    collision_positions = []
    current_time = start_time
    while current_time <= end_time:
        sat = Satrec.twoline2rv(sat_tle['TLE_LINE1'], sat_tle['TLE_LINE2'])
        pos1, vel1 = get_satellite_position_velocity(sat, current_time, True)
        positions.append(pos1)
        for other_sat_id, tle in all_sats.items():
            other_sat = Satrec.twoline2rv(tle['TLE_LINE1'], tle['TLE_LINE2'])
            pos2, vel2 = get_satellite_position_velocity(other_sat, current_time, False)
            if check_collision(pos1, pos2):
                collision_positions.append((current_time, other_sat_id, pos1, pos2))
        current_time += time_step
    return positions, collision_positions

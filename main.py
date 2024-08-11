import helpers
import orbital as ob
import json
import sgp4_function as sg
from datetime import datetime, timedelta, timezone



def individual_satellites():
    # Example orbital elements
    epoch = datetime(2022, 1, 31, 12, 0, 0)
    sat1 = ob.OrbitalElements(a=7000e3, e=0.0012345, i=98.7654, raan=123.4567, arg_pe=45.6789, true_anomaly=321.6543, mean_motion=15.12345678, epoch = epoch)
    sat2 = ob.OrbitalElements(a=7500e3, e=0.0032345, i=94.7654, raan=129.4567, arg_pe=49.6789, true_anomaly=327.6543, mean_motion=11.12345678, epoch = epoch)
    sat1.mean_anomaly = 0
    sat2.mean_anomaly = 0

    # Simulation parameters
    time_span = 86400  # One day in seconds
    time_step = 60  # One minute time steps
    collision_threshold = 100  # Collision threshold in meters

    # Simulate orbits and check for collisions
    satellites = [sat1, sat2]
    all_positions, collisions = helpers.simulate_orbits(satellites, time_span, time_step, collision_threshold)
    if collisions:
        if len(collisions) <= 9:
            for collision in collisions:
                print(f"Potential collision between satellite {collision[0]} and satellite {collision[1]} at step {collision[2]}")
        else:
            for val in range(4):
                print(f"Potential collision between satellite {collisions[val][0]} and satellite {collisions[val][1]} at step {collisions[val][2]}")
            print("...")
            for val in range(len(collisions) - 4, len(collisions)):
                print(f"Potential collision between satellite {collisions[val][0]} and satellite {collisions[val][1]} at step {collisions[val][2]}")
    else:
        print("No collisions")
    print(all_positions)


    # Plotting the orbit with Earth
    helpers.plot_orbit_with_earth(all_positions)

def total_sats():
    # Load the data from the JSON file
    with open('satellites.json', 'r') as json_file:
        tle_data = json.load(json_file)

    # Create the dictionary
    tle_dict = {item['OBJECT_ID']: {'TLE_LINE1': item['TLE_LINE1'], 'TLE_LINE2': item['TLE_LINE2']} for item in tle_data}
    epoch = datetime(2022, 1, 31, 12, 0, 0)
    sat1 = ob.OrbitalElements(a=7000e3, e=0.0012345, i=98.7654, raan=123.4567, arg_pe=45.6789, true_anomaly=321.6543, mean_motion=15.12345678, epoch = epoch)
    sat1_tle = sat1.to_tle()
    epoch_start = datetime(2024, 8, 6, 18, 35, 0)
    epoch_end = datetime(2024, 8, 7, 18, 35, 0)
    time_step = timedelta(seconds=100)
    sat_positions, collision_positions = sg.check_specific_satellites_collision(sat1_tle, tle_dict, epoch_start, epoch_end, time_step)

    
    helpers.plot_points_and_orbit(tle_dict, sat_positions)

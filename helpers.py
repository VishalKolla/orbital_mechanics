import numpy as np
import matplotlib.pyplot as plt
import orbital as ob
import sgp4_function as sg
from sgp4.api import Satrec, jday
import datetime

MU = 3.986e14  # Gravitational parameter for Earth (m^3/s^2)
R_EARTH = 6371e3  # Radius of Earth (m)

def solve_kepler(M, e, tol=1e-6):
    E = M  # Initial guess: mean anomaly
    while True:
        delta_E = (M - (E - e * np.sin(E))) / (1 - e * np.cos(E))
        E += delta_E
        if abs(delta_E) < tol:
            break
    return E

def true_anomaly(E, e):
    return 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

def simulate_step(satellites, time_step):
    positions = []
    for sat in satellites:
        position, _ = sat.calculate_position()
        positions.append(position)
        # Update mean anomaly for the next step
        sat.mean_anomaly += np.sqrt(MU / sat.a**3) * time_step
        sat.mean_anomaly = sat.mean_anomaly % (2 * np.pi)  # Normalize to 0-2pi
    return np.array(positions)

def simulate_orbits(satellites, time_span, time_step, collision_threshold):
    num_steps = int(time_span / time_step)
    all_positions = {i: [] for i in range(len(satellites))}
    collisions = []
    
    for step in range(num_steps):
        positions = simulate_step(satellites, time_step)
        for i, pos in enumerate(positions):
            all_positions[i].append(pos)
        collision = check_for_collisions(positions, collision_threshold)
        if collision:
            collision.append(step)
            collisions.append(collision)

    for i in all_positions:
        all_positions[i] = np.array(all_positions[i])
    
    return all_positions, collisions

def calculate_distances(positions1, positions2):
    return np.linalg.norm(positions1 - positions2, axis=1)

def check_for_collisions(positions, collision_threshold):
    num_sats = len(positions)
    
    for i in range(num_sats):
        for j in range(i + 1, num_sats):
            distance = np.linalg.norm(positions[i] - positions[j])
            if distance < collision_threshold:
                return [i+1, j+1]  # Return True if a collision is detected
    return False  # Return False if no collisions are detected

def plot_orbit_with_earth(positions_dict):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot Earth
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = R_EARTH / 1000 * np.outer(np.cos(u), np.sin(v))  # Convert to km
    y = R_EARTH / 1000 * np.outer(np.sin(u), np.sin(v))  # Convert to km
    z = R_EARTH / 1000 * np.outer(np.ones(np.size(u)), np.cos(v))  # Convert to km
    
    ax.plot_surface(x, y, z, color='b', alpha=0.3)
    
    # Plot orbits
    for sat_id, positions in positions_dict.items():
        positions = np.array(positions)  # Ensure positions is a numpy array
        ax.plot(positions[:, 0] / 1000, positions[:, 1] / 1000, positions[:, 2] / 1000, label=f'Satellite {sat_id+1}')  # Convert to km
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    plt.show()

def plot_points_and_orbit(all_sats, input_positions):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    current_time = datetime.datetime.now()
    
    # Plot Earth
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = R_EARTH / 1000 * np.outer(np.cos(u), np.sin(v))  # Convert to km
    y = R_EARTH / 1000 * np.outer(np.sin(u), np.sin(v))  # Convert to km
    z = R_EARTH / 1000 * np.outer(np.ones(np.size(u)), np.cos(v))  # Convert to km
    ax.plot_surface(x, y, z, color='b', alpha=0.5)

    for other_sat_id, tle in all_sats.items():
        other_sat = Satrec.twoline2rv(tle['TLE_LINE1'], tle['TLE_LINE2'])
        position, _ = sg.get_satellite_position_velocity(other_sat, current_time, False)
        ax.scatter(position[0], position[1], position[2], color='g', s=7)

    x_line = [point[0] for point in input_positions]
    y_line = [point[1] for point in input_positions]
    z_line = [point[2] for point in input_positions]

    ax.plot(x_line, y_line, z_line, color="r", linewidth = 2, linestyle="dotted")
    
    
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    ax.legend()
    plt.show()

def parse_tle(line2, id):
    # Parse line 2
    inclination = float(line2[8:16])
    raan = float(line2[17:25])
    eccentricity = float('0.' + line2[26:33].strip())
    arg_perigee = float(line2[34:42])
    mean_anomaly = float(line2[43:51])
    mean_motion = float(line2[52:63])
    
    # Calculate semi-major axis from mean motion
    mean_motion_rad = mean_motion * 2 * np.pi / 86400  # Convert to radians per second
    a = (MU / (mean_motion_rad ** 2)) ** (1/3)  # Semi-major axis in meters
    
    # Create OrbitalElements object
    elements = ob.OrbitalElements(a, eccentricity, inclination, raan, arg_perigee, mean_anomaly, mean_motion, id)

    return elements
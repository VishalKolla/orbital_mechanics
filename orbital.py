import math
import numpy as np
from datetime import datetime, timedelta

# Constants
G = 6.67430e-11  # Gravitational constant (m^3 kg^-1 s^-2)
M_EARTH = 5.972e24  # Mass of Earth (kg)
MU = 3.986e14  # Gravitational parameter for Earth (m^3/s^2)
R_EARTH = 6371e3  # Radius of Earth (m)

class OrbitalElements:
    def __init__(self, a, e, i, raan, arg_pe, true_anomaly, mean_motion, epoch, mean_motion_dot=0, mean_motion_ddot="00000-0", b_star = "00000-0", revolution_number = 56353):
        self.a = a  # Semi-major axis in meters
        self.e = e  # Eccentricity
        self.i = np.radians(i)  # Inclination in radians
        self.raan = np.radians(raan)  # Right ascension of ascending node in radians
        self.arg_pe = np.radians(arg_pe)  # Argument of periapsis in radians
        self.true_anomaly = np.radians(true_anomaly)  # True anomaly in radians
        self.mean_motion = mean_motion  # Mean motion in revolutions per day
        self.mean_anomaly = self.calculate_mean_anomaly()
        self.epoch = epoch  # Epoch time as a datetime object
        self.mean_motion_dot = mean_motion_dot  # First time derivative of mean motion
        self.mean_motion_ddot = mean_motion_ddot  # Second time derivative of mean motion
        self.b_star = b_star
        self.revolution_number = revolution_number

    def calculate_mean_motion(self):
        # Calculate mean motion in radians per second
        n_rad_per_sec = math.sqrt(MU / self.a**3)

        # Convert to revolutions per day
        n_rev_per_day = (n_rad_per_sec * 86400) / (2 * math.pi)

        return n_rev_per_day

    def calculate_position(self):
        M = self.mean_anomaly
        E = solve_kepler(M, self.e)
        v = true_anomaly(E, self.e)
        
        # Calculate distance
        r = self.a * (1 - self.e * np.cos(E))
        
        # Position in the orbital plane
        x_orb = r * np.cos(v)
        y_orb = r * np.sin(v)

        cos_raan = np.cos(self.raan)
        sin_raan = np.sin(self.raan)
        cos_arg_pe_nu = np.cos(self.arg_pe + v)
        sin_arg_pe_nu = np.sin(self.arg_pe + v)
        cos_i = np.cos(self.i)
        sin_i = np.sin(self.i)
        
        # Rotate to 3D space using orbital elements
        x = (cos_raan * cos_arg_pe_nu - sin_raan * sin_arg_pe_nu * cos_i) * x_orb + (-cos_raan * sin_arg_pe_nu - sin_raan * cos_arg_pe_nu * cos_i) * y_orb
        y = (sin_raan * cos_arg_pe_nu + cos_raan * sin_arg_pe_nu * cos_i) * x_orb + (-sin_raan * sin_arg_pe_nu + cos_raan * cos_arg_pe_nu * cos_i) * y_orb
        z = (sin_i * sin_arg_pe_nu) * x_orb + (sin_i * cos_arg_pe_nu) * y_orb
            
        return np.array([x, y, z]), v
    
    def calculate_mean_anomaly(self):
        # Mean anomaly from true anomaly and eccentricity
        E = 2 * np.arctan(np.sqrt((1 - self.e) / (1 + self.e)) * np.tan(self.true_anomaly / 2))
        return E - self.e * np.sin(E)


    def format_tle_line1(self, catalog_number, classification, intl_designator, ephemeris_type, element_set):
        # Ensure catalog_number is an integer
        catalog_number = int(catalog_number)
        while len(self.mean_motion_ddot) < 8:
            self.mean_motion_ddot =  " " + self.mean_motion_ddot

        while len(self.b_star) < 8:
            self.b_star =  " " + self.b_star
        
        epoch_year = self.epoch.year % 100
        epoch_day_of_year = (self.epoch - datetime(self.epoch.year, 1, 1)).days + 1 + self.epoch.hour / 24 + self.epoch.minute / 1440 + self.epoch.second / 86400
        
        line1 = (f"1 {catalog_number:5d}{classification} {intl_designator:<8} "
                f"{epoch_year:02d}{epoch_day_of_year:012.8f} "
                f"{self.mean_motion_dot:8.8f} "
                + self.mean_motion_ddot + " " +
                self.b_star + " " +
                f"{ephemeris_type} "
                f"{element_set:3d}")
        
        # Calculate checksum
        checksum = calculate_checksum(line1)
        
        return line1 + str(checksum)

    def format_tle_line2(self, catalog_number, inclination, raan, eccentricity, arg_pe, mean_anomaly, mean_motion):
        # Ensure all values are correctly formatted
        catalog_number = int(catalog_number)
        eccentricity = float(eccentricity)  # Ensure eccentricity is a float

        
        # Convert eccentricity to TLE format (without decimal point)
        eccentricity_str = f"{eccentricity:.7f}".replace("0.", "")
        
        # Format TLE line 2
        line2 = (f"2 {catalog_number:5d} "
                f"{inclination:8.4f} "
                f"{raan:8.4f} " +
                eccentricity_str + " " +  # Eccentricity without decimal point
                f"{arg_pe:8.4f} "
                f"{mean_anomaly:8.4f} "
                f"{mean_motion:11.8f}" + 
                str(self.revolution_number))
        
        # Calculate checksum
        checksum = calculate_checksum(line2)
        
        return line2 + str(checksum)


    def to_tle(self, catalog_number="25544", classification="U", intl_designator="98067A", ephemeris_type=0, element_set=9999):
        line1 = self.format_tle_line1(catalog_number, classification, intl_designator, ephemeris_type, element_set)
        line2 = self.format_tle_line2(catalog_number, np.degrees(self.i), np.degrees(self.raan), self.e, np.degrees(self.arg_pe), np.degrees(self.mean_anomaly), self.mean_motion)
        tle = {
            'TLE_LINE1': line1,
            'TLE_LINE2': line2
        }
        return tle
    
def calculate_checksum(tle_line):
        checksum = 0
        for char in tle_line:
            if char.isdigit():
                checksum += int(char)
            elif char == '-':
                checksum += 1
        checksum %= 10
        return checksum

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

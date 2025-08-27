"""
Base class for magnetic field calculations and simulations.
Provides common functionality for processing measurement files and calculating fields.
"""
import os
import re
import math
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from .constants import (
    WAVE_NUMBER, PERMEABILITY_VACUUM, ANTENNA_SURFACE, 
    ANGULAR_FREQUENCY, DEFAULT_RESOLUTION
)
from .utils import (
    cartesian_to_spherical_coordinates,
    calculate_magnetic_field_components,
    calculate_current_from_radial_field,
    calculate_current_from_tangential_field,
    filter_points_excluding_origin,
    select_closest_points_by_axis
)


class MagneticFieldCalculator:
    """
    Base class for magnetic field calculations and analysis.
    
    This class provides common functionality for:
    - Processing measurement files
    - Calculating magnetic field components
    - Interpolating points
    - Visualizing results
    """
    
    def __init__(self, data_directory, resolution=DEFAULT_RESOLUTION):
        """
        Initialize the magnetic field calculator.
        
        Args:
            data_directory: Directory containing measurement files
            resolution: Resolution for interpolation operations
        """
        self.data_directory = data_directory
        self.resolution = resolution
        self.measurement_results = []
        self.high_resolution_points = []
        self.average_current = None
    
    def extract_amplitude_and_position_from_file(self, file_content, filename):
        """
        Extract amplitude measurements and position from file content.
        
        Args:
            file_content: Content of the measurement file
            filename: Name of the file for position extraction
            
        Returns:
            tuple: (amplitude_x, amplitude_y, amplitude_z, x, y, z)
        """
        # Extract amplitudes using regex
        matches = re.findall(r"C(\d+)_pkpk:(\d+\.\d+)", file_content)
        amplitudes = {f"C{axis}": float(amp) for axis, amp in matches}
        
        # Extract position from filename
        pos_match = re.search(
            r"logMeasureOscillo_\((-?\d+\.?\d*) (-?\d+\.?\d*) (-?\d+\.?\d*)\)", 
            filename
        )
        
        if pos_match:
            x, y, z = map(float, pos_match.groups())
        else:
            x = y = z = 0.0
        
        # Return amplitudes and position
        amplitude_x = amplitudes.get('C1', 0)
        amplitude_y = amplitudes.get('C2', 0)  # Usually 0 for single-channel measurements
        amplitude_z = amplitudes.get('C3', 0)  # Usually 0 for single-channel measurements
        
        return amplitude_x, amplitude_y, amplitude_z, x, y, z
    
    def process_measurement_files(self):
        """
        Process all measurement files in the data directory.
        
        Returns:
            list: List of measurement results with calculated field components
        """
        self.measurement_results = []
        
        for filename in os.listdir(self.data_directory):
            if filename.startswith("logMeasureOscillo"):
                filepath = os.path.join(self.data_directory, filename)
                
                with open(filepath, 'r') as file:
                    file_content = file.read()
                
                # Extract amplitudes and position
                amp_x, amp_y, amp_z, x, y, z = self.extract_amplitude_and_position_from_file(
                    file_content, filename
                )
                
                # Calculate magnetic field components
                (Bx, By, Bz), (Hx, Hy, Hz), H_magnitude, theta, phi = calculate_magnetic_field_components(
                    amp_x, amp_y, amp_z
                )
                
                # Store results
                self.measurement_results.append({
                    'filename': filename,
                    'x': x, 'y': y, 'z': z,
                    'Bx': Bx, 'By': By, 'Bz': Bz,
                    'Hx': Hx, 'Hy': Hy, 'Hz': Hz,
                    '|H|': H_magnitude,
                    'Htheta': theta,
                    'Hphi': phi
                })
        
        return self.measurement_results
    
    def select_nearest_points(self):
        """
        Select the nearest points to the origin along each axis.
        
        Returns:
            list: Selected points closest to origin
        """
        points_excluding_origin = filter_points_excluding_origin(self.measurement_results)
        
        selected_points = []
        for axis in ['x', 'y', 'z']:
            axis_points = select_closest_points_by_axis(points_excluding_origin, axis)
            selected_points.extend(axis_points)
        
        return selected_points
    
    def calculate_average_current(self, selected_points):
        """
        Calculate average current from selected measurement points.
        
        Args:
            selected_points: List of selected measurement points
            
        Returns:
            float: Average current value or None if calculation fails
        """
        current_values = []
        
        for point in selected_points:
            r, theta, _ = cartesian_to_spherical_coordinates(point['x'], point['y'], point['z'])
            H_radial = point['|H|']
            H_tangential = point['Htheta']
            
            # Store spherical coordinates in point
            point['r'] = r
            point['theta'] = math.degrees(theta)
            
            # Calculate current from both radial and tangential components
            current_radial = calculate_current_from_radial_field(H_radial, r, math.degrees(theta))
            current_tangential = calculate_current_from_tangential_field(H_tangential, r, math.degrees(theta))
            
            # Use average if both values exist, otherwise use the available one
            if current_radial is not None and current_tangential is not None:
                current_average = (current_radial + current_tangential) / 2
                current_values.append(current_average)
            elif current_radial is not None:
                current_values.append(current_radial)
            elif current_tangential is not None:
                current_values.append(current_tangential)
        
        # Calculate global average
        if current_values:
            self.average_current = sum(current_values) / len(current_values)
            return self.average_current
        else:
            return None
    
    def calculate_radial_field_component(self, x, y, z, current):
        """
        Calculate radial component of magnetic field for given position and current.
        
        Args:
            x, y, z: Position coordinates
            current: Current value
            
        Returns:
            float: Radial field component magnitude
        """
        r, theta, _ = cartesian_to_spherical_coordinates(x, y, z)
        
        if r == 0:
            return 0
        
        factor = (1j / (WAVE_NUMBER**2 * r**2)) + (1 / (WAVE_NUMBER**3 * r**3))
        H_radial = (current * ANTENNA_SURFACE * WAVE_NUMBER**3 / (2 * math.pi)) * \
                   factor * math.cos(theta) * np.exp(-1j * WAVE_NUMBER * r)
        
        return abs(H_radial)
    
    def calculate_tangential_field_component(self, x, y, z, current):
        """
        Calculate tangential component of magnetic field for given position and current.
        
        Args:
            x, y, z: Position coordinates
            current: Current value
            
        Returns:
            float: Tangential field component magnitude
        """
        r, theta, _ = cartesian_to_spherical_coordinates(x, y, z)
        
        if r == 0:
            return 0
        
        sin_theta = np.sin(theta)
        if np.isclose(sin_theta, 0):
            return 0
        
        factor = (-1 / (WAVE_NUMBER * r) + 1j / (WAVE_NUMBER**2 * r**2) + 1 / (WAVE_NUMBER**3 * r**3))
        H_tangential = (current * ANTENNA_SURFACE * WAVE_NUMBER**3 / (4 * math.pi)) * \
                       factor * sin_theta * np.exp(-1j * WAVE_NUMBER * r)
        
        return abs(H_tangential)
    
    def calculate_total_field_magnitude(self, x, y, z, current):
        """
        Calculate total magnetic field magnitude at given position.
        
        Args:
            x, y, z: Position coordinates
            current: Current value
            
        Returns:
            float: Total field magnitude
        """
        H_radial = self.calculate_radial_field_component(x, y, z, current)
        H_tangential = self.calculate_tangential_field_component(x, y, z, current)
        
        return math.sqrt(H_radial**2 + H_tangential**2)
    
    def visualize_3d_vectors(self, points=None):
        """
        Create 3D visualization of magnetic field vectors.
        
        Args:
            points: List of points to visualize (defaults to measurement results)
        """
        if points is None:
            points = self.measurement_results
        
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        for point in points:
            x, y, z = point['x'], point['y'], point['z']
            Hx, Hy, Hz = point['Hx'], point['Hy'], point['Hz']
            H_magnitude = point['|H|']
            
            # Plot vector arrows
            ax.quiver(x, y, z, Hx, Hy, Hz, length=H_magnitude/100, normalize=True, color='blue')
            
            # Plot position points
            ax.scatter(x, y, z, color='red', s=20)
        
        ax.set_xlabel('X Position (mm)')
        ax.set_ylabel('Y Position (mm)')
        ax.set_zlabel('Z Position (mm)')
        ax.set_title('Magnetic Field Vectors in 3D Space')
        
        plt.show()
    
    def execute_analysis_pipeline(self):
        """
        Execute the complete magnetic field analysis pipeline.
        
        Returns:
            dict: Analysis results including average current and processed points
        """
        # Process measurement files
        self.process_measurement_files()
        
        # Select nearest points for current calculation
        nearest_points = self.select_nearest_points()
        
        # Calculate average current
        average_current = self.calculate_average_current(nearest_points)
        
        # Return results
        return {
            'measurement_results': self.measurement_results,
            'nearest_points': nearest_points,
            'average_current': average_current
        }
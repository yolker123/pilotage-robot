"""
Magnetic field interpolation and resolution enhancement functionality.
Extends the base magnetic field calculator with interpolation capabilities.
"""
import itertools
import math
import numpy as np

from .magnetic_field_base import MagneticFieldCalculator
from .utils import interpolate_linear, cartesian_to_spherical_coordinates


class MagneticFieldInterpolator(MagneticFieldCalculator):
    """
    Extended magnetic field calculator with interpolation capabilities.
    
    Provides functionality for:
    - Linear interpolation between measurement points
    - Trilinear interpolation for 3D grids
    - Resolution enhancement
    """
    
    def interpolate_points_between_two(self, point1, point2, resolution, current):
        """
        Interpolate points between two measurement points.
        
        Args:
            point1, point2: Two measurement points to interpolate between
            resolution: Number of interpolated points to create
            current: Current value for field calculations
            
        Returns:
            list: List of interpolated points with calculated field values
        """
        interpolated_points = []
        
        for i in range(1, resolution):
            fraction = i / resolution
            
            # Linear interpolation of coordinates
            interpolated_position = interpolate_linear(point1, point2, fraction)
            x, y, z = interpolated_position['x'], interpolated_position['y'], interpolated_position['z']
            
            # Calculate field magnitude for interpolated position
            r, theta, _ = cartesian_to_spherical_coordinates(x, y, z)
            H_magnitude = self.calculate_radial_field_component(x, y, z, current)
            
            # Create interpolated point with all necessary attributes
            interpolated_point = {
                'x': x, 'y': y, 'z': z,
                '|H|': H_magnitude,
                'r': r,
                'theta': math.degrees(theta),
                # Copy other field components (could be improved with better interpolation)
                'Hx': point1['Hx'],
                'Hy': point1['Hy'], 
                'Hz': point1['Hz'],
                'Bx': point1['Bx'],
                'By': point1['By'],
                'Bz': point1['Bz']
            }
            
            interpolated_points.append(interpolated_point)
        
        return interpolated_points
    
    def increase_resolution_linear(self, points, resolution, current):
        """
        Increase resolution using linear interpolation between consecutive points.
        
        Args:
            points: List of measurement points
            resolution: Number of interpolated points between each pair
            current: Current value for field calculations
            
        Returns:
            list: High resolution point list
        """
        high_resolution_points = []
        
        for i in range(len(points) - 1):
            point1 = points[i]
            point2 = points[i + 1]
            
            # Add the first point
            high_resolution_points.append(point1)
            
            # Add interpolated points
            interpolated = self.interpolate_points_between_two(point1, point2, resolution, current)
            high_resolution_points.extend(interpolated)
        
        # Add the last point
        if points:
            high_resolution_points.append(points[-1])
        
        self.high_resolution_points = high_resolution_points
        return high_resolution_points
    
    def trilinear_interpolation(self, cube_vertices, u, v, w):
        """
        Perform trilinear interpolation between 8 cube vertices.
        
        Args:
            cube_vertices: List of 8 vertices defining a cube
            u, v, w: Interpolation parameters (0 to 1)
            
        Returns:
            dict: Interpolated point with field components
        """
        # Extract field components from vertices
        Hx_values = [vertex['Hx'] for vertex in cube_vertices]
        Hy_values = [vertex['Hy'] for vertex in cube_vertices]
        Hz_values = [vertex['Hz'] for vertex in cube_vertices]
        
        # Trilinear interpolation formula
        def interpolate_component(values):
            return (
                values[0] * (1 - u) * (1 - v) * (1 - w) +
                values[1] * u * (1 - v) * (1 - w) +
                values[2] * (1 - u) * v * (1 - w) +
                values[3] * u * v * (1 - w) +
                values[4] * (1 - u) * (1 - v) * w +
                values[5] * u * (1 - v) * w +
                values[6] * (1 - u) * v * w +
                values[7] * u * v * w
            )
        
        # Interpolate field components
        Hx_interpolated = interpolate_component(Hx_values)
        Hy_interpolated = interpolate_component(Hy_values)
        Hz_interpolated = interpolate_component(Hz_values)
        
        # Interpolate coordinates
        x_coords = [vertex['x'] for vertex in cube_vertices]
        y_coords = [vertex['y'] for vertex in cube_vertices]
        z_coords = [vertex['z'] for vertex in cube_vertices]
        
        x_interpolated = interpolate_component(x_coords)
        y_interpolated = interpolate_component(y_coords)
        z_interpolated = interpolate_component(z_coords)
        
        return {
            'x': x_interpolated,
            'y': y_interpolated,
            'z': z_interpolated,
            'Hx': Hx_interpolated,
            'Hy': Hy_interpolated,
            'Hz': Hz_interpolated
        }
    
    def increase_resolution_trilinear(self, points, algorithm="linear", current=0):
        """
        Increase resolution using trilinear interpolation for 3D grid data.
        
        Args:
            points: List of measurement points arranged in a 3D grid
            algorithm: Interpolation algorithm ("linear" or others)
            current: Current value for field calculations
            
        Returns:
            list: High resolution interpolated points
        """
        interpolated_points = []
        grid_points = np.array([[point['x'], point['y'], point['z']] for point in points])
        
        # Extract unique coordinates for grid definition
        x_unique = np.sort(np.unique(grid_points[:, 0]))
        y_unique = np.sort(np.unique(grid_points[:, 1]))
        z_unique = np.sort(np.unique(grid_points[:, 2]))
        
        # Create point lookup dictionary for fast access
        points_dict = {
            (round(p['x'], 5), round(p['y'], 5), round(p['z'], 5)): p 
            for p in points
        }
        
        def find_point_by_coordinates(x, y, z):
            """Find point in dictionary with tolerance for floating point errors."""
            rounded_key = (round(x, 5), round(y, 5), round(z, 5))
            return points_dict.get(rounded_key)
        
        # Process each cube in the grid
        for i in range(len(x_unique) - 1):
            for j in range(len(y_unique) - 1):
                for k in range(len(z_unique) - 1):
                    # Get the 8 vertices of the current cube
                    cube_vertices = []
                    for dx, dy, dz in itertools.product([0, 1], repeat=3):
                        x = x_unique[i + dx]
                        y = y_unique[j + dy]
                        z = z_unique[k + dz]
                        vertex = find_point_by_coordinates(x, y, z)
                        
                        if vertex is None:
                            raise ValueError(f"Point not found for coordinates x={x}, y={y}, z={z}")
                        
                        cube_vertices.append(vertex)
                    
                    # Generate interpolated points within the cube
                    steps = self.resolution + 1
                    for u, v, w in itertools.product(np.linspace(0, 1, steps), repeat=3):
                        # Skip corner vertices to avoid duplication
                        if (u == 0 and v == 0 and w == 0) or (u == 1 and v == 1 and w == 1):
                            continue
                        
                        interpolated_point = self.trilinear_interpolation(cube_vertices, u, v, w)
                        
                        if algorithm == "linear":
                            interpolated_points.append(interpolated_point)
        
        self.high_resolution_points = interpolated_points
        return interpolated_points
    
    def execute_interpolation_pipeline(self, interpolation_method="linear"):
        """
        Execute the complete interpolation pipeline.
        
        Args:
            interpolation_method: Method to use ("linear" or "trilinear")
            
        Returns:
            dict: Results including high resolution points
        """
        # Execute base analysis first
        base_results = self.execute_analysis_pipeline()
        
        if interpolation_method == "linear":
            # Use linear interpolation between nearest points
            nearest_points = base_results['nearest_points']
            if self.average_current and nearest_points:
                high_res_points = self.increase_resolution_linear(
                    nearest_points, self.resolution, self.average_current
                )
        elif interpolation_method == "trilinear":
            # Use trilinear interpolation on full grid
            measurement_points = base_results['measurement_results']
            if measurement_points:
                high_res_points = self.increase_resolution_trilinear(
                    measurement_points, "linear", self.average_current or 0
                )
        else:
            high_res_points = []
        
        return {
            **base_results,
            'high_resolution_points': high_res_points,
            'interpolation_method': interpolation_method
        }
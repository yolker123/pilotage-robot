"""
Utility functions for magnetic field calculations.
Contains common coordinate transformations and mathematical operations.
"""
import math
import numpy as np
from .constants import WAVE_NUMBER, PERMEABILITY_VACUUM, ANTENNA_SURFACE, ANGULAR_FREQUENCY


def cartesian_to_spherical_coordinates(x, y, z):
    """
    Convert cartesian coordinates to spherical coordinates.
    
    Args:
        x, y, z: Cartesian coordinates
        
    Returns:
        tuple: (r, theta, phi) in spherical coordinates
               r: radial distance
               theta: polar angle (0 to pi)
               phi: azimuthal angle (0 to 2pi)
    """
    r = math.sqrt(x**2 + y**2 + z**2)
    
    if r == 0:
        return 0, 0, 0
        
    theta = math.acos(z / r) if r != 0 else 0.0
    phi = math.atan2(y, x)
    
    return r, theta, phi


def calculate_magnetic_field_components(amplitude_x, amplitude_y, amplitude_z):
    """
    Calculate magnetic field components from measured amplitudes.
    
    Args:
        amplitude_x, amplitude_y, amplitude_z: Measured amplitudes
        
    Returns:
        tuple: ((Bx, By, Bz), (Hx, Hy, Hz), H_magnitude, theta, phi)
    """
    # Calculate B field components
    Bx = amplitude_x / (ANTENNA_SURFACE * ANGULAR_FREQUENCY)
    By = amplitude_y / (ANTENNA_SURFACE * ANGULAR_FREQUENCY)
    Bz = amplitude_z / (ANTENNA_SURFACE * ANGULAR_FREQUENCY)
    
    # Calculate H field components
    Hx = Bx / PERMEABILITY_VACUUM
    Hy = By / PERMEABILITY_VACUUM
    Hz = Bz / PERMEABILITY_VACUUM
    
    # Calculate magnitude and angles
    H_magnitude = math.sqrt(Hx**2 + Hy**2 + Hz**2)
    
    theta = math.degrees(math.acos(Hz / H_magnitude)) if H_magnitude != 0 else 0.0
    phi = math.degrees(math.atan2(Hy, Hx))
    
    return (Bx, By, Bz), (Hx, Hy, Hz), H_magnitude, theta, phi


def calculate_current_from_radial_field(H_radial, r, theta):
    """
    Calculate current I from radial magnetic field component.
    
    Args:
        H_radial: Radial component of magnetic field
        r: Radial distance
        theta: Polar angle in degrees
        
    Returns:
        float: Current value or None if calculation fails
    """
    if r == 0:
        return None
        
    numerator = 2 * math.pi * H_radial
    factor1 = 1j / (WAVE_NUMBER**2 * r**2)
    factor2 = 1 / (WAVE_NUMBER**3 * r**3)
    denominator = ANTENNA_SURFACE * (WAVE_NUMBER**3) * (factor1 + factor2) * math.cos(math.radians(theta))
    
    if denominator == 0:
        return None
        
    current = numerator / denominator
    return abs(current)


def calculate_current_from_tangential_field(H_tangential, r, theta):
    """
    Calculate current I from tangential magnetic field component.
    
    Args:
        H_tangential: Tangential component of magnetic field
        r: Radial distance
        theta: Polar angle in degrees
        
    Returns:
        float: Current value or None if calculation fails
    """
    if np.isclose(theta, 0) or r == 0:
        return None
        
    sin_theta = np.sin(np.radians(theta))
    if np.isclose(sin_theta, 0):
        return None
        
    factor = -1 / (WAVE_NUMBER * r) + 1j / (WAVE_NUMBER**2 * r**2) + 1 / (WAVE_NUMBER**3 * r**3)
    current = (4 * np.pi * H_tangential) / (ANTENNA_SURFACE * WAVE_NUMBER**3 * factor * sin_theta)
    
    return abs(current)


def interpolate_linear(point1, point2, fraction):
    """
    Linear interpolation between two points.
    
    Args:
        point1, point2: Dictionaries containing point data with 'x', 'y', 'z' keys
        fraction: Interpolation fraction (0 to 1)
        
    Returns:
        dict: Interpolated point
    """
    x = point1['x'] + (point2['x'] - point1['x']) * fraction
    y = point1['y'] + (point2['y'] - point1['y']) * fraction
    z = point1['z'] + (point2['z'] - point1['z']) * fraction
    
    return {'x': x, 'y': y, 'z': z}


def filter_points_excluding_origin(points):
    """
    Filter out points at the origin (0, 0, 0).
    
    Args:
        points: List of point dictionaries
        
    Returns:
        list: Filtered points excluding origin
    """
    return [point for point in points if not (point['x'] == 0 and point['y'] == 0 and point['z'] == 0)]


def select_closest_points_by_axis(points, axis):
    """
    Select the closest points along a specific axis.
    
    Args:
        points: List of point dictionaries
        axis: Axis to consider ('x', 'y', or 'z')
        
    Returns:
        list: Up to 2 closest points along the specified axis
    """
    axis_points = sorted(
        points,
        key=lambda point: abs(point[axis]) if all(point[a] == 0 for a in 'xyz' if a != axis) else float('inf')
    )[:2]
    
    # Filter to keep only one point if they have the same sign
    if len(axis_points) > 1 and (axis_points[0][axis] * axis_points[1][axis] > 0):
        axis_points = [axis_points[0]]
        
    return axis_points
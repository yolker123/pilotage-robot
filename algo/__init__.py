"""
Magnetic field analysis package.

This package provides refactored and improved functionality for:
- Magnetic field calculations and analysis
- File processing for measurement data
- Point interpolation and resolution enhancement
- 3D visualization of field vectors

Main classes:
- MagneticFieldCalculator: Base class for field calculations
- MagneticFieldInterpolator: Extended class with interpolation capabilities
- MeasurementFileProcessor: Utility for processing measurement files

Usage:
    from algo.magnetic_field_interpolation import MagneticFieldInterpolator
    
    calculator = MagneticFieldInterpolator('path/to/data', resolution=5)
    results = calculator.execute_interpolation_pipeline()
"""

from .constants import *
from .utils import *
from .magnetic_field_base import MagneticFieldCalculator
from .magnetic_field_interpolation import MagneticFieldInterpolator
from .file_processor import MeasurementFileProcessor

__version__ = "1.0.0"
__author__ = "Refactored magnetic field analysis tools"

__all__ = [
    # Constants
    'SPEED_OF_LIGHT', 'PERMEABILITY_VACUUM', 'FREQUENCY', 'ANGULAR_FREQUENCY',
    'ANTENNA_SURFACE', 'WAVE_NUMBER', 'DEFAULT_RESOLUTION',
    
    # Utility functions
    'cartesian_to_spherical_coordinates', 'calculate_magnetic_field_components',
    'calculate_current_from_radial_field', 'calculate_current_from_tangential_field',
    'interpolate_linear', 'filter_points_excluding_origin', 'select_closest_points_by_axis',
    
    # Main classes
    'MagneticFieldCalculator', 'MagneticFieldInterpolator', 'MeasurementFileProcessor'
]
"""
Shared constants for magnetic field calculations.
Centralizes physical constants and configuration values used across the project.
"""
import math

# Physical constants
SPEED_OF_LIGHT = 3e8  # Speed of light in m/s
PERMEABILITY_VACUUM = 4 * math.pi * 1e-7  # Magnetic permeability of vacuum in T·m/A
FREQUENCY = 13.56e6  # Frequency in Hz
ANGULAR_FREQUENCY = 2 * math.pi * FREQUENCY  # Angular frequency in rad/s
ANTENNA_SURFACE = 1e-4  # Antenna surface in m² (1 cm²)
WAVE_NUMBER = ANGULAR_FREQUENCY / SPEED_OF_LIGHT  # Wave number in rad/m

# Default calculation parameters
DEFAULT_RESOLUTION = 5  # Default interpolation resolution
DEFAULT_CPU_LIMIT = 80  # Default CPU usage limit for processing
DEFAULT_RAM_LIMIT = 80  # Default RAM usage limit for processing
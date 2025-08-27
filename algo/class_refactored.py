"""
Refactored magnetic field calculation class.
This replaces the original class.py with improved structure and function names.

DEPRECATED: This file is provided for compatibility. 
Use MagneticFieldInterpolator from magnetic_field_interpolation module instead.
"""
import os
import warnings
from .magnetic_field_interpolation import MagneticFieldInterpolator


class ChampMagnetique(MagneticFieldInterpolator):
    """
    DEPRECATED: Legacy compatibility class.
    
    This class provides backward compatibility with the original ChampMagnetique class
    while using the new refactored implementation underneath.
    
    New code should use MagneticFieldInterpolator directly.
    """
    
    def __init__(self, dossier, resolution=5):
        warnings.warn(
            "ChampMagnetique is deprecated. Use MagneticFieldInterpolator instead.",
            DeprecationWarning,
            stacklevel=2
        )
        
        # Initialize with new base class
        super().__init__(dossier, resolution)
        
        # Legacy attribute names for compatibility
        from .constants import (
            SPEED_OF_LIGHT, PERMEABILITY_VACUUM, FREQUENCY, 
            ANGULAR_FREQUENCY, ANTENNA_SURFACE, WAVE_NUMBER
        )
        
        self.dossier = dossier
        self.c = SPEED_OF_LIGHT
        self.mu0 = PERMEABILITY_VACUUM
        self.F = FREQUENCY
        self.omega = ANGULAR_FREQUENCY
        self.S = ANTENNA_SURFACE
        self.k = WAVE_NUMBER
        self.resultats = self.measurement_results
        self.points_haute_resolution = self.high_resolution_points
    
    # Legacy method names for backward compatibility
    def extraire_amplitude_et_position(self, file_content, filename):
        """Legacy method name. Use extract_amplitude_and_position_from_file instead."""
        return self.extract_amplitude_and_position_from_file(file_content, filename)
    
    def calcul_champ_magnetique(self, amp_x, amp_y, amp_z):
        """Legacy method name. Use calculate_magnetic_field_components instead."""
        from .utils import calculate_magnetic_field_components
        return calculate_magnetic_field_components(amp_x, amp_y, amp_z)
    
    def traiter_fichiers(self):
        """Legacy method name. Use process_measurement_files instead."""
        results = self.process_measurement_files()
        self.resultats = results  # Update legacy attribute
        return results
    
    def calcul_r_teta(self, x, y, z):
        """Legacy method name. Use cartesian_to_spherical_coordinates instead."""
        from .utils import cartesian_to_spherical_coordinates
        r, theta, _ = cartesian_to_spherical_coordinates(x, y, z)
        return r, theta
    
    def selectionner_points_proches(self):
        """Legacy method name. Use select_nearest_points instead."""
        return self.select_nearest_points()
    
    def moyenne_I(self, points_proches):
        """Legacy method name. Use calculate_average_current instead."""
        return self.calculate_average_current(points_proches)
    
    def calculer_I(self, H_r, r, theta):
        """Legacy method name. Use calculate_current_from_radial_field instead."""
        from .utils import calculate_current_from_radial_field
        return calculate_current_from_radial_field(H_r, r, theta)
    
    def calculer_I_Htetha(self, H_theta, r, theta):
        """Legacy method name. Use calculate_current_from_tangential_field instead."""
        from .utils import calculate_current_from_tangential_field
        return calculate_current_from_tangential_field(H_theta, r, theta)
    
    def augmenter_resolution(self, points_proches):
        """Legacy method name. Use increase_resolution_linear instead."""
        if self.average_current:
            self.points_haute_resolution = self.increase_resolution_linear(
                points_proches, self.resolution, self.average_current
            )
        else:
            self.points_haute_resolution = []
    
    def interpoler_points(self, p1, p2, resolution):
        """Legacy method name. Use interpolate_points_between_two instead."""
        if self.average_current:
            return self.interpolate_points_between_two(p1, p2, resolution, self.average_current)
        return []
    
    def calculer_H_total(self, x, y, z, I):
        """Legacy method name. Use calculate_total_field_magnitude instead."""
        return self.calculate_total_field_magnitude(x, y, z, I)
    
    def calculer_Hr(self, x, y, z, I):
        """Legacy method name. Use calculate_radial_field_component instead."""
        return self.calculate_radial_field_component(x, y, z, I)
    
    def calculer_Htheta(self, x, y, z, I):
        """Legacy method name. Use calculate_tangential_field_component instead."""
        return self.calculate_tangential_field_component(x, y, z, I)
    
    def afficher_vecteurs_3D(self, points=None):
        """Legacy method name. Use visualize_3d_vectors instead."""
        self.visualize_3d_vectors(points)
    
    def execute_pipeline(self):
        """Legacy method name. Use execute_interpolation_pipeline instead."""
        results = self.execute_interpolation_pipeline()
        
        # Update legacy attributes for compatibility
        self.resultats = self.measurement_results
        self.points_haute_resolution = self.high_resolution_points
        
        return results


# Example usage with backward compatibility
if __name__ == "__main__":
    # This maintains the original usage pattern
    dossier_mesures = "logMeasure"
    resolution_interpolation = 10
    
    # Create instance using legacy class name
    champ_magnetique = ChampMagnetique(dossier=dossier_mesures, resolution=resolution_interpolation)
    
    # Execute using legacy method name
    champ_magnetique.execute_pipeline()
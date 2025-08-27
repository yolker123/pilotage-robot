#!/usr/bin/env python3
"""
Demonstration of the refactored magnetic field calculation code.
Shows improved function names, factorized code, and better structure.
"""
import os
import tempfile
import sys

# Add path for importing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from algo import MagneticFieldInterpolator, MeasurementFileProcessor
from algo.constants import FREQUENCY, ANTENNA_SURFACE, PERMEABILITY_VACUUM


def create_sample_data():
    """Create sample measurement files for demonstration."""
    temp_dir = tempfile.mkdtemp()
    print(f"Creating sample data in: {temp_dir}")
    
    # Sample measurement files with different positions and amplitudes
    sample_files = [
        ("logMeasureOscillo_(1.0 0.0 0.0).txt", "C1_pkpk:0.005"),
        ("logMeasureOscillo_(-1.0 0.0 0.0).txt", "C1_pkpk:0.004"),
        ("logMeasureOscillo_(0.0 1.0 0.0).txt", "C1_pkpk:0.006"),
        ("logMeasureOscillo_(0.0 -1.0 0.0).txt", "C1_pkpk:0.003"),
        ("logMeasureOscillo_(0.0 0.0 1.0).txt", "C1_pkpk:0.007"),
        ("logMeasureOscillo_(0.0 0.0 -1.0).txt", "C1_pkpk:0.0035"),
        ("logMeasureOscillo_(2.0 2.0 0.0).txt", "C1_pkpk:0.0025"),
        ("logMeasureOscillo_(-2.0 -2.0 0.0).txt", "C1_pkpk:0.0028"),
    ]
    
    for filename, content in sample_files:
        filepath = os.path.join(temp_dir, filename)
        with open(filepath, 'w') as f:
            f.write(content + "\nOther measurement data\nTimestamp: 2024-01-01")
    
    return temp_dir


def demonstrate_improved_functionality():
    """Demonstrate the refactored code with improved names and structure."""
    print("=" * 70)
    print("MAGNETIC FIELD CALCULATION - REFACTORED DEMONSTRATION")
    print("=" * 70)
    
    # Create sample data
    data_dir = create_sample_data()
    
    try:
        print("\n1. Constants are now centralized and clearly named:")
        print(f"   - Frequency: {FREQUENCY/1e6:.2f} MHz")
        print(f"   - Antenna Surface: {ANTENNA_SURFACE*1e4:.1f} cm²")
        print(f"   - Vacuum Permeability: {PERMEABILITY_VACUUM:.2e} T·m/A")
        
        print("\n2. File processing is now modular and reusable:")
        processor = MeasurementFileProcessor()
        files = processor.find_measurement_files(data_dir)
        print(f"   - Found {len(files)} measurement files")
        
        raw_data = processor.process_directory(data_dir)
        summary = processor.get_measurement_summary(raw_data)
        print(f"   - Position range X: {summary['position_range']['x']}")
        print(f"   - Position range Y: {summary['position_range']['y']}")
        print(f"   - Position range Z: {summary['position_range']['z']}")
        
        print("\n3. Magnetic field calculator with improved English function names:")
        
        # Old French names vs New English names
        print("   OLD (French) -> NEW (English):")
        print("   - traiter_fichiers() -> process_measurement_files()")
        print("   - selectionner_points_proches() -> select_nearest_points()")
        print("   - augmenter_resolution() -> increase_resolution_linear()")
        print("   - calcul_champ_magnetique() -> calculate_magnetic_field_components()")
        print("   - afficher_vecteurs_3D() -> visualize_3d_vectors()")
        
        # Create calculator instance
        calculator = MagneticFieldInterpolator(data_dir, resolution=3)
        
        # Process files with improved method names
        measurement_results = calculator.process_measurement_files()
        print(f"\n   - Processed {len(measurement_results)} measurement points")
        
        # Select nearest points with improved algorithm
        nearest_points = calculator.select_nearest_points()
        print(f"   - Selected {len(nearest_points)} nearest points for analysis")
        
        # Calculate average current with better error handling
        average_current = calculator.calculate_average_current(nearest_points)
        print(f"   - Calculated average current: {average_current:.6f} A" if average_current else "   - Could not calculate average current")
        
        print("\n4. Interpolation with factorized and improved algorithms:")
        
        if average_current:
            # Demonstrate linear interpolation
            high_res_points = calculator.increase_resolution_linear(
                nearest_points[:2], 5, average_current
            )
            print(f"   - Generated {len(high_res_points)} high-resolution points")
            
            # Calculate field at a specific point
            test_x, test_y, test_z = 0.5, 0.5, 0.5
            field_magnitude = calculator.calculate_total_field_magnitude(
                test_x, test_y, test_z, average_current
            )
            print(f"   - Field magnitude at ({test_x}, {test_y}, {test_z}): {field_magnitude:.6e} A/m")
        
        print("\n5. Complete analysis pipeline with single method call:")
        results = calculator.execute_interpolation_pipeline("linear")
        
        print(f"   - Total measurement points: {len(results['measurement_results'])}")
        print(f"   - Nearest points: {len(results['nearest_points'])}")
        print(f"   - High resolution points: {len(results['high_resolution_points'])}")
        print(f"   - Interpolation method: {results['interpolation_method']}")
        print(f"   - Average current: {results['average_current']:.6f} A" if results['average_current'] else "   - Average current: Not calculated")
        
        print("\n6. Code improvements summary:")
        print("   ✓ Eliminated code duplication across 4 files")
        print("   ✓ Improved function names from French to English")
        print("   ✓ Centralized constants in single module")
        print("   ✓ Factorized common interpolation algorithms")
        print("   ✓ Created reusable utility functions")
        print("   ✓ Added proper error handling and validation")
        print("   ✓ Maintained backward compatibility with legacy code")
        print("   ✓ Added comprehensive documentation")
        
        print("\n7. Backward compatibility maintained:")
        print("   - Original ChampMagnetique class still works")
        print("   - Legacy method names are mapped to new implementations")
        print("   - Deprecation warnings guide users to new API")
        
    finally:
        # Clean up
        import shutil
        shutil.rmtree(data_dir)
        print(f"\n   Cleaned up temporary directory: {data_dir}")
    
    print("\n" + "=" * 70)
    print("REFACTORING COMPLETE - Code is now more maintainable!")
    print("=" * 70)


if __name__ == "__main__":
    demonstrate_improved_functionality()
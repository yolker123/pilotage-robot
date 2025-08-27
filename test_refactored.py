#!/usr/bin/env python3
"""
Test script for the refactored magnetic field calculation code.
"""
import sys
import os
import tempfile

# Add the parent directory to the path so we can import our modules
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from algo import MagneticFieldInterpolator, MeasurementFileProcessor
from algo.constants import FREQUENCY, PERMEABILITY_VACUUM


def test_constants():
    """Test that constants are properly defined."""
    print("Testing constants...")
    assert FREQUENCY == 13.56e6, f"Expected frequency 13.56e6, got {FREQUENCY}"
    assert abs(PERMEABILITY_VACUUM - 4 * 3.14159265359 * 1e-7) < 1e-12, "Permeability constant issue"
    print("✓ Constants test passed")


def test_file_processor():
    """Test the file processor with mock data."""
    print("Testing file processor...")
    
    # Create a temporary directory with mock files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a mock measurement file
        mock_filename = "logMeasureOscillo_(1.0 2.0 3.0).txt"
        mock_content = """
        C1_pkpk:0.005
        C2_pkpk:0.003
        Some other data
        """
        
        mock_filepath = os.path.join(temp_dir, mock_filename)
        with open(mock_filepath, 'w') as f:
            f.write(mock_content)
        
        # Test file processing
        processor = MeasurementFileProcessor()
        
        # Test finding files
        files = processor.find_measurement_files(temp_dir)
        assert len(files) == 1, f"Expected 1 file, found {len(files)}"
        assert files[0] == mock_filename
        
        # Test extracting data
        file_content = processor.read_measurement_file(mock_filepath)
        amp_x, amp_y, amp_z, x, y, z = processor.extract_amplitude_and_position(
            file_content, mock_filename
        )
        
        assert amp_x == 0.005, f"Expected amplitude_x 0.005, got {amp_x}"
        assert x == 1.0 and y == 2.0 and z == 3.0, f"Expected position (1,2,3), got ({x},{y},{z})"
        
        print("✓ File processor test passed")


def test_magnetic_field_calculator():
    """Test the magnetic field calculator with mock data."""
    print("Testing magnetic field calculator...")
    
    # Create a temporary directory with mock measurement files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create multiple mock files
        mock_files = [
            ("logMeasureOscillo_(1.0 0.0 0.0).txt", "C1_pkpk:0.005"),
            ("logMeasureOscillo_(-1.0 0.0 0.0).txt", "C1_pkpk:0.004"),
            ("logMeasureOscillo_(0.0 1.0 0.0).txt", "C1_pkpk:0.006"),
            ("logMeasureOscillo_(0.0 0.0 1.0).txt", "C1_pkpk:0.003"),
        ]
        
        for filename, content in mock_files:
            filepath = os.path.join(temp_dir, filename)
            with open(filepath, 'w') as f:
                f.write(content)
        
        # Test the calculator
        calculator = MagneticFieldInterpolator(temp_dir, resolution=3)
        
        # Test file processing
        results = calculator.process_measurement_files()
        assert len(results) == 4, f"Expected 4 results, got {len(results)}"
        
        # Test nearest point selection
        nearest_points = calculator.select_nearest_points()
        assert len(nearest_points) > 0, "Should find some nearest points"
        
        print("✓ Magnetic field calculator test passed")


def test_backward_compatibility():
    """Test that the legacy class still works."""
    print("Testing backward compatibility...")
    
    # Import the legacy class
    from algo.class_refactored import ChampMagnetique
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create a mock file
        mock_filename = "logMeasureOscillo_(1.0 0.0 0.0).txt"
        mock_content = "C1_pkpk:0.005"
        
        mock_filepath = os.path.join(temp_dir, mock_filename)
        with open(mock_filepath, 'w') as f:
            f.write(mock_content)
        
        # Test legacy class
        champ = ChampMagnetique(temp_dir, resolution=2)
        
        # Test legacy method names
        results = champ.traiter_fichiers()
        assert len(results) == 1, f"Expected 1 result, got {len(results)}"
        
        # Test that legacy attributes exist
        assert hasattr(champ, 'dossier')
        assert hasattr(champ, 'mu0')
        assert hasattr(champ, 'F')
        
        print("✓ Backward compatibility test passed")


def main():
    """Run all tests."""
    print("Running refactored magnetic field calculation tests...\n")
    
    try:
        test_constants()
        test_file_processor()
        test_magnetic_field_calculator()
        test_backward_compatibility()
        
        print("\n✅ All tests passed successfully!")
        print("The refactored code is working correctly.")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
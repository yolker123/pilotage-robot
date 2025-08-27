"""
File processing utilities for magnetic field measurement data.
Handles reading, parsing, and extracting data from measurement files.
"""
import os
import re
from typing import List, Dict, Tuple, Optional


class MeasurementFileProcessor:
    """
    Utility class for processing measurement files and extracting data.
    
    Handles various file formats and provides standardized data extraction.
    """
    
    @staticmethod
    def extract_amplitude_and_position(file_content: str, filename: str) -> Tuple[float, float, float, float, float, float]:
        """
        Extract amplitude measurements and position coordinates from file content.
        
        Args:
            file_content: Content of the measurement file
            filename: Name of the file containing position information
            
        Returns:
            tuple: (amplitude_x, amplitude_y, amplitude_z, x, y, z)
        """
        # Extract amplitudes using regex pattern
        matches = re.findall(r"C(\d+)_pkpk:(\d+\.\d+)", file_content)
        amplitudes = {f"C{axis}": float(amp) for axis, amp in matches}
        
        # Extract position coordinates from filename
        pos_match = re.search(
            r"logMeasureOscillo_\((-?\d+\.?\d*) (-?\d+\.?\d*) (-?\d+\.?\d*)\)", 
            filename
        )
        
        if pos_match:
            x, y, z = map(float, pos_match.groups())
        else:
            x = y = z = 0.0  # Default values if position is not found
        
        # Get amplitudes (default to 0 if not found)
        amplitude_x = amplitudes.get('C1', 0.0)
        amplitude_y = amplitudes.get('C2', 0.0)
        amplitude_z = amplitudes.get('C3', 0.0)
        
        return amplitude_x, amplitude_y, amplitude_z, x, y, z
    
    @staticmethod
    def find_measurement_files(directory: str, file_pattern: str = "logMeasureOscillo") -> List[str]:
        """
        Find all measurement files in a directory matching the pattern.
        
        Args:
            directory: Directory to search in
            file_pattern: Pattern to match filenames against
            
        Returns:
            list: List of matching filenames
        """
        if not os.path.exists(directory):
            raise FileNotFoundError(f"Directory not found: {directory}")
        
        matching_files = []
        for filename in os.listdir(directory):
            if filename.startswith(file_pattern):
                matching_files.append(filename)
        
        return sorted(matching_files)  # Sort for consistent ordering
    
    @staticmethod
    def read_measurement_file(filepath: str) -> str:
        """
        Read content from a measurement file.
        
        Args:
            filepath: Full path to the file
            
        Returns:
            str: File content
            
        Raises:
            FileNotFoundError: If file doesn't exist
            IOError: If file can't be read
        """
        try:
            with open(filepath, 'r', encoding='utf-8') as file:
                return file.read()
        except UnicodeDecodeError:
            # Try with different encoding if UTF-8 fails
            with open(filepath, 'r', encoding='latin-1') as file:
                return file.read()
    
    @classmethod
    def process_directory(cls, directory: str, file_pattern: str = "logMeasureOscillo") -> List[Dict]:
        """
        Process all measurement files in a directory and extract data.
        
        Args:
            directory: Directory containing measurement files
            file_pattern: Pattern to match filenames against
            
        Returns:
            list: List of dictionaries containing extracted data
        """
        results = []
        
        # Find all matching files
        measurement_files = cls.find_measurement_files(directory, file_pattern)
        
        for filename in measurement_files:
            filepath = os.path.join(directory, filename)
            
            try:
                # Read file content
                file_content = cls.read_measurement_file(filepath)
                
                # Extract amplitude and position data
                amp_x, amp_y, amp_z, x, y, z = cls.extract_amplitude_and_position(
                    file_content, filename
                )
                
                # Store extracted data
                results.append({
                    'filename': filename,
                    'filepath': filepath,
                    'amplitude_x': amp_x,
                    'amplitude_y': amp_y,
                    'amplitude_z': amp_z,
                    'x': x,
                    'y': y,
                    'z': z,
                    'file_content': file_content
                })
                
            except Exception as e:
                print(f"Warning: Failed to process file {filename}: {e}")
                continue
        
        return results
    
    @staticmethod
    def validate_measurement_data(data: Dict) -> bool:
        """
        Validate that measurement data contains required fields.
        
        Args:
            data: Dictionary containing measurement data
            
        Returns:
            bool: True if data is valid, False otherwise
        """
        required_fields = ['filename', 'amplitude_x', 'x', 'y', 'z']
        
        for field in required_fields:
            if field not in data:
                return False
            
            # Check for numeric fields
            if field in ['amplitude_x', 'x', 'y', 'z']:
                try:
                    float(data[field])
                except (ValueError, TypeError):
                    return False
        
        return True
    
    @staticmethod
    def filter_valid_measurements(measurements: List[Dict]) -> List[Dict]:
        """
        Filter out invalid measurement data entries.
        
        Args:
            measurements: List of measurement dictionaries
            
        Returns:
            list: List of valid measurement dictionaries
        """
        valid_measurements = []
        
        for measurement in measurements:
            if MeasurementFileProcessor.validate_measurement_data(measurement):
                valid_measurements.append(measurement)
            else:
                print(f"Warning: Invalid measurement data in {measurement.get('filename', 'unknown')}")
        
        return valid_measurements
    
    @staticmethod
    def get_measurement_summary(measurements: List[Dict]) -> Dict:
        """
        Generate a summary of processed measurements.
        
        Args:
            measurements: List of measurement dictionaries
            
        Returns:
            dict: Summary statistics and information
        """
        if not measurements:
            return {
                'total_files': 0,
                'valid_files': 0,
                'position_range': None,
                'amplitude_range': None
            }
        
        # Extract numeric data
        positions_x = [m['x'] for m in measurements]
        positions_y = [m['y'] for m in measurements]
        positions_z = [m['z'] for m in measurements]
        amplitudes = [m['amplitude_x'] for m in measurements]
        
        return {
            'total_files': len(measurements),
            'valid_files': len(measurements),
            'position_range': {
                'x': (min(positions_x), max(positions_x)),
                'y': (min(positions_y), max(positions_y)),
                'z': (min(positions_z), max(positions_z))
            },
            'amplitude_range': (min(amplitudes), max(amplitudes)),
            'files_processed': [m['filename'] for m in measurements]
        }
import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def validate_json_structure(json_data):
    """Validate the basic structure of the FDTD JSON file"""
    required_sections = ['format', 'general', 'boundary', 'materials', 
                        'mesh', 'probes', 'sources', 'materialAssociations']
    
    for section in required_sections:
        if section not in json_data:
            print(f"Error: Missing required section '{section}'")
            return False
    return True

def validate_parameters(json_data):
    """Validate specific parameters in the simulation"""
    # Check time step and number of steps
    time_step = json_data['general']['timeStep']
    num_steps = json_data['general']['numberOfSteps']
    total_time = time_step * num_steps
    
    # Calculate expected propagation delay
    wire_length = 30  # meters (from mesh configuration)
    speed_of_light = 3e8  # m/s
    expected_delay = wire_length / speed_of_light
    
    print(f"\nSimulation Parameters:")
    print(f"Time step: {time_step:.2e} s")
    print(f"Number of steps: {num_steps}")
    print(f"Total simulation time: {total_time:.2e} s")
    print(f"Expected propagation delay: {expected_delay:.2e} s")
    
    if total_time < expected_delay:
        print(f"Warning: Total simulation time ({total_time:.2e} s) is less than expected propagation delay ({expected_delay:.2e} s)")
        return False
    return True

def visualize_circuit(json_data):
    """Create a simple visualization of the circuit layout"""
    # Extract coordinates
    coords = json_data['mesh']['coordinates']
    elements = json_data['mesh']['elements']
    
    # Create figure
    plt.figure(figsize=(12, 4))
    
    # Plot wires
    for element in elements:
        if element['type'] == 'polyline':
            coord_ids = element['coordinateIds']
            x = [coords[coord_ids[0]-1]['relativePosition'][0], 
                 coords[coord_ids[1]-1]['relativePosition'][0]]
            z = [coords[coord_ids[0]-1]['relativePosition'][2], 
                 coords[coord_ids[1]-1]['relativePosition'][2]]
            plt.plot(z, x, 'b-', linewidth=2)
    
    # Plot nodes
    for element in elements:
        if element['type'] == 'node':
            coord_id = element['coordinateIds'][0]
            pos = coords[coord_id-1]['relativePosition']
            plt.plot(pos[2], pos[0], 'ro', markersize=10)
    
    plt.title('Circuit Layout')
    plt.xlabel('Z position (cells)')
    plt.ylabel('X position (cells)')
    plt.grid(True)
    plt.savefig('circuit_layout.png')
    plt.close()

def main():
    # Load JSON file
    json_path = Path(__file__).parent / 'closedCircuit.fdtd.json'
    try:
        with open(json_path, 'r') as f:
            json_data = json.load(f)
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return
    
    print("Validating FDTD simulation file...")
    
    # Validate structure
    if not validate_json_structure(json_data):
        print("JSON structure validation failed!")
        return
    
    # Validate parameters
    if not validate_parameters(json_data):
        print("Parameter validation failed!")
        return
    
    # Create visualization
    print("\nCreating circuit visualization...")
    visualize_circuit(json_data)
    print("Circuit layout saved as 'circuit_layout.png'")
    
    print("\nValidation complete! The simulation file appears to be correctly configured.")
    print("You can now run the FDTD simulation.")

if __name__ == "__main__":
    main() 
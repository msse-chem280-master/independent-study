"""
Simulation Loop
"""

import monte_carlo as mcsim
import random

if __name__ == "__main__":
    # Set number of steps
    num_steps = 10000

    # Generate random coordinates in a box
    coordinates, box_length = mcsim.generate_random_coordinates(num_atoms=500, density=0.9)

    # set random seed for reproducible results.
    random.seed(0)
    sim_coords = mcsim.run_simulation(
        coordinates, box_length, 3, 0.9, num_steps, freq=1000
    )




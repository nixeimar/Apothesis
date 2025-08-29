import os
from pathlib import Path
import subprocess
import shutil
import time
import threading
import pygame
import numpy as np


class ApothesisJob:
    def __init__(self, config, lattice, processes=[],colors = None, cell_size = None , fps = None, workdir="."):
        self.config = config
        self.lattice = lattice
        self.processes = processes
        self.workdir = os.getcwd() if workdir == "." else Path(workdir).resolve()
        self.output_dir = self.workdir
        
        # Simulation parameters
        if cell_size is not None:
            self.cell_size = cell_size
        else:
            self.cell_size = 6
        if fps is not None:
            self.fps = fps
        else :
            self.fps = 100
        if colors is not None:
            self.colors = colors
        else:
            self.colors = {
                "X": (0, 0, 0),      # Black for "X"
                "O*": (0, 0, 255),    # Blue for "O*"
                "CO*": (255, 0, 0),    # Red for "CO*"
                "CO2*": (0, 255, 0)   # Green for CO2
            }
        
        # Process control
        self.apothesis_process = None
        self.simulation_running = False
        

    def write_inputs(self):
        # Pass the shared config object to each component
        for process in self.processes:
            process.save(self.config)
        self.lattice.save(self.config)
        
        # Now, save the fully modified config object once at the end
        self.config.save()

    def get_latest_surface_file(self):
        """
        Look in output_dir for the latest .dat file.
        Returns the full file path (or None if not found).
        """
        try:
            files = [f for f in os.listdir(self.output_dir) if f.endswith(".dat")]
        except FileNotFoundError:
            print("Output directory not found:", self.output_dir)
            return None

        if not files:
            return None

        # Sort files by modification time (newest first)
        files.sort(key=lambda f: os.path.getmtime(os.path.join(self.output_dir, f)), reverse=True)
        return os.path.join(self.output_dir, files[0])

    def read_surface_file(self, path):
        """
        Read and parse the surface file.
        """
        if not os.path.exists(path):
            print("File not found:", path)
            return None, None

        surface_data = []
        time_value = "Time: Unknown"
        try:
            with open(path, "r") as file:
                for line in file:
                    line = line.strip()
                    if line.startswith("Time (s):"):
                        time_value = f"Time (s): {line.split(':')[1].strip()}"
                    else:
                        # Assume grid rows are separated by spaces.
                        row = line.split()
                        if row:  # Only add non-empty rows
                            surface_data.append(row)
        except Exception as e:
            print(f"Error reading surface file: {e}")
            return None, None
            
        # Convert list to NumPy array for ease of processing (if needed)
        if surface_data:
            grid = np.array(surface_data)
        else:
            grid = None
        return time_value, grid

    def run_simulation_display(self):
        """
        Runs the Pygame simulation: periodically reload the latest surface file,
        update the display (draw the grid as circles), and overlay the time stamp.
        """
        pygame.init()
        clock = pygame.time.Clock()

        # Wait for initial surface file to be generated
        print("Waiting for surface files to be generated...")
        initial_file = None
        timeout = 10  # 10 seconds timeout
        start_time = time.time()
        
        while initial_file is None and (time.time() - start_time) < timeout:
            initial_file = self.get_latest_surface_file()
            if not initial_file:
                time.sleep(0.5)
                continue

        if not initial_file:
            print("No surface files found within timeout. Exiting simulation.")
            return

        time_str, grid = self.read_surface_file(initial_file)
        if grid is None:
            print("Failed to read grid data from file.")
            return

        grid_height, grid_width = grid.shape
        window_height = grid_height * self.cell_size + 40  # extra space for time text
        window_width = grid_width * self.cell_size
        screen = pygame.display.set_mode((window_width, window_height))
        pygame.display.set_caption("Real-Time Apothesis Surface Simulation")

        # Pygame font for displaying the time
        font = pygame.font.Font(None, 30)

        print("Starting real-time visualization...")
        self.simulation_running = True
        
        while self.simulation_running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    self.simulation_running = False
                    break

            # Load the latest surface file on each loop iteration
            latest_file = self.get_latest_surface_file()
            if latest_file:
                time_str, grid = self.read_surface_file(latest_file)
            else:
                grid = None

            screen.fill((0, 0, 0))  # Clear screen with black background

            # If grid data is available, draw it using circles
            if grid is not None:
                rows, cols = grid.shape
                radius = self.cell_size // 2
                for y in range(rows):
                    for x in range(cols):
                        species = grid[y, x]
                        # Get color based on species; default to gray if not found
                        color = self.colors.get(species, (128, 128, 128))
                        # Center position for the circle in each cell
                        pos_x = x * self.cell_size + radius
                        pos_y = y * self.cell_size + radius
                        pygame.draw.circle(screen, color, (pos_x, pos_y), radius)

            # Draw time information at the bottom of the window
            text_surface = font.render(time_str, True, (255, 255, 255))
            screen.blit(text_surface, (10, grid_height * self.cell_size + 5))

            pygame.display.flip()
            clock.tick(self.fps)

        pygame.quit()
        print("Simulation display closed.")

    def run(self, executable="apothesis", simulation=False):
        """
        Run the Apothesis simulation.
        
        Args:
            executable (str): Name or path of the executable
            simulation (bool): If True, start real-time visualization
        """
        self.write_inputs()
        exe_path = os.path.join(self.workdir, executable)
        
        if exe_path is None:
            raise FileNotFoundError(
                f"Could not find '{executable}' in your PATH. "
                "Please install it or provide the full path to the executable."
            )
        
        if simulation:
            # Start Apothesis as a background process
            print("Starting Apothesis simulation...")
            self.apothesis_process = subprocess.Popen([exe_path], cwd=str(self.workdir))
            
            try:
                # Start the visualization in the main thread
                self.run_simulation_display()
            finally:
                # Clean up: terminate the Apothesis process
                if self.apothesis_process and self.apothesis_process.poll() is None:
                    self.apothesis_process.terminate()
                    print("Terminated Apothesis process.")
        else:
            # Run normally (blocking)
            subprocess.run([exe_path], cwd=str(self.workdir), check=True)

    def stop_simulation(self):
        """
        Stop the running simulation.
        """
        self.simulation_running = False
        if self.apothesis_process and self.apothesis_process.poll() is None:
            self.apothesis_process.terminate()
            print("Simulation stopped.")

    def set_simulation_params(self, cell_size=None, fps=None, colors=None):
        """
        Configure simulation visualization parameters.
        
        Args:
            cell_size (int): Size of each cell in pixels
            fps (int): Frames per second for the display
            colors (dict): Color mapping for different species
        """
        if cell_size is not None:
            self.cell_size = cell_size
        if fps is not None:
            self.fps = fps
        if colors is not None:
            self.colors.update(colors)
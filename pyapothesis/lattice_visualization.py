"""
Enhanced Lattice Visualization Module for PyApothesis
Advanced 3D visualization, interactive plots, and real-time monitoring
"""

import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from matplotlib.colors import ListedColormap
from matplotlib.animation import FuncAnimation
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
import pandas as pd
import seaborn as sns

class LatticeVisualizer:
    """
    Advanced lattice visualization for PyApothesis simulations.
    """
    
    def __init__(self, lattice=None):
        """
        Initialize the lattice visualizer.
        
        Args:
            lattice: Lattice object from pyapothesis
        """
        self.lattice = lattice
        self.default_colors = {
            "X": "#2C2C2C",      # Dark gray for empty sites
            "O*": "#E74C3C",     # Red for oxygen
            "CO*": "#3498DB",    # Blue for CO
            "CO2": "#27AE60",    # Green for CO2
            "CO2*": "#9B59B6",   # Purple for CO2*
            "H*": "#F39C12",     # Orange for hydrogen
            "OH*": "#1ABC9C",    # Teal for OH
        }
        
    def plot_3d_lattice(self, surface_file: str = None, 
                       height_data: np.ndarray = None,
                       save_path: Optional[str] = None) -> go.Figure:
        """
        Create a 3D interactive visualization of the lattice.
        
        Args:
            surface_file: Path to surface data file
            height_data: Optional height data for 3D surface
            save_path: Optional path to save HTML file
            
        Returns:
            Plotly 3D figure
        """
        if surface_file:
            # Read surface data
            surface_data = self._read_surface_file(surface_file)
        elif self.lattice:
            surface_data = self._extract_lattice_data()
        else:
            raise ValueError("Either surface_file or lattice object must be provided")
            
        nx, ny = surface_data.shape
        
        # Create coordinate meshes
        x = np.arange(nx)
        y = np.arange(ny)
        X, Y = np.meshgrid(x, y, indexing='ij')
        
        # If no height data, create flat surface with small random variations
        if height_data is None:
            Z = np.random.normal(0, 0.01, (nx, ny))
        else:
            Z = height_data
            
        # Map species to colors
        unique_species = np.unique(surface_data)
        color_map = {species: i for i, species in enumerate(unique_species)}
        color_array = np.array([[color_map[surface_data[i, j]] for j in range(ny)] for i in range(nx)])
        
        # Create 3D surface plot
        fig = go.Figure()
        
        # Add surface
        fig.add_trace(go.Surface(
            x=X, y=Y, z=Z,
            surfacecolor=color_array,
            colorscale=[[i/(len(unique_species)-1), self.default_colors.get(species, '#CCCCCC')] 
                       for i, species in enumerate(unique_species)],
            showscale=False,
            name='Lattice Surface'
        ))
        
        # Add scatter points for better species visualization
        for i, species in enumerate(unique_species):
            mask = surface_data == species
            if np.any(mask):
                x_points, y_points = np.where(mask)
                z_points = Z[mask] + 0.1  # Slightly above surface
                
                fig.add_trace(go.Scatter3d(
                    x=x_points, y=y_points, z=z_points,
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=self.default_colors.get(species, '#CCCCCC'),
                        symbol='circle'
                    ),
                    name=species,
                    hovertemplate=f'<b>{species}</b><br>' +
                                 'X: %{x}<br>' +
                                 'Y: %{y}<br>' +
                                 'Z: %{z:.3f}<br>' +
                                 '<extra></extra>'
                ))
        
        # Update layout
        fig.update_layout(
            title=dict(
                text="3D Lattice Visualization",
                x=0.5,
                font=dict(size=16)
            ),
            scene=dict(
                xaxis_title="X Direction",
                yaxis_title="Y Direction", 
                zaxis_title="Height",
                camera=dict(
                    eye=dict(x=1.5, y=1.5, z=1.5)
                )
            ),
            width=800,
            height=600
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"3D lattice visualization saved to {save_path}")
            
        return fig
        
    def plot_coverage_map(self, surface_file: str = None, 
                         save_path: Optional[str] = None) -> plt.Figure:
        """
        Create a detailed coverage map with site information.
        
        Args:
            surface_file: Path to surface data file
            save_path: Optional path to save figure
            
        Returns:
            Matplotlib figure
        """
        if surface_file:
            surface_data = self._read_surface_file(surface_file)
        elif self.lattice:
            surface_data = self._extract_lattice_data()
        else:
            raise ValueError("Either surface_file or lattice object must be provided")
            
        nx, ny = surface_data.shape
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Main coverage map
        unique_species = np.unique(surface_data)
        species_to_int = {species: i for i, species in enumerate(unique_species)}
        int_surface = np.array([[species_to_int[surface_data[i, j]] for j in range(ny)] for i in range(nx)])
        
        colors = [self.default_colors.get(species, '#CCCCCC') for species in unique_species]
        cmap = ListedColormap(colors)
        
        im = axes[0, 0].imshow(int_surface.T, origin='lower', cmap=cmap, interpolation='nearest')
        axes[0, 0].set_title("Coverage Map", fontsize=14, fontweight='bold')
        axes[0, 0].set_xlabel("X Direction")
        axes[0, 0].set_ylabel("Y Direction")
        
        # Add grid
        axes[0, 0].set_xticks(np.arange(-0.5, nx, 1), minor=True)
        axes[0, 0].set_yticks(np.arange(-0.5, ny, 1), minor=True)
        axes[0, 0].grid(which='minor', color='white', linestyle='-', linewidth=0.5, alpha=0.5)
        
        # Species distribution pie chart
        species_counts = {species: np.sum(surface_data == species) for species in unique_species}
        species_fractions = [count / (nx * ny) for count in species_counts.values()]
        
        wedges, texts, autotexts = axes[0, 1].pie(
            species_fractions, 
            labels=unique_species, 
            colors=[self.default_colors.get(species, '#CCCCCC') for species in unique_species],
            autopct='%1.1f%%',
            startangle=90
        )
        axes[0, 1].set_title("Species Distribution", fontsize=14, fontweight='bold')
        
        # Coverage along X direction (averaged over Y)
        x_coverage = {}
        for species in unique_species:
            coverage_x = np.mean(surface_data == species, axis=1)
            axes[1, 0].plot(range(nx), coverage_x, 
                           label=species, color=self.default_colors.get(species, '#CCCCCC'),
                           marker='o', linewidth=2)
        
        axes[1, 0].set_xlabel("X Position")
        axes[1, 0].set_ylabel("Coverage Fraction")
        axes[1, 0].set_title("Coverage Profile (X Direction)", fontsize=14, fontweight='bold')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Coverage along Y direction (averaged over X)
        for species in unique_species:
            coverage_y = np.mean(surface_data == species, axis=0)
            axes[1, 1].plot(range(ny), coverage_y,
                           label=species, color=self.default_colors.get(species, '#CCCCCC'),
                           marker='s', linewidth=2)
        
        axes[1, 1].set_xlabel("Y Position")
        axes[1, 1].set_ylabel("Coverage Fraction")
        axes[1, 1].set_title("Coverage Profile (Y Direction)", fontsize=14, fontweight='bold')
        axes[1, 1].legend()
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Coverage map saved to {save_path}")
            
        return fig
        
    def plot_site_coordination(self, surface_file: str = None,
                             neighbor_distance: float = 1.5,
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Analyze and visualize site coordination and local environment.
        
        Args:
            surface_file: Path to surface data file
            neighbor_distance: Distance threshold for neighbor detection
            save_path: Optional path to save figure
            
        Returns:
            Matplotlib figure
        """
        if surface_file:
            surface_data = self._read_surface_file(surface_file)
        elif self.lattice:
            surface_data = self._extract_lattice_data()
        else:
            raise ValueError("Either surface_file or lattice object must be provided")
            
        nx, ny = surface_data.shape
        
        # Calculate coordination numbers and local environments
        coordination_data = np.zeros((nx, ny))
        local_env_data = {}
        
        for i in range(nx):
            for j in range(ny):
                # Find neighbors within distance threshold
                neighbors = []
                for di in range(-2, 3):
                    for dj in range(-2, 3):
                        if di == 0 and dj == 0:
                            continue
                        ni, nj = i + di, j + dj
                        if 0 <= ni < nx and 0 <= nj < ny:
                            distance = np.sqrt(di**2 + dj**2)
                            if distance <= neighbor_distance:
                                neighbors.append(surface_data[ni, nj])
                
                coordination_data[i, j] = len(neighbors)
                
                # Count neighbor species
                center_species = surface_data[i, j]
                if center_species not in local_env_data:
                    local_env_data[center_species] = {}
                
                for neighbor_species in neighbors:
                    if neighbor_species not in local_env_data[center_species]:
                        local_env_data[center_species][neighbor_species] = 0
                    local_env_data[center_species][neighbor_species] += 1
        
        # Create visualization
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        
        # Coordination number map
        im1 = axes[0, 0].imshow(coordination_data.T, origin='lower', cmap='viridis')
        axes[0, 0].set_title("Coordination Number Map", fontsize=14, fontweight='bold')
        axes[0, 0].set_xlabel("X Direction")
        axes[0, 0].set_ylabel("Y Direction")
        plt.colorbar(im1, ax=axes[0, 0], label='Coordination Number')
        
        # Species map for reference
        unique_species = np.unique(surface_data)
        species_to_int = {species: i for i, species in enumerate(unique_species)}
        int_surface = np.array([[species_to_int[surface_data[i, j]] for j in range(ny)] for i in range(nx)])
        
        colors = [self.default_colors.get(species, '#CCCCCC') for species in unique_species]
        cmap = ListedColormap(colors)
        
        axes[0, 1].imshow(int_surface.T, origin='lower', cmap=cmap, interpolation='nearest')
        axes[0, 1].set_title("Species Map (Reference)", fontsize=14, fontweight='bold')
        axes[0, 1].set_xlabel("X Direction")
        axes[0, 1].set_ylabel("Y Direction")
        
        # Coordination number distribution
        coord_counts = np.bincount(coordination_data.astype(int).flatten())
        axes[1, 0].bar(range(len(coord_counts)), coord_counts, alpha=0.7, color='skyblue', edgecolor='black')
        axes[1, 0].set_xlabel("Coordination Number")
        axes[1, 0].set_ylabel("Number of Sites")
        axes[1, 0].set_title("Coordination Number Distribution", fontsize=14, fontweight='bold')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Local environment analysis (heatmap)
        if local_env_data:
            # Create matrix for heatmap
            all_species = list(unique_species)
            env_matrix = np.zeros((len(all_species), len(all_species)))
            
            for i, center_species in enumerate(all_species):
                for j, neighbor_species in enumerate(all_species):
                    if (center_species in local_env_data and 
                        neighbor_species in local_env_data[center_species]):
                        env_matrix[i, j] = local_env_data[center_species][neighbor_species]
            
            # Normalize by row (center species)
            row_sums = env_matrix.sum(axis=1, keepdims=True)
            env_matrix_norm = np.divide(env_matrix, row_sums, 
                                      out=np.zeros_like(env_matrix), where=row_sums!=0)
            
            sns.heatmap(env_matrix_norm, 
                       xticklabels=all_species, 
                       yticklabels=all_species,
                       annot=True, fmt='.2f', cmap='Blues',
                       ax=axes[1, 1])
            axes[1, 1].set_title("Local Environment Matrix\n(Normalized by Center Species)", 
                                fontsize=14, fontweight='bold')
            axes[1, 1].set_xlabel("Neighbor Species")
            axes[1, 1].set_ylabel("Center Species")
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Site coordination analysis saved to {save_path}")
            
        return fig
        
    def create_animated_lattice(self, surface_files: List[str], 
                              save_path: Optional[str] = None,
                              fps: int = 5) -> FuncAnimation:
        """
        Create an animated visualization of lattice evolution.
        
        Args:
            surface_files: List of surface data files in chronological order
            save_path: Optional path to save animation (mp4 or gif)
            fps: Frames per second for animation
            
        Returns:
            Matplotlib animation object
        """
        if not surface_files:
            raise ValueError("No surface files provided")
            
        # Read all surface data
        surface_data_list = []
        time_points = []
        
        for surface_file in surface_files:
            data = self._read_surface_file(surface_file)
            surface_data_list.append(data)
            
            # Extract time from filename if possible
            file_path = Path(surface_file)
            try:
                # Assuming filename format like "SurfaceSpecies_0.123456.dat"
                time_str = file_path.stem.split('_')[-1]
                time_points.append(float(time_str))
            except:
                time_points.append(len(time_points))
        
        # Set up the figure and axis
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Get unique species across all files
        all_species = set()
        for data in surface_data_list:
            all_species.update(np.unique(data))
        all_species = sorted(list(all_species))
        
        species_to_int = {species: i for i, species in enumerate(all_species)}
        colors = [self.default_colors.get(species, '#CCCCCC') for species in all_species]
        cmap = ListedColormap(colors)
        
        # Initialize the plot
        int_surface = np.array([[species_to_int[surface_data_list[0][i, j]] 
                               for j in range(surface_data_list[0].shape[1])] 
                              for i in range(surface_data_list[0].shape[0])])
        
        im = ax.imshow(int_surface.T, origin='lower', cmap=cmap, 
                      interpolation='nearest', animated=True)
        
        title = ax.set_title(f"Lattice Evolution - Time: {time_points[0]:.4f}", 
                           fontsize=14, fontweight='bold')
        ax.set_xlabel("X Direction")
        ax.set_ylabel("Y Direction")
        
        # Add legend
        legend_elements = [plt.Rectangle((0,0),1,1, color=self.default_colors.get(species, '#CCCCCC'), 
                                       label=species) for species in all_species]
        ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
        
        def animate(frame):
            """Animation function"""
            data = surface_data_list[frame]
            int_surface = np.array([[species_to_int[data[i, j]] 
                                   for j in range(data.shape[1])] 
                                  for i in range(data.shape[0])])
            
            im.set_array(int_surface.T)
            title.set_text(f"Lattice Evolution - Time: {time_points[frame]:.4f}")
            return [im, title]
        
        # Create animation
        anim = FuncAnimation(fig, animate, frames=len(surface_data_list),
                           interval=1000//fps, blit=True, repeat=True)
        
        plt.tight_layout()
        
        if save_path:
            save_path = Path(save_path)
            if save_path.suffix.lower() == '.mp4':
                anim.save(save_path, writer='ffmpeg', fps=fps, dpi=150)
                print(f"Animation saved as MP4 to {save_path}")
            elif save_path.suffix.lower() == '.gif':
                anim.save(save_path, writer='pillow', fps=fps, dpi=150)
                print(f"Animation saved as GIF to {save_path}")
            else:
                print("Unsupported format. Use .mp4 or .gif")
                
        return anim
        
    def _read_surface_file(self, surface_file: str) -> np.ndarray:
        """Read surface data from file."""
        with open(surface_file, 'r') as f:
            lines = f.readlines()
        
        surface_data = []
        for line in lines:
            line = line.strip()
            if not line.startswith("Time"):
                row = line.split()
                if row:
                    surface_data.append(row)
        
        return np.array(surface_data)
        
    def _extract_lattice_data(self) -> np.ndarray:
        """Extract surface data from lattice object."""
        if not self.lattice:
            raise ValueError("No lattice object available")
            
        nx, ny = self.lattice.nx, self.lattice.ny
        surface_data = np.empty((nx, ny), dtype=object)
        
        for i in range(nx):
            for j in range(ny):
                surface_data[i, j] = self.lattice.nodes[i, j, 0].species
                
        return surface_data

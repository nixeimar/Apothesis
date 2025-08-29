"""
Enhanced Visualization Module for PyApothesis
Inspired by pyZacros visualization capabilities with additional ML features
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from matplotlib.colors import ListedColormap
from matplotlib.animation import FuncAnimation
import networkx as nx
from typing import Union, Dict, List, Optional, Tuple
from pathlib import Path
import re
from scipy import signal
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import warnings
warnings.filterwarnings('ignore')

class EnhancedVisualizer:
    """
    Enhanced visualization class for PyApothesis simulations.
    Provides advanced plotting, interactive visualizations, and ML-based analysis.
    """
    
    def __init__(self, results_obj=None):
        """
        Initialize the enhanced visualizer.
        
        Args:
            results_obj: Results object from pyapothesis.results
        """
        self.results = results_obj
        self.default_colors = {
            "X": "#2C2C2C",      # Dark gray for empty sites
            "O*": "#E74C3C",     # Red for oxygen
            "CO*": "#3498DB",    # Blue for CO
            "CO2": "#27AE60",    # Green for CO2
            "CO2*": "#9B59B6",   # Purple for CO2*
            "H*": "#F39C12",     # Orange for hydrogen
            "OH*": "#1ABC9C",    # Teal for OH
        }
        
    def plot_interactive_coverages(self, save_path: Optional[str] = None) -> go.Figure:
        """
        Create an interactive coverage plot using Plotly.
        
        Args:
            save_path: Optional path to save the HTML file
            
        Returns:
            Plotly figure object
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        
        fig = go.Figure()
        
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
                
            fig.add_trace(go.Scatter(
                x=coverages_df['Time (s)'],
                y=coverages_df[species],
                mode='lines',
                name=species,
                line=dict(color=self.default_colors.get(species, 'black'), width=2),
                hovertemplate=f'<b>{species}</b><br>' +
                             'Time: %{x:.3f} s<br>' +
                             'Coverage: %{y:.4f}<br>' +
                             '<extra></extra>'
            ))
        
        fig.update_layout(
            title=dict(
                text="Interactive Species Coverage vs. Time",
                x=0.5,
                font=dict(size=16, family="Arial")
            ),
            xaxis_title="Time (s)",
            yaxis_title="Coverage",
            template="plotly_white",
            hovermode='x unified',
            legend=dict(
                orientation="v",
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.01
            )
        )
        
        if save_path:
            fig.write_html(save_path)
            print(f"Interactive coverage plot saved to {save_path}")
            
        return fig
        
    def plot_phase_space(self, species_x: str, species_y: str, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create a phase space plot between two species coverages.
        
        Args:
            species_x: First species for x-axis
            species_y: Second species for y-axis
            save_path: Optional path to save the figure
            
        Returns:
            Matplotlib figure object
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        
        if species_x not in coverages_df.columns or species_y not in coverages_df.columns:
            raise ValueError(f"Species {species_x} or {species_y} not found in data")
            
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create scatter plot with time-based coloring
        scatter = ax.scatter(
            coverages_df[species_x],
            coverages_df[species_y],
            c=coverages_df['Time (s)'],
            cmap='viridis',
            alpha=0.7,
            s=20
        )
        
        # Add arrow to show trajectory direction
        n_arrows = 10
        step = len(coverages_df) // n_arrows
        for i in range(0, len(coverages_df) - step, step):
            ax.annotate('', 
                       xy=(coverages_df[species_x].iloc[i + step], 
                           coverages_df[species_y].iloc[i + step]),
                       xytext=(coverages_df[species_x].iloc[i], 
                              coverages_df[species_y].iloc[i]),
                       arrowprops=dict(arrowstyle='->', color='red', alpha=0.6))
        
        ax.set_xlabel(f"{species_x} Coverage", fontsize=12)
        ax.set_ylabel(f"{species_y} Coverage", fontsize=12)
        ax.set_title(f"Phase Space: {species_x} vs {species_y}", fontsize=14, fontweight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Time (s)', fontsize=12)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Phase space plot saved to {save_path}")
            
        return fig
        
    def plot_correlation_heatmap(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create a correlation heatmap between different species coverages.
        
        Args:
            save_path: Optional path to save the figure
            
        Returns:
            Matplotlib figure object
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        
        # Remove time column for correlation
        species_data = coverages_df.drop('Time (s)', axis=1)
        
        # Calculate correlation matrix
        corr_matrix = species_data.corr()
        
        fig, ax = plt.subplots(figsize=(10, 8))
        
        # Create heatmap
        mask = np.triu(np.ones_like(corr_matrix, dtype=bool))
        sns.heatmap(corr_matrix, 
                   mask=mask,
                   annot=True, 
                   cmap='RdBu_r', 
                   center=0,
                   square=True,
                   linewidths=0.5,
                   cbar_kws={"shrink": .8},
                   ax=ax)
        
        ax.set_title("Species Coverage Correlation Matrix", fontsize=14, fontweight='bold')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Correlation heatmap saved to {save_path}")
            
        return fig
        
    def plot_lattice_snapshots(self, surface_files: List[str], 
                             titles: Optional[List[str]] = None,
                             save_path: Optional[str] = None) -> plt.Figure:
        """
        Plot multiple lattice snapshots side by side.
        
        Args:
            surface_files: List of surface file paths
            titles: Optional titles for each snapshot
            save_path: Optional path to save the figure
            
        Returns:
            Matplotlib figure object
        """
        n_snapshots = len(surface_files)
        fig, axes = plt.subplots(1, n_snapshots, figsize=(5*n_snapshots, 5))
        
        if n_snapshots == 1:
            axes = [axes]
            
        for i, surface_file in enumerate(surface_files):
            try:
                # Read surface data
                with open(surface_file, 'r') as f:
                    lines = f.readlines()
                
                surface_data = []
                for line in lines:
                    line = line.strip()
                    if not line.startswith("Time"):
                        row = line.split()
                        if row:
                            surface_data.append(row)
                
                if surface_data:
                    grid = np.array(surface_data)
                    
                    # Map species to integers
                    unique_species = sorted(set(grid.flatten()))
                    species_to_int = {sp: idx for idx, sp in enumerate(unique_species)}
                    int_grid = np.array([[species_to_int[cell] for cell in row] for row in grid])
                    
                    # Create colormap
                    colors = [self.default_colors.get(sp, 'gray') for sp in unique_species]
                    cmap = ListedColormap(colors)
                    
                    # Plot
                    im = axes[i].imshow(int_grid, cmap=cmap, interpolation='nearest')
                    
                    if titles and i < len(titles):
                        axes[i].set_title(titles[i], fontsize=12, fontweight='bold')
                    else:
                        axes[i].set_title(f"Snapshot {i+1}", fontsize=12, fontweight='bold')
                        
                    axes[i].set_xticks([])
                    axes[i].set_yticks([])
                    
                    # Add species legend for first plot
                    if i == 0:
                        legend_elements = [plt.Rectangle((0,0),1,1, color=self.default_colors.get(sp, 'gray'), label=sp) 
                                         for sp in unique_species]
                        axes[i].legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5))
                        
            except Exception as e:
                print(f"Error reading {surface_file}: {e}")
                axes[i].text(0.5, 0.5, f"Error loading\n{surface_file}", 
                           ha='center', va='center', transform=axes[i].transAxes)
                
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Lattice snapshots saved to {save_path}")
            
        return fig
        
    def plot_steady_state_analysis(self, window_size: int = 100, 
                                 save_path: Optional[str] = None) -> plt.Figure:
        """
        Analyze and visualize approach to steady state.
        
        Args:
            window_size: Window size for rolling statistics
            save_path: Optional path to save the figure
            
        Returns:
            Matplotlib figure object
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        
        fig, axes = plt.subplots(2, 2, figsize=(15, 10))
        
        # Plot 1: Coverage vs time with rolling mean
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
                
            # Calculate rolling mean and std
            rolling_mean = coverages_df[species].rolling(window=window_size, center=True).mean()
            rolling_std = coverages_df[species].rolling(window=window_size, center=True).std()
            
            axes[0, 0].plot(coverages_df['Time (s)'], coverages_df[species], 
                           alpha=0.3, color=self.default_colors.get(species, 'black'))
            axes[0, 0].plot(coverages_df['Time (s)'], rolling_mean, 
                           label=f"{species} (mean)", 
                           color=self.default_colors.get(species, 'black'), linewidth=2)
            
        axes[0, 0].set_xlabel("Time (s)")
        axes[0, 0].set_ylabel("Coverage")
        axes[0, 0].set_title("Coverage vs Time (with rolling mean)")
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # Plot 2: Variance analysis
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
                
            rolling_var = coverages_df[species].rolling(window=window_size, center=True).var()
            axes[0, 1].plot(coverages_df['Time (s)'], rolling_var, 
                           label=species, color=self.default_colors.get(species, 'black'))
            
        axes[0, 1].set_xlabel("Time (s)")
        axes[0, 1].set_ylabel("Rolling Variance")
        axes[0, 1].set_title("Coverage Variance vs Time")
        axes[0, 1].legend()
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_yscale('log')
        
        # Plot 3: Derivative analysis (rate of change)
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
                
            # Calculate derivative using gradient
            dt = np.gradient(coverages_df['Time (s)'])
            dcov_dt = np.gradient(coverages_df[species]) / dt
            
            axes[1, 0].plot(coverages_df['Time (s)'], np.abs(dcov_dt), 
                           label=species, color=self.default_colors.get(species, 'black'))
            
        axes[1, 0].set_xlabel("Time (s)")
        axes[1, 0].set_ylabel("|dCoverage/dt|")
        axes[1, 0].set_title("Rate of Coverage Change")
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].set_yscale('log')
        
        # Plot 4: Steady state indicator
        total_variance = sum(coverages_df[species].rolling(window=window_size, center=True).var().fillna(0)
                           for species in coverages_df.columns if species != 'Time (s)')
        
        axes[1, 1].plot(coverages_df['Time (s)'], total_variance, 'k-', linewidth=2)
        axes[1, 1].set_xlabel("Time (s)")
        axes[1, 1].set_ylabel("Total System Variance")
        axes[1, 1].set_title("Steady State Indicator")
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_yscale('log')
        
        # Add horizontal line for steady state threshold
        threshold = np.nanpercentile(total_variance, 10)  # 10th percentile as threshold
        axes[1, 1].axhline(y=threshold, color='red', linestyle='--', 
                          label=f'Steady state threshold: {threshold:.2e}')
        axes[1, 1].legend()
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Steady state analysis saved to {save_path}")
            
        return fig
        
    def create_reaction_network_graph(self, processes: List, save_path: Optional[str] = None) -> plt.Figure:
        """
        Create a network graph visualization of the reaction mechanism.
        
        Args:
            processes: List of Process objects
            save_path: Optional path to save the figure
            
        Returns:
            Matplotlib figure object
        """
        G = nx.DiGraph()
        
        # Add nodes for all species
        species = set()
        for process in processes:
            species.update(process.reactants)
            species.update(process.products)
            
        for sp in species:
            G.add_node(sp)
            
        # Add edges for reactions
        for i, process in enumerate(processes):
            # Create a reaction node
            reaction_node = f"R{i+1}"
            G.add_node(reaction_node, type='reaction')
            
            # Connect reactants to reaction
            for reactant in process.reactants:
                G.add_edge(reactant, reaction_node)
                
            # Connect reaction to products
            for product in process.products:
                G.add_edge(reaction_node, product)
                
        fig, ax = plt.subplots(figsize=(12, 10))
        
        # Separate species and reaction nodes
        species_nodes = [n for n in G.nodes() if not n.startswith('R')]
        reaction_nodes = [n for n in G.nodes() if n.startswith('R')]
        
        # Use hierarchical layout
        pos = nx.spring_layout(G, k=3, iterations=50)
        
        # Draw species nodes
        nx.draw_networkx_nodes(G, pos, nodelist=species_nodes, 
                              node_color=[self.default_colors.get(n, 'lightblue') for n in species_nodes],
                              node_size=1000, alpha=0.8, ax=ax)
        
        # Draw reaction nodes
        nx.draw_networkx_nodes(G, pos, nodelist=reaction_nodes,
                              node_color='red', node_shape='s',
                              node_size=500, alpha=0.8, ax=ax)
        
        # Draw edges
        nx.draw_networkx_edges(G, pos, edge_color='gray', arrows=True, 
                              arrowsize=20, alpha=0.6, ax=ax)
        
        # Draw labels
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_weight='bold')
        
        ax.set_title("Reaction Network Graph", fontsize=16, fontweight='bold')
        ax.axis('off')
        
        # Add legend
        legend_elements = [
            plt.Rectangle((0,0),1,1, color='lightblue', label='Species'),
            plt.Rectangle((0,0),1,1, color='red', label='Reactions')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Reaction network graph saved to {save_path}")
            
        return fig

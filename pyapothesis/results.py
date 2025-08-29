import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import re
import glob
from matplotlib.animation import FuncAnimation
import seaborn as sns

class ApothesisResults:
    """
    Results parser and visualizer for Apothesis KMC simulations
    Handles the actual output format from Apothesis
    """
    
    def __init__(self, workdir=".", log_file="Output.log"):
        self.workdir = Path(workdir)
        self.log_file = self.workdir / log_file
        
        # Parsed data containers
        self.time_data = None
        self.coverage_data = {}
        self.process_data = {}
        self.diagnostic_data = {}
        self.simulation_params = {}
        self.processes_info = []
        
        # Parse the log file
        self._parse_apothesis_log()
        
    def _parse_apothesis_log(self):
        """Parse the Apothesis output log file"""
        if not self.log_file.exists():
            print(f"Log file {self.log_file} not found")
            return
            
        with open(self.log_file, 'r') as f:
            content = f.read()
        
        lines = content.strip().split('\n')
        
        # Parse header information
        self._parse_header(lines)
        
        # Find the data table
        data_start_idx = None
        for i, line in enumerate(lines):
            if line.startswith('Time (s)'):
                data_start_idx = i
                break
                
        if data_start_idx is None:
            print("Could not find data table in log file")
            return
            
        # Parse column headers
        header_line = lines[data_start_idx]
        self.column_names = [col.strip() for col in header_line.split('\t')]
        
        # Parse data rows
        self._parse_data_table(lines[data_start_idx + 1:])
        
    def _parse_header(self, lines):
        """Parse simulation parameters from header"""
        for line in lines:
            if 'End time' in line:
                self.simulation_params['end_time'] = float(re.search(r'(\d+\.?\d*)', line).group(1))
            elif 'Temperature' in line:
                self.simulation_params['temperature'] = float(re.search(r'(\d+\.?\d*)', line).group(1))
            elif 'Pressure' in line:
                self.simulation_params['pressure'] = float(re.search(r'(\d+\.?\d*)', line).group(1))
            elif 'Random init num' in line:
                self.simulation_params['seed'] = int(re.search(r'(\d+)', line).group(1))
            elif 'Lattice' in line:
                numbers = re.findall(r'(\d+)', line)
                if len(numbers) >= 2:
                    self.simulation_params['lattice_size'] = (int(numbers[0]), int(numbers[1]))
            elif line.strip().startswith(('CO +', 'O2 +', 'CO* +', 'CO2* ->')):
                self.processes_info.append(line.strip())
                
    def _parse_data_table(self, data_lines):
        """Parse the main data table"""
        # Initialize data containers
        data_dict = {col: [] for col in self.column_names}
        
        for line in data_lines:
            if not line.strip():
                continue
                
            values = line.split('\t')
            if len(values) != len(self.column_names):
                continue
                
            for i, (col, val) in enumerate(zip(self.column_names, values)):
                try:
                    data_dict[col].append(float(val))
                except ValueError:
                    data_dict[col].append(0.0)
        
        # Convert to numpy arrays
        for col in data_dict:
            data_dict[col] = np.array(data_dict[col])
            
        # Extract time data
        self.time_data = data_dict['Time (s)']
        
        # Separate data by type
        self._categorize_data(data_dict)
        
    def _categorize_data(self, data_dict):
        """Categorize data into coverages, processes, and diagnostics"""
        
        for col, data in data_dict.items():
            col_lower = col.lower()
            
            # Coverage data (ends with 'coverage')
            if 'coverage' in col_lower:
                species_name = col.replace('(coverage)', '').replace(' (coverage)', '').strip()
                self.coverage_data[species_name] = data
                
            # Diagnostic data
            elif any(keyword in col_lower for keyword in ['growth rate', 'rms', 'roughness']):
                self.diagnostic_data[col] = data
                
            # Process data (reaction event counts)
            elif any(symbol in col for symbol in ['->', '+', '*']) and 'class size' not in col_lower:
                # Clean up process name
                process_name = col.replace(' (1 N)', '').replace(' (2 N)', '').replace(' (3 N)', '').replace(' (4 N)', '').replace(' (5 N)', '')
                process_name = process_name.replace(' (0 V)', '').replace(' (1 V)', '').replace(' (2 V)', '').replace(' (3 V)', '').replace(' (4 V)', '')
                process_name = process_name.strip()
                
                if process_name not in self.process_data:
                    self.process_data[process_name] = data
                else:
                    # Sum different neighbor variants
                    self.process_data[process_name] = self.process_data[process_name] + data
    
    def get_final_coverages(self):
        """Get final coverage values"""
        if not self.coverage_data:
            return {}
            
        final_coverages = {}
        for species, data in self.coverage_data.items():
            final_coverages[species] = data[-1] if len(data) > 0 else 0.0
            
        return final_coverages
    
    def plot_coverages(self, save_path=None, figsize=(15, 10)):
        """Plot coverage evolution with multiple views"""
        if not self.coverage_data:
            print("No coverage data available")
            return None
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        
        # Colors for different species
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        # 1. Linear coverage evolution
        for i, (species, data) in enumerate(self.coverage_data.items()):
            color = colors[i % len(colors)]
            ax1.plot(self.time_data, data, linewidth=2, label=species, 
                    color=color, marker='o', markersize=2, alpha=0.8)
            
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Coverage')
        ax1.set_title('Surface Coverage Evolution')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1)
        
        # 2. Stacked area plot
        species_names = list(self.coverage_data.keys())
        coverage_matrix = np.array([self.coverage_data[species] for species in species_names])
        
        ax2.stackplot(self.time_data, *coverage_matrix, 
                     labels=species_names, colors=colors[:len(species_names)], alpha=0.7)
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Coverage')
        ax2.set_title('Stacked Coverage (Total Surface Composition)')
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        ax2.grid(True, alpha=0.3)
        ax2.set_ylim(0, 1)
        
        # 3. Coverage rates (derivatives)
        for i, (species, data) in enumerate(self.coverage_data.items()):
            if len(data) > 1:
                rate = np.gradient(data, self.time_data)
                color = colors[i % len(colors)]
                ax3.plot(self.time_data, rate, linewidth=2, label=f'd{species}/dt', 
                        color=color, alpha=0.8)
                
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Coverage Rate (1/s)')
        ax3.set_title('Coverage Change Rates')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        ax3.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        
        # 4. Final coverage pie chart
        final_cov = self.get_final_coverages()
        # Filter out zero coverages
        nonzero_cov = {k: v for k, v in final_cov.items() if v > 0.001}
        
        if nonzero_cov:
            wedges, texts, autotexts = ax4.pie(nonzero_cov.values(), 
                                              labels=nonzero_cov.keys(),
                                              autopct='%1.1f%%',
                                              colors=colors[:len(nonzero_cov)])
            ax4.set_title('Final Surface Composition')
        else:
            ax4.text(0.5, 0.5, 'No significant\ncoverages', 
                    ha='center', va='center', transform=ax4.transAxes)
            
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_process_rates(self, save_path=None, figsize=(15, 10)):
        """Plot process rates and event counts"""
        if not self.process_data:
            print("No process data available")
            return None
            
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=figsize)
        
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        # 1. Event counts over time
        for i, (process, data) in enumerate(self.process_data.items()):
            color = colors[i % len(colors)]
            ax1.plot(self.time_data, data, linewidth=2, label=process, 
                    color=color, marker='o', markersize=2, alpha=0.8)
            
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Cumulative Event Count')
        ax1.set_title('Process Event Counts')
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax1.grid(True, alpha=0.3)
        
        # 2. Event rates (derivatives of event counts)
        for i, (process, data) in enumerate(self.process_data.items()):
            if len(data) > 1:
                rate = np.gradient(data, self.time_data)
                color = colors[i % len(colors)]
                ax2.plot(self.time_data, rate, linewidth=2, label=process, 
                        color=color, alpha=0.8)
                
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('Event Rate (events/s)')
        ax2.set_title('Process Rates')
        ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax2.grid(True, alpha=0.3)
        ax2.set_yscale('log')
        
        # 3. Final event counts bar chart
        final_counts = {process: data[-1] for process, data in self.process_data.items()}
        
        processes = list(final_counts.keys())
        counts = list(final_counts.values())
        
        bars = ax3.bar(range(len(processes)), counts, 
                      color=colors[:len(processes)], alpha=0.7)
        ax3.set_xlabel('Process')
        ax3.set_ylabel('Total Events')
        ax3.set_title('Total Events by Process')
        ax3.set_xticks(range(len(processes)))
        ax3.set_xticklabels([p[:20] + '...' if len(p) > 20 else p for p in processes], 
                           rotation=45, ha='right')
        
        # Add value labels on bars
        for bar, count in zip(bars, counts):
            height = bar.get_height()
            ax3.annotate(f'{int(count)}', xy=(bar.get_x() + bar.get_width()/2, height),
                        xytext=(0, 3), textcoords="offset points", 
                        ha='center', va='bottom', fontsize=8)
        
        # 4. Process efficiency (events per unit time)
        if self.time_data[-1] > 0:
            efficiency = {process: data[-1] / self.time_data[-1] 
                         for process, data in self.process_data.items()}
            
            processes = list(efficiency.keys())
            eff_values = list(efficiency.values())
            
            bars2 = ax4.bar(range(len(processes)), eff_values, 
                           color=colors[:len(processes)], alpha=0.7)
            ax4.set_xlabel('Process')
            ax4.set_ylabel('Average Rate (events/s)')
            ax4.set_title('Average Process Rates')
            ax4.set_xticks(range(len(processes)))
            ax4.set_xticklabels([p[:20] + '...' if len(p) > 20 else p for p in processes], 
                               rotation=45, ha='right')
            ax4.set_yscale('log')
            
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_diagnostics(self, save_path=None, figsize=(12, 8)):
        """Plot surface diagnostics (roughness, growth rate, etc.)"""
        if not self.diagnostic_data:
            print("No diagnostic data available")
            return None
            
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        axes = axes.flatten()
        
        colors = ['red', 'blue', 'green', 'orange']
        
        for i, (diagnostic, data) in enumerate(self.diagnostic_data.items()):
            if i < len(axes):
                ax = axes[i]
                color = colors[i % len(colors)]
                ax.plot(self.time_data, data, linewidth=2, color=color, 
                       marker='o', markersize=3, alpha=0.8)
                ax.set_xlabel('Time (s)')
                ax.set_ylabel(diagnostic)
                ax.set_title(f'{diagnostic} Evolution')
                ax.grid(True, alpha=0.3)
                
        # Hide unused subplots
        for i in range(len(self.diagnostic_data), len(axes)):
            axes[i].set_visible(False)
            
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def plot_comprehensive_overview(self, save_path=None, figsize=(20, 12)):
        """Create a comprehensive overview dashboard"""
        fig = plt.figure(figsize=figsize)
        
        # Create a grid layout
        gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
        
        colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        # 1. Coverage evolution (top-left, spanning 2 columns)
        ax1 = fig.add_subplot(gs[0, :2])
        for i, (species, data) in enumerate(self.coverage_data.items()):
            color = colors[i % len(colors)]
            ax1.plot(self.time_data, data, linewidth=2, label=species, 
                    color=color, alpha=0.8)
        ax1.set_xlabel('Time (s)')
        ax1.set_ylabel('Coverage')
        ax1.set_title('Surface Coverage Evolution')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        ax1.set_ylim(0, 1)
        
        # 2. Final composition pie chart (top-right)
        ax2 = fig.add_subplot(gs[0, 2])
        final_cov = self.get_final_coverages()
        nonzero_cov = {k: v for k, v in final_cov.items() if v > 0.001}
        if nonzero_cov:
            ax2.pie(nonzero_cov.values(), labels=nonzero_cov.keys(), 
                   autopct='%1.1f%%', colors=colors[:len(nonzero_cov)])
            ax2.set_title('Final Composition')
            
        # 3. Simulation parameters text (top-right corner)
        ax3 = fig.add_subplot(gs[0, 3])
        ax3.axis('off')
        param_text = "Simulation Parameters:\n"
        for key, value in self.simulation_params.items():
            param_text += f"{key}: {value}\n"
        ax3.text(0.05, 0.95, param_text, transform=ax3.transAxes, 
                verticalalignment='top', fontfamily='monospace', fontsize=10)
        
        # 4. Process rates (middle-left, spanning 2 columns)
        ax4 = fig.add_subplot(gs[1, :2])
        for i, (process, data) in enumerate(self.process_data.items()):
            if len(data) > 1:
                rate = np.gradient(data, self.time_data)
                color = colors[i % len(colors)]
                ax4.plot(self.time_data, rate, linewidth=2, 
                        label=process[:20] + '...' if len(process) > 20 else process, 
                        color=color, alpha=0.8)
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Event Rate (events/s)')
        ax4.set_title('Process Rates')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
        ax4.set_yscale('log')
        
        # 5. Coverage rates (middle-right)
        ax5 = fig.add_subplot(gs[1, 2:])
        for i, (species, data) in enumerate(self.coverage_data.items()):
            if len(data) > 1:
                rate = np.gradient(data, self.time_data)
                color = colors[i % len(colors)]
                ax5.plot(self.time_data, rate, linewidth=2, label=species, 
                        color=color, alpha=0.8)
        ax5.set_xlabel('Time (s)')
        ax5.set_ylabel('Coverage Rate (1/s)')
        ax5.set_title('Coverage Change Rates')
        ax5.legend()
        ax5.grid(True, alpha=0.3)
        ax5.axhline(y=0, color='black', linestyle='--', alpha=0.5)
        
        # 6. Total events bar chart (bottom-left)
        ax6 = fig.add_subplot(gs[2, :2])
        final_counts = {process: data[-1] for process, data in self.process_data.items()}
        processes = list(final_counts.keys())
        counts = list(final_counts.values())
        
        bars = ax6.bar(range(len(processes)), counts, 
                      color=colors[:len(processes)], alpha=0.7)
        ax6.set_xlabel('Process')
        ax6.set_ylabel('Total Events')
        ax6.set_title('Total Events by Process')
        ax6.set_xticks(range(len(processes)))
        ax6.set_xticklabels([p[:15] + '...' if len(p) > 15 else p for p in processes], 
                           rotation=45, ha='right')
        
        # 7. Diagnostics (bottom-right)
        ax7 = fig.add_subplot(gs[2, 2:])
        for i, (diagnostic, data) in enumerate(self.diagnostic_data.items()):
            if i < 2:  # Show only first 2 diagnostics
                color = colors[i]
                ax7.plot(self.time_data, data, linewidth=2, label=diagnostic, 
                        color=color, alpha=0.8)
        ax7.set_xlabel('Time (s)')
        ax7.set_ylabel('Value')
        ax7.set_title('Surface Diagnostics')
        ax7.legend()
        ax7.grid(True, alpha=0.3)
        
        plt.suptitle('Apothesis KMC Simulation Analysis', fontsize=16, y=0.98)
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            
        return fig
    
    def print_summary(self):
        """Print a summary of the simulation results"""
        print("="*60)
        print("APOTHESIS SIMULATION SUMMARY")
        print("="*60)
        
        # Simulation parameters
        print("\nSimulation Parameters:")
        print("-"*25)
        for key, value in self.simulation_params.items():
            print(f"{key:15}: {value}")
        
        # Processes
        print(f"\nProcesses ({len(self.processes_info)}):")
        print("-"*15)
        for process in self.processes_info:
            print(f"  {process}")
        
        # Final coverages
        print("\nFinal Surface Coverages:")
        print("-"*25)
        final_cov = self.get_final_coverages()
        total_coverage = sum(final_cov.values())
        for species, coverage in sorted(final_cov.items(), key=lambda x: x[1], reverse=True):
            percentage = (coverage / total_coverage * 100) if total_coverage > 0 else 0
            print(f"{species:>8}: {coverage:7.4f} ({percentage:5.1f}%)")
        
        # Process statistics
        print(f"\nTotal Events by Process:")
        print("-"*25)
        final_counts = {process: data[-1] for process, data in self.process_data.items()}
        total_events = sum(final_counts.values())
        for process, count in sorted(final_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / total_events * 100) if total_events > 0 else 0
            print(f"{process[:20]:20}: {int(count):8d} ({percentage:5.1f}%)")
        
        print(f"\nTotal Events: {int(total_events):,}")
        if self.time_data is not None and len(self.time_data) > 0:
            print(f"Average Rate: {total_events/self.time_data[-1]:,.1f} events/s")
        print("="*60)
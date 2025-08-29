import pandas as pd
import matplotlib.pyplot as plt
import re
from pathlib import Path
from typing import Union, Dict, List

# It's good practice to have a default color map for consistency
DEFAULT_SPECIES_COLORS = {
    "X": "gray",
    "O*": "red",
    "CO*": "blue",
    "CO2": "green",
    "CO2*": "purple",
}

class Results:
    """
    Parses, analyzes, and visualizes the output from an Apothesis simulation.
    """
    def __init__(self, workdir: Union[str, Path]):
        self.workdir = Path(workdir)
        # Assuming the log file is named 'Output.log' as per your job script
        self.log_file = self.workdir / "Output.log" 
        self._df = None
        
        if not self.log_file.exists():
            raise FileNotFoundError(f"'{self.log_file.name}' not found in directory: {self.workdir}")

    def get_time_series(self) -> pd.DataFrame:
        """
        Parses the main data table from Output.log into a pandas DataFrame.
        
        This method uses pandas.read_csv for robust parsing of the log file.
        """
        if self._df is not None:
            return self._df

        with open(self.log_file, 'r') as f:
            lines = f.readlines()

        # Find the line number of the header row, which starts with "Time (s)"
        header_line_index = -1
        for i, line in enumerate(lines):
            if line.strip().startswith("Time (s)"):
                header_line_index = i
                break
        
        if header_line_index == -1:
            raise ValueError("Could not find the data header row in Output.log")

        # Use pandas to read the data directly from the file, skipping the preamble.
        # The separator is a tab ('\t').
        df = pd.read_csv(self.log_file, sep='\t', skiprows=header_line_index)
        
        # The last column might be empty due to a trailing tab, so we drop it if it's all NaNs
        df = df.dropna(axis='columns', how='all')

        # Clean up column names by stripping any extra whitespace
        df.columns = df.columns.str.strip()

        # Clean up coverage column names (e.g., "CO* (coverage)" -> "CO*")
        df = df.rename(columns=lambda c: c.replace(' (coverage)', '').strip() if 'coverage' in c else c)
        
        self._df = df
        return self._df

    def get_coverages_df(self) -> pd.DataFrame:
        """
        Returns a DataFrame containing only the time and species coverages.
        """
        full_df = self.get_time_series()
        # Select the 'Time (s)' column and all columns that are in our color map
        coverage_cols = ['Time (s)'] + [col for col in full_df.columns if col in DEFAULT_SPECIES_COLORS]
        
        # Check which of the desired columns actually exist in the DataFrame
        existing_cols = [col for col in coverage_cols if col in full_df.columns]
        return full_df[existing_cols]

    def get_final_coverages(self) -> Dict[str, float]:
        """
        Returns a dictionary of the final surface coverages for each species.
        """
        df = self.get_coverages_df()
        if df.empty:
            return {}
        # Select the last row and drop the time column
        final_vals = df.iloc[-1].drop('Time (s)').to_dict()
        return final_vals

    def plot_coverages(self, save_path: Union[str, Path] = None, **kwargs):
        """
        Plots the species coverages as a function of time.

        Args:
            save_path (str, optional): Path to save the figure. If None, displays the plot.
            **kwargs: Additional keyword arguments passed to plt.figure().
        """
        coverages_df = self.get_coverages_df()
        
        plt.style.use('seaborn-v0_8-whitegrid') # Use a nice style
        fig, ax = plt.subplots(**kwargs)

        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
            
            ax.plot(
                coverages_df['Time (s)'], 
                coverages_df[species], 
                label=species,
                color=DEFAULT_SPECIES_COLORS.get(species, 'black') # Use default color
            )

        ax.set_xlabel("Time (s)", fontsize=12)
        ax.set_ylabel("Coverage", fontsize=12)
        ax.set_title("Species Coverage vs. Time", fontsize=14, fontweight='bold')
        ax.legend(title="Species")
        ax.tick_params(axis='both', which='major', labelsize=10)
        fig.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"Coverage plot saved to {save_path}")
        else:
            plt.show()

    ### --- New Suggested Methods --- ###

    def get_process_rates_df(self) -> pd.DataFrame:
        """
        Returns a DataFrame of the process rates over time.
        """
        full_df = self.get_time_series()
        # Select time and all columns that are process rates (contain '->')
        rate_cols = ['Time (s)'] + [col for col in full_df.columns if '->' in col]
        
        # Check which of the desired columns actually exist in the DataFrame
        existing_cols = [col for col in rate_cols if col in full_df.columns]
        return full_df[existing_cols]

    def plot_process_rates(self, save_path: Union[str, Path] = None, **kwargs):
        """
        Plots the process rates (turnover frequencies) as a function of time.
        """
        rates_df = self.get_process_rates_df()
        
        plt.style.use('seaborn-v0_8-whitegrid')
        fig, ax = plt.subplots(**kwargs)

        for process in rates_df.columns:
            if process == 'Time (s)':
                continue
            ax.plot(rates_df['Time (s)'], rates_df[process], label=process)

        ax.set_xlabel("Time (s)", fontsize=12)
        ax.set_ylabel("Process Rate (events/site/s)", fontsize=12)
        ax.set_title("Process Rates vs. Time", fontsize=14, fontweight='bold')
        ax.legend(title="Process", fontsize='small')
        ax.tick_params(axis='both', which='major', labelsize=10)
        fig.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            print(f"Process rate plot saved to {save_path}")
        else:
            plt.show()

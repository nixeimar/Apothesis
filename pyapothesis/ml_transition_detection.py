"""
Machine Learning Module for Transition Point Detection in KMC Simulations
Inspired by advanced analysis techniques used in catalysis research
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import IsolationForest, RandomForestClassifier
from sklearn.cluster import DBSCAN, KMeans
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.svm import OneClassSVM
from sklearn.neural_network import MLPClassifier
from scipy import signal, stats
from scipy.signal import find_peaks, savgol_filter
from scipy.stats import zscore
import warnings
from typing import List, Dict, Tuple, Optional, Union
from pathlib import Path
import joblib

warnings.filterwarnings('ignore')

class TransitionPointDetector:
    """
    Advanced ML-based detector for transition points in KMC simulations.
    
    Capabilities:
    - Automatic detection of phase transitions
    - Steady-state identification
    - Poisoning detection
    - Oscillatory behavior analysis
    - Catalyst deactivation prediction
    """
    
    def __init__(self, results_obj=None):
        """
        Initialize the transition point detector.
        
        Args:
            results_obj: Results object from pyapothesis.results
        """
        self.results = results_obj
        self.scaler = StandardScaler()
        self.models = {}
        self.transition_points = {}
        
    def extract_features(self, data: pd.DataFrame, window_size: int = 50) -> pd.DataFrame:
        """
        Extract statistical and dynamic features from time series data.
        
        Args:
            data: DataFrame with time series data
            window_size: Window size for rolling statistics
            
        Returns:
            DataFrame with extracted features
        """
        features = pd.DataFrame()
        
        # Remove time column for feature extraction
        if 'Time (s)' in data.columns:
            time_col = data['Time (s)']
            data_no_time = data.drop('Time (s)', axis=1)
        else:
            data_no_time = data
            time_col = np.arange(len(data))
            
        features['time'] = time_col
        
        # Ensure we have enough data for window operations
        if len(data_no_time) < window_size:
            window_size = max(3, len(data_no_time) // 3)
            
        for col in data_no_time.columns:
            series = data_no_time[col]
            
            # Basic rolling statistics
            try:
                features[f'{col}_mean'] = series.rolling(window=window_size, center=True, min_periods=1).mean()
                features[f'{col}_std'] = series.rolling(window=window_size, center=True, min_periods=1).std()
                features[f'{col}_var'] = series.rolling(window=window_size, center=True, min_periods=1).var()
            except Exception as e:
                print(f"Warning: Error calculating rolling statistics for {col}: {e}")
                features[f'{col}_mean'] = series
                features[f'{col}_std'] = 0
                features[f'{col}_var'] = 0
            
            # Derivatives (rate of change)
            try:
                features[f'{col}_derivative'] = np.gradient(series)
                features[f'{col}_derivative2'] = np.gradient(features[f'{col}_derivative'])
            except Exception as e:
                print(f"Warning: Error calculating derivatives for {col}: {e}")
                features[f'{col}_derivative'] = 0
                features[f'{col}_derivative2'] = 0
            
        # Cross-correlation features (simplified)
        species_cols = list(data_no_time.columns)
        if len(species_cols) > 1:
            for i, col1 in enumerate(species_cols):
                for col2 in species_cols[i+1:]:
                    try:
                        corr = data_no_time[col1].rolling(window=window_size, center=True, min_periods=1).corr(data_no_time[col2])
                        features[f'{col1}_{col2}_corr'] = corr
                    except Exception as e:
                        print(f"Warning: Error calculating correlation between {col1} and {col2}: {e}")
                        features[f'{col1}_{col2}_corr'] = 0
                        
        # Fill NaN values
        features = features.fillna(method='ffill').fillna(method='bfill').fillna(0)
        
        return features
        
    def detect_phase_transitions(self, contamination: float = 0.1) -> Dict[str, List[int]]:
        """
        Detect phase transitions using anomaly detection.
        
        Args:
            contamination: Fraction of outliers in the data
            
        Returns:
            Dictionary with detected transition points
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        features = self.extract_features(coverages_df)
        
        # Remove non-numeric columns
        numeric_features = features.select_dtypes(include=[np.number])
        numeric_features = numeric_features.dropna()
        
        if len(numeric_features) == 0:
            raise ValueError("No numeric features available for analysis")
            
        # Scale features
        features_scaled = self.scaler.fit_transform(numeric_features)
        
        # Apply multiple anomaly detection methods
        detectors = {
            'isolation_forest': IsolationForest(contamination=contamination, random_state=42),
            'one_class_svm': OneClassSVM(nu=contamination),
            'dbscan': DBSCAN(eps=0.5, min_samples=5)
        }
        
        transition_points = {}
        
        for name, detector in detectors.items():
            if name == 'dbscan':
                labels = detector.fit_predict(features_scaled)
                # Consider points in small clusters as anomalies
                unique_labels, counts = np.unique(labels, return_counts=True)
                small_clusters = unique_labels[counts < len(features_scaled) * 0.05]  # Less than 5% of data
                anomalies = np.isin(labels, small_clusters)
                outlier_indices = np.where(anomalies)[0]
            else:
                outliers = detector.fit_predict(features_scaled)
                outlier_indices = np.where(outliers == -1)[0]
                
            transition_points[name] = outlier_indices.tolist()
            self.models[name] = detector
            
        # Consensus-based transitions (detected by multiple methods)
        all_transitions = set()
        for transitions in transition_points.values():
            all_transitions.update(transitions)
            
        consensus_transitions = []
        for transition in all_transitions:
            count = sum(1 for transitions in transition_points.values() if transition in transitions)
            if count >= 2:  # Detected by at least 2 methods
                consensus_transitions.append(transition)
                
        transition_points['consensus'] = consensus_transitions
        self.transition_points = transition_points
        
        return transition_points
        
    def detect_steady_state(self, variance_threshold: float = 1e-6) -> Dict[str, Union[int, float]]:
        """
        Detect when the system reaches steady state.
        
        Args:
            variance_threshold: Threshold for considering system at steady state
            
        Returns:
            Dictionary with steady state information
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        
        # Calculate rolling variance for each species
        window_size = min(100, len(coverages_df) // 10)
        
        steady_state_info = {}
        
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
                
            rolling_var = coverages_df[species].rolling(window=window_size, center=True).var()
            
            # Find first point where variance drops below threshold and stays there
            below_threshold = rolling_var < variance_threshold
            
            if below_threshold.any():
                # Find the first index where it goes below threshold and stays there for at least 50 points
                for i in range(len(below_threshold) - 50):
                    if below_threshold.iloc[i:i+50].all():
                        steady_state_info[f'{species}_steady_state_index'] = i
                        steady_state_info[f'{species}_steady_state_time'] = coverages_df['Time (s)'].iloc[i]
                        break
                        
        # Overall system steady state (when all species are at steady state)
        species_steady_times = [v for k, v in steady_state_info.items() if k.endswith('_steady_state_time')]
        if species_steady_times:
            steady_state_info['system_steady_state_time'] = max(species_steady_times)
            
        return steady_state_info
        
    def detect_oscillations(self, min_period: int = 10) -> Dict[str, Dict]:
        """
        Detect oscillatory behavior in species coverages.
        
        Args:
            min_period: Minimum period for oscillation detection
            
        Returns:
            Dictionary with oscillation information for each species
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        oscillation_info = {}
        
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
                
            series = coverages_df[species].values
            time = coverages_df['Time (s)'].values
            
            # Detrend the signal
            detrended = signal.detrend(series)
            
            # Apply FFT to find dominant frequencies
            fft = np.fft.fft(detrended)
            freqs = np.fft.fftfreq(len(series), d=np.mean(np.diff(time)))
            
            # Find peaks in the power spectrum
            power_spectrum = np.abs(fft) ** 2
            peaks, properties = find_peaks(power_spectrum[:len(power_spectrum)//2], 
                                         height=np.max(power_spectrum) * 0.1)
            
            if len(peaks) > 0:
                # Get the dominant frequency
                dominant_freq_idx = peaks[np.argmax(power_spectrum[peaks])]
                dominant_freq = freqs[dominant_freq_idx]
                
                if dominant_freq > 0:
                    period = 1 / dominant_freq
                    
                    # Check if the period is reasonable
                    if period >= min_period and period <= len(series) / 3:
                        # Calculate autocorrelation to confirm periodicity
                        autocorr = np.correlate(detrended, detrended, mode='full')
                        autocorr = autocorr[autocorr.size // 2:]
                        autocorr = autocorr / autocorr[0]  # Normalize
                        
                        # Find peaks in autocorrelation
                        autocorr_peaks, _ = find_peaks(autocorr[1:], height=0.3)
                        
                        if len(autocorr_peaks) > 0:
                            estimated_period = autocorr_peaks[0] + 1
                            
                            oscillation_info[species] = {
                                'has_oscillations': True,
                                'dominant_frequency': dominant_freq,
                                'period': period,
                                'estimated_period_points': estimated_period,
                                'autocorr_peak_height': autocorr[autocorr_peaks[0] + 1] if len(autocorr_peaks) > 0 else 0,
                                'amplitude': np.std(detrended)
                            }
                        else:
                            oscillation_info[species] = {'has_oscillations': False}
                    else:
                        oscillation_info[species] = {'has_oscillations': False}
                else:
                    oscillation_info[species] = {'has_oscillations': False}
            else:
                oscillation_info[species] = {'has_oscillations': False}
                
        return oscillation_info
        
    def predict_poisoning(self, poisoning_species: List[str], 
                         training_fraction: float = 0.7) -> Dict[str, Union[float, np.ndarray]]:
        """
        Train a model to predict catalyst poisoning based on coverage patterns.
        
        Args:
            poisoning_species: List of species that cause poisoning
            training_fraction: Fraction of data to use for training
            
        Returns:
            Dictionary with model performance and predictions
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        coverages_df = self.results.get_coverages_df()
        features = self.extract_features(coverages_df)
        
        # Create labels based on poisoning species coverage
        poisoning_threshold = 0.8  # Threshold for considering the surface poisoned
        
        poisoning_coverages = []
        for species in poisoning_species:
            if species in coverages_df.columns:
                poisoning_coverages.append(coverages_df[species])
                
        if not poisoning_coverages:
            raise ValueError("None of the specified poisoning species found in data")
            
        # Label as poisoned if any poisoning species exceeds threshold
        total_poisoning_coverage = sum(poisoning_coverages)
        labels = (total_poisoning_coverage > poisoning_threshold).astype(int)
        
        # Prepare features
        numeric_features = features.select_dtypes(include=[np.number]).fillna(0)
        
        # Split data
        split_index = int(len(numeric_features) * training_fraction)
        X_train, X_test = numeric_features.iloc[:split_index], numeric_features.iloc[split_index:]
        y_train, y_test = labels.iloc[:split_index], labels.iloc[split_index:]
        
        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)
        
        # Train multiple models
        models = {
            'random_forest': RandomForestClassifier(n_estimators=100, random_state=42),
            'mlp': MLPClassifier(hidden_layer_sizes=(100, 50), max_iter=500, random_state=42)
        }
        
        results = {}
        
        for name, model in models.items():
            # Train model
            model.fit(X_train_scaled, y_train)
            
            # Make predictions
            y_pred = model.predict(X_test_scaled)
            y_pred_proba = model.predict_proba(X_test_scaled)[:, 1] if hasattr(model, 'predict_proba') else None
            
            # Calculate metrics
            from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
            
            results[name] = {
                'model': model,
                'accuracy': accuracy_score(y_test, y_pred),
                'precision': precision_score(y_test, y_pred, zero_division=0),
                'recall': recall_score(y_test, y_pred, zero_division=0),
                'f1_score': f1_score(y_test, y_pred, zero_division=0),
                'predictions': y_pred,
                'probabilities': y_pred_proba,
                'test_labels': y_test.values
            }
            
        self.models.update({f'poisoning_{k}': v['model'] for k, v in results.items()})
        
        return results
        
    def visualize_transitions(self, save_path: Optional[str] = None) -> plt.Figure:
        """
        Visualize detected transition points on coverage plots.
        
        Args:
            save_path: Optional path to save the figure
            
        Returns:
            Matplotlib figure object
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        if not self.transition_points:
            self.detect_phase_transitions()
            
        coverages_df = self.results.get_coverages_df()
        
        fig, axes = plt.subplots(2, 1, figsize=(15, 10))
        
        # Plot 1: Coverage with transition points
        for species in coverages_df.columns:
            if species == 'Time (s)':
                continue
            axes[0].plot(coverages_df['Time (s)'], coverages_df[species], 
                        label=species, alpha=0.8, linewidth=2)
            
        # Mark transition points
        consensus_transitions = self.transition_points.get('consensus', [])
        for transition in consensus_transitions:
            if transition < len(coverages_df):
                axes[0].axvline(x=coverages_df['Time (s)'].iloc[transition], 
                               color='red', linestyle='--', alpha=0.7)
                
        axes[0].set_xlabel("Time (s)")
        axes[0].set_ylabel("Coverage")
        axes[0].set_title("Species Coverage with Detected Transition Points")
        axes[0].legend()
        axes[0].grid(True, alpha=0.3)
        
        # Plot 2: Transition detection summary
        method_names = [name for name in self.transition_points.keys() if name != 'consensus']
        transition_counts = [len(self.transition_points[name]) for name in method_names]
        
        bars = axes[1].bar(method_names, transition_counts, alpha=0.7)
        axes[1].set_ylabel("Number of Transitions Detected")
        axes[1].set_title("Transition Detection Summary by Method")
        axes[1].grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, count in zip(bars, transition_counts):
            axes[1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                        str(count), ha='center', va='bottom')
            
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"Transition visualization saved to {save_path}")
            
        return fig
        
    def save_models(self, directory: Union[str, Path]):
        """
        Save trained models to disk.
        
        Args:
            directory: Directory to save models
        """
        directory = Path(directory)
        directory.mkdir(exist_ok=True)
        
        for name, model in self.models.items():
            model_path = directory / f"{name}_model.pkl"
            joblib.dump(model, model_path)
            print(f"Model {name} saved to {model_path}")
            
        # Save scaler
        scaler_path = directory / "scaler.pkl"
        joblib.dump(self.scaler, scaler_path)
        print(f"Scaler saved to {scaler_path}")
        
    def load_models(self, directory: Union[str, Path]):
        """
        Load trained models from disk.
        
        Args:
            directory: Directory containing saved models
        """
        directory = Path(directory)
        
        # Load scaler
        scaler_path = directory / "scaler.pkl"
        if scaler_path.exists():
            self.scaler = joblib.load(scaler_path)
            print(f"Scaler loaded from {scaler_path}")
            
        # Load models
        for model_file in directory.glob("*_model.pkl"):
            model_name = model_file.stem.replace("_model", "")
            self.models[model_name] = joblib.load(model_file)
            print(f"Model {model_name} loaded from {model_file}")
            
    def generate_report(self, save_path: Optional[str] = None) -> str:
        """
        Generate a comprehensive analysis report.
        
        Args:
            save_path: Optional path to save the report
            
        Returns:
            Report as string
        """
        if not self.results:
            raise ValueError("No results object provided")
            
        # Run all analyses
        transitions = self.detect_phase_transitions()
        steady_state = self.detect_steady_state()
        oscillations = self.detect_oscillations()
        
        report_lines = [
            "="*60,
            "PyApothesis KMC Simulation Analysis Report",
            "="*60,
            "",
            "1. PHASE TRANSITION ANALYSIS",
            "-"*30,
        ]
        
        for method, points in transitions.items():
            report_lines.append(f"{method.upper()}: {len(points)} transitions detected")
            if points:
                report_lines.append(f"  Transition indices: {points}")
                
        report_lines.extend([
            "",
            "2. STEADY STATE ANALYSIS", 
            "-"*30,
        ])
        
        if steady_state:
            for key, value in steady_state.items():
                if isinstance(value, float):
                    report_lines.append(f"{key}: {value:.4f}")
                else:
                    report_lines.append(f"{key}: {value}")
        else:
            report_lines.append("No steady state detected")
            
        report_lines.extend([
            "",
            "3. OSCILLATION ANALYSIS",
            "-"*30,
        ])
        
        for species, info in oscillations.items():
            if info.get('has_oscillations', False):
                report_lines.append(f"{species}: OSCILLATING")
                report_lines.append(f"  Period: {info['period']:.2f} time units")
                report_lines.append(f"  Amplitude: {info['amplitude']:.4f}")
            else:
                report_lines.append(f"{species}: No oscillations detected")
                
        report_lines.extend([
            "",
            "="*60,
            "End of Report",
            "="*60
        ])
        
        report = "\n".join(report_lines)
        
        if save_path:
            with open(save_path, 'w') as f:
                f.write(report)
            print(f"Report saved to {save_path}")
            
        return report

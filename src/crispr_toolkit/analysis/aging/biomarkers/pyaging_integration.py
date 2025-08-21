"""
PyAging Integration Module

Provides a comprehensive Python wrapper for the pyaging library, enabling
easy access to 100+ aging clocks with specific focus on ReHMGB1/RAGE signaling
and senescence pathways.

PyAging provides GPU-optimized aging clocks covering:
- DNA methylation (Horvath, Hannum, PhenoAge, GrimAge, etc.)
- Transcriptomics (GTEx-based clocks, tissue-specific clocks)
- Histone mark ChIP-Seq (H3K4me3, H3K27me3, etc.)
- ATAC-Seq chromatin accessibility
- Multi-modal aging clocks

Author: CRISPR Toolkit Development Team
"""

import os
import warnings
from typing import Any, Dict, List, Optional, Union

import pandas as pd

# Handle optional pyaging import
try:
    import pyaging as pya
    PYAGING_AVAILABLE = True
except ImportError:
    PYAGING_AVAILABLE = False
    warnings.warn("pyaging not installed. Install with: pip install pyaging")

try:
    import torch
    TORCH_AVAILABLE = True
except ImportError:
    TORCH_AVAILABLE = False


class PyAgingAnalyzer:
    """
    Comprehensive wrapper for pyaging library with ReHMGB1/senescence focus.

    Provides easy access to 100+ aging clocks with specialized analysis for
    aging and senescence pathways, particularly ReHMGB1/RAGE signaling.
    """

    def __init__(self, device: Optional[str] = None, verbose: bool = True):
        """
        Initialize PyAging analyzer.

        Args:
            device: PyTorch device ('cpu', 'cuda', 'mps'). Auto-detected if None.
            verbose: Enable verbose logging
        """
        if not PYAGING_AVAILABLE:
            raise ImportError(
                "pyaging is required but not installed. "
                "Install with: pip install pyaging"
            )

        self.verbose = verbose
        self.device = self._setup_device(device)
        self.available_clocks = self._load_available_clocks()
        self.senescence_clocks = self._identify_senescence_clocks()

        if self.verbose:
            print(f"PyAging initialized with {len(self.available_clocks)} "
                  f"available clocks on device: {self.device}")

    def _setup_device(self, device: Optional[str]) -> str:
        """Setup PyTorch device for computation."""
        if device is not None:
            return device

        if TORCH_AVAILABLE:
            if torch.cuda.is_available():
                return 'cuda'
            elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
                return 'mps'

        return 'cpu'

    def _load_available_clocks(self) -> Dict[str, Dict]:
        """Load metadata for all available aging clocks."""
        try:
            # Get all available clocks from pyaging
            clock_list = pya.utils.get_all_clock_names()

            clocks_metadata = {}
            for clock_name in clock_list:
                try:
                    metadata = pya.utils.get_clock_metadata(clock_name)
                    clocks_metadata[clock_name] = {
                        'name': clock_name,
                        'species': metadata.get('species', 'unknown'),
                        'tissue': metadata.get('tissue', 'unknown'),
                        'data_type': metadata.get('data_type', 'unknown'),
                        'year': metadata.get('year', 'unknown'),
                        'pmid': metadata.get('pmid', 'unknown'),
                        'doi': metadata.get('doi', 'unknown')
                    }
                except Exception as e:
                    if self.verbose:
                        print(f"Warning: Could not load metadata for {clock_name}: {e}")
                    clocks_metadata[clock_name] = {
                        'name': clock_name,
                        'species': 'unknown',
                        'tissue': 'unknown',
                        'data_type': 'unknown'
                    }

            return clocks_metadata

        except Exception as e:
            if self.verbose:
                print(f"Warning: Could not load clock list: {e}")
            return {}

    def _identify_senescence_clocks(self) -> Dict[str, List[str]]:
        """Identify aging clocks most relevant to senescence research."""

        senescence_relevant = {
            'dna_methylation_senescence': [
                'horvath2013',
                'hannum2013',
                'phenoage',
                'grimage',
                'dunedinpoam38',
                'dunedinpace',
                'mpoa',
                'zhang2019',
                'lin2016'
            ],
            'transcriptomic_senescence': [
                'mammalianlifespan',
                'mammalianfemalelifespan',
                'mammalianmalelifespan',
                'gtexaging',
                'gtexagingmale',
                'gtexagingfemale'
            ],
            'composite_senescence': [
                'phenoage',
                'grimage',
                'grimage2',
                'dunedinpace'
            ],
            'inflammatory_aging': [
                'inflammage',
                'grimage',
                'phenoage'
            ]
        }

        # Filter to only include clocks that are actually available
        available_senescence = {}
        for category, clock_list in senescence_relevant.items():
            available_clocks = [
                clock for clock in clock_list
                if clock in self.available_clocks
            ]
            if available_clocks:
                available_senescence[category] = available_clocks

        return available_senescence

    def get_clock_info(self, clock_name: str) -> Dict:
        """Get detailed information about a specific aging clock."""
        if clock_name not in self.available_clocks:
            raise ValueError(f"Clock '{clock_name}' not found in available clocks")

        return self.available_clocks[clock_name]

    def list_clocks_by_datatype(self, data_type: str) -> List[str]:
        """List all clocks that work with a specific data type."""
        return [
            name for name, info in self.available_clocks.items()
            if info.get('data_type', '').lower() == data_type.lower()
        ]

    def list_senescence_clocks(self, category: Optional[str] = None) -> Union[Dict, List]:
        """
        List aging clocks most relevant to senescence research.

        Args:
            category: Specific category ('dna_methylation_senescence',
                     'transcriptomic_senescence', 'composite_senescence',
                     'inflammatory_aging'). If None, returns all categories.

        Returns:
            Dictionary of categories or list of clocks if category specified
        """
        if category is None:
            return self.senescence_clocks

        if category not in self.senescence_clocks:
            available_cats = list(self.senescence_clocks.keys())
            raise ValueError(f"Category '{category}' not found. "
                           f"Available: {available_cats}")

        return self.senescence_clocks[category]

    def predict_age(
        self,
        data: pd.DataFrame,
        clock_names: Union[str, List[str]],
        sample_column: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Predict biological age using specified aging clocks.

        Args:
            data: Input data (samples x features)
            clock_names: Clock name(s) to use for prediction
            sample_column: Column name for sample IDs. If None, uses index.

        Returns:
            DataFrame with age predictions for each clock and sample
        """
        if isinstance(clock_names, str):
            clock_names = [clock_names]

        # Validate clocks
        for clock_name in clock_names:
            if clock_name not in self.available_clocks:
                available = list(self.available_clocks.keys())[:10]
                raise ValueError(f"Clock '{clock_name}' not available. "
                               f"Available clocks (first 10): {available}")

        results = []

        for clock_name in clock_names:
            try:
                if self.verbose:
                    print(f"Predicting age with {clock_name}...")

                # Load the specific clock
                clock = pya.age(data, clock_name, device=self.device)

                # Extract predictions
                predictions = clock.predicted_age

                # Create result dataframe
                if sample_column and sample_column in data.columns:
                    sample_ids = data[sample_column]
                else:
                    sample_ids = data.index

                result_df = pd.DataFrame({
                    'sample_id': sample_ids,
                    'clock': clock_name,
                    'predicted_age': predictions,
                    'clock_species': self.available_clocks[clock_name].get('species', 'unknown'),
                    'clock_tissue': self.available_clocks[clock_name].get('tissue', 'unknown'),
                    'data_type': self.available_clocks[clock_name].get('data_type', 'unknown')
                })

                results.append(result_df)

            except Exception as e:
                if self.verbose:
                    print(f"Warning: Failed to predict with {clock_name}: {e}")
                continue

        if not results:
            raise RuntimeError("No successful age predictions generated")

        # Combine results
        combined_results = pd.concat(results, ignore_index=True)

        return combined_results

    def analyze_senescence_aging(
        self,
        data: pd.DataFrame,
        chronological_age: Optional[pd.Series] = None,
        sample_column: Optional[str] = None
    ) -> Dict[str, Any]:
        """
        Comprehensive senescence-focused aging analysis.

        Args:
            data: Input omics data
            chronological_age: Chronological ages for comparison
            sample_column: Column name for sample IDs

        Returns:
            Comprehensive analysis results including ReHMGB1-relevant metrics
        """
        results = {
            'senescence_clocks_used': [],
            'age_predictions': None,
            'age_acceleration': None,
            'senescence_summary': {},
            'rehmgb1_relevant_metrics': {}
        }

        # Get all available senescence clocks
        all_senescence_clocks = []
        for category, clocks in self.senescence_clocks.items():
            all_senescence_clocks.extend(clocks)

        # Remove duplicates
        all_senescence_clocks = list(set(all_senescence_clocks))

        if not all_senescence_clocks:
            raise ValueError("No senescence-relevant clocks available")

        # Predict ages with senescence clocks
        try:
            age_predictions = self.predict_age(
                data, all_senescence_clocks, sample_column
            )
            results['age_predictions'] = age_predictions
            results['senescence_clocks_used'] = all_senescence_clocks
        except Exception as e:
            raise RuntimeError(f"Failed to predict ages: {e}")

        # Calculate age acceleration if chronological age provided
        if chronological_age is not None:
            age_acceleration = self._calculate_age_acceleration(
                age_predictions, chronological_age, sample_column
            )
            results['age_acceleration'] = age_acceleration

        # Generate senescence summary statistics
        results['senescence_summary'] = self._generate_senescence_summary(
            age_predictions, chronological_age
        )

        # Calculate ReHMGB1-relevant metrics
        results['rehmgb1_relevant_metrics'] = self._calculate_rehmgb1_metrics(
            age_predictions, chronological_age
        )

        return results

    def _calculate_age_acceleration(
        self,
        age_predictions: pd.DataFrame,
        chronological_age: pd.Series,
        sample_column: Optional[str]
    ) -> pd.DataFrame:
        """Calculate age acceleration (predicted - chronological age)."""

        # Prepare chronological age mapping
        if sample_column and sample_column in chronological_age.index.names:
            chron_age_map = chronological_age.to_dict()
        else:
            chron_age_map = dict(zip(chronological_age.index, chronological_age.values))

        # Calculate age acceleration
        age_predictions['chronological_age'] = age_predictions['sample_id'].map(chron_age_map)
        age_predictions['age_acceleration'] = (
            age_predictions['predicted_age'] - age_predictions['chronological_age']
        )

        return age_predictions

    def _generate_senescence_summary(
        self,
        age_predictions: pd.DataFrame,
        chronological_age: Optional[pd.Series]
    ) -> Dict:
        """Generate summary statistics for senescence analysis."""

        summary = {
            'total_samples': age_predictions['sample_id'].nunique(),
            'clocks_used': age_predictions['clock'].nunique(),
            'clock_list': age_predictions['clock'].unique().tolist()
        }

        # Clock-specific summaries
        clock_summaries = {}
        for clock in age_predictions['clock'].unique():
            clock_data = age_predictions[age_predictions['clock'] == clock]
            clock_summaries[clock] = {
                'mean_predicted_age': clock_data['predicted_age'].mean(),
                'std_predicted_age': clock_data['predicted_age'].std(),
                'min_predicted_age': clock_data['predicted_age'].min(),
                'max_predicted_age': clock_data['predicted_age'].max()
            }

            if 'age_acceleration' in clock_data.columns:
                clock_summaries[clock].update({
                    'mean_age_acceleration': clock_data['age_acceleration'].mean(),
                    'std_age_acceleration': clock_data['age_acceleration'].std()
                })

        summary['clock_summaries'] = clock_summaries

        return summary

    def _calculate_rehmgb1_metrics(
        self,
        age_predictions: pd.DataFrame,
        chronological_age: Optional[pd.Series]
    ) -> Dict:
        """Calculate metrics specifically relevant to ReHMGB1 research."""

        metrics = {
            'inflammatory_aging_clocks': [],
            'dna_methylation_senescence_score': None,
            'senescence_acceleration_risk': None
        }

        # Identify inflammatory aging clocks used
        inflammatory_clocks = self.senescence_clocks.get('inflammatory_aging', [])
        used_inflammatory = [
            clock for clock in age_predictions['clock'].unique()
            if clock in inflammatory_clocks
        ]
        metrics['inflammatory_aging_clocks'] = used_inflammatory

        # Calculate DNA methylation senescence score
        dna_meth_clocks = self.senescence_clocks.get('dna_methylation_senescence', [])
        dna_meth_data = age_predictions[
            age_predictions['clock'].isin(dna_meth_clocks)
        ]

        if not dna_meth_data.empty and 'age_acceleration' in dna_meth_data.columns:
            # Average age acceleration across DNA methylation clocks
            sample_scores = dna_meth_data.groupby('sample_id')['age_acceleration'].mean()
            metrics['dna_methylation_senescence_score'] = {
                'mean': sample_scores.mean(),
                'std': sample_scores.std(),
                'samples_with_acceleration': (sample_scores > 0).sum(),
                'samples_with_deceleration': (sample_scores < 0).sum()
            }

        # Calculate senescence acceleration risk
        if 'age_acceleration' in age_predictions.columns:
            # Samples with consistent positive age acceleration across clocks
            sample_acceleration = age_predictions.groupby('sample_id')['age_acceleration'].agg([
                'mean', 'std', 'count'
            ])

            # High risk: positive mean acceleration with multiple clocks
            high_risk_samples = sample_acceleration[
                (sample_acceleration['mean'] > 5) &  # >5 years acceleration
                (sample_acceleration['count'] >= 2)   # At least 2 clocks
            ]

            metrics['senescence_acceleration_risk'] = {
                'high_risk_samples': len(high_risk_samples),
                'total_samples': len(sample_acceleration),
                'high_risk_percentage': (len(high_risk_samples) / len(sample_acceleration)) * 100
            }

        return metrics

    def export_results(
        self,
        results: Dict,
        output_dir: str,
        prefix: str = 'pyaging_analysis'
    ) -> Dict[str, str]:
        """
        Export analysis results to files.

        Args:
            results: Results from analyze_senescence_aging()
            output_dir: Output directory
            prefix: File prefix

        Returns:
            Dictionary mapping result types to file paths
        """
        os.makedirs(output_dir, exist_ok=True)

        exported_files = {}

        # Export age predictions
        if results.get('age_predictions') is not None:
            age_pred_file = os.path.join(output_dir, f'{prefix}_age_predictions.csv')
            results['age_predictions'].to_csv(age_pred_file, index=False)
            exported_files['age_predictions'] = age_pred_file

        # Export summary as JSON
        import json
        summary_file = os.path.join(output_dir, f'{prefix}_summary.json')

        # Prepare JSON-serializable summary
        json_summary = {
            'senescence_clocks_used': results.get('senescence_clocks_used', []),
            'senescence_summary': results.get('senescence_summary', {}),
            'rehmgb1_relevant_metrics': results.get('rehmgb1_relevant_metrics', {})
        }

        with open(summary_file, 'w') as f:
            json.dump(json_summary, f, indent=2, default=str)
        exported_files['summary'] = summary_file

        if self.verbose:
            print(f"Results exported to {output_dir}")
            for result_type, file_path in exported_files.items():
                print(f"  {result_type}: {file_path}")

        return exported_files


def create_demo_pyaging_analysis() -> str:
    """Create a demonstration script for PyAging integration."""

    demo_script = '''#!/usr/bin/env python3
"""
PyAging Integration Demo for ReHMGB1 Senescence Research

This script demonstrates how to use the PyAging integration for analyzing
aging biomarkers with focus on ReHMGB1/RAGE signaling pathways.
"""

import pandas as pd
import numpy as np
from crispr_toolkit.analysis.aging.biomarkers import PyAgingAnalyzer

def main():
    print("🧬 PyAging Integration Demo - ReHMGB1 Senescence Research")
    print("=" * 60)

    # Initialize PyAging analyzer
    print("\\n1. Initializing PyAging analyzer...")
    try:
        pyaging = PyAgingAnalyzer(verbose=True)
        print(f"   ✅ Successfully initialized with {len(pyaging.available_clocks)} clocks")
    except ImportError as e:
        print(f"   ❌ PyAging not available: {e}")
        print("   📦 Install with: pip install pyaging")
        return

    # List senescence-relevant clocks
    print("\\n2. Available senescence-relevant aging clocks:")
    senescence_clocks = pyaging.list_senescence_clocks()
    for category, clocks in senescence_clocks.items():
        print(f"   📊 {category}: {len(clocks)} clocks")
        for clock in clocks[:3]:  # Show first 3
            info = pyaging.get_clock_info(clock)
            print(f"      - {clock} ({info.get('data_type', 'unknown')} data)")

    # Generate synthetic methylation data for demo
    print("\\n3. Generating synthetic DNA methylation data...")
    n_samples = 50
    n_cpgs = 1000

    # Create synthetic data that mimics aging patterns
    np.random.seed(42)
    age_range = np.random.uniform(20, 80, n_samples)

    # Generate methylation data with age-related patterns
    methylation_data = np.random.beta(2, 2, (n_samples, n_cpgs))

    # Add age-related trend to some CpGs
    for i in range(0, n_cpgs, 10):  # Every 10th CpG shows aging pattern
        age_effect = (age_range - 50) / 50 * 0.2  # Small age effect
        methylation_data[:, i] += age_effect
        methylation_data[:, i] = np.clip(methylation_data[:, i], 0, 1)

    # Create DataFrame
    cpg_names = [f"cg{i:07d}" for i in range(n_cpgs)]
    sample_names = [f"Sample_{i:03d}" for i in range(n_samples)]

    data = pd.DataFrame(
        methylation_data,
        index=sample_names,
        columns=cpg_names
    )

    chronological_age = pd.Series(age_range, index=sample_names)

    print(f"   📊 Generated data: {data.shape[0]} samples x {data.shape[1]} CpGs")
    print(f"   📅 Age range: {chronological_age.min():.1f} - {chronological_age.max():.1f} years")

    # Perform senescence aging analysis
    print("\\n4. Performing senescence-focused aging analysis...")
    try:
        results = pyaging.analyze_senescence_aging(
            data=data,
            chronological_age=chronological_age
        )

        print("   ✅ Analysis completed successfully!")

        # Display results
        print(f"   📊 Clocks used: {len(results['senescence_clocks_used'])}")
        print(f"   📈 Samples analyzed: {results['senescence_summary']['total_samples']}")

        # Show ReHMGB1-relevant metrics
        rehmgb1_metrics = results['rehmgb1_relevant_metrics']
        print("\\n5. ReHMGB1-Relevant Aging Metrics:")

        if rehmgb1_metrics.get('inflammatory_aging_clocks'):
            print(f"   🔥 Inflammatory aging clocks: {rehmgb1_metrics['inflammatory_aging_clocks']}")

        if rehmgb1_metrics.get('dna_methylation_senescence_score'):
            score = rehmgb1_metrics['dna_methylation_senescence_score']
            print(f"   🧬 DNA methylation senescence score: {score['mean']:.2f} ± {score['std']:.2f}")
            print(f"   📈 Samples with age acceleration: {score['samples_with_acceleration']}")
            print(f"   📉 Samples with age deceleration: {score['samples_with_deceleration']}")

        if rehmgb1_metrics.get('senescence_acceleration_risk'):
            risk = rehmgb1_metrics['senescence_acceleration_risk']
            print(f"   ⚠️  High senescence risk samples: {risk['high_risk_samples']} "
                  f"({risk['high_risk_percentage']:.1f}%)")

        # Export results
        print("\\n6. Exporting results...")
        exported = pyaging.export_results(results, './pyaging_demo_results')
        print("   ✅ Results exported successfully!")

    except Exception as e:
        print(f"   ❌ Analysis failed: {e}")
        print("   💡 This is expected with synthetic data - real data needed for actual clocks")

    print("\\n🎉 PyAging integration demo completed!")
    print("\\n💡 Next steps for ReHMGB1 research:")
    print("   - Load real DNA methylation data (450K, EPIC arrays)")
    print("   - Use senescence-specific clocks (PhenoAge, GrimAge)")
    print("   - Correlate age acceleration with ReHMGB1 expression")
    print("   - Analyze RAGE pathway methylation patterns")

if __name__ == "__main__":
    main()
'''

    return demo_script
    return demo_script

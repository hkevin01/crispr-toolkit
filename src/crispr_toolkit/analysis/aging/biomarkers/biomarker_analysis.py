"""
Combined Aging Biomarker Analysis Module

Integrates pyaging and biolearn libraries to provide comprehensive aging
biomarker analysis with specific focus on ReHMGB1/RAGE signaling and
senescence research applications.

This module combines the strengths of both libraries:
- PyAging: 100+ GPU-optimized aging clocks across multiple data types
- Biolearn: Standardized implementations and clinical validation

Author: CRISPR Toolkit Development Team
"""

import os
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from .biolearn_integration import BIOLEARN_AVAILABLE, BiolearnAnalyzer
from .pyaging_integration import PYAGING_AVAILABLE, PyAgingAnalyzer


class AgingBiomarkerAnalyzer:
    """
    Comprehensive aging biomarker analysis combining pyaging and biolearn.

    Provides unified interface for multi-library aging analysis with
    specialized focus on ReHMGB1/RAGE signaling and senescence pathways.
    """

    def __init__(
        self,
        use_pyaging: bool = True,
        use_biolearn: bool = True,
        device: Optional[str] = None,
        cache_dir: Optional[str] = None,
        verbose: bool = True
    ):
        """
        Initialize combined aging biomarker analyzer.

        Args:
            use_pyaging: Enable pyaging integration
            use_biolearn: Enable biolearn integration
            device: PyTorch device for pyaging
            cache_dir: Cache directory for biolearn
            verbose: Enable verbose logging
        """
        self.verbose = verbose
        self.pyaging_analyzer = None
        self.biolearn_analyzer = None

        # Initialize pyaging if requested and available
        if use_pyaging:
            if PYAGING_AVAILABLE:
                try:
                    self.pyaging_analyzer = PyAgingAnalyzer(
                        device=device, verbose=verbose
                    )
                    if verbose:
                        print("✅ PyAging integration enabled")
                except Exception as e:
                    if verbose:
                        print(f"⚠️  PyAging initialization failed: {e}")
            else:
                if verbose:
                    print("⚠️  PyAging not available - install with: pip install pyaging")

        # Initialize biolearn if requested and available
        if use_biolearn:
            if BIOLEARN_AVAILABLE:
                try:
                    self.biolearn_analyzer = BiolearnAnalyzer(
                        cache_dir=cache_dir, verbose=verbose
                    )
                    if verbose:
                        print("✅ Biolearn integration enabled")
                except Exception as e:
                    if verbose:
                        print(f"⚠️  Biolearn initialization failed: {e}")
            else:
                if verbose:
                    print("⚠️  Biolearn not available - install with: pip install biolearn")

        # Check that at least one library is available
        if self.pyaging_analyzer is None and self.biolearn_analyzer is None:
            raise RuntimeError(
                "Neither pyaging nor biolearn could be initialized. "
                "Install at least one: pip install pyaging biolearn"
            )

        # Load combined clock information
        self.available_clocks = self._combine_clock_info()
        self.senescence_clocks = self._combine_senescence_clocks()

        if verbose:
            total_clocks = len(self.available_clocks)
            print(f"🧬 Combined analyzer initialized with {total_clocks} total clocks")

    def _combine_clock_info(self) -> Dict[str, Dict]:
        """Combine clock information from both libraries."""
        combined_clocks = {}

        # Add pyaging clocks
        if self.pyaging_analyzer:
            for clock_name, info in self.pyaging_analyzer.available_clocks.items():
                combined_clocks[f"pyaging_{clock_name}"] = {
                    **info,
                    'library': 'pyaging'
                }

        # Add biolearn clocks
        if self.biolearn_analyzer:
            for clock_name, info in self.biolearn_analyzer.available_clocks.items():
                combined_clocks[f"biolearn_{clock_name}"] = {
                    **info,
                    'library': 'biolearn'
                }

        return combined_clocks

    def _combine_senescence_clocks(self) -> Dict[str, Dict]:
        """Combine senescence-relevant clocks from both libraries."""
        combined_senescence = {
            'pyaging_clocks': {},
            'biolearn_clocks': {},
            'consensus_clocks': []
        }

        # Add pyaging senescence clocks
        if self.pyaging_analyzer:
            combined_senescence['pyaging_clocks'] = (
                self.pyaging_analyzer.senescence_clocks
            )

        # Add biolearn senescence clocks
        if self.biolearn_analyzer:
            combined_senescence['biolearn_clocks'] = (
                self.biolearn_analyzer.senescence_clocks
            )

        # Identify consensus clocks (available in both libraries)
        if self.pyaging_analyzer and self.biolearn_analyzer:
            pyaging_all = set()
            for clocks in self.pyaging_analyzer.senescence_clocks.values():
                pyaging_all.update(clocks)

            biolearn_all = set()
            for clocks in self.biolearn_analyzer.senescence_clocks.values():
                biolearn_all.update(clocks)

            # Find common clock names (normalize case and naming)
            consensus = []
            for py_clock in pyaging_all:
                for bio_clock in biolearn_all:
                    if self._clocks_are_similar(py_clock, bio_clock):
                        consensus.append({
                            'pyaging_name': py_clock,
                            'biolearn_name': bio_clock,
                            'consensus_name': py_clock
                        })

            combined_senescence['consensus_clocks'] = consensus

        return combined_senescence

    def _clocks_are_similar(self, clock1: str, clock2: str) -> bool:
        """Check if two clock names refer to the same clock."""
        # Normalize names for comparison
        norm1 = clock1.lower().replace('_', '').replace('-', '')
        norm2 = clock2.lower().replace('_', '').replace('-', '')

        # Common name mappings
        name_mappings = {
            'horvath2013': ['horvath', 'horvath2013'],
            'hannum2013': ['hannum', 'hannum2013'],
            'phenoage': ['phenoage', 'pheno'],
            'grimage': ['grimage', 'grim'],
            'dunedinpace': ['dunedinpace', 'pace'],
            'dunedinpoam38': ['dunedinpoam38', 'poam']
        }

        for canonical, variants in name_mappings.items():
            if norm1 in variants and norm2 in variants:
                return True

        return norm1 == norm2

    def get_available_libraries(self) -> List[str]:
        """Get list of available aging analysis libraries."""
        libraries = []
        if self.pyaging_analyzer:
            libraries.append('pyaging')
        if self.biolearn_analyzer:
            libraries.append('biolearn')
        return libraries

    def list_all_clocks(self, library: Optional[str] = None) -> Dict[str, List[str]]:
        """
        List all available aging clocks.

        Args:
            library: Filter by library ('pyaging', 'biolearn'). If None, show all.

        Returns:
            Dictionary organized by library and categories
        """
        if library == 'pyaging' and self.pyaging_analyzer:
            return {
                'pyaging': {
                    'senescence_clocks': self.pyaging_analyzer.senescence_clocks,
                    'all_clocks': list(self.pyaging_analyzer.available_clocks.keys())
                }
            }
        elif library == 'biolearn' and self.biolearn_analyzer:
            return {
                'biolearn': {
                    'senescence_clocks': self.biolearn_analyzer.senescence_clocks,
                    'all_clocks': list(self.biolearn_analyzer.available_clocks.keys())
                }
            }
        else:
            result = {}
            if self.pyaging_analyzer:
                result['pyaging'] = {
                    'senescence_clocks': self.pyaging_analyzer.senescence_clocks,
                    'all_clocks': list(self.pyaging_analyzer.available_clocks.keys())
                }
            if self.biolearn_analyzer:
                result['biolearn'] = {
                    'senescence_clocks': self.biolearn_analyzer.senescence_clocks,
                    'all_clocks': list(self.biolearn_analyzer.available_clocks.keys())
                }
            return result

    def comprehensive_senescence_analysis(
        self,
        data: pd.DataFrame,
        chronological_age: Optional[pd.Series] = None,
        metadata: Optional[pd.DataFrame] = None,
        sample_column: Optional[str] = None,
        use_consensus_clocks: bool = True
    ) -> Dict[str, Any]:
        """
        Comprehensive senescence analysis using both libraries.

        Args:
            data: Input omics data (samples x features)
            chronological_age: Chronological ages for comparison
            metadata: Additional sample metadata
            sample_column: Column name for sample IDs
            use_consensus_clocks: Prioritize clocks available in both libraries

        Returns:
            Comprehensive analysis combining both libraries
        """
        results = {
            'pyaging_results': None,
            'biolearn_results': None,
            'consensus_analysis': None,
            'integrated_summary': {},
            'rehmgb1_integrated_metrics': {}
        }

        # Run pyaging analysis
        if self.pyaging_analyzer:
            try:
                if self.verbose:
                    print("🔬 Running PyAging senescence analysis...")

                pyaging_results = self.pyaging_analyzer.analyze_senescence_aging(
                    data=data,
                    chronological_age=chronological_age,
                    sample_column=sample_column
                )
                results['pyaging_results'] = pyaging_results

                if self.verbose:
                    print("   ✅ PyAging analysis completed")

            except Exception as e:
                if self.verbose:
                    print(f"   ⚠️  PyAging analysis failed: {e}")

        # Run biolearn analysis
        if self.biolearn_analyzer:
            try:
                if self.verbose:
                    print("🔬 Running Biolearn senescence analysis...")

                # Prepare chronological age column for biolearn
                if chronological_age is not None and metadata is None:
                    metadata = pd.DataFrame({'age': chronological_age})

                biolearn_results = self.biolearn_analyzer.analyze_rehmgb1_biomarkers(
                    data=data,
                    metadata=metadata,
                    chronological_age_column='age'
                )
                results['biolearn_results'] = biolearn_results

                if self.verbose:
                    print("   ✅ Biolearn analysis completed")

            except Exception as e:
                if self.verbose:
                    print(f"   ⚠️  Biolearn analysis failed: {e}")

        # Perform consensus analysis if both libraries successful
        if (results['pyaging_results'] and results['biolearn_results'] and
            use_consensus_clocks):
            try:
                if self.verbose:
                    print("🔬 Performing consensus analysis...")

                results['consensus_analysis'] = self._perform_consensus_analysis(
                    results['pyaging_results'],
                    results['biolearn_results'],
                    chronological_age
                )

                if self.verbose:
                    print("   ✅ Consensus analysis completed")

            except Exception as e:
                if self.verbose:
                    print(f"   ⚠️  Consensus analysis failed: {e}")

        # Generate integrated summary
        results['integrated_summary'] = self._generate_integrated_summary(results)

        # Calculate integrated ReHMGB1 metrics
        results['rehmgb1_integrated_metrics'] = self._calculate_integrated_rehmgb1_metrics(
            results
        )

        return results

    def _perform_consensus_analysis(
        self,
        pyaging_results: Dict,
        biolearn_results: Dict,
        chronological_age: Optional[pd.Series]
    ) -> Dict[str, Any]:
        """Perform consensus analysis between libraries."""

        consensus = {
            'consensus_clocks_used': [],
            'correlation_analysis': {},
            'agreement_metrics': {},
            'combined_risk_scores': {}
        }

        # Find consensus clocks
        consensus_clocks = self.senescence_clocks.get('consensus_clocks', [])
        if not consensus_clocks:
            return consensus

        consensus['consensus_clocks_used'] = consensus_clocks

        # Compare predictions for consensus clocks
        pyaging_preds = pyaging_results.get('age_predictions')
        biolearn_preds = biolearn_results.get('senescence_biomarkers')

        if pyaging_preds is not None and biolearn_preds is not None:
            correlations = {}

            for consensus_clock in consensus_clocks:
                py_name = consensus_clock['pyaging_name']
                bio_name = consensus_clock['biolearn_name']

                # Get predictions from each library
                py_data = pyaging_preds[pyaging_preds['clock'] == py_name]
                bio_data = biolearn_preds[biolearn_preds['clock'] == bio_name]

                if not py_data.empty and not bio_data.empty:
                    # Merge on sample ID
                    merged = py_data.merge(
                        bio_data,
                        on='sample_id',
                        suffixes=('_pyaging', '_biolearn')
                    )

                    if len(merged) > 0:
                        # Calculate correlation
                        corr = merged['predicted_age'].corr(
                            merged['predicted_value']
                        )
                        correlations[consensus_clock['consensus_name']] = {
                            'correlation': corr,
                            'n_samples': len(merged),
                            'pyaging_mean': merged['predicted_age'].mean(),
                            'biolearn_mean': merged['predicted_value'].mean()
                        }

            consensus['correlation_analysis'] = correlations

        return consensus

    def _generate_integrated_summary(self, results: Dict) -> Dict:
        """Generate integrated summary across both libraries."""

        summary = {
            'libraries_used': self.get_available_libraries(),
            'total_analyses_completed': 0,
            'pyaging_summary': {},
            'biolearn_summary': {},
            'consensus_summary': {}
        }

        # PyAging summary
        if results.get('pyaging_results'):
            summary['total_analyses_completed'] += 1
            py_results = results['pyaging_results']
            summary['pyaging_summary'] = {
                'clocks_used': len(py_results.get('senescence_clocks_used', [])),
                'samples_analyzed': py_results.get('senescence_summary', {}).get('total_samples', 0),
                'rehmgb1_metrics_available': bool(py_results.get('rehmgb1_relevant_metrics'))
            }

        # Biolearn summary
        if results.get('biolearn_results'):
            summary['total_analyses_completed'] += 1
            bio_results = results['biolearn_results']
            summary['biolearn_summary'] = {
                'biomarkers_calculated': bio_results.get('rehmgb1_summary', {}).get('biomarkers_calculated', 0),
                'samples_analyzed': bio_results.get('rehmgb1_summary', {}).get('total_samples', 0),
                'mortality_risk_available': bool(bio_results.get('mortality_risk_scores')),
                'aging_pace_available': bool(bio_results.get('aging_pace_metrics'))
            }

        # Consensus summary
        if results.get('consensus_analysis'):
            consensus = results['consensus_analysis']
            summary['consensus_summary'] = {
                'consensus_clocks_found': len(consensus.get('consensus_clocks_used', [])),
                'correlations_calculated': len(consensus.get('correlation_analysis', {}))
            }

        return summary

    def _calculate_integrated_rehmgb1_metrics(self, results: Dict) -> Dict:
        """Calculate integrated ReHMGB1-relevant metrics."""

        integrated_metrics = {
            'senescence_acceleration_score': None,
            'inflammatory_aging_index': None,
            'mortality_risk_composite': None,
            'rehmgb1_pathway_relevance': {}
        }

        # Combine senescence acceleration from both libraries
        py_metrics = results.get('pyaging_results', {}).get('rehmgb1_relevant_metrics', {})
        bio_metrics = results.get('biolearn_results', {}).get('rehmgb1_summary', {})

        # Senescence acceleration score
        accelerations = []

        if py_metrics.get('dna_methylation_senescence_score'):
            py_score = py_metrics['dna_methylation_senescence_score']['mean']
            accelerations.append(py_score)

        if bio_metrics.get('age_acceleration_analysis'):
            bio_score = bio_metrics['age_acceleration_analysis']['mean_acceleration']
            accelerations.append(bio_score)

        if accelerations:
            integrated_metrics['senescence_acceleration_score'] = {
                'combined_score': np.mean(accelerations),
                'score_std': np.std(accelerations) if len(accelerations) > 1 else 0,
                'n_libraries': len(accelerations)
            }

        # Inflammatory aging index from pyaging
        if py_metrics.get('inflammatory_aging_clocks'):
            integrated_metrics['inflammatory_aging_index'] = {
                'clocks_used': py_metrics['inflammatory_aging_clocks'],
                'n_inflammatory_clocks': len(py_metrics['inflammatory_aging_clocks'])
            }

        # Mortality risk composite from biolearn
        bio_results = results.get('biolearn_results', {})
        if bio_results.get('mortality_risk_scores'):
            risk_summary = bio_results['mortality_risk_scores']['summary']
            integrated_metrics['mortality_risk_composite'] = {
                'high_risk_percentage': (
                    risk_summary['high_risk_samples'] /
                    risk_summary['total_samples'] * 100
                ),
                'mean_risk_score': risk_summary['mean_mortality_risk']
            }

        # ReHMGB1 pathway relevance scoring
        integrated_metrics['rehmgb1_pathway_relevance'] = {
            'rage_signaling_clocks': 0,
            'inflammatory_markers': 0,
            'senescence_indicators': 0
        }

        # Count relevant biomarkers across libraries
        if py_metrics.get('inflammatory_aging_clocks'):
            integrated_metrics['rehmgb1_pathway_relevance']['inflammatory_markers'] += len(
                py_metrics['inflammatory_aging_clocks']
            )

        if bio_results.get('senescence_biomarkers') is not None:
            biomarker_data = bio_results['senescence_biomarkers']
            senescence_clocks = ['PhenoAge', 'GrimAge', 'DunedinPACE']
            senescence_count = biomarker_data[
                biomarker_data['clock'].isin(senescence_clocks)
            ]['clock'].nunique()
            integrated_metrics['rehmgb1_pathway_relevance']['senescence_indicators'] = senescence_count

        return integrated_metrics

    def export_integrated_results(
        self,
        results: Dict,
        output_dir: str,
        prefix: str = 'integrated_aging_analysis'
    ) -> Dict[str, str]:
        """
        Export comprehensive integrated analysis results.

        Args:
            results: Results from comprehensive_senescence_analysis()
            output_dir: Output directory
            prefix: File prefix

        Returns:
            Dictionary mapping result types to file paths
        """
        os.makedirs(output_dir, exist_ok=True)

        exported_files = {}

        # Export pyaging results
        if results.get('pyaging_results'):
            py_files = self.pyaging_analyzer.export_results(
                results['pyaging_results'],
                output_dir,
                f'{prefix}_pyaging'
            )
            exported_files.update({f'pyaging_{k}': v for k, v in py_files.items()})

        # Export biolearn results
        if results.get('biolearn_results'):
            bio_files = self.biolearn_analyzer.export_results(
                results['biolearn_results'],
                output_dir,
                f'{prefix}_biolearn'
            )
            exported_files.update({f'biolearn_{k}': v for k, v in bio_files.items()})

        # Export integrated summary
        import json
        summary_file = os.path.join(output_dir, f'{prefix}_integrated_summary.json')

        json_summary = {
            'integrated_summary': results.get('integrated_summary', {}),
            'rehmgb1_integrated_metrics': results.get('rehmgb1_integrated_metrics', {}),
            'consensus_analysis': results.get('consensus_analysis', {})
        }

        with open(summary_file, 'w') as f:
            json.dump(json_summary, f, indent=2, default=str)
        exported_files['integrated_summary'] = summary_file

        if self.verbose:
            print(f"Integrated results exported to {output_dir}")
            for result_type, file_path in exported_files.items():
                print(f"  {result_type}: {file_path}")

        return exported_files


def create_demo_integrated_analysis() -> str:
    """Create demonstration script for integrated aging biomarker analysis."""

    demo_script = '''#!/usr/bin/env python3
"""
Integrated Aging Biomarker Analysis Demo

Demonstrates comprehensive aging biomarker analysis combining pyaging and
biolearn libraries for ReHMGB1 senescence research.
"""

import pandas as pd
import numpy as np
from crispr_toolkit.analysis.aging.biomarkers import AgingBiomarkerAnalyzer

def main():
    print("🧬 Integrated Aging Biomarker Analysis Demo")
    print("=" * 60)

    # Initialize integrated analyzer
    print("\\n1. Initializing integrated aging biomarker analyzer...")
    try:
        analyzer = AgingBiomarkerAnalyzer(verbose=True)
        libraries = analyzer.get_available_libraries()
        print(f"   ✅ Initialized with libraries: {libraries}")
    except Exception as e:
        print(f"   ❌ Initialization failed: {e}")
        return

    # Show available clocks
    print("\\n2. Available aging clocks across libraries:")
    all_clocks = analyzer.list_all_clocks()
    for library, clock_info in all_clocks.items():
        print(f"   📚 {library.upper()}:")
        if 'senescence_clocks' in clock_info:
            for category, clocks in clock_info['senescence_clocks'].items():
                print(f"      {category}: {len(clocks)} clocks")

    # Generate demo data
    print("\\n3. Generating synthetic omics data...")
    n_samples = 40
    n_features = 800

    np.random.seed(42)
    age_range = np.random.uniform(30, 80, n_samples)

    # Generate realistic omics data
    data = pd.DataFrame(
        np.random.beta(2, 2, (n_samples, n_features)),
        index=[f"Sample_{i:03d}" for i in range(n_samples)],
        columns=[f"Feature_{i:06d}" for i in range(n_features)]
    )

    chronological_age = pd.Series(age_range, index=data.index)

    metadata = pd.DataFrame({
        'age': age_range,
        'sex': np.random.choice(['M', 'F'], n_samples),
        'bmi': np.random.normal(25, 5, n_samples)
    }, index=data.index)

    print(f"   📊 Generated: {data.shape[0]} samples x {data.shape[1]} features")
    print(f"   📅 Age range: {age_range.min():.1f} - {age_range.max():.1f} years")

    # Perform comprehensive analysis
    print("\\n4. Performing comprehensive senescence analysis...")
    try:
        results = analyzer.comprehensive_senescence_analysis(
            data=data,
            chronological_age=chronological_age,
            metadata=metadata,
            use_consensus_clocks=True
        )

        print("   ✅ Comprehensive analysis completed!")

        # Display integrated summary
        print("\\n5. Integrated Analysis Summary:")
        summary = results['integrated_summary']
        print(f"   📚 Libraries used: {summary['libraries_used']}")
        print(f"   ✅ Analyses completed: {summary['total_analyses_completed']}")

        if summary.get('pyaging_summary'):
            py_sum = summary['pyaging_summary']
            print(f"   🔬 PyAging: {py_sum['clocks_used']} clocks, {py_sum['samples_analyzed']} samples")

        if summary.get('biolearn_summary'):
            bio_sum = summary['biolearn_summary']
            print(f"   🔬 Biolearn: {bio_sum['biomarkers_calculated']} biomarkers, {bio_sum['samples_analyzed']} samples")

        # Show ReHMGB1 integrated metrics
        print("\\n6. ReHMGB1 Integrated Metrics:")
        rehmgb1_metrics = results['rehmgb1_integrated_metrics']

        if rehmgb1_metrics.get('senescence_acceleration_score'):
            acc_score = rehmgb1_metrics['senescence_acceleration_score']
            print(f"   📈 Senescence acceleration: {acc_score['combined_score']:.2f} years")
            print(f"      (combined from {acc_score['n_libraries']} libraries)")

        if rehmgb1_metrics.get('inflammatory_aging_index'):
            inflam_idx = rehmgb1_metrics['inflammatory_aging_index']
            print(f"   🔥 Inflammatory aging clocks: {inflam_idx['n_inflammatory_clocks']}")

        if rehmgb1_metrics.get('mortality_risk_composite'):
            mort_risk = rehmgb1_metrics['mortality_risk_composite']
            print(f"   ⚠️  High mortality risk samples: {mort_risk['high_risk_percentage']:.1f}%")

        # Show consensus analysis if available
        if results.get('consensus_analysis'):
            print("\\n7. Consensus Analysis:")
            consensus = results['consensus_analysis']
            correlations = consensus.get('correlation_analysis', {})

            if correlations:
                print("   🤝 Cross-library clock correlations:")
                for clock, corr_data in correlations.items():
                    print(f"      {clock}: r = {corr_data['correlation']:.3f} (n={corr_data['n_samples']})")
            else:
                print("   ℹ️  No consensus clocks found (expected with synthetic data)")

        # Export results
        print("\\n8. Exporting integrated results...")
        exported = analyzer.export_integrated_results(
            results, './integrated_aging_demo_results'
        )
        print("   ✅ All results exported successfully!")

    except Exception as e:
        print(f"   ❌ Analysis failed: {e}")
        print("   💡 This is expected with synthetic data - real omics data needed")

    print("\\n🎉 Integrated aging biomarker analysis demo completed!")
    print("\\n💡 Next steps for ReHMGB1 research:")
    print("   - Load real DNA methylation and/or gene expression data")
    print("   - Compare aging acceleration across treatment conditions")
    print("   - Correlate biomarker scores with ReHMGB1 expression levels")
    print("   - Integrate with CRISPR screen results from Phase 1")
    print("   - Use consensus clocks for robust aging measurements")

if __name__ == "__main__":
    main()
'''

    return demo_script
    return demo_script

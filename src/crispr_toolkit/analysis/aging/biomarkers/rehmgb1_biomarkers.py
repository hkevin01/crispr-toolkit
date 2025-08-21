"""
ReHMGB1-Specific Biomarker Analysis Module

Specialized aging biomarker analysis focused on ReHMGB1/RAGE signaling
pathways and their role in senescence. Integrates with both pyaging and
biolearn to provide pathway-specific aging metrics.

ReHMGB1 Research Focus:
- RAGE-mediated JAK/STAT and NF-κB pathway aging
- Senescence-associated secretory phenotype (SASP) biomarkers
- Inflammatory aging signatures
- DNA methylation patterns in senescence genes
- Clinical aging metrics relevant to ReHMGB1 function

Author: CRISPR Toolkit Development Team
"""

import warnings
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

try:
    from ..screens.senescence_analysis import SenescenceScreenAnalyzer
    SENESCENCE_ANALYZER_AVAILABLE = True
except ImportError:
    SENESCENCE_ANALYZER_AVAILABLE = False
    warnings.warn("SenescenceScreenAnalyzer not available")


class ReHMGB1BiomarkerAnalyzer:
    """
    Specialized analyzer for ReHMGB1-focused aging biomarker analysis.

    Combines aging clock predictions with ReHMGB1 pathway knowledge to
    provide targeted analysis for senescence research.
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize ReHMGB1 biomarker analyzer.

        Args:
            verbose: Enable verbose logging
        """
        self.verbose = verbose

        # Load ReHMGB1 pathway gene sets
        self.rehmgb1_pathways = self._load_rehmgb1_pathways()
        self.aging_biomarker_genes = self._load_aging_biomarker_genes()
        self.senescence_signatures = self._load_senescence_signatures()

        # Initialize senescence analyzer if available
        self.senescence_analyzer = None
        if SENESCENCE_ANALYZER_AVAILABLE:
            try:
                self.senescence_analyzer = SenescenceScreenAnalyzer()
                if verbose:
                    print("✅ Senescence analyzer integration enabled")
            except Exception as e:
                if verbose:
                    print(f"⚠️  Senescence analyzer failed to initialize: {e}")

        if verbose:
            print(f"🧬 ReHMGB1 analyzer initialized with "
                  f"{len(self.rehmgb1_pathways)} pathway modules")

    def _load_rehmgb1_pathways(self) -> Dict[str, List[str]]:
        """Load ReHMGB1/RAGE signaling pathway gene sets."""
        return {
            'hmgb1_core': [
                'HMGB1',  # High mobility group box 1
                'HMGB2',  # High mobility group box 2
                'HMGB3'   # High mobility group box 3
            ],
            'rage_signaling': [
                'AGER',    # RAGE receptor
                'MYD88',   # Myeloid differentiation primary response 88
                'IRAK1',   # Interleukin-1 receptor-associated kinase 1
                'IRAK4',   # Interleukin-1 receptor-associated kinase 4
                'TRAF6',   # TNF receptor associated factor 6
                'MAP3K7',  # TAK1 - Mitogen-activated protein kinase kinase kinase 7
                'TAB1',    # TGF-beta activated kinase 1 binding protein 1
                'TAB2'     # TGF-beta activated kinase 1 binding protein 2
            ],
            'mapk_cascade': [
                'MAPK14',  # p38 MAPK
                'MAPK8',   # JNK1
                'MAPK9',   # JNK2
                'MAPK10',  # JNK3
                'MAPK1',   # ERK2
                'MAPK3'    # ERK1
            ],
            'jak_stat_pathway': [
                'JAK1', 'JAK2', 'JAK3', 'TYK2',
                'STAT1', 'STAT2', 'STAT3', 'STAT4', 'STAT5A', 'STAT5B', 'STAT6'
            ],
            'nfkb_signaling': [
                'NFKB1', 'NFKB2', 'RELA', 'RELB', 'REL',
                'IKBKA', 'IKBKB', 'IKBKG', 'NFKBIA', 'NFKBIB', 'NFKBIE'
            ],
            'oxidative_stress_response': [
                'SOD1', 'SOD2', 'SOD3',
                'CAT',  # Catalase
                'GPX1', 'GPX4',
                'PRDX1', 'PRDX3', 'PRDX4',
                'NFE2L2',  # NRF2
                'KEAP1'
            ],
            'sasp_factors': [
                'IL6', 'IL1A', 'IL1B', 'TNF', 'IFNA1', 'IFNB1',
                'CXCL8', 'CCL2', 'CCL20', 'CXCL1', 'CXCL2',
                'MMP1', 'MMP3', 'MMP9', 'MMP12',
                'TIMP1', 'TIMP2'
            ]
        }

    def _load_aging_biomarker_genes(self) -> Dict[str, List[str]]:
        """Load genes commonly used in aging biomarker development."""
        return {
            'dna_methylation_clocks': [
                # Genes frequently used in DNA methylation aging clocks
                'ELOVL2', 'FHL2', 'PENK', 'KLF14', 'TRIM59',
                'CCDC102B', 'HOXC4', 'NPTX2', 'MIR29B2CHG'
            ],
            'cellular_senescence_markers': [
                'CDKN1A',  # p21
                'CDKN2A',  # p16
                'CDKN1B',  # p27
                'RB1', 'TP53', 'TP21',
                'LMNB1',   # Lamin B1 (senescence marker)
                'GLB1'     # SA-β-gal marker
            ],
            'telomere_maintenance': [
                'TERT', 'TERC', 'DKC1', 'TINF2', 'TPP1',
                'POT1', 'RTEL1', 'CTC1', 'STN1', 'TEN1'
            ],
            'autophagy_aging': [
                'ATG5', 'ATG7', 'ATG12', 'BECN1', 'LC3A', 'LC3B',
                'SQSTM1', 'ULK1', 'MTOR', 'TFEB'
            ],
            'mitochondrial_aging': [
                'PGC1A', 'TFAM', 'NRF1', 'NRF2', 'SIRT1', 'SIRT3',
                'PINK1', 'PARKIN', 'MFN1', 'MFN2', 'OPA1'
            ]
        }

    def _load_senescence_signatures(self) -> Dict[str, Dict]:
        """Load validated senescence gene expression signatures."""
        return {
            'fridman_senescence_up': {
                'description': 'Genes upregulated in senescence (Fridman et al.)',
                'genes': [
                    'CDKN1A', 'IL6', 'IL1B', 'CXCL8', 'CCL2', 'MMP1',
                    'SERPINE1', 'IGFBP3', 'GADD45A', 'DDIT3'
                ],
                'direction': 'up'
            },
            'fridman_senescence_down': {
                'description': 'Genes downregulated in senescence (Fridman et al.)',
                'genes': [
                    'PCNA', 'MCM2', 'MCM3', 'MCM7', 'CDC20', 'CCNA2',
                    'CCNB1', 'CCNE1', 'CDK1', 'CDK2'
                ],
                'direction': 'down'
            },
            'rehmgb1_induced_senescence': {
                'description': 'Genes associated with ReHMGB1-induced senescence',
                'genes': [
                    'HMGB1', 'AGER', 'IL6', 'TNF', 'NFKB1', 'STAT3',
                    'CDKN1A', 'CDKN2A', 'MMP1', 'CXCL8'
                ],
                'direction': 'up'
            }
        }

    def analyze_rehmgb1_aging_profile(
        self,
        aging_predictions: pd.DataFrame,
        expression_data: Optional[pd.DataFrame] = None,
        methylation_data: Optional[pd.DataFrame] = None,
        chronological_age: Optional[pd.Series] = None
    ) -> Dict[str, Any]:
        """
        Comprehensive ReHMGB1-focused aging profile analysis.

        Args:
            aging_predictions: Aging clock predictions from pyaging/biolearn
            expression_data: Gene expression data (optional)
            methylation_data: DNA methylation data (optional)
            chronological_age: Chronological ages for comparison

        Returns:
            Comprehensive ReHMGB1 aging analysis results
        """
        results = {
            'rehmgb1_pathway_aging': {},
            'senescence_signature_scores': {},
            'aging_acceleration_by_pathway': {},
            'clinical_rehmgb1_metrics': {},
            'integrated_senescence_score': None
        }

        # Analyze aging patterns by ReHMGB1 pathway modules
        results['rehmgb1_pathway_aging'] = self._analyze_pathway_aging_patterns(
            aging_predictions, chronological_age
        )

        # Calculate senescence signature scores if expression data available
        if expression_data is not None:
            results['senescence_signature_scores'] = self._calculate_senescence_signatures(
                expression_data
            )

        # Analyze aging acceleration by pathway involvement
        if chronological_age is not None:
            results['aging_acceleration_by_pathway'] = self._analyze_pathway_age_acceleration(
                aging_predictions, chronological_age, expression_data
            )

        # Generate clinical ReHMGB1 metrics
        results['clinical_rehmgb1_metrics'] = self._generate_clinical_rehmgb1_metrics(
            aging_predictions, expression_data, chronological_age
        )

        # Calculate integrated senescence score
        results['integrated_senescence_score'] = self._calculate_integrated_senescence_score(
            results
        )

        return results

    def _analyze_pathway_aging_patterns(
        self,
        aging_predictions: pd.DataFrame,
        chronological_age: Optional[pd.Series]
    ) -> Dict[str, Any]:
        """Analyze aging patterns specific to ReHMGB1 pathways."""

        pathway_analysis = {}

        # Categorize aging clocks by relevance to ReHMGB1 pathways
        clock_pathway_relevance = {
            'inflammatory_aging': ['GrimAge', 'PhenoAge', 'inflammage'],
            'dna_methylation': ['Horvath2013', 'Hannum2013', 'zhang2019'],
            'mortality_prediction': ['GrimAge', 'PhenoAge'],
            'senescence_pace': ['DunedinPACE', 'dunedinpace']
        }

        for pathway_category, relevant_clocks in clock_pathway_relevance.items():
            # Find predictions from relevant clocks
            pathway_preds = aging_predictions[
                aging_predictions['clock'].isin(relevant_clocks)
            ]

            if not pathway_preds.empty:
                # Calculate pathway-specific aging metrics
                sample_scores = pathway_preds.groupby('sample_id').agg({
                    'predicted_age': 'mean' if 'predicted_age' in pathway_preds.columns
                                    else 'first',
                    'predicted_value': 'mean' if 'predicted_value' in pathway_preds.columns
                                      else 'first'
                }).fillna(0)

                # Use the available prediction column
                pred_col = 'predicted_age' if 'predicted_age' in sample_scores.columns else 'predicted_value'

                pathway_analysis[pathway_category] = {
                    'n_clocks_used': pathway_preds['clock'].nunique(),
                    'clocks_used': pathway_preds['clock'].unique().tolist(),
                    'mean_predicted_age': sample_scores[pred_col].mean(),
                    'std_predicted_age': sample_scores[pred_col].std(),
                    'n_samples': len(sample_scores)
                }

                # Calculate age acceleration if chronological age available
                if chronological_age is not None:
                    sample_scores['chronological_age'] = sample_scores.index.map(
                        chronological_age.to_dict()
                    )
                    sample_scores['age_acceleration'] = (
                        sample_scores[pred_col] - sample_scores['chronological_age']
                    )

                    pathway_analysis[pathway_category].update({
                        'mean_age_acceleration': sample_scores['age_acceleration'].mean(),
                        'std_age_acceleration': sample_scores['age_acceleration'].std(),
                        'samples_with_acceleration': (sample_scores['age_acceleration'] > 0).sum(),
                        'percent_accelerated': (sample_scores['age_acceleration'] > 0).mean() * 100
                    })

        return pathway_analysis

    def _calculate_senescence_signatures(
        self,
        expression_data: pd.DataFrame
    ) -> Dict[str, Any]:
        """Calculate senescence signature scores from expression data."""

        signature_scores = {}

        for sig_name, sig_info in self.senescence_signatures.items():
            sig_genes = sig_info['genes']
            direction = sig_info['direction']

            # Find genes present in the expression data
            available_genes = [gene for gene in sig_genes if gene in expression_data.columns]

            if available_genes:
                # Calculate signature score
                sig_expression = expression_data[available_genes]

                if direction == 'up':
                    # For upregulated signatures, higher expression = higher score
                    signature_score = sig_expression.mean(axis=1)
                else:
                    # For downregulated signatures, lower expression = higher score
                    signature_score = -sig_expression.mean(axis=1)

                signature_scores[sig_name] = {
                    'scores': signature_score,
                    'genes_used': available_genes,
                    'genes_available': len(available_genes),
                    'genes_total': len(sig_genes),
                    'coverage': len(available_genes) / len(sig_genes),
                    'mean_score': signature_score.mean(),
                    'std_score': signature_score.std(),
                    'description': sig_info['description']
                }
            else:
                signature_scores[sig_name] = {
                    'scores': None,
                    'genes_used': [],
                    'genes_available': 0,
                    'genes_total': len(sig_genes),
                    'coverage': 0,
                    'message': 'No signature genes found in expression data'
                }

        return signature_scores

    def _analyze_pathway_age_acceleration(
        self,
        aging_predictions: pd.DataFrame,
        chronological_age: pd.Series,
        expression_data: Optional[pd.DataFrame]
    ) -> Dict[str, Any]:
        """Analyze age acceleration patterns by pathway involvement."""

        acceleration_analysis = {}

        # Calculate overall age acceleration
        if 'predicted_age' in aging_predictions.columns:
            pred_col = 'predicted_age'
        else:
            pred_col = 'predicted_value'

        # Merge with chronological age
        age_data = aging_predictions.copy()
        age_data['chronological_age'] = age_data['sample_id'].map(
            chronological_age.to_dict()
        )
        age_data['age_acceleration'] = (
            age_data[pred_col] - age_data['chronological_age']
        )

        # Group by sample and calculate mean acceleration
        sample_acceleration = age_data.groupby('sample_id')['age_acceleration'].mean()

        acceleration_analysis['overall'] = {
            'mean_acceleration': sample_acceleration.mean(),
            'std_acceleration': sample_acceleration.std(),
            'high_accelerators': (sample_acceleration > 5).sum(),  # >5 years acceleration
            'normal_aging': ((sample_acceleration >= -2) & (sample_acceleration <= 5)).sum(),
            'slow_agers': (sample_acceleration < -2).sum()
        }

        # Analyze acceleration by pathway expression if available
        if expression_data is not None:
            for pathway_name, pathway_genes in self.rehmgb1_pathways.items():
                available_genes = [g for g in pathway_genes if g in expression_data.columns]

                if available_genes:
                    # Calculate pathway expression score
                    pathway_expression = expression_data[available_genes].mean(axis=1)

                    # Correlate with age acceleration
                    merged_data = pd.DataFrame({
                        'sample_id': expression_data.index,
                        'pathway_expression': pathway_expression
                    }).merge(
                        sample_acceleration.reset_index(),
                        on='sample_id'
                    )

                    if len(merged_data) > 3:  # Need at least 3 samples for correlation
                        correlation = merged_data['pathway_expression'].corr(
                            merged_data['age_acceleration']
                        )

                        acceleration_analysis[pathway_name] = {
                            'genes_available': len(available_genes),
                            'correlation_with_acceleration': correlation,
                            'mean_pathway_expression': pathway_expression.mean(),
                            'std_pathway_expression': pathway_expression.std(),
                            'interpretation': self._interpret_pathway_acceleration_correlation(
                                correlation, pathway_name
                            )
                        }

        return acceleration_analysis

    def _interpret_pathway_acceleration_correlation(
        self,
        correlation: float,
        pathway_name: str
    ) -> str:
        """Interpret pathway-age acceleration correlation in ReHMGB1 context."""

        if pd.isna(correlation):
            return "Insufficient data for correlation analysis"

        abs_corr = abs(correlation)

        if abs_corr < 0.2:
            strength = "weak"
        elif abs_corr < 0.5:
            strength = "moderate"
        else:
            strength = "strong"

        direction = "positive" if correlation > 0 else "negative"

        # Pathway-specific interpretations
        interpretations = {
            'hmgb1_core': f"{strength} {direction} correlation suggests ReHMGB1 levels {'promote' if correlation > 0 else 'protect against'} aging acceleration",
            'rage_signaling': f"{strength} {direction} correlation indicates RAGE signaling {'contributes to' if correlation > 0 else 'opposes'} accelerated aging",
            'sasp_factors': f"{strength} {direction} correlation shows SASP {'enhances' if correlation > 0 else 'reduces'} age acceleration",
            'jak_stat_pathway': f"{strength} {direction} correlation suggests JAK/STAT activity {'promotes' if correlation > 0 else 'inhibits'} aging",
            'nfkb_signaling': f"{strength} {direction} correlation indicates NF-κB {'drives' if correlation > 0 else 'protects against'} aging acceleration"
        }

        return interpretations.get(
            pathway_name,
            f"{strength} {direction} correlation with age acceleration"
        )

    def _generate_clinical_rehmgb1_metrics(
        self,
        aging_predictions: pd.DataFrame,
        expression_data: Optional[pd.DataFrame],
        chronological_age: Optional[pd.Series]
    ) -> Dict[str, Any]:
        """Generate clinically relevant ReHMGB1 aging metrics."""

        clinical_metrics = {
            'senescence_risk_score': None,
            'inflammatory_aging_index': None,
            'rehmgb1_pathway_dysregulation': None,
            'clinical_recommendations': []
        }

        # Calculate senescence risk score
        senescence_clocks = ['PhenoAge', 'GrimAge', 'DunedinPACE', 'phenoage', 'grimage', 'dunedinpace']
        senescence_preds = aging_predictions[
            aging_predictions['clock'].isin(senescence_clocks)
        ]

        if not senescence_preds.empty and chronological_age is not None:
            pred_col = 'predicted_age' if 'predicted_age' in senescence_preds.columns else 'predicted_value'

            # Calculate age acceleration for senescence clocks
            senescence_preds['chronological_age'] = senescence_preds['sample_id'].map(
                chronological_age.to_dict()
            )
            senescence_preds['age_acceleration'] = (
                senescence_preds[pred_col] - senescence_preds['chronological_age']
            )

            sample_risk = senescence_preds.groupby('sample_id')['age_acceleration'].mean()

            # Risk categorization
            def categorize_risk(acceleration):
                if acceleration > 10:
                    return 'very_high'
                elif acceleration > 5:
                    return 'high'
                elif acceleration > 0:
                    return 'moderate'
                elif acceleration > -5:
                    return 'low'
                else:
                    return 'very_low'

            risk_categories = sample_risk.apply(categorize_risk)

            clinical_metrics['senescence_risk_score'] = {
                'sample_scores': sample_risk,
                'risk_categories': risk_categories,
                'high_risk_samples': (risk_categories.isin(['high', 'very_high'])).sum(),
                'mean_acceleration': sample_risk.mean(),
                'interpretation': f"Average {sample_risk.mean():.1f} years of senescence acceleration"
            }

        # Calculate inflammatory aging index
        inflammatory_clocks = ['GrimAge', 'grimage', 'inflammage']
        inflammatory_preds = aging_predictions[
            aging_predictions['clock'].isin(inflammatory_clocks)
        ]

        if not inflammatory_preds.empty:
            pred_col = 'predicted_age' if 'predicted_age' in inflammatory_preds.columns else 'predicted_value'
            inflammatory_scores = inflammatory_preds.groupby('sample_id')[pred_col].mean()

            clinical_metrics['inflammatory_aging_index'] = {
                'sample_scores': inflammatory_scores,
                'mean_score': inflammatory_scores.mean(),
                'std_score': inflammatory_scores.std(),
                'interpretation': "Higher scores indicate greater inflammatory aging burden"
            }

        # Assess ReHMGB1 pathway dysregulation if expression data available
        if expression_data is not None:
            dysregulation_scores = {}

            for pathway_name, pathway_genes in self.rehmgb1_pathways.items():
                available_genes = [g for g in pathway_genes if g in expression_data.columns]

                if available_genes:
                    pathway_expr = expression_data[available_genes]

                    # Calculate pathway dysregulation as variance in expression
                    pathway_cv = pathway_expr.std(axis=1) / (pathway_expr.mean(axis=1) + 1e-6)
                    dysregulation_scores[pathway_name] = pathway_cv.mean()

            if dysregulation_scores:
                clinical_metrics['rehmgb1_pathway_dysregulation'] = {
                    'pathway_scores': dysregulation_scores,
                    'overall_dysregulation': np.mean(list(dysregulation_scores.values())),
                    'most_dysregulated': max(dysregulation_scores.items(), key=lambda x: x[1]),
                    'interpretation': "Higher scores indicate greater pathway dysregulation"
                }

        # Generate clinical recommendations
        recommendations = []

        if clinical_metrics.get('senescence_risk_score'):
            risk_data = clinical_metrics['senescence_risk_score']
            if risk_data['high_risk_samples'] > 0:
                recommendations.append(
                    f"{risk_data['high_risk_samples']} samples show high senescence risk - "
                    "consider anti-aging interventions"
                )

        if clinical_metrics.get('inflammatory_aging_index'):
            recommendations.append(
                "Monitor inflammatory markers and consider anti-inflammatory strategies"
            )

        if clinical_metrics.get('rehmgb1_pathway_dysregulation'):
            dysreg = clinical_metrics['rehmgb1_pathway_dysregulation']
            most_dysreg = dysreg['most_dysregulated']
            recommendations.append(
                f"Focus on {most_dysreg[0]} pathway - shows highest dysregulation"
            )

        clinical_metrics['clinical_recommendations'] = recommendations

        return clinical_metrics

    def _calculate_integrated_senescence_score(
        self,
        results: Dict[str, Any]
    ) -> Optional[Dict[str, Any]]:
        """Calculate an integrated senescence score combining multiple metrics."""

        scores = []
        components = []

        # Add pathway aging patterns
        pathway_aging = results.get('rehmgb1_pathway_aging', {})
        for pathway, data in pathway_aging.items():
            if data.get('mean_age_acceleration') is not None:
                scores.append(data['mean_age_acceleration'])
                components.append(f"{pathway}_acceleration")

        # Add senescence signature scores
        sig_scores = results.get('senescence_signature_scores', {})
        for sig_name, sig_data in sig_scores.items():
            if sig_data.get('mean_score') is not None:
                scores.append(sig_data['mean_score'])
                components.append(f"{sig_name}_signature")

        # Add clinical risk scores
        clinical = results.get('clinical_rehmgb1_metrics', {})
        if clinical.get('senescence_risk_score', {}).get('mean_acceleration') is not None:
            scores.append(clinical['senescence_risk_score']['mean_acceleration'])
            components.append('clinical_senescence_risk')

        if not scores:
            return None

        # Calculate integrated score
        integrated_score = np.mean(scores)
        score_std = np.std(scores) if len(scores) > 1 else 0

        # Interpret score
        def interpret_integrated_score(score):
            if score > 8:
                return "Very high senescence burden - immediate intervention recommended"
            elif score > 4:
                return "High senescence burden - intervention recommended"
            elif score > 0:
                return "Moderate senescence burden - monitoring recommended"
            elif score > -2:
                return "Normal aging profile"
            else:
                return "Slow aging profile - protective factors present"

        return {
            'integrated_score': integrated_score,
            'score_std': score_std,
            'n_components': len(scores),
            'components_used': components,
            'interpretation': interpret_integrated_score(integrated_score),
            'percentile_rank': None  # Would need population data to calculate
        }

    def export_rehmgb1_analysis(
        self,
        results: Dict[str, Any],
        output_dir: str,
        prefix: str = 'rehmgb1_biomarker_analysis'
    ) -> Dict[str, str]:
        """
        Export ReHMGB1 biomarker analysis results.

        Args:
            results: Results from analyze_rehmgb1_aging_profile()
            output_dir: Output directory
            prefix: File prefix

        Returns:
            Dictionary mapping result types to file paths
        """
        import json
        import os

        os.makedirs(output_dir, exist_ok=True)
        exported_files = {}

        # Export pathway aging analysis
        if results.get('rehmgb1_pathway_aging'):
            pathway_file = os.path.join(output_dir, f'{prefix}_pathway_aging.json')
            with open(pathway_file, 'w') as f:
                json.dump(results['rehmgb1_pathway_aging'], f, indent=2, default=str)
            exported_files['pathway_aging'] = pathway_file

        # Export senescence signatures
        if results.get('senescence_signature_scores'):
            sig_file = os.path.join(output_dir, f'{prefix}_senescence_signatures.json')

            # Prepare data for JSON export (handle pandas Series)
            sig_data_export = {}
            for sig_name, sig_info in results['senescence_signature_scores'].items():
                sig_export = {k: v for k, v in sig_info.items() if k != 'scores'}
                if sig_info.get('scores') is not None:
                    sig_export['sample_scores'] = sig_info['scores'].to_dict()
                sig_data_export[sig_name] = sig_export

            with open(sig_file, 'w') as f:
                json.dump(sig_data_export, f, indent=2, default=str)
            exported_files['senescence_signatures'] = sig_file

        # Export clinical metrics
        if results.get('clinical_rehmgb1_metrics'):
            clinical_file = os.path.join(output_dir, f'{prefix}_clinical_metrics.json')

            # Prepare clinical data for export
            clinical_export = {}
            for metric_name, metric_data in results['clinical_rehmgb1_metrics'].items():
                if isinstance(metric_data, dict):
                    metric_export = {}
                    for k, v in metric_data.items():
                        if hasattr(v, 'to_dict'):  # pandas Series
                            metric_export[k] = v.to_dict()
                        else:
                            metric_export[k] = v
                    clinical_export[metric_name] = metric_export
                else:
                    clinical_export[metric_name] = metric_data

            with open(clinical_file, 'w') as f:
                json.dump(clinical_export, f, indent=2, default=str)
            exported_files['clinical_metrics'] = clinical_file

        # Export integrated score
        if results.get('integrated_senescence_score'):
            score_file = os.path.join(output_dir, f'{prefix}_integrated_score.json')
            with open(score_file, 'w') as f:
                json.dump(results['integrated_senescence_score'], f, indent=2, default=str)
            exported_files['integrated_score'] = score_file

        if self.verbose:
            print(f"ReHMGB1 analysis results exported to {output_dir}")
            for result_type, file_path in exported_files.items():
                print(f"  {result_type}: {file_path}")

        return exported_files
        return exported_files

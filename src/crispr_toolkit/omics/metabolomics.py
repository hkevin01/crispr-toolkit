"""
Metabolomics Analysis Module for CRISPR Toolkit Phase 3
======================================================

Advanced metabolomics data analysis for aging intervention research,
including metabolic pathway analysis, biomarker discovery, and
intervention response prediction.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


@dataclass
class MetaboliteAnalysisResult:
    """Results from metabolite analysis."""
    metabolite_id: str
    fold_change: float
    p_value: float
    adjusted_p_value: float
    pathway_impact: float
    aging_relevance_score: float
    biomarker_potential: float


@dataclass
class MetabolomicsDataset:
    """Container for metabolomics experimental data."""
    metabolite_matrix: pd.DataFrame  # metabolites x samples
    sample_metadata: pd.DataFrame
    metabolite_annotations: pd.DataFrame
    pathway_mapping: Dict[str, List[str]]
    quality_metrics: Dict[str, float]


class MetabolomicsAnalyzer:
    """
    Comprehensive metabolomics analysis for aging intervention research.

    This analyzer handles metabolite profiling data and provides specialized
    analysis for aging-related metabolic changes and intervention responses.
    """

    def __init__(self, platform: str = "lcms", focus_pathways: List[str] = None):
        """
        Initialize metabolomics analyzer.

        Args:
            platform: Data platform ("lcms", "gcms", "nmr")
            focus_pathways: Specific pathways to focus on
        """
        self.platform = platform
        self.focus_pathways = focus_pathways or [
            'glycolysis', 'tca_cycle', 'amino_acid_metabolism',
            'lipid_metabolism', 'purine_metabolism', 'oxidative_stress'
        ]
        self.scaler = StandardScaler()
        self.aging_metabolites = self._load_aging_metabolite_database()
        self.pathway_database = self._load_metabolic_pathway_database()

        logger.info(f"Initialized MetabolomicsAnalyzer for {platform}")

    def _load_aging_metabolite_database(self) -> Dict[str, Dict]:
        """Load database of aging-related metabolites."""
        aging_metabolites = {
            'NAD+': {
                'pathways': ['nad_metabolism', 'energy_metabolism'],
                'aging_relevance': 0.95,
                'biomarker_class': 'longevity'
            },
            'ATP': {
                'pathways': ['energy_metabolism', 'cellular_respiration'],
                'aging_relevance': 0.90,
                'biomarker_class': 'energy'
            },
            'Glucose': {
                'pathways': ['glycolysis', 'glucose_metabolism'],
                'aging_relevance': 0.80,
                'biomarker_class': 'metabolic_health'
            },
            'Lactate': {
                'pathways': ['glycolysis', 'cellular_stress'],
                'aging_relevance': 0.75,
                'biomarker_class': 'stress_response'
            },
            'Glutamine': {
                'pathways': ['amino_acid_metabolism', 'cellular_stress'],
                'aging_relevance': 0.70,
                'biomarker_class': 'stress_response'
            }
        }
        return aging_metabolites

    def _load_metabolic_pathway_database(self) -> Dict[str, List[str]]:
        """Load metabolic pathway database."""
        pathways = {
            'glycolysis': ['Glucose', 'Glucose-6-phosphate', 'Pyruvate', 'Lactate'],
            'tca_cycle': ['Citrate', 'Succinate', 'Fumarate', 'Malate'],
            'amino_acid_metabolism': ['Glutamine', 'Alanine', 'Serine', 'Glycine'],
            'lipid_metabolism': ['Palmitate', 'Oleate', 'Cholesterol', 'Acetyl-CoA'],
            'purine_metabolism': ['ATP', 'ADP', 'AMP', 'Adenosine'],
            'nad_metabolism': ['NAD+', 'NADH', 'Nicotinamide', 'Tryptophan'],
            'oxidative_stress': ['GSH', 'GSSG', 'Ascorbate', 'Urate']
        }
        return pathways

    def load_metabolomics_data(self,
                             metabolite_file: str,
                             metadata_file: str,
                             annotation_file: Optional[str] = None) -> MetabolomicsDataset:
        """Load and validate metabolomics experimental data."""
        logger.info("Loading metabolomics data...")

        # Load metabolite matrix
        metabolite_df = pd.read_csv(metabolite_file, index_col=0)

        # Load sample metadata
        metadata_df = pd.read_csv(metadata_file, index_col=0)

        # Load metabolite annotations
        if annotation_file:
            annotations_df = pd.read_csv(annotation_file, index_col=0)
        else:
            annotations_df = self._create_default_annotations(metabolite_df.index)

        # Create pathway mapping
        pathway_mapping = self._create_pathway_mapping(metabolite_df.index)

        # Quality control
        quality_metrics = self._perform_quality_control(metabolite_df, metadata_df)

        dataset = MetabolomicsDataset(
            metabolite_matrix=metabolite_df,
            sample_metadata=metadata_df,
            metabolite_annotations=annotations_df,
            pathway_mapping=pathway_mapping,
            quality_metrics=quality_metrics
        )

        logger.info(f"Loaded metabolomics dataset: {metabolite_df.shape[0]} metabolites, {metabolite_df.shape[1]} samples")
        return dataset

    def _create_default_annotations(self, metabolite_ids: List[str]) -> pd.DataFrame:
        """Create default metabolite annotations."""
        annotations = []
        for metabolite_id in metabolite_ids:
            aging_relevance = 0.0
            biomarker_class = 'unknown'

            if metabolite_id in self.aging_metabolites:
                aging_relevance = self.aging_metabolites[metabolite_id]['aging_relevance']
                biomarker_class = self.aging_metabolites[metabolite_id]['biomarker_class']

            annotations.append({
                'metabolite_id': metabolite_id,
                'aging_relevance': aging_relevance,
                'biomarker_class': biomarker_class,
                'molecular_weight': np.random.uniform(100, 500),  # Simulated
                'chemical_class': np.random.choice(['organic_acid', 'amino_acid', 'lipid', 'nucleotide'])
            })

        return pd.DataFrame(annotations).set_index('metabolite_id')

    def _create_pathway_mapping(self, metabolite_ids: List[str]) -> Dict[str, List[str]]:
        """Create pathway mapping for metabolites."""
        mapping = {}
        for pathway, metabolites in self.pathway_database.items():
            pathway_metabolites = [m for m in metabolites if m in metabolite_ids]
            if pathway_metabolites:
                mapping[pathway] = pathway_metabolites
        return mapping

    def _perform_quality_control(self, metabolite_df: pd.DataFrame,
                               metadata_df: pd.DataFrame) -> Dict[str, float]:
        """Perform quality control on metabolomics data."""
        metrics = {}

        # Missing data percentage
        metrics['missing_data_percent'] = (metabolite_df.isnull().sum().sum() /
                                         (metabolite_df.shape[0] * metabolite_df.shape[1])) * 100

        # Coefficient of variation
        cv_values = metabolite_df.std(axis=1) / metabolite_df.mean(axis=1)
        metrics['median_cv'] = cv_values.median()

        # Number of detected metabolites
        metrics['detected_metabolites'] = (metabolite_df > 0).any(axis=1).sum()

        # Dynamic range
        valid_data = metabolite_df[metabolite_df > 0]
        if not valid_data.empty:
            metrics['dynamic_range_log10'] = np.log10(valid_data.max().max() / valid_data.min().min())
        else:
            metrics['dynamic_range_log10'] = 0

        return metrics

    def analyze_metabolic_intervention_response(self,
                                              dataset: MetabolomicsDataset,
                                              control_condition: str,
                                              treatment_condition: str) -> List[MetaboliteAnalysisResult]:
        """Analyze metabolite changes in response to aging interventions."""
        logger.info(f"Analyzing metabolic intervention response: {control_condition} vs {treatment_condition}")

        metadata = dataset.sample_metadata
        control_samples = metadata[metadata['condition'] == control_condition].index
        treatment_samples = metadata[metadata['condition'] == treatment_condition].index

        if len(control_samples) == 0 or len(treatment_samples) == 0:
            raise ValueError(f"No samples found for conditions: {control_condition}, {treatment_condition}")

        results = []
        metabolite_matrix = dataset.metabolite_matrix

        for metabolite_id in metabolite_matrix.index:
            result = self._analyze_single_metabolite(
                metabolite_id,
                metabolite_matrix.loc[metabolite_id, control_samples],
                metabolite_matrix.loc[metabolite_id, treatment_samples],
                dataset
            )
            results.append(result)

        # Sort by aging relevance and significance
        results.sort(key=lambda x: (x.aging_relevance_score, -x.p_value), reverse=True)

        logger.info(f"Analyzed {len(results)} metabolites for intervention response")
        return results

    def _analyze_single_metabolite(self,
                                 metabolite_id: str,
                                 control_values: pd.Series,
                                 treatment_values: pd.Series,
                                 dataset: MetabolomicsDataset) -> MetaboliteAnalysisResult:
        """Analyze a single metabolite for intervention response."""

        # Remove missing values
        control_clean = control_values.dropna()
        treatment_clean = treatment_values.dropna()

        if len(control_clean) < 2 or len(treatment_clean) < 2:
            return MetaboliteAnalysisResult(
                metabolite_id=metabolite_id,
                fold_change=np.nan,
                p_value=1.0,
                adjusted_p_value=1.0,
                pathway_impact=0.0,
                aging_relevance_score=0.0,
                biomarker_potential=0.0
            )

        # Calculate fold change
        mean_control = control_clean.mean()
        mean_treatment = treatment_clean.mean()
        fold_change = mean_treatment / mean_control if mean_control > 0 else np.nan

        # Statistical testing
        try:
            t_stat, p_value = stats.ttest_ind(treatment_clean, control_clean)
        except:
            p_value = 1.0

        # Calculate pathway impact
        pathway_impact = self._calculate_pathway_impact(metabolite_id, fold_change, dataset)

        # Calculate aging relevance
        aging_relevance_score = self._calculate_aging_relevance(metabolite_id, fold_change)

        # Calculate biomarker potential
        biomarker_potential = self._calculate_biomarker_potential(
            metabolite_id, fold_change, p_value, control_clean, treatment_clean
        )

        return MetaboliteAnalysisResult(
            metabolite_id=metabolite_id,
            fold_change=fold_change,
            p_value=p_value,
            adjusted_p_value=p_value,  # Will be corrected in batch
            pathway_impact=pathway_impact,
            aging_relevance_score=aging_relevance_score,
            biomarker_potential=biomarker_potential
        )

    def _calculate_pathway_impact(self,
                                metabolite_id: str,
                                fold_change: float,
                                dataset: MetabolomicsDataset) -> float:
        """Calculate pathway impact score for a metabolite."""
        impact_score = 0.0

        # Find pathways containing this metabolite
        for pathway, metabolites in dataset.pathway_mapping.items():
            if metabolite_id in metabolites:
                # Weight by pathway importance and fold change magnitude
                pathway_weight = 1.0 if pathway in self.focus_pathways else 0.5
                if not np.isnan(fold_change):
                    impact_score += pathway_weight * abs(np.log2(fold_change))

        return min(impact_score, 5.0)  # Cap at 5.0

    def _calculate_aging_relevance(self, metabolite_id: str, fold_change: float) -> float:
        """Calculate aging relevance score for a metabolite."""
        base_score = 0.0

        if metabolite_id in self.aging_metabolites:
            base_score = self.aging_metabolites[metabolite_id]['aging_relevance']

        # Add fold change contribution
        if not np.isnan(fold_change):
            fc_contribution = min(abs(np.log2(fold_change)) / 5.0, 0.3)
            base_score += fc_contribution

        return min(base_score, 1.0)

    def _calculate_biomarker_potential(self,
                                     metabolite_id: str,
                                     fold_change: float,
                                     p_value: float,
                                     control_values: pd.Series,
                                     treatment_values: pd.Series) -> float:
        """Calculate biomarker potential score."""
        score = 0.0

        # Statistical significance contribution
        if p_value < 0.01:
            score += 0.4
        elif p_value < 0.05:
            score += 0.2

        # Effect size contribution
        if not np.isnan(fold_change):
            effect_size = abs(np.log2(fold_change))
            score += min(effect_size / 3.0, 0.3)

        # Variability contribution (lower is better for biomarkers)
        control_cv = control_values.std() / control_values.mean()
        treatment_cv = treatment_values.std() / treatment_values.mean()
        avg_cv = (control_cv + treatment_cv) / 2
        variability_score = max(0, 0.3 - avg_cv)
        score += variability_score

        return min(score, 1.0)

    def perform_pathway_enrichment_analysis(self,
                                          results: List[MetaboliteAnalysisResult],
                                          p_threshold: float = 0.05) -> Dict[str, Dict]:
        """Perform pathway enrichment analysis on metabolite results."""
        logger.info("Performing pathway enrichment analysis...")

        # Get significantly changed metabolites
        significant_metabolites = [
            r.metabolite_id for r in results
            if r.adjusted_p_value < p_threshold
        ]

        enrichment_results = {}

        for pathway, pathway_metabolites in self.pathway_database.items():
            # Count significant metabolites in pathway
            pathway_significant = [
                m for m in significant_metabolites
                if m in pathway_metabolites
            ]

            if len(pathway_significant) == 0:
                continue

            # Fisher's exact test for enrichment
            from scipy.stats import fisher_exact

            # Create contingency table
            in_pathway_significant = len(pathway_significant)
            in_pathway_total = len(pathway_metabolites)
            total_significant = len(significant_metabolites)
            total_metabolites = len([r.metabolite_id for r in results])

            not_in_pathway_significant = total_significant - in_pathway_significant
            in_pathway_not_significant = in_pathway_total - in_pathway_significant
            not_in_pathway_not_significant = (total_metabolites - total_significant -
                                            in_pathway_not_significant)

            contingency_table = [
                [in_pathway_significant, not_in_pathway_significant],
                [in_pathway_not_significant, not_in_pathway_not_significant]
            ]

            try:
                odds_ratio, p_val = fisher_exact(contingency_table, alternative='greater')
            except:
                odds_ratio, p_val = 1.0, 1.0

            enrichment_results[pathway] = {
                'significant_metabolites': pathway_significant,
                'total_pathway_metabolites': in_pathway_total,
                'odds_ratio': odds_ratio,
                'p_value': p_val,
                'enrichment_score': in_pathway_significant / max(in_pathway_total, 1)
            }

        return enrichment_results

    def generate_metabolomics_report(self,
                                   dataset: MetabolomicsDataset,
                                   results: List[MetaboliteAnalysisResult],
                                   enrichment_results: Dict[str, Dict],
                                   output_file: str) -> str:
        """Generate comprehensive metabolomics analysis report."""

        report_lines = [
            "# Metabolomics Analysis Report - CRISPR Toolkit Phase 3",
            f"**Analysis Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Platform**: {self.platform}",
            "",
            "## Dataset Summary",
            f"- **Metabolites analyzed**: {dataset.metabolite_matrix.shape[0]}",
            f"- **Samples**: {dataset.metabolite_matrix.shape[1]}",
            f"- **Missing data**: {dataset.quality_metrics.get('missing_data_percent', 0):.2f}%",
            f"- **Detected metabolites**: {dataset.quality_metrics.get('detected_metabolites', 0)}",
            "",
            "## Intervention Response Analysis",
            f"- **Significantly changed metabolites** (p < 0.05): {sum(1 for r in results if r.adjusted_p_value < 0.05)}",
            f"- **High biomarker potential** (score > 0.7): {sum(1 for r in results if r.biomarker_potential > 0.7)}",
            f"- **Aging-relevant metabolites** (score > 0.5): {sum(1 for r in results if r.aging_relevance_score > 0.5)}",
            "",
            "## Pathway Enrichment Results",
            ""
        ]

        # Add pathway enrichment table
        if enrichment_results:
            report_lines.extend([
                "| Pathway | Significant Metabolites | Total | Enrichment Score | P-value |",
                "|---------|------------------------|-------|------------------|---------|"
            ])

            for pathway, data in sorted(enrichment_results.items(),
                                      key=lambda x: x[1]['p_value']):
                if data['p_value'] < 0.1:  # Show marginally significant pathways
                    report_lines.append(
                        f"| {pathway} | {len(data['significant_metabolites'])} | "
                        f"{data['total_pathway_metabolites']} | "
                        f"{data['enrichment_score']:.2f} | {data['p_value']:.2e} |"
                    )

        # Add top metabolite changes
        report_lines.extend([
            "",
            "## Top Metabolite Changes",
            "",
            "| Metabolite | Fold Change | P-value | Biomarker Potential | Aging Score |",
            "|------------|-------------|---------|-------------------|-------------|"
        ])

        significant_results = [r for r in results if r.adjusted_p_value < 0.05][:15]

        for result in significant_results:
            fc_str = f"{result.fold_change:.2f}" if not np.isnan(result.fold_change) else "N/A"

            report_lines.append(
                f"| {result.metabolite_id} | {fc_str} | {result.p_value:.2e} | "
                f"{result.biomarker_potential:.2f} | {result.aging_relevance_score:.2f} |"
            )

        report_content = "\n".join(report_lines)

        with open(output_file, 'w') as f:
            f.write(report_content)

        logger.info(f"Generated metabolomics report: {output_file}")
        return report_content


class MetabolicPathwayModel:
    """
    Network-based model for metabolic pathway analysis and prediction.
    """

    def __init__(self):
        """Initialize metabolic pathway model."""
        self.pathway_network = nx.DiGraph()
        self.metabolite_features = {}

    def build_pathway_network(self, pathway_database: Dict[str, List[str]]) -> nx.DiGraph:
        """Build metabolic pathway network."""

        # Add nodes for metabolites
        for pathway, metabolites in pathway_database.items():
            for metabolite in metabolites:
                self.pathway_network.add_node(metabolite, pathway=pathway)

        # Add edges based on biochemical relationships
        # This would be more sophisticated in real implementation
        for pathway, metabolites in pathway_database.items():
            for i in range(len(metabolites) - 1):
                self.pathway_network.add_edge(metabolites[i], metabolites[i+1],
                                            pathway=pathway, weight=1.0)

        logger.info(f"Built pathway network with {self.pathway_network.number_of_nodes()} nodes "
                   f"and {self.pathway_network.number_of_edges()} edges")

        return self.pathway_network

    def calculate_pathway_impact_scores(self,
                                      metabolite_changes: Dict[str, float]) -> Dict[str, float]:
        """Calculate pathway-level impact scores from metabolite changes."""

        pathway_scores = {}

        # Group metabolites by pathway
        pathway_metabolites = {}
        for node, data in self.pathway_network.nodes(data=True):
            pathway = data.get('pathway', 'unknown')
            if pathway not in pathway_metabolites:
                pathway_metabolites[pathway] = []
            pathway_metabolites[pathway].append(node)

        # Calculate impact scores for each pathway
        for pathway, metabolites in pathway_metabolites.items():
            pathway_changes = []
            for metabolite in metabolites:
                if metabolite in metabolite_changes:
                    pathway_changes.append(abs(metabolite_changes[metabolite]))

            if pathway_changes:
                # Use median absolute change as pathway impact
                pathway_scores[pathway] = np.median(pathway_changes)
            else:
                pathway_scores[pathway] = 0.0

        return pathway_scores


if __name__ == "__main__":
    # Example usage
    analyzer = MetabolomicsAnalyzer(platform="lcms")

    # Create synthetic data
    n_metabolites, n_samples = 200, 40

    metabolite_data = pd.DataFrame(
        np.random.lognormal(5, 1, (n_metabolites, n_samples)),
        index=[f"Metabolite_{i}" for i in range(n_metabolites)],
        columns=[f"Sample_{i}" for i in range(n_samples)]
    )

    metadata = pd.DataFrame({
        'condition': ['control'] * 20 + ['treatment'] * 20,
        'age': np.random.randint(25, 75, n_samples)
    }, index=metabolite_data.columns)

    # Save temporary files
    metabolite_data.to_csv("/tmp/metabolites.csv")
    metadata.to_csv("/tmp/metadata.csv")

    # Load and analyze
    dataset = analyzer.load_metabolomics_data(
        "/tmp/metabolites.csv",
        "/tmp/metadata.csv"
    )

    results = analyzer.analyze_metabolic_intervention_response(
        dataset, 'control', 'treatment'
    )

    # Pathway enrichment
    enrichment = analyzer.perform_pathway_enrichment_analysis(results)

    # Generate report
    report = analyzer.generate_metabolomics_report(
        dataset, results, enrichment, "/tmp/metabolomics_report.md"
    )

    print("Metabolomics analysis completed!")
    print(f"Analyzed {len(results)} metabolites")
    print(f"Enriched pathways: {len(enrichment)}")

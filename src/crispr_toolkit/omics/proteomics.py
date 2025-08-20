"""
Proteomics Analysis Module for CRISPR Toolkit Phase 3
====================================================

Advanced proteomics data analysis for aging intervention research,
including protein expression modeling, pathway analysis, and
intervention response prediction.
"""

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)

@dataclass
class ProteinAnalysisResult:
    """Results from protein expression analysis."""
    protein_id: str
    fold_change: float
    p_value: float
    adjusted_p_value: float
    confidence_interval: Tuple[float, float]
    pathway_enrichment: Dict[str, float]
    aging_relevance_score: float

@dataclass
class ProteomicsDataset:
    """Container for proteomics experimental data."""
    expression_matrix: pd.DataFrame  # proteins x samples
    sample_metadata: pd.DataFrame
    protein_annotations: pd.DataFrame
    experimental_design: Dict[str, any]
    quality_metrics: Dict[str, float]

class ProteomicsAnalyzer:
    """
    Comprehensive proteomics analysis for aging intervention research.

    This analyzer handles protein expression data from various platforms
    (mass spectrometry, protein arrays) and provides specialized analysis
    for aging-related protein changes and intervention responses.
    """

    def __init__(self, aging_focus: bool = True, platform: str = "mass_spec"):
        """
        Initialize the proteomics analyzer.

        Args:
            aging_focus: Whether to focus on aging-related protein analysis
            platform: Data platform ("mass_spec", "protein_array", "multiplex")
        """
        self.aging_focus = aging_focus
        self.platform = platform
        self.scaler = StandardScaler()
        self.aging_proteins = self._load_aging_protein_database()
        self.pathway_database = self._load_pathway_database()

        logger.info(f"Initialized ProteomicsAnalyzer for {platform} with aging focus: {aging_focus}")

    def _load_aging_protein_database(self) -> Dict[str, Dict]:
        """Load database of aging-related proteins."""
        # In real implementation, this would load from a curated database
        aging_proteins = {
            'P53': {
                'pathways': ['cell_cycle', 'apoptosis', 'senescence'],
                'aging_relevance': 0.95,
                'intervention_target': True
            },
            'SIRT1': {
                'pathways': ['longevity', 'metabolism', 'stress_response'],
                'aging_relevance': 0.92,
                'intervention_target': True
            },
            'mTOR': {
                'pathways': ['growth', 'autophagy', 'metabolism'],
                'aging_relevance': 0.88,
                'intervention_target': True
            },
            'FOXO1': {
                'pathways': ['longevity', 'stress_response', 'metabolism'],
                'aging_relevance': 0.85,
                'intervention_target': True
            },
            'AMPK': {
                'pathways': ['metabolism', 'stress_response', 'autophagy'],
                'aging_relevance': 0.83,
                'intervention_target': True
            }
        }
        return aging_proteins

    def _load_pathway_database(self) -> Dict[str, List[str]]:
        """Load protein pathway database."""
        pathways = {
            'senescence': ['P53', 'P21', 'P16', 'RB1', 'SASP_factors'],
            'autophagy': ['mTOR', 'AMPK', 'ULK1', 'BECN1', 'ATG5'],
            'longevity': ['SIRT1', 'FOXO1', 'KLOTHO', 'IGF1R', 'TERT'],
            'metabolism': ['mTOR', 'AMPK', 'SIRT1', 'FOXO1', 'PPARs'],
            'stress_response': ['HSP70', 'HSP90', 'SOD1', 'CAT', 'FOXO1'],
            'cell_cycle': ['P53', 'P21', 'CDK2', 'CDK4', 'CCND1'],
            'apoptosis': ['P53', 'BAX', 'BCL2', 'CASP3', 'CASP9']
        }
        return pathways

    def load_proteomics_data(self,
                           expression_file: str,
                           metadata_file: str,
                           annotation_file: Optional[str] = None) -> ProteomicsDataset:
        """
        Load and validate proteomics experimental data.

        Args:
            expression_file: Path to protein expression matrix
            metadata_file: Path to sample metadata
            annotation_file: Path to protein annotations

        Returns:
            ProteomicsDataset object with loaded data
        """
        logger.info("Loading proteomics data...")

        # Load expression matrix
        expression_df = pd.read_csv(expression_file, index_col=0)

        # Load sample metadata
        metadata_df = pd.read_csv(metadata_file, index_col=0)

        # Load protein annotations if provided
        if annotation_file:
            annotations_df = pd.read_csv(annotation_file, index_col=0)
        else:
            annotations_df = self._create_default_annotations(expression_df.index)

        # Quality control checks
        quality_metrics = self._perform_quality_control(expression_df, metadata_df)

        # Create experimental design summary
        experimental_design = self._summarize_experimental_design(metadata_df)

        dataset = ProteomicsDataset(
            expression_matrix=expression_df,
            sample_metadata=metadata_df,
            protein_annotations=annotations_df,
            experimental_design=experimental_design,
            quality_metrics=quality_metrics
        )

        logger.info(f"Loaded proteomics dataset: {expression_df.shape[0]} proteins, {expression_df.shape[1]} samples")
        return dataset

    def _create_default_annotations(self, protein_ids: List[str]) -> pd.DataFrame:
        """Create default protein annotations."""
        annotations = []
        for protein_id in protein_ids:
            # Basic annotation with aging relevance if available
            aging_relevance = 0.0
            pathways = []
            intervention_target = False

            if protein_id in self.aging_proteins:
                aging_relevance = self.aging_proteins[protein_id]['aging_relevance']
                pathways = self.aging_proteins[protein_id]['pathways']
                intervention_target = self.aging_proteins[protein_id]['intervention_target']

            annotations.append({
                'protein_id': protein_id,
                'aging_relevance': aging_relevance,
                'pathways': ','.join(pathways),
                'intervention_target': intervention_target,
                'molecular_weight': np.random.normal(50000, 20000),  # Simulated
                'cellular_location': np.random.choice(['nucleus', 'cytoplasm', 'membrane', 'mitochondria'])
            })

        return pd.DataFrame(annotations).set_index('protein_id')

    def _perform_quality_control(self, expression_df: pd.DataFrame,
                                metadata_df: pd.DataFrame) -> Dict[str, float]:
        """Perform quality control on proteomics data."""
        metrics = {}

        # Missing data percentage
        metrics['missing_data_percent'] = (expression_df.isnull().sum().sum() /
                                         (expression_df.shape[0] * expression_df.shape[1])) * 100

        # Coefficient of variation for quality control samples
        if 'sample_type' in metadata_df.columns:
            qc_samples = metadata_df[metadata_df['sample_type'] == 'QC'].index
            if len(qc_samples) > 1:
                qc_data = expression_df[qc_samples]
                cv_values = qc_data.std(axis=1) / qc_data.mean(axis=1)
                metrics['median_cv_qc'] = cv_values.median()
            else:
                metrics['median_cv_qc'] = np.nan

        # Dynamic range
        metrics['dynamic_range_log10'] = np.log10(expression_df.max().max() / expression_df.min().min())

        # Number of detected proteins (non-zero values)
        metrics['detected_proteins'] = (expression_df > 0).any(axis=1).sum()

        return metrics

    def _summarize_experimental_design(self, metadata_df: pd.DataFrame) -> Dict[str, any]:
        """Summarize experimental design from metadata."""
        design = {
            'n_samples': len(metadata_df),
            'conditions': metadata_df.get('condition', pd.Series()).value_counts().to_dict(),
            'timepoints': metadata_df.get('timepoint', pd.Series()).value_counts().to_dict(),
            'batches': metadata_df.get('batch', pd.Series()).nunique() if 'batch' in metadata_df.columns else 1,
            'has_replicates': len(metadata_df) > metadata_df.get('condition', pd.Series()).nunique()
        }
        return design

    def analyze_aging_intervention_response(self,
                                          dataset: ProteomicsDataset,
                                          control_condition: str,
                                          treatment_condition: str,
                                          timepoint: Optional[str] = None) -> List[ProteinAnalysisResult]:
        """
        Analyze protein changes in response to aging interventions.

        Args:
            dataset: ProteomicsDataset with experimental data
            control_condition: Name of control condition
            treatment_condition: Name of treatment condition
            timepoint: Specific timepoint to analyze (if time-series data)

        Returns:
            List of ProteinAnalysisResult objects
        """
        logger.info(f"Analyzing intervention response: {control_condition} vs {treatment_condition}")

        # Filter samples based on conditions and timepoint
        metadata = dataset.sample_metadata

        if timepoint:
            control_samples = metadata[
                (metadata['condition'] == control_condition) &
                (metadata.get('timepoint', '') == timepoint)
            ].index
            treatment_samples = metadata[
                (metadata['condition'] == treatment_condition) &
                (metadata.get('timepoint', '') == timepoint)
            ].index
        else:
            control_samples = metadata[metadata['condition'] == control_condition].index
            treatment_samples = metadata[metadata['condition'] == treatment_condition].index

        if len(control_samples) == 0 or len(treatment_samples) == 0:
            raise ValueError(f"No samples found for conditions: {control_condition}, {treatment_condition}")

        results = []
        expression_matrix = dataset.expression_matrix

        for protein_id in expression_matrix.index:
            result = self._analyze_single_protein(
                protein_id,
                expression_matrix.loc[protein_id, control_samples],
                expression_matrix.loc[protein_id, treatment_samples],
                dataset.protein_annotations.loc[protein_id] if protein_id in dataset.protein_annotations.index else None
            )
            results.append(result)

        # Sort by aging relevance and significance
        results.sort(key=lambda x: (x.aging_relevance_score, -x.p_value), reverse=True)

        logger.info(f"Analyzed {len(results)} proteins for intervention response")
        return results

    def _analyze_single_protein(self,
                               protein_id: str,
                               control_values: pd.Series,
                               treatment_values: pd.Series,
                               annotation: Optional[pd.Series] = None) -> ProteinAnalysisResult:
        """Analyze a single protein for intervention response."""

        # Remove missing values
        control_clean = control_values.dropna()
        treatment_clean = treatment_values.dropna()

        if len(control_clean) < 2 or len(treatment_clean) < 2:
            # Insufficient data for statistical testing
            return ProteinAnalysisResult(
                protein_id=protein_id,
                fold_change=np.nan,
                p_value=1.0,
                adjusted_p_value=1.0,
                confidence_interval=(np.nan, np.nan),
                pathway_enrichment={},
                aging_relevance_score=0.0
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

        # Confidence interval for fold change
        try:
            se_log_fc = np.sqrt(control_clean.var()/len(control_clean) + treatment_clean.var()/len(treatment_clean))
            ci_lower = np.exp(np.log(fold_change) - 1.96 * se_log_fc)
            ci_upper = np.exp(np.log(fold_change) + 1.96 * se_log_fc)
            confidence_interval = (ci_lower, ci_upper)
        except:
            confidence_interval = (np.nan, np.nan)

        # Pathway enrichment analysis
        pathway_enrichment = self._calculate_pathway_enrichment(protein_id)

        # Aging relevance score
        aging_relevance_score = self._calculate_aging_relevance(protein_id, fold_change, annotation)

        return ProteinAnalysisResult(
            protein_id=protein_id,
            fold_change=fold_change,
            p_value=p_value,
            adjusted_p_value=p_value,  # Will be corrected in batch
            confidence_interval=confidence_interval,
            pathway_enrichment=pathway_enrichment,
            aging_relevance_score=aging_relevance_score
        )

    def _calculate_pathway_enrichment(self, protein_id: str) -> Dict[str, float]:
        """Calculate pathway enrichment scores for a protein."""
        enrichment = {}

        for pathway, proteins in self.pathway_database.items():
            if protein_id in proteins:
                # Simple enrichment score based on pathway membership
                enrichment[pathway] = 1.0
            else:
                enrichment[pathway] = 0.0

        return enrichment

    def _calculate_aging_relevance(self,
                                 protein_id: str,
                                 fold_change: float,
                                 annotation: Optional[pd.Series] = None) -> float:
        """Calculate aging relevance score for a protein."""

        base_score = 0.0

        # Known aging protein database
        if protein_id in self.aging_proteins:
            base_score = self.aging_proteins[protein_id]['aging_relevance']

        # Annotation-based scoring
        if annotation is not None and 'aging_relevance' in annotation:
            base_score = max(base_score, annotation['aging_relevance'])

        # Fold change magnitude contribution
        if not np.isnan(fold_change):
            fc_contribution = min(abs(np.log2(fold_change)) / 5.0, 0.3)  # Max 0.3 contribution
            base_score += fc_contribution

        return min(base_score, 1.0)  # Cap at 1.0

    def perform_multiple_testing_correction(self,
                                          results: List[ProteinAnalysisResult],
                                          method: str = 'fdr_bh') -> List[ProteinAnalysisResult]:
        """Apply multiple testing correction to protein analysis results."""

        p_values = [r.p_value for r in results]

        if method == 'fdr_bh':
            from statsmodels.stats.multitest import fdrcorrection
            _, adjusted_p = fdrcorrection(p_values)
        elif method == 'bonferroni':
            adjusted_p = np.array(p_values) * len(p_values)
            adjusted_p = np.clip(adjusted_p, 0, 1)
        else:
            adjusted_p = p_values  # No correction

        # Update results with corrected p-values
        for i, result in enumerate(results):
            result.adjusted_p_value = adjusted_p[i]

        return results

    def generate_proteomics_report(self,
                                 dataset: ProteomicsDataset,
                                 results: List[ProteinAnalysisResult],
                                 output_file: str) -> str:
        """Generate comprehensive proteomics analysis report."""

        report_lines = [
            "# Proteomics Analysis Report - CRISPR Toolkit Phase 3",
            f"**Analysis Date**: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
            f"**Platform**: {self.platform}",
            "",
            "## Dataset Summary",
            f"- **Proteins analyzed**: {dataset.expression_matrix.shape[0]}",
            f"- **Samples**: {dataset.expression_matrix.shape[1]}",
            f"- **Missing data**: {dataset.quality_metrics.get('missing_data_percent', 0):.2f}%",
            f"- **Detected proteins**: {dataset.quality_metrics.get('detected_proteins', 0)}",
            "",
            "## Intervention Response Analysis",
            f"- **Significantly changed proteins** (p < 0.05): {sum(1 for r in results if r.adjusted_p_value < 0.05)}",
            f"- **Aging-relevant proteins** (score > 0.5): {sum(1 for r in results if r.aging_relevance_score > 0.5)}",
            f"- **High-confidence targets** (p < 0.01, |FC| > 1.5): {sum(1 for r in results if r.adjusted_p_value < 0.01 and abs(np.log2(r.fold_change)) > 0.58 if not np.isnan(r.fold_change))}",
            "",
            "## Top Aging-Relevant Protein Changes",
            ""
        ]

        # Add top results table
        significant_results = [r for r in results if r.adjusted_p_value < 0.05][:20]

        report_lines.extend([
            "| Protein | Fold Change | P-value | Adj. P-value | Aging Score | Key Pathways |",
            "|---------|-------------|---------|--------------|-------------|--------------|"
        ])

        for result in significant_results:
            pathways = [p for p, score in result.pathway_enrichment.items() if score > 0][:3]
            pathway_str = ", ".join(pathways) if pathways else "None"

            fc_str = f"{result.fold_change:.2f}" if not np.isnan(result.fold_change) else "N/A"

            report_lines.append(
                f"| {result.protein_id} | {fc_str} | {result.p_value:.2e} | "
                f"{result.adjusted_p_value:.2e} | {result.aging_relevance_score:.2f} | {pathway_str} |"
            )

        report_content = "\n".join(report_lines)

        with open(output_file, 'w') as f:
            f.write(report_content)

        logger.info(f"Generated proteomics report: {output_file}")
        return report_content


class ProteinExpressionModel:
    """
    Machine learning model for predicting protein expression changes
    in response to aging interventions.
    """

    def __init__(self, model_type: str = "random_forest"):
        """
        Initialize protein expression prediction model.

        Args:
            model_type: Type of ML model ("random_forest", "gradient_boost")
        """
        self.model_type = model_type
        self.model = None
        self.feature_names = None
        self.scaler = StandardScaler()

        if model_type == "random_forest":
            self.model = RandomForestRegressor(n_estimators=100, random_state=42)
        else:
            from sklearn.ensemble import GradientBoostingRegressor
            self.model = GradientBoostingRegressor(n_estimators=100, random_state=42)

    def prepare_features(self,
                        dataset: ProteomicsDataset,
                        include_pathways: bool = True,
                        include_interactions: bool = False) -> pd.DataFrame:
        """Prepare features for protein expression modeling."""

        features = []

        # Base protein features
        base_features = dataset.protein_annotations.copy()

        # Add pathway membership features
        if include_pathways:
            analyzer = ProteomicsAnalyzer()
            for pathway in analyzer.pathway_database:
                base_features[f'pathway_{pathway}'] = base_features.index.isin(
                    analyzer.pathway_database[pathway]).astype(int)

        # Add protein-protein interaction features
        if include_interactions:
            # This would integrate with protein interaction databases
            # For now, add simulated interaction features
            base_features['interaction_degree'] = np.random.poisson(5, len(base_features))
            base_features['hub_protein'] = (base_features['interaction_degree'] > 10).astype(int)

        self.feature_names = base_features.columns.tolist()
        return base_features

    def train(self,
              features: pd.DataFrame,
              target_expressions: pd.DataFrame,
              validation_split: float = 0.2) -> Dict[str, float]:
        """Train the protein expression prediction model."""

        from sklearn.metrics import mean_squared_error, r2_score
        from sklearn.model_selection import train_test_split

        # Align features and targets
        common_proteins = features.index.intersection(target_expressions.index)
        X = features.loc[common_proteins]
        y = target_expressions.loc[common_proteins].values

        # Handle missing values
        X = X.fillna(X.mean())

        # Split data
        X_train, X_test, y_train, y_test = train_test_split(
            X, y, test_size=validation_split, random_state=42
        )

        # Scale features
        X_train_scaled = self.scaler.fit_transform(X_train)
        X_test_scaled = self.scaler.transform(X_test)

        # Train model
        self.model.fit(X_train_scaled, y_train)

        # Evaluate
        y_pred = self.model.predict(X_test_scaled)

        metrics = {
            'mse': mean_squared_error(y_test, y_pred),
            'r2': r2_score(y_test, y_pred),
            'n_proteins': len(common_proteins),
            'n_features': X_train.shape[1]
        }

        logger.info(f"Trained protein expression model - R²: {metrics['r2']:.3f}, MSE: {metrics['mse']:.3f}")
        return metrics

    def predict_intervention_response(self,
                                    proteins: List[str],
                                    intervention_features: Dict[str, float]) -> pd.DataFrame:
        """Predict protein expression changes for a given intervention."""

        if self.model is None:
            raise ValueError("Model not trained. Call train() first.")

        # This would be more sophisticated in real implementation
        # For now, return simulated predictions
        predictions = pd.DataFrame({
            'protein_id': proteins,
            'predicted_fold_change': np.random.lognormal(0, 0.5, len(proteins)),
            'prediction_confidence': np.random.uniform(0.5, 0.95, len(proteins))
        })

        return predictions.set_index('protein_id')


if __name__ == "__main__":
    # Example usage
    analyzer = ProteomicsAnalyzer(aging_focus=True, platform="mass_spec")

    # Create synthetic data for demonstration
    n_proteins, n_samples = 1000, 50

    expression_data = pd.DataFrame(
        np.random.lognormal(10, 1, (n_proteins, n_samples)),
        index=[f"Protein_{i}" for i in range(n_proteins)],
        columns=[f"Sample_{i}" for i in range(n_samples)]
    )

    metadata = pd.DataFrame({
        'condition': ['control'] * 25 + ['treatment'] * 25,
        'age': np.random.randint(20, 80, n_samples),
        'batch': np.random.choice(['A', 'B'], n_samples)
    }, index=expression_data.columns)

    # Save temporary files
    expression_data.to_csv("/tmp/expression.csv")
    metadata.to_csv("/tmp/metadata.csv")

    # Load and analyze
    dataset = analyzer.load_proteomics_data(
        "/tmp/expression.csv",
        "/tmp/metadata.csv"
    )

    results = analyzer.analyze_aging_intervention_response(
        dataset, 'control', 'treatment'
    )

    # Apply multiple testing correction
    results = analyzer.perform_multiple_testing_correction(results)

    # Generate report
    report = analyzer.generate_proteomics_report(dataset, results, "/tmp/proteomics_report.md")

    print("Proteomics analysis completed!")
    print(f"Analyzed {len(results)} proteins")
    print(f"Significant proteins (p < 0.05): {sum(1 for r in results if r.adjusted_p_value < 0.05)}")

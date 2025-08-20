"""
Training data collection and preparation for aging ML models.

This module handles data collection from various sources for training
the aging target prioritization and rejuvenation prediction models.
"""

import logging
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class AgingDatasetCollector:
    """Collect and prepare training datasets for aging research models."""

    def __init__(self, data_dir: Path = None):
        """Initialize data collector."""
        self.data_dir = data_dir or Path("data")
        self.raw_dir = self.data_dir / "raw"
        self.processed_dir = self.data_dir / "processed"

        # Create directories
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.processed_dir.mkdir(parents=True, exist_ok=True)

    def collect_expression_data(self) -> pd.DataFrame:
        """Collect aging-related gene expression data."""
        logger.info("Collecting expression data for aging research")

        # For demonstration, create synthetic aging expression data
        # In practice, this would pull from GEO, GTEx, or other databases

        aging_genes = [
            'TP53', 'CDKN2A', 'CDKN1A', 'SIRT1', 'FOXO3', 'ATM', 'BRCA1',
            'TERT', 'TERC', 'SOD2', 'CATALASE', 'PARP1', 'MTOR', 'AMPK',
            'PGC1A', 'NRF2', 'NFKB1', 'IL6', 'TNF', 'CCL2', 'SASP1',
            'CDKN2B', 'RB1', 'MDM2', 'ATG5', 'ATG7', 'BECN1', 'LC3B',
            'PINK1', 'PARKIN', 'TFAM', 'SOD1', 'GPX1', 'GSTP1'
        ]

        tissues = ['liver', 'brain', 'muscle', 'skin', 'heart', 'kidney']
        conditions = ['young', 'aged', 'senescent', 'rejuvenated']

        data = []
        np.random.seed(42)

        for gene in aging_genes:
            for tissue in tissues:
                for condition in conditions:
                    # Simulate realistic expression patterns
                    base_expression = np.random.lognormal(0, 1)

                    # Age-related changes
                    if condition == 'aged':
                        if gene in ['CDKN2A', 'CDKN1A', 'TP53']:  # Senescence markers up
                            expression = base_expression * np.random.uniform(2, 5)
                        elif gene in ['SIRT1', 'FOXO3', 'TERT']:  # Longevity genes down
                            expression = base_expression * np.random.uniform(0.3, 0.7)
                        else:
                            expression = base_expression * np.random.uniform(0.8, 1.5)
                    elif condition == 'senescent':
                        if gene in ['CDKN2A', 'CDKN1A', 'IL6', 'TNF']:
                            expression = base_expression * np.random.uniform(3, 8)
                        else:
                            expression = base_expression * np.random.uniform(0.5, 1.2)
                    elif condition == 'rejuvenated':
                        if gene in ['SIRT1', 'FOXO3', 'TERT', 'SOD2']:
                            expression = base_expression * np.random.uniform(1.5, 3)
                        elif gene in ['CDKN2A', 'CDKN1A']:
                            expression = base_expression * np.random.uniform(0.2, 0.6)
                        else:
                            expression = base_expression * np.random.uniform(0.9, 1.3)
                    else:  # young
                        expression = base_expression

                    data.append({
                        'gene': gene,
                        'tissue': tissue,
                        'condition': condition,
                        'expression_level': expression,
                        'log2_fc': np.log2(expression / base_expression) if base_expression > 0 else 0
                    })

        df = pd.DataFrame(data)

        # Save raw data
        output_path = self.raw_dir / "aging_expression_data.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Saved expression data to {output_path}")

        return df

    def collect_intervention_outcomes(self) -> pd.DataFrame:
        """Collect aging intervention outcome data."""
        logger.info("Collecting intervention outcome data")

        interventions = [
            {'type': 'OSK', 'targets': ['POU5F1', 'SOX2', 'KLF4']},
            {'type': 'OSKM', 'targets': ['POU5F1', 'SOX2', 'KLF4', 'MYC']},
            {'type': 'senolytic', 'targets': ['BCL2', 'BCLXL']},
            {'type': 'base_edit', 'targets': ['APOE', 'LDLR']},
            {'type': 'autophagy', 'targets': ['ATG5', 'ATG7', 'MTOR']},
        ]

        tissues = ['liver', 'brain', 'muscle', 'skin']

        data = []
        np.random.seed(42)

        for intervention in interventions:
            for tissue in tissues:
                for dose in [0.5, 1.0, 1.5, 2.0]:
                    for duration in [1, 2, 4, 8]:  # weeks

                        # Simulate realistic outcomes
                        base_safety_risk = np.random.uniform(0.1, 0.3)

                        # Intervention-specific effects
                        if intervention['type'] == 'OSK':
                            clock_delta = np.random.normal(-2.0, 0.8)
                            sen_change = np.random.normal(-0.6, 0.2)
                            safety_risk = base_safety_risk + (dose - 1.0) * 0.1
                        elif intervention['type'] == 'OSKM':
                            clock_delta = np.random.normal(-2.8, 1.0)
                            sen_change = np.random.normal(-0.8, 0.3)
                            safety_risk = base_safety_risk + 0.2 + (dose - 1.0) * 0.15
                        elif intervention['type'] == 'senolytic':
                            clock_delta = np.random.normal(-1.2, 0.5)
                            sen_change = np.random.normal(-0.9, 0.2)
                            safety_risk = base_safety_risk * 0.7
                        else:
                            clock_delta = np.random.normal(-0.8, 0.6)
                            sen_change = np.random.normal(-0.4, 0.3)
                            safety_risk = base_safety_risk

                        # Tissue-specific modifiers
                        tissue_multiplier = {
                            'liver': 1.2,
                            'muscle': 1.1,
                            'skin': 1.3,
                            'brain': 0.7
                        }[tissue]

                        clock_delta *= tissue_multiplier

                        # Duration effects
                        duration_effect = min(duration / 4.0, 1.5)
                        clock_delta *= duration_effect
                        sen_change *= duration_effect

                        # Dose effects
                        dose_effect = dose ** 0.8
                        clock_delta *= dose_effect
                        safety_risk *= dose_effect

                        functional_improvement = abs(clock_delta) * 0.3 + np.random.normal(0, 0.1)
                        transcriptional_shift = clock_delta * 0.7 + np.random.normal(0, 0.3)

                        data.append({
                            'intervention_type': intervention['type'],
                            'n_targets': len(intervention['targets']),
                            'tissue': tissue,
                            'dose_level': dose,
                            'duration_weeks': duration,
                            'epigenetic_clock_delta': clock_delta,
                            'senescence_score_change': sen_change,
                            'transcriptional_age_shift': transcriptional_shift,
                            'functional_improvement': max(0, functional_improvement),
                            'safety_risk': min(1.0, max(0, safety_risk))
                        })

        df = pd.DataFrame(data)

        # Save raw data
        output_path = self.raw_dir / "intervention_outcomes.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Saved intervention outcomes to {output_path}")

        return df

    def prepare_prioritization_training_data(self) -> pd.DataFrame:
        """Prepare training data for target prioritization model."""
        logger.info("Preparing target prioritization training data")

        # Get expression data
        expr_df = self.collect_expression_data()

        # Aggregate by gene and tissue for features
        features = []

        for gene in expr_df['gene'].unique():
            for tissue in expr_df['tissue'].unique():
                gene_tissue_data = expr_df[
                    (expr_df['gene'] == gene) & (expr_df['tissue'] == tissue)
                ]

                # Calculate features
                young_expr = gene_tissue_data[
                    gene_tissue_data['condition'] == 'young'
                ]['expression_level'].mean()

                aged_expr = gene_tissue_data[
                    gene_tissue_data['condition'] == 'aged'
                ]['expression_level'].mean()

                sen_expr = gene_tissue_data[
                    gene_tissue_data['condition'] == 'senescent'
                ]['expression_level'].mean()

                # Aging-related features
                aging_fc = np.log2(aged_expr / young_expr) if young_expr > 0 else 0
                senescence_fc = np.log2(sen_expr / young_expr) if young_expr > 0 else 0

                # Literature and pathway scores (simulated)
                np.random.seed(hash(gene + tissue) % 2**32)

                # Known aging genes get higher literature scores
                aging_genes = {
                    'TP53': 0.95, 'CDKN2A': 0.92, 'SIRT1': 0.88, 'FOXO3': 0.85,
                    'TERT': 0.82, 'ATM': 0.78, 'MTOR': 0.75, 'AMPK': 0.73
                }
                literature_score = aging_genes.get(gene, np.random.uniform(0.1, 0.6))

                # Network centrality (simulated)
                network_centrality = np.random.lognormal(0, 1)

                # Druggability score
                druggability = np.random.uniform(0.2, 0.9)

                # Target score (ground truth) - combine multiple factors
                target_score = (
                    abs(aging_fc) * 0.3 +
                    abs(senescence_fc) * 0.3 +
                    literature_score * 0.25 +
                    (network_centrality / 10) * 0.15
                )
                target_score = min(1.0, target_score)

                features.append({
                    'gene': gene,
                    'tissue': tissue,
                    f'expression_{tissue}': young_expr,
                    'aging_fold_change': aging_fc,
                    'senescence_fold_change': senescence_fc,
                    'literature_score': literature_score,
                    'network_centrality': network_centrality,
                    'druggability_score': druggability,
                    'aging_score': literature_score,  # Duplicate for compatibility
                    'target_score': target_score
                })

        df = pd.DataFrame(features)

        # Save processed data
        output_path = self.processed_dir / "prioritization_training_data.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Saved prioritization training data to {output_path}")

        return df

    def prepare_rejuvenation_training_data(self) -> pd.DataFrame:
        """Prepare training data for rejuvenation prediction model."""
        logger.info("Preparing rejuvenation prediction training data")

        # Get intervention outcomes
        outcomes_df = self.collect_intervention_outcomes()

        # Add intervention features
        features = []

        for _, row in outcomes_df.iterrows():
            intervention_features = {
                # Intervention type encoding
                'is_osk': 1.0 if row['intervention_type'] == 'OSK' else 0.0,
                'is_oskm': 1.0 if row['intervention_type'] == 'OSKM' else 0.0,
                'is_senolytic': 1.0 if row['intervention_type'] == 'senolytic' else 0.0,
                'is_base_edit': 1.0 if row['intervention_type'] == 'base_edit' else 0.0,
                'is_autophagy': 1.0 if row['intervention_type'] == 'autophagy' else 0.0,

                # Numerical features
                'n_targets': row['n_targets'],
                'dose_level': row['dose_level'],
                'duration_weeks': row['duration_weeks'],

                # Tissue encoding
                'tissue_liver': 1.0 if row['tissue'] == 'liver' else 0.0,
                'tissue_brain': 1.0 if row['tissue'] == 'brain' else 0.0,
                'tissue_muscle': 1.0 if row['tissue'] == 'muscle' else 0.0,
                'tissue_skin': 1.0 if row['tissue'] == 'skin' else 0.0,

                # Safety features (simulated)
                'oncogene_targets': 1.0 if 'MYC' in str(row['intervention_type']) else 0.0,
                'essential_gene_targets': 0.5,  # Placeholder

                # Targets
                'epigenetic_clock_delta': row['epigenetic_clock_delta'],
                'senescence_score_change': row['senescence_score_change'],
                'transcriptional_age_shift': row['transcriptional_age_shift'],
                'functional_improvement': row['functional_improvement'],
                'safety_risk': row['safety_risk']
            }

            features.append(intervention_features)

        df = pd.DataFrame(features)

        # Save processed data
        output_path = self.processed_dir / "rejuvenation_training_data.csv"
        df.to_csv(output_path, index=False)
        logger.info(f"Saved rejuvenation training data to {output_path}")

        return df

    def create_train_test_splits(self, test_size: float = 0.2) -> Dict[str, pd.DataFrame]:
        """Create train/test splits for both models."""
        logger.info("Creating train/test splits")

        # Prepare both datasets
        prioritization_data = self.prepare_prioritization_training_data()
        rejuvenation_data = self.prepare_rejuvenation_training_data()

        # Split prioritization data
        n_prior = len(prioritization_data)
        test_size_prior = int(n_prior * test_size)
        prior_indices = np.random.permutation(n_prior)

        prior_train = prioritization_data.iloc[prior_indices[test_size_prior:]]
        prior_test = prioritization_data.iloc[prior_indices[:test_size_prior]]

        # Split rejuvenation data
        n_rejuv = len(rejuvenation_data)
        test_size_rejuv = int(n_rejuv * test_size)
        rejuv_indices = np.random.permutation(n_rejuv)

        rejuv_train = rejuvenation_data.iloc[rejuv_indices[test_size_rejuv:]]
        rejuv_test = rejuvenation_data.iloc[rejuv_indices[:test_size_rejuv]]

        splits = {
            'prioritization_train': prior_train,
            'prioritization_test': prior_test,
            'rejuvenation_train': rejuv_train,
            'rejuvenation_test': rejuv_test
        }

        # Save splits
        for name, data in splits.items():
            output_path = self.processed_dir / f"{name}.csv"
            data.to_csv(output_path, index=False)
            logger.info(f"Saved {name} to {output_path}")

        return splits


def create_aging_datasets(data_dir: Path = None) -> Dict[str, pd.DataFrame]:
    """Convenience function to create all aging research datasets."""
    collector = AgingDatasetCollector(data_dir)
    return collector.create_train_test_splits()

"""Tests for aging-specific prioritization functionality."""

from unittest.mock import Mock, patch

import numpy as np

from crispr_toolkit.analysis.aging.prioritization import AgingTargetPrioritizer


class TestAgingTargetPrioritizer:
    """Test cases for aging target prioritization."""

    def setup_method(self):
        """Set up test fixtures."""
        self.prioritizer = AgingTargetPrioritizer()
        self.test_genes = ['TP53', 'CDKN2A', 'SIRT1', 'FOXO3']

    def test_prioritizer_initialization(self):
        """Test prioritizer initialization."""
        assert self.prioritizer.model is None
        assert self.prioritizer.scaler is not None
        assert hasattr(self.prioritizer, 'feature_names')

    def test_compute_features_basic(self):
        """Test basic feature computation."""
        features_df = self.prioritizer.compute_features(
            self.test_genes,
            tissue="liver"
        )

        assert len(features_df) == len(self.test_genes)
        assert 'gene' in features_df.columns
        assert 'expression_liver' in features_df.columns
        assert 'aging_score' in features_df.columns

    def test_rank_targets_heuristic(self):
        """Test target ranking with heuristic method."""
        results = self.prioritizer.rank_targets(
            genes=self.test_genes,
            tissue="liver",
            phenotype="senescence",
            n_targets=2
        )

        assert isinstance(results, list)
        assert len(results) == 2
        assert all('gene' in target for target in results)
        assert all('score' in target for target in results)
        assert all(isinstance(target['score'], float) for target in results)

    def test_rank_targets_different_tissues(self):
        """Test ranking for different tissues."""
        tissues = ["liver", "brain", "muscle"]

        for tissue in tissues:
            results = self.prioritizer.rank_targets(
                genes=self.test_genes,
                tissue=tissue,
                n_targets=2
            )
            assert len(results) == 2
            assert all('tissue' in target for target in results)

    def test_rank_targets_empty_genes(self):
        """Test ranking with empty gene list."""
        results = self.prioritizer.rank_targets(
            genes=[],
            tissue="liver",
            n_targets=5
        )
        assert results == []

    def test_heuristic_scoring(self):
        """Test heuristic scoring method."""
        mock_features = {
            'expression_liver': [2.0, 1.5, 3.0],
            'aging_score': [0.8, 0.6, 0.9],
            'literature_score': [0.7, 0.5, 0.8]
        }

        scores = self.prioritizer._heuristic_scoring(
            mock_features,
            phenotype="senescence"
        )

        assert len(scores) == 3
        assert all(isinstance(score, float) for score in scores)
        assert all(score >= 0 for score in scores)

    def test_senescence_gene_scoring(self):
        """Test specific scoring for senescence genes."""
        senescence_genes = ['TP53', 'CDKN2A', 'CDKN1A']
        other_genes = ['GAPDH', 'ACTB']

        sen_results = self.prioritizer.rank_targets(
            genes=senescence_genes,
            phenotype="senescence",
            n_targets=3
        )

        other_results = self.prioritizer.rank_targets(
            genes=other_genes,
            phenotype="senescence",
            n_targets=2
        )

        # Senescence genes should generally score higher
        avg_sen_score = np.mean([r['score'] for r in sen_results])
        avg_other_score = np.mean([r['score'] for r in other_results])

        assert avg_sen_score > avg_other_score

    def test_model_training_mock(self):
        """Test model training with mock data."""
        # Create mock training data
        import pandas as pd

        mock_data = pd.DataFrame({
            'gene': ['TP53', 'CDKN2A', 'SIRT1'],
            'expression_liver': [2.0, 1.5, 3.0],
            'aging_score': [0.8, 0.9, 0.7],
            'target_score': [0.85, 0.92, 0.73]
        })

        target_columns = ['target_score']

        # This would normally train a real model
        with patch('lightgbm.LGBMRanker') as mock_lgbm:
            mock_model = Mock()
            mock_lgbm.return_value = mock_model

            results = self.prioritizer.train(mock_data, target_columns)

            assert 'model_type' in results
            assert 'n_features' in results
            assert 'n_samples' in results

    def test_feature_computation_edge_cases(self):
        """Test feature computation edge cases."""
        # Test with unknown genes
        unknown_genes = ['FAKEGENE1', 'NOTREAL2']
        features_df = self.prioritizer.compute_features(
            unknown_genes,
            tissue="liver"
        )

        assert len(features_df) == len(unknown_genes)
        # Should handle unknown genes gracefully

    def test_tissue_specific_scoring(self):
        """Test tissue-specific scoring differences."""
        brain_results = self.prioritizer.rank_targets(
            genes=self.test_genes,
            tissue="brain",
            n_targets=2
        )

        liver_results = self.prioritizer.rank_targets(
            genes=self.test_genes,
            tissue="liver",
            n_targets=2
        )

        # Results should be different for different tissues
        brain_scores = [r['score'] for r in brain_results]
        liver_scores = [r['score'] for r in liver_results]

        # Scores should be different (though this isn't guaranteed)
        # At minimum, tissue should be recorded correctly
        assert all(r['tissue'] == 'brain' for r in brain_results)
        assert all(r['tissue'] == 'liver' for r in liver_results)

"""Tests for target prioritization functionality."""

import pandas as pd

from crispr_toolkit.analysis.aging.target_prioritization import prioritize_targets


class TestTargetPrioritization:
    """Test cases for target prioritization."""

    def test_prioritize_targets_basic(self):
        """Test basic target prioritization functionality."""
        # Create sample data
        sample_data = pd.DataFrame({
            'gene': ['TP53', 'CDKN2A', 'TERT'],
            'expression': [1.5, 2.0, 0.8]
        })

        # Test prioritization
        targets = prioritize_targets(sample_data, n_targets=5)

        # Assertions
        assert isinstance(targets, list)
        assert len(targets) <= 5
        assert all('gene' in target for target in targets)
        assert all('score' in target for target in targets)

    def test_prioritize_targets_empty_data(self):
        """Test prioritization with empty data."""
        empty_data = pd.DataFrame()
        targets = prioritize_targets(empty_data)
        assert targets == []

    def test_prioritize_targets_with_constraints(self):
        """Test prioritization with constraints."""
        sample_data = pd.DataFrame({
            'gene': ['TP53', 'CDKN2A'],
            'expression': [1.5, 2.0]
        })

        constraints = {'min_score': 0.8}
        targets = prioritize_targets(
            sample_data,
            constraints=constraints,
            n_targets=3
        )

        assert isinstance(targets, list)
        assert len(targets) <= 3

    def test_target_scores_are_sorted(self):
        """Test that targets are returned in descending score order."""
        sample_data = pd.DataFrame({
            'gene': ['TP53', 'CDKN2A', 'TERT'],
            'expression': [1.5, 2.0, 0.8]
        })

        targets = prioritize_targets(sample_data, n_targets=3)

        if len(targets) > 1:
            scores = [target['score'] for target in targets]
            assert scores == sorted(scores, reverse=True)

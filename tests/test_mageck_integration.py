"""
Tests for MAGeCK Integration

Test suite for the CRISPR screen analysis functionality including
MAGeCK wrapper, quality control, and senescence analysis.
"""

import tempfile
import unittest
from unittest.mock import MagicMock, patch

import numpy as np
import pandas as pd

from crispr_toolkit.analysis.screens import (
    MAGeCKAnalyzer,
    PathwayEnrichment,
    ScreenQualityControl,
    SenescenceScreenAnalyzer,
)


class TestMAGeCKAnalyzer(unittest.TestCase):
    """Test MAGeCK analyzer functionality."""

    def setUp(self):
        """Set up test fixtures."""
        self.temp_dir = tempfile.mkdtemp()

    def test_mageck_initialization(self):
        """Test MAGeCK analyzer initialization."""
        with patch('subprocess.run') as mock_run:
            mock_run.return_value.stdout = "MAGeCK v0.5.9"
            analyzer = MAGeCKAnalyzer(temp_dir=self.temp_dir)
            self.assertIsInstance(analyzer, MAGeCKAnalyzer)

    def test_senescence_pathway_analysis(self):
        """Test senescence pathway analysis."""
        # Create mock MAGeCK results
        gene_summary = pd.DataFrame({
            'id': ['HMGB1', 'AGER', 'JAK1', 'CONTROL1', 'CONTROL2'],
            'pos|lfc': [-1.5, -1.2, -0.8, 0.1, -0.1],
            'pos|fdr': [0.001, 0.005, 0.05, 0.5, 0.8]
        })

        mock_results = MagicMock()
        mock_results.gene_summary = gene_summary

        with patch('subprocess.run'):
            analyzer = MAGeCKAnalyzer(temp_dir=self.temp_dir)

            # Test senescence analysis
            senescence_results = analyzer.senescence_pathway_analysis(mock_results)

            self.assertIn('senescence_hits', senescence_results)
            self.assertIn('enrichment_stats', senescence_results)
            self.assertGreater(
                senescence_results['enrichment_stats']['significant_hits'], 0
            )


class TestScreenQualityControl(unittest.TestCase):
    """Test screen quality control functionality."""

    def setUp(self):
        """Set up test fixtures."""
        np.random.seed(42)
        # Create mock count matrix
        self.count_matrix = pd.DataFrame(
            np.random.poisson(100, (100, 5)),
            columns=['Sample1', 'Sample2', 'Sample3', 'Sample4', 'Sample5'],
            index=[f'Gene_{i}' for i in range(100)]
        )

    def test_qc_analysis(self):
        """Test quality control analysis."""
        qc = ScreenQualityControl()
        results = qc.analyze_count_matrix(self.count_matrix)

        self.assertIn('read_distribution', results)
        self.assertIn('library_representation', results)
        self.assertIn('sample_correlation', results)
        self.assertIn('recommendations', results)

    def test_read_distribution_analysis(self):
        """Test read distribution analysis."""
        qc = ScreenQualityControl()
        results = qc._analyze_read_distribution(self.count_matrix)

        self.assertIn('sample_read_counts', results)
        self.assertIn('gene_read_counts', results)
        self.assertIn('mean', results['sample_read_counts'])
        self.assertIn('cv', results['sample_read_counts'])

    def test_library_representation(self):
        """Test library representation analysis."""
        qc = ScreenQualityControl()
        results = qc._analyze_library_representation(self.count_matrix)

        self.assertIn('total_genes', results)
        self.assertIn('total_samples', results)
        self.assertIn('mean_detection_rate', results)
        self.assertEqual(results['total_genes'], 100)
        self.assertEqual(results['total_samples'], 5)


class TestPathwayEnrichment(unittest.TestCase):
    """Test pathway enrichment analysis."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock screen results
        self.screen_results = pd.DataFrame({
            'gene': ['HMGB1', 'AGER', 'JAK1', 'STAT1', 'CONTROL1', 'CONTROL2'],
            'log2fc': [-1.5, -1.2, -0.8, -1.0, 0.1, -0.1],
            'fdr': [0.001, 0.005, 0.05, 0.02, 0.5, 0.8]
        })

    def test_enrichment_analysis(self):
        """Test pathway enrichment analysis."""
        pathway_analyzer = PathwayEnrichment()
        results = pathway_analyzer.run_enrichment_analysis(
            self.screen_results,
            fdr_threshold=0.1
        )

        self.assertIsInstance(results, dict)
        self.assertIn('kegg_aging', results)

    def test_enrichment_calculation(self):
        """Test enrichment calculation for specific pathway."""
        pathway_analyzer = PathwayEnrichment()

        significant_genes = ['HMGB1', 'AGER', 'JAK1', 'STAT1']
        pathway_genes = ['HMGB1', 'AGER', 'OTHER1', 'OTHER2']

        enrichment = pathway_analyzer._calculate_enrichment(
            significant_genes=significant_genes,
            pathway_genes=pathway_genes,
            total_genes=100,
            screen_results=self.screen_results
        )

        self.assertIn('num_hits', enrichment)
        self.assertIn('enrichment_ratio', enrichment)
        self.assertEqual(enrichment['num_hits'], 2)  # HMGB1 and AGER

    def test_create_enrichment_summary(self):
        """Test enrichment summary creation."""
        pathway_analyzer = PathwayEnrichment()

        # Mock enrichment results
        mock_results = {
            'database1': {
                'pathway1': {
                    'pathway_size': 10,
                    'num_hits': 3,
                    'enrichment_ratio': 2.0,
                    'p_value': 0.01,
                    'mean_log2fc': -1.0,
                    'hit_genes': ['GENE1', 'GENE2', 'GENE3']
                }
            }
        }

        summary_df = pathway_analyzer.create_enrichment_summary(mock_results)

        self.assertIsInstance(summary_df, pd.DataFrame)
        if len(summary_df) > 0:
            self.assertIn('database', summary_df.columns)
            self.assertIn('pathway', summary_df.columns)
            self.assertIn('p_value', summary_df.columns)


class TestSenescenceScreenAnalyzer(unittest.TestCase):
    """Test senescence-specific screen analysis."""

    def setUp(self):
        """Set up test fixtures."""
        # Create mock screen results with senescence genes
        self.screen_results = pd.DataFrame({
            'gene': [
                'HMGB1', 'AGER', 'JAK1', 'STAT1', 'NFKB1', 'CDKN1A',
                'CONTROL1', 'CONTROL2', 'CONTROL3'
            ],
            'log2fc': [-1.5, -1.2, -0.8, -1.0, -0.9, -1.1, 0.1, -0.1, 0.05],
            'fdr': [0.001, 0.005, 0.05, 0.02, 0.03, 0.01, 0.5, 0.8, 0.9],
            'pvalue': [0.0001, 0.0005, 0.005, 0.002, 0.003, 0.001, 0.05, 0.08, 0.09]
        })

    def test_senescence_analyzer_initialization(self):
        """Test senescence analyzer initialization."""
        analyzer = SenescenceScreenAnalyzer()
        self.assertIsInstance(analyzer.senescence_pathways, dict)
        self.assertIn('rehmgb1_rage_signaling', analyzer.senescence_pathways)
        self.assertIn('jak_stat_pathway', analyzer.senescence_pathways)

    def test_senescence_screen_analysis(self):
        """Test comprehensive senescence screen analysis."""
        analyzer = SenescenceScreenAnalyzer()
        results = analyzer.analyze_senescence_screen(
            self.screen_results,
            fdr_threshold=0.1,
            lfc_threshold=0.5
        )

        self.assertIn('pathway_enrichment', results)
        self.assertIn('rehmgb1_pathway_hits', results)
        self.assertIn('senescence_modulators', results)
        self.assertIn('therapeutic_targets', results)

    def test_rehmgb1_pathway_analysis(self):
        """Test ReHMGB1 pathway specific analysis."""
        analyzer = SenescenceScreenAnalyzer()
        results = analyzer._analyze_rehmgb1_pathway(
            self.screen_results,
            fdr_threshold=0.1
        )

        self.assertIn('total_rehmgb1_genes', results)
        self.assertIn('significant_hits', results)
        self.assertIn('hit_details', results)
        self.assertGreater(results['significant_hits'], 0)

    def test_senescence_modulators_identification(self):
        """Test identification of senescence modulators."""
        analyzer = SenescenceScreenAnalyzer()
        results = analyzer._identify_senescence_modulators(
            self.screen_results,
            fdr_threshold=0.1,
            lfc_threshold=0.5
        )

        self.assertIn('pro_senescence_hits', results)
        self.assertIn('anti_senescence_hits', results)

        # All our mock hits should be anti-senescence (negative LFC)
        self.assertEqual(results['pro_senescence_hits']['total'], 0)
        self.assertGreater(results['anti_senescence_hits']['total'], 0)

    def test_therapeutic_target_identification(self):
        """Test therapeutic target identification."""
        analyzer = SenescenceScreenAnalyzer()
        results = analyzer._identify_therapeutic_targets(
            self.screen_results,
            fdr_threshold=0.1,
            lfc_threshold=0.5
        )

        self.assertIn('senescence_suppressors', results)
        self.assertIn('druggable_targets', results)
        self.assertIn('pathway_targets', results)


class TestIntegration(unittest.TestCase):
    """Integration tests for the complete screen analysis workflow."""

    def test_complete_workflow(self):
        """Test complete analysis workflow."""
        # Create comprehensive mock data
        np.random.seed(42)

        # Screen results
        genes = ['HMGB1', 'AGER', 'JAK1', 'STAT1'] + [f'CTRL_{i}' for i in range(96)]
        screen_data = pd.DataFrame({
            'gene': genes,
            'log2fc': np.random.normal(-0.5, 1.0, 100),
            'fdr': np.random.beta(2, 10, 100),
            'pvalue': np.random.beta(1, 20, 100)
        })

        # Count matrix
        count_matrix = pd.DataFrame(
            np.random.poisson(50, (100, 4)),
            columns=['Ctrl1', 'Ctrl2', 'Treat1', 'Treat2'],
            index=genes
        )

        # Run QC analysis
        qc = ScreenQualityControl()
        qc_results = qc.analyze_count_matrix(count_matrix)
        self.assertIn('recommendations', qc_results)

        # Run pathway enrichment
        pathway_analyzer = PathwayEnrichment()
        enrichment_results = pathway_analyzer.run_enrichment_analysis(
            screen_data, fdr_threshold=0.5
        )
        self.assertIsInstance(enrichment_results, dict)

        # Run senescence analysis
        senescence_analyzer = SenescenceScreenAnalyzer()
        senescence_results = senescence_analyzer.analyze_senescence_screen(
            screen_data, fdr_threshold=0.5
        )
        self.assertIn('rehmgb1_pathway_hits', senescence_results)

        # All analyses should complete without error
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
    unittest.main()

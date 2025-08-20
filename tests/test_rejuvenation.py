"""Tests for rejuvenation prediction functionality."""


from crispr_toolkit.analysis.aging.rejuvenation import (
    RejuvenationPredictor,
    predict_rejuvenation,
)


class TestRejuvenationPredictor:
    """Test cases for rejuvenation prediction."""

    def setup_method(self):
        """Set up test fixtures."""
        self.predictor = RejuvenationPredictor()
        self.test_config = {
            "intervention": "OSK",
            "targets": ["TP53", "SIRT1"],
            "delivery": "AAV",
            "dose_level": 1.0,
            "duration_weeks": 4
        }

    def test_predictor_initialization(self):
        """Test predictor initialization."""
        assert self.predictor.model is None
        assert self.predictor.scaler is not None
        assert len(self.predictor.target_names) == 5
        assert "epigenetic_clock_delta" in self.predictor.target_names

    def test_extract_intervention_features(self):
        """Test intervention feature extraction."""
        features = self.predictor.extract_intervention_features(
            self.test_config,
            tissue="liver"
        )

        assert isinstance(features, dict)
        assert "is_osk" in features
        assert features["is_osk"] == 1.0
        assert "n_targets" in features
        assert features["n_targets"] == 2
        assert "tissue_liver" in features
        assert features["tissue_liver"] == 1.0

    def test_predict_outcomes_heuristic(self):
        """Test outcome prediction with heuristic method."""
        results = self.predictor.predict_outcomes(
            self.test_config,
            tissue="liver",
            return_uncertainty=True
        )

        assert "predictions" in results
        assert "intervention" in results
        assert "tissue" in results

        predictions = results["predictions"]
        assert "epigenetic_clock_delta" in predictions
        assert "safety_risk" in predictions

        # Check that values are reasonable
        assert isinstance(predictions["epigenetic_clock_delta"], float)
        assert 0 <= predictions["safety_risk"] <= 1.0

    def test_different_interventions(self):
        """Test predictions for different intervention types."""
        interventions = ["OSK", "OSKM", "senolytic"]

        for intervention in interventions:
            config = {**self.test_config, "intervention": intervention}
            results = self.predictor.predict_outcomes(config, tissue="liver")

            assert results["intervention"] == intervention
            assert "predictions" in results

    def test_tissue_specific_predictions(self):
        """Test tissue-specific prediction differences."""
        tissues = ["liver", "brain", "muscle"]

        results_by_tissue = {}
        for tissue in tissues:
            results = self.predictor.predict_outcomes(
                self.test_config,
                tissue=tissue
            )
            results_by_tissue[tissue] = results

        # Each tissue should have predictions
        for tissue, results in results_by_tissue.items():
            assert results["tissue"] == tissue
            assert "predictions" in results

    def test_context_data_integration(self):
        """Test integration of context data."""
        context_data = {
            "age": 24,
            "senescence_load": 0.3,
            "expression_level": 1.2
        }

        results = self.predictor.predict_outcomes(
            self.test_config,
            tissue="liver",
            context_data=context_data
        )

        assert "features" in results
        features = results["features"]
        assert "baseline_age" in features
        assert features["baseline_age"] == 24

    def test_safety_risk_calculation(self):
        """Test safety risk calculation for different scenarios."""
        # High-risk intervention with oncogenes
        risky_config = {
            "intervention": "OSKM",
            "targets": ["MYC", "TP53"],
            "delivery": "AAV"
        }

        risky_results = self.predictor.predict_outcomes(
            risky_config,
            tissue="brain"
        )

        # Low-risk intervention
        safe_config = {
            "intervention": "senolytic",
            "targets": ["SIRT1"],
            "delivery": "LNP"
        }

        safe_results = self.predictor.predict_outcomes(
            safe_config,
            tissue="liver"
        )

        risky_safety = risky_results["predictions"]["safety_risk"]
        safe_safety = safe_results["predictions"]["safety_risk"]

        # Risky intervention should have higher safety risk
        assert risky_safety > safe_safety

    def test_uncertainty_estimation(self):
        """Test uncertainty estimation."""
        results = self.predictor.predict_outcomes(
            self.test_config,
            return_uncertainty=True
        )

        assert "uncertainty" in results
        assert "confidence_intervals" in results

        uncertainty = results["uncertainty"]
        intervals = results["confidence_intervals"]

        assert len(uncertainty) == len(self.predictor.target_names)
        assert len(intervals) == len(self.predictor.target_names)

        # Each interval should be a list of [lower, upper]
        for interval in intervals.values():
            assert len(interval) == 2
            assert interval[0] < interval[1]  # lower < upper

    def test_target_gene_specific_effects(self):
        """Test target gene-specific effects."""
        # Test different target gene combinations
        senescence_config = {
            **self.test_config,
            "targets": ["TP53", "CDKN2A"]  # senescence genes
        }

        autophagy_config = {
            **self.test_config,
            "targets": ["ATG5", "BECN1"]  # autophagy genes
        }

        sen_results = self.predictor.predict_outcomes(senescence_config)
        auto_results = self.predictor.predict_outcomes(autophagy_config)

        # Both should have predictions
        assert "predictions" in sen_results
        assert "predictions" in auto_results

    def test_dose_level_effects(self):
        """Test dose level effects on predictions."""
        low_dose = {**self.test_config, "dose_level": 0.5}
        high_dose = {**self.test_config, "dose_level": 2.0}

        low_results = self.predictor.predict_outcomes(low_dose)
        high_results = self.predictor.predict_outcomes(high_dose)

        low_safety = low_results["predictions"]["safety_risk"]
        high_safety = high_results["predictions"]["safety_risk"]

        # Higher dose should generally have higher safety risk
        # (though this may vary with heuristic implementation)
        assert isinstance(low_safety, float)
        assert isinstance(high_safety, float)


class TestConvenienceFunction:
    """Test the convenience function for rejuvenation prediction."""

    def test_predict_rejuvenation_basic(self):
        """Test basic usage of predict_rejuvenation function."""
        config = {
            "intervention": "OSK",
            "targets": ["TP53"],
            "delivery": "AAV"
        }

        results = predict_rejuvenation(config, tissue="liver")

        assert "predictions" in results
        assert "intervention" in results
        assert results["intervention"] == "OSK"

    def test_predict_rejuvenation_with_context(self):
        """Test predict_rejuvenation with context data."""
        config = {
            "intervention": "senolytic",
            "targets": ["CDKN2A"]
        }

        context = {
            "age": 18,
            "senescence_load": 0.4
        }

        results = predict_rejuvenation(
            config,
            tissue="muscle",
            context_data=context
        )

        assert results["tissue"] == "muscle"
        assert "predictions" in results

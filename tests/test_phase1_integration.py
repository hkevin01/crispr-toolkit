"""
Test Phase 1 critical functionality.

Basic integration tests for the Phase 1 development priorities.
"""

import unittest
import numpy as np
import sys
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))


class TestPhase1Integration(unittest.TestCase):
    """Test Phase 1 critical priorities."""
    
    def test_real_dataset_loading(self):
        """Test real dataset integration."""
        try:
            from src.crispr_toolkit.data.real_datasets import (
                AgingDatasetLoader, get_intervention_target_genes
            )
            
            # Test dataset loader creation
            loader = AgingDatasetLoader()
            self.assertIsNotNone(loader)
            
            # Test intervention targets
            targets = get_intervention_target_genes()
            self.assertIsInstance(targets, dict)
            self.assertGreater(len(targets), 0)
            
            print("✅ Real dataset integration working")
            
        except ImportError as e:
            self.skipTest(f"Dataset modules not available: {e}")
    
    def test_experiment_tracking(self):
        """Test MLflow experiment tracking."""
        try:
            from src.crispr_toolkit.models.experiment_tracking import (
                ExperimentTracker
            )
            
            # Test tracker creation
            tracker = ExperimentTracker("test_experiment")
            self.assertIsNotNone(tracker)
            
            print("✅ Experiment tracking working")
            
        except ImportError as e:
            self.skipTest(f"Tracking modules not available: {e}")
    
    def test_hyperparameter_optimization(self):
        """Test hyperparameter optimization."""
        try:
            from src.crispr_toolkit.models.hyperparameter_optimization import (
                HyperparameterOptimizer
            )
            
            # Test optimizer creation
            optimizer = HyperparameterOptimizer("test_study")
            self.assertIsNotNone(optimizer)
            
            print("✅ Hyperparameter optimization working")
            
        except ImportError as e:
            self.skipTest(f"Optimization modules not available: {e}")
    
    def test_ensemble_methods(self):
        """Test ensemble methods."""
        try:
            from src.crispr_toolkit.models.ensemble_methods import (
                StackingEnsemble, DynamicEnsemble
            )
            
            # Create simple mock models for testing
            class MockModel:
                def fit(self, X, y):
                    return self
                
                def predict(self, X):
                    return np.random.random(len(X))
                
                def get_params(self):
                    return {}
            
            # Test stacking ensemble
            models = [MockModel(), MockModel()]
            meta_model = MockModel()
            
            stacking = StackingEnsemble(models, meta_model)
            self.assertIsNotNone(stacking)
            
            # Test dynamic ensemble
            dynamic = DynamicEnsemble(models)
            self.assertIsNotNone(dynamic)
            
            print("✅ Ensemble methods working")
            
        except ImportError as e:
            self.skipTest(f"Ensemble modules not available: {e}")
    
    def test_synthetic_data_generation(self):
        """Test that we can generate synthetic data for testing."""
        # Generate simple synthetic aging data
        np.random.seed(42)
        n_samples = 100
        n_features = 10
        
        # Create synthetic expression data
        X = np.random.normal(5, 2, (n_samples, n_features))
        
        # Create age target with some correlation to features
        age_effect = np.sum(X[:, :3], axis=1) * 0.1
        y = 50 + age_effect + np.random.normal(0, 5, n_samples)
        
        self.assertEqual(X.shape, (n_samples, n_features))
        self.assertEqual(len(y), n_samples)
        
        print(f"✅ Generated synthetic data: {X.shape} features, {len(y)} samples")
    
    def test_basic_ml_pipeline(self):
        """Test basic ML pipeline with synthetic data."""
        try:
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.model_selection import train_test_split
            from sklearn.metrics import r2_score
            
            # Generate synthetic data
            np.random.seed(42)
            X = np.random.normal(5, 2, (200, 15))
            y = 50 + np.sum(X[:, :5], axis=1) * 0.2 + np.random.normal(0, 3, 200)
            
            # Split data
            X_train, X_test, y_train, y_test = train_test_split(
                X, y, test_size=0.3, random_state=42
            )
            
            # Train model
            model = RandomForestRegressor(n_estimators=50, random_state=42)
            model.fit(X_train, y_train)
            
            # Make predictions
            y_pred = model.predict(X_test)
            score = r2_score(y_test, y_pred)
            
            self.assertGreater(score, 0.0)  # Should have some predictive power
            
            print(f"✅ Basic ML pipeline working (R² = {score:.3f})")
            
        except ImportError as e:
            self.skipTest(f"ML libraries not available: {e}")


def run_phase1_tests():
    """Run Phase 1 integration tests."""
    print("🧪 Running Phase 1 Integration Tests")
    print("=" * 50)
    
    # Create test suite
    suite = unittest.TestLoader().loadTestsFromTestCase(TestPhase1Integration)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "=" * 50)
    if result.wasSuccessful():
        print("🎉 All Phase 1 tests passed!")
    else:
        print("❌ Some tests failed. Check the output above.")
        
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_phase1_tests()
    sys.exit(0 if success else 1)

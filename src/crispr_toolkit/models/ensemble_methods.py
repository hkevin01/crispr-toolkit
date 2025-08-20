"""
Ensemble methods for aging intervention prediction models.

This module implements advanced ensemble techniques to improve
prediction accuracy for aging and longevity interventions.
"""

import logging
import numpy as np
from typing import Dict, List, Any
from sklearn.ensemble import (
    VotingRegressor, RandomForestRegressor, GradientBoostingRegressor
)
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
from sklearn.base import BaseEstimator, RegressorMixin, clone

logger = logging.getLogger(__name__)


class StackingEnsemble(BaseEstimator, RegressorMixin):
    """Advanced stacking ensemble for aging prediction models."""
    
    def __init__(self, base_models: List[Any], meta_model: Any,
                 cv_folds: int = 5, random_state: int = 42):
        """Initialize stacking ensemble."""
        self.base_models = base_models
        self.meta_model = meta_model
        self.cv_folds = cv_folds
        self.random_state = random_state
        self.fitted_base_models = []
        self.fitted_meta_model = None
        
    def fit(self, X: np.ndarray, y: np.ndarray) -> 'StackingEnsemble':
        """Fit the stacking ensemble."""
        X = np.array(X)
        y = np.array(y)
        
        # Create meta-features using cross-validation
        meta_features = np.zeros((X.shape[0], len(self.base_models)))
        
        from sklearn.model_selection import KFold
        kf = KFold(n_splits=self.cv_folds, shuffle=True, 
                  random_state=self.random_state)
        
        for i, model in enumerate(self.base_models):
            model_predictions = np.zeros(X.shape[0])
            
            for train_idx, val_idx in kf.split(X):
                X_train, X_val = X[train_idx], X[val_idx]
                y_train = y[train_idx]
                
                # Clone and fit model
                model_clone = clone(model)
                model_clone.fit(X_train, y_train)
                
                # Predict validation set
                val_pred = model_clone.predict(X_val)
                model_predictions[val_idx] = val_pred
                
            meta_features[:, i] = model_predictions
            
        # Fit base models on full data
        self.fitted_base_models = []
        for model in self.base_models:
            fitted_model = clone(model)
            fitted_model.fit(X, y)
            self.fitted_base_models.append(fitted_model)
            
        # Fit meta-model
        self.fitted_meta_model = clone(self.meta_model)
        self.fitted_meta_model.fit(meta_features, y)
        
        return self
        
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions using the stacking ensemble."""
        X = np.array(X)
        
        # Get base model predictions
        base_predictions = np.zeros((X.shape[0], len(self.fitted_base_models)))
        
        for i, model in enumerate(self.fitted_base_models):
            base_predictions[:, i] = model.predict(X)
            
        # Meta-model prediction
        final_predictions = self.fitted_meta_model.predict(base_predictions)
        
        return final_predictions


class DynamicEnsemble(BaseEstimator, RegressorMixin):
    """Dynamic ensemble that adapts weights based on prediction confidence."""
    
    def __init__(self, models: List[Any], confidence_threshold: float = 0.8):
        """Initialize dynamic ensemble."""
        self.models = models
        self.confidence_threshold = confidence_threshold
        self.fitted_models = []
        self.model_weights = None
        
    def _calculate_prediction_confidence(self, predictions: np.ndarray
                                       ) -> np.ndarray:
        """Calculate prediction confidence based on model agreement."""
        # Use standard deviation as inverse confidence measure
        std_dev = np.std(predictions, axis=1)
        max_std = np.max(std_dev)
        
        if max_std > 0:
            confidence = 1 - (std_dev / max_std)
        else:
            confidence = np.ones(len(predictions))
            
        return confidence
        
    def fit(self, X: np.ndarray, y: np.ndarray) -> 'DynamicEnsemble':
        """Fit the dynamic ensemble."""
        X = np.array(X)
        y = np.array(y)
        
        # Fit all models
        self.fitted_models = []
        for model in self.models:
            fitted_model = clone(model)
            fitted_model.fit(X, y)
            self.fitted_models.append(fitted_model)
            
        # Calculate model weights based on performance
        self.model_weights = self._calculate_model_weights(X, y)
        
        return self
        
    def _calculate_model_weights(self, X: np.ndarray, y: np.ndarray
                               ) -> np.ndarray:
        """Calculate model weights based on cross-validation performance."""
        from sklearn.model_selection import cross_val_score
        
        weights = []
        for model in self.fitted_models:
            try:
                cv_scores = cross_val_score(
                    model, X, y, cv=3, scoring='r2', n_jobs=-1
                )
                weight = max(0, cv_scores.mean())  # Ensure non-negative
                weights.append(weight)
            except Exception as e:
                logger.warning(f"Failed to calculate weight for model: {e}")
                weights.append(0.1)  # Default small weight
                
        weights = np.array(weights)
        # Normalize weights
        if weights.sum() > 0:
            weights = weights / weights.sum()
        else:
            weights = np.ones(len(weights)) / len(weights)
            
        return weights
        
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions using dynamic weighting."""
        X = np.array(X)
        
        # Get predictions from all models
        all_predictions = np.zeros((X.shape[0], len(self.fitted_models)))
        
        for i, model in enumerate(self.fitted_models):
            all_predictions[:, i] = model.predict(X)
            
        # Calculate confidence for each prediction
        confidence = self._calculate_prediction_confidence(all_predictions)
        
        # Apply dynamic weighting
        final_predictions = np.zeros(X.shape[0])
        
        for i in range(X.shape[0]):
            if confidence[i] >= self.confidence_threshold:
                # High confidence: use weighted average
                pred = np.average(all_predictions[i], weights=self.model_weights)
            else:
                # Low confidence: use simple average
                pred = np.mean(all_predictions[i])
                
            final_predictions[i] = pred
            
        return final_predictions


class AdaptiveEnsemble(BaseEstimator, RegressorMixin):
    """Adaptive ensemble that learns optimal model combinations."""
    
    def __init__(self, models: List[Any], adaptation_rate: float = 0.1):
        """Initialize adaptive ensemble."""
        self.models = models
        self.adaptation_rate = adaptation_rate
        self.fitted_models = []
        self.combination_weights = None
        
    def fit(self, X: np.ndarray, y: np.ndarray) -> 'AdaptiveEnsemble':
        """Fit the adaptive ensemble."""
        X = np.array(X)
        y = np.array(y)
        
        # Fit all models
        self.fitted_models = []
        for model in self.models:
            fitted_model = clone(model)
            fitted_model.fit(X, y)
            self.fitted_models.append(fitted_model)
            
        # Learn optimal combination weights
        self.combination_weights = self._learn_combination_weights(X, y)
        
        return self
        
    def _learn_combination_weights(self, X: np.ndarray, y: np.ndarray
                                 ) -> np.ndarray:
        """Learn optimal combination weights using gradient descent."""
        n_models = len(self.fitted_models)
        weights = np.ones(n_models) / n_models  # Initialize equal weights
        
        # Get base predictions
        base_predictions = np.zeros((X.shape[0], n_models))
        for i, model in enumerate(self.fitted_models):
            base_predictions[:, i] = model.predict(X)
            
        # Simple optimization loop
        for iteration in range(100):
            # Calculate current ensemble prediction
            ensemble_pred = np.dot(base_predictions, weights)
            
            # Calculate gradient
            error = ensemble_pred - y
            gradient = np.dot(base_predictions.T, error) / len(y)
            
            # Update weights
            weights = weights - self.adaptation_rate * gradient
            
            # Apply constraints (non-negative, sum to 1)
            weights = np.maximum(weights, 0)
            if weights.sum() > 0:
                weights = weights / weights.sum()
            else:
                weights = np.ones(n_models) / n_models
                
        return weights
        
    def predict(self, X: np.ndarray) -> np.ndarray:
        """Make predictions using learned combination weights."""
        X = np.array(X)
        
        # Get predictions from all models
        predictions = np.zeros((X.shape[0], len(self.fitted_models)))
        
        for i, model in enumerate(self.fitted_models):
            predictions[:, i] = model.predict(X)
            
        # Combine using learned weights
        final_predictions = np.dot(predictions, self.combination_weights)
        
        return final_predictions


class EnsembleOptimizer:
    """Optimize ensemble configurations for aging prediction."""
    
    def __init__(self, random_state: int = 42):
        """Initialize ensemble optimizer."""
        self.random_state = random_state
        self.best_ensemble = None
        self.best_score = -np.inf
        
    def create_base_models(self) -> List[Any]:
        """Create diverse base models for ensembling."""
        from sklearn.linear_model import Ridge, Lasso
        from sklearn.svm import SVR
        from sklearn.neighbors import KNeighborsRegressor
        
        models = [
            RandomForestRegressor(
                n_estimators=100, random_state=self.random_state
            ),
            GradientBoostingRegressor(
                n_estimators=100, random_state=self.random_state
            ),
            Ridge(alpha=1.0),
            Lasso(alpha=0.1),
            SVR(kernel='rbf', C=1.0),
            KNeighborsRegressor(n_neighbors=5)
        ]
        
        return models
        
    def optimize_voting_ensemble(self, X: np.ndarray, y: np.ndarray
                               ) -> Dict[str, Any]:
        """Optimize voting ensemble configuration."""
        base_models = self.create_base_models()
        
        best_score = -np.inf
        best_models = None
        
        # Try different combinations of models
        from itertools import combinations
        
        for r in range(3, len(base_models) + 1):
            for model_combo in combinations(base_models, r):
                try:
                    # Create voting ensemble
                    ensemble = VotingRegressor([
                        (f'model_{i}', model) 
                        for i, model in enumerate(model_combo)
                    ])
                    
                    # Evaluate with cross-validation
                    cv_scores = cross_val_score(
                        ensemble, X, y, cv=5, scoring='r2', n_jobs=-1
                    )
                    mean_score = cv_scores.mean()
                    
                    if mean_score > best_score:
                        best_score = mean_score
                        best_models = model_combo
                        
                except Exception as e:
                    logger.warning(f"Failed to evaluate model combo: {e}")
                    continue
                    
        if best_models:
            self.best_ensemble = VotingRegressor([
                (f'model_{i}', model) 
                for i, model in enumerate(best_models)
            ])
            self.best_score = best_score
            
        return {
            'best_score': best_score,
            'best_models': [str(model) for model in best_models] if best_models else [],
            'n_models': len(best_models) if best_models else 0
        }
        
    def optimize_stacking_ensemble(self, X: np.ndarray, y: np.ndarray
                                 ) -> Dict[str, Any]:
        """Optimize stacking ensemble configuration."""
        base_models = self.create_base_models()
        
        # Try different meta-models
        from sklearn.linear_model import Ridge, LinearRegression
        meta_models = [
            LinearRegression(),
            Ridge(alpha=1.0),
            RandomForestRegressor(n_estimators=50, random_state=self.random_state)
        ]
        
        best_score = -np.inf
        best_config = None
        
        for meta_model in meta_models:
            try:
                ensemble = StackingEnsemble(
                    base_models=base_models[:4],  # Use subset for efficiency
                    meta_model=meta_model,
                    random_state=self.random_state
                )
                
                # Evaluate with cross-validation
                cv_scores = cross_val_score(
                    ensemble, X, y, cv=3, scoring='r2', n_jobs=-1
                )
                mean_score = cv_scores.mean()
                
                if mean_score > best_score:
                    best_score = mean_score
                    best_config = {
                        'base_models': base_models[:4],
                        'meta_model': meta_model
                    }
                    
            except Exception as e:
                logger.warning(f"Failed to evaluate stacking config: {e}")
                continue
                
        if best_config:
            self.best_ensemble = StackingEnsemble(
                base_models=best_config['base_models'],
                meta_model=best_config['meta_model'],
                random_state=self.random_state
            )
            self.best_score = best_score
            
        return {
            'best_score': best_score,
            'meta_model': str(best_config['meta_model']) if best_config else None
        }


def create_aging_ensemble(X: np.ndarray, y: np.ndarray,
                         ensemble_type: str = 'voting') -> Any:
    """Create optimized ensemble for aging prediction."""
    
    optimizer = EnsembleOptimizer()
    
    if ensemble_type == 'voting':
        config = optimizer.optimize_voting_ensemble(X, y)
        logger.info(f"Best voting ensemble score: {config['best_score']:.4f}")
        return optimizer.best_ensemble
        
    elif ensemble_type == 'stacking':
        config = optimizer.optimize_stacking_ensemble(X, y)
        logger.info(f"Best stacking ensemble score: {config['best_score']:.4f}")
        return optimizer.best_ensemble
        
    elif ensemble_type == 'dynamic':
        base_models = optimizer.create_base_models()
        ensemble = DynamicEnsemble(models=base_models[:5])
        ensemble.fit(X, y)
        logger.info("Created dynamic ensemble")
        return ensemble
        
    elif ensemble_type == 'adaptive':
        base_models = optimizer.create_base_models()
        ensemble = AdaptiveEnsemble(models=base_models[:4])
        ensemble.fit(X, y)
        logger.info("Created adaptive ensemble")
        return ensemble
        
    else:
        raise ValueError(f"Unknown ensemble type: {ensemble_type}")


def evaluate_ensemble_performance(ensemble: Any, X_test: np.ndarray,
                                y_test: np.ndarray) -> Dict[str, float]:
    """Evaluate ensemble performance on test data."""
    
    predictions = ensemble.predict(X_test)
    
    metrics = {
        'r2_score': r2_score(y_test, predictions),
        'rmse': np.sqrt(mean_squared_error(y_test, predictions)),
        'mae': mean_absolute_error(y_test, predictions)
    }
    
    return metrics

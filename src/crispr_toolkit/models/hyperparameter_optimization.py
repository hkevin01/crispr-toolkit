"""
Hyperparameter optimization for aging intervention models using Optuna.
"""

import logging
import optuna
import numpy as np
import pandas as pd
from typing import Dict, Any, Optional, List, Union
from sklearn.model_selection import cross_val_score, StratifiedKFold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import lightgbm as lgb
import xgboost as xgb

logger = logging.getLogger(__name__)


class HyperparameterOptimizer:
    """Automated hyperparameter optimization for aging intervention models."""
    
    def __init__(self, study_name: str = "aging_intervention_optimization",
                 storage_url: Optional[str] = None):
        """Initialize hyperparameter optimizer."""
        self.study_name = study_name
        self.storage_url = storage_url or "sqlite:///optuna_studies.db"
        self.study: Optional[optuna.Study] = None
        self.best_params: Optional[Dict[str, Any]] = None
        self.best_score: Optional[float] = None
        
    def create_study(self, direction: str = "maximize",
                     sampler_name: str = "TPE") -> optuna.Study:
        """Create or load Optuna study."""
        # Choose sampler
        if sampler_name == "TPE":
            sampler = optuna.samplers.TPESampler(seed=42)
        elif sampler_name == "Random":
            sampler = optuna.samplers.RandomSampler(seed=42)
        elif sampler_name == "CmaEs":
            sampler = optuna.samplers.CmaEsSampler(seed=42)
        else:
            sampler = optuna.samplers.TPESampler(seed=42)
            
        try:
            # Load existing study or create new one
            self.study = optuna.load_study(
                study_name=self.study_name,
                storage=self.storage_url
            )
            logger.info(f"Loaded existing study: {self.study_name}")
        except KeyError:
            # Create new study
            self.study = optuna.create_study(
                study_name=self.study_name,
                storage=self.storage_url,
                direction=direction,
                sampler=sampler
            )
            logger.info(f"Created new study: {self.study_name}")
            
        return self.study
        
    def optimize_random_forest(self, X: np.ndarray, y: np.ndarray,
                               n_trials: int = 100,
                               cv_folds: int = 5) -> Dict[str, Any]:
        """Optimize Random Forest hyperparameters."""
        
        def objective(trial):
            # Suggest hyperparameters
            params = {
                'n_estimators': trial.suggest_int('n_estimators', 50, 500),
                'max_depth': trial.suggest_int('max_depth', 3, 20),
                'min_samples_split': trial.suggest_int(
                    'min_samples_split', 2, 20),
                'min_samples_leaf': trial.suggest_int(
                    'min_samples_leaf', 1, 10),
                'max_features': trial.suggest_categorical(
                    'max_features', ['sqrt', 'log2', None]),
                'bootstrap': trial.suggest_categorical(
                    'bootstrap', [True, False]),
                'random_state': 42
            }
            
            # Create model
            model = RandomForestRegressor(**params)
            
            # Cross-validation
            cv_scores = cross_val_score(
                model, X, y, cv=cv_folds, 
                scoring='r2', n_jobs=-1
            )
            
            return cv_scores.mean()
            
        # Run optimization
        if self.study is None:
            self.create_study(direction="maximize")
            
        if self.study is not None:
            self.study.optimize(objective, n_trials=n_trials)
            
            # Store results
            self.best_params = self.study.best_params
            self.best_score = self.study.best_value
            
            logger.info(f"Best RF score: {self.best_score:.4f}")
            logger.info(f"Best RF params: {self.best_params}")
            
            return {
                'best_params': self.best_params,
                'best_score': self.best_score,
                'n_trials': len(self.study.trials)
            }
        
        return {}
        
    def optimize_lightgbm(self, X: np.ndarray, y: np.ndarray,
                          n_trials: int = 100,
                          cv_folds: int = 5) -> Dict[str, Any]:
        """Optimize LightGBM hyperparameters."""
        
        def objective(trial):
            # Suggest hyperparameters
            params = {
                'objective': 'regression',
                'metric': 'rmse',
                'boosting_type': 'gbdt',
                'num_leaves': trial.suggest_int('num_leaves', 10, 300),
                'learning_rate': trial.suggest_float(
                    'learning_rate', 0.01, 0.3),
                'feature_fraction': trial.suggest_float(
                    'feature_fraction', 0.4, 1.0),
                'bagging_fraction': trial.suggest_float(
                    'bagging_fraction', 0.4, 1.0),
                'bagging_freq': trial.suggest_int('bagging_freq', 1, 7),
                'min_child_samples': trial.suggest_int(
                    'min_child_samples', 5, 100),
                'reg_alpha': trial.suggest_float('reg_alpha', 0, 10),
                'reg_lambda': trial.suggest_float('reg_lambda', 0, 10),
                'random_state': 42,
                'verbosity': -1
            }
            
            # Cross-validation with LightGBM
            dtrain = lgb.Dataset(X, label=y)
            cv_results = lgb.cv(
                params,
                dtrain,
                num_boost_round=1000,
                nfold=cv_folds,
                stratified=False,
                shuffle=True,
                seed=42,
                return_cvbooster=True,
                callbacks=[lgb.early_stopping(50), lgb.log_evaluation(0)]
            )
            
            # Return negative RMSE (since we want to minimize RMSE)
            return -cv_results['valid rmse-mean'][-1]
            
        # Run optimization
        if self.study is None:
            self.create_study(direction="maximize")
            
        if self.study is not None:
            self.study.optimize(objective, n_trials=n_trials)
            
            # Store results
            self.best_params = self.study.best_params
            self.best_score = self.study.best_value
            
            logger.info(f"Best LGB score: {self.best_score:.4f}")
            logger.info(f"Best LGB params: {self.best_params}")
            
            return {
                'best_params': self.best_params,
                'best_score': self.best_score,
                'n_trials': len(self.study.trials)
            }
        
        return {}
        
    def compare_models(self, X: np.ndarray, y: np.ndarray,
                       models_to_test: Optional[List[str]] = None,
                       n_trials: int = 50) -> pd.DataFrame:
        """Compare multiple models with optimized hyperparameters."""
        
        if models_to_test is None:
            models_to_test = ['random_forest', 'lightgbm']
            
        results = []
        
        for model_name in models_to_test:
            logger.info(f"Optimizing {model_name}...")
            
            # Create separate study for each model
            self.study_name = f"aging_intervention_{model_name}"
            self.study = None
            
            try:
                if model_name == 'random_forest':
                    result = self.optimize_random_forest(X, y, n_trials)
                elif model_name == 'lightgbm':
                    result = self.optimize_lightgbm(X, y, n_trials)
                else:
                    logger.warning(f"Unknown model: {model_name}")
                    continue
                    
                if result:
                    results.append({
                        'model': model_name,
                        'best_score': result['best_score'],
                        'n_trials': result['n_trials'],
                        'best_params': result['best_params']
                    })
                
            except Exception as e:
                logger.error(f"Failed to optimize {model_name}: {e}")
                continue
                
        # Create comparison DataFrame
        comparison_df = pd.DataFrame(results)
        if not comparison_df.empty:
            comparison_df = comparison_df.sort_values(
                'best_score', ascending=False)
            
            logger.info("Model comparison results:")
            logger.info(comparison_df[['model', 'best_score']].to_string())
        
        return comparison_df


class CrossValidationEvaluator:
    """Enhanced cross-validation for aging intervention models."""
    
    def __init__(self, cv_folds: int = 5, random_state: int = 42):
        """Initialize cross-validation evaluator."""
        self.cv_folds = cv_folds
        self.random_state = random_state
        
    def evaluate_model(self, model, X: np.ndarray, y: np.ndarray,
                       scoring_metrics: Optional[List[str]] = None
                       ) -> Dict[str, Any]:
        """Comprehensive model evaluation with cross-validation."""
        
        if scoring_metrics is None:
            scoring_metrics = [
                'r2', 'neg_mean_squared_error', 'neg_mean_absolute_error']
            
        results = {}
        
        # Perform cross-validation for each metric
        for metric in scoring_metrics:
            cv_scores = cross_val_score(
                model, X, y, cv=self.cv_folds,
                scoring=metric, n_jobs=-1
            )
            
            results[metric] = {
                'mean': cv_scores.mean(),
                'std': cv_scores.std(),
                'scores': cv_scores.tolist()
            }
            
        return results


def optimize_aging_models(X: np.ndarray, y: np.ndarray,
                          models: Optional[List[str]] = None,
                          n_trials: int = 100) -> Dict[str, Any]:
    """Convenience function for optimizing aging intervention models."""
    
    if models is None:
        models = ['random_forest', 'lightgbm']
        
    optimizer = HyperparameterOptimizer()
    results = optimizer.compare_models(X, y, models, n_trials)
    
    if not results.empty:
        return {
            'comparison_results': results,
            'best_model': results.iloc[0]['model'],
            'best_score': results.iloc[0]['best_score']
        }
    
    return {'comparison_results': results}
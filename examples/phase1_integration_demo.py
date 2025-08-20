"""
Comprehensive integration example for CRISPR Toolkit Phase 1 enhancements.

This example demonstrates the complete workflow using:
- Real dataset integration
- MLflow experiment tracking
- Hyperparameter optimization
- Ensemble methods
"""

import logging
import sys
from pathlib import Path

import pandas as pd

# Add project root to path
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from src.crispr_toolkit.data.real_datasets import (
    get_intervention_target_genes,
    load_comprehensive_aging_dataset,
)
from src.crispr_toolkit.models.ensemble_methods import (
    create_aging_ensemble,
    evaluate_ensemble_performance,
)
from src.crispr_toolkit.models.experiment_tracking import ExperimentTracker
from src.crispr_toolkit.models.hyperparameter_optimization import optimize_aging_models

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def main():
    """Run comprehensive aging intervention prediction pipeline."""

    logger.info("🧬 Starting CRISPR Toolkit Phase 1 Integration Demo")

    # =====================================================
    # STEP 1: Load Real Aging Datasets
    # =====================================================
    logger.info("📊 Step 1: Loading real aging datasets...")

    try:
        X, y, feature_names = load_comprehensive_aging_dataset()
        logger.info(f"✅ Loaded dataset: {X.shape[0]} samples, {X.shape[1]} features")

        # Get intervention targets
        targets = get_intervention_target_genes()
        logger.info(f"✅ Retrieved {len(targets)} intervention categories")

    except Exception as e:
        logger.error(f"❌ Failed to load datasets: {e}")
        return

    # =====================================================
    # STEP 2: Initialize Experiment Tracking
    # =====================================================
    logger.info("🔬 Step 2: Setting up MLflow experiment tracking...")

    try:
        tracker = ExperimentTracker(
            experiment_name="aging_intervention_phase1",
            tracking_uri="sqlite:///mlruns.db"
        )

        run_id = tracker.start_run(
            run_name=f"comprehensive_pipeline_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}"
        )

        # Log dataset info
        tracker.log_params({
            'dataset_samples': X.shape[0],
            'dataset_features': X.shape[1],
            'target_variable': 'age',
            'intervention_categories': len(targets)
        })

        logger.info(f"✅ Started experiment run: {run_id}")

    except Exception as e:
        logger.error(f"❌ Failed to setup tracking: {e}")
        return

    # =====================================================
    # STEP 3: Split Data for Training/Testing
    # =====================================================
    logger.info("🔀 Step 3: Splitting data for training and testing...")

    from sklearn.model_selection import train_test_split

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=0.2, random_state=42, stratify=None
    )

    logger.info(f"✅ Train: {X_train.shape}, Test: {X_test.shape}")

    # =====================================================
    # STEP 4: Hyperparameter Optimization
    # =====================================================
    logger.info("⚙️ Step 4: Running hyperparameter optimization...")

    try:
        # Run limited optimization for demo (fewer trials)
        optimization_results = optimize_aging_models(
            X_train, y_train,
            models=['random_forest', 'lightgbm'],
            n_trials=10  # Reduced for demo speed
        )

        if 'best_model' in optimization_results:
            best_model_name = optimization_results['best_model']
            best_score = optimization_results['best_score']

            logger.info(f"✅ Best model: {best_model_name} (R² = {best_score:.4f})")

            # Log optimization results
            tracker.log_metrics({
                'best_cv_score': best_score,
                'optimization_trials': 10
            })

            tracker.log_params({
                'best_model_type': best_model_name
            })

        else:
            logger.warning("⚠️ No optimization results obtained")

    except Exception as e:
        logger.error(f"❌ Hyperparameter optimization failed: {e}")
        # Continue with default models

    # =====================================================
    # STEP 5: Create and Evaluate Ensemble Models
    # =====================================================
    logger.info("🤝 Step 5: Creating ensemble models...")

    ensemble_results = {}

    for ensemble_type in ['voting', 'dynamic']:
        try:
            logger.info(f"  Creating {ensemble_type} ensemble...")

            # Create ensemble
            ensemble = create_aging_ensemble(
                X_train, y_train,
                ensemble_type=ensemble_type
            )

            if ensemble is not None:
                # Train ensemble
                ensemble.fit(X_train, y_train)

                # Evaluate on test set
                performance = evaluate_ensemble_performance(
                    ensemble, X_test, y_test
                )

                ensemble_results[ensemble_type] = performance

                logger.info(f"  ✅ {ensemble_type.title()} ensemble - R²: {performance['r2_score']:.4f}")

                # Log ensemble performance
                for metric, value in performance.items():
                    tracker.log_metrics({
                        f'{ensemble_type}_ensemble_{metric}': value
                    })

            else:
                logger.warning(f"  ⚠️ Failed to create {ensemble_type} ensemble")

        except Exception as e:
            logger.error(f"  ❌ {ensemble_type} ensemble failed: {e}")
            continue

    # =====================================================
    # STEP 6: Feature Importance Analysis
    # =====================================================
    logger.info("🔍 Step 6: Analyzing feature importance...")

    try:
        # Create a simple random forest for feature analysis
        from sklearn.ensemble import RandomForestRegressor

        rf_model = RandomForestRegressor(n_estimators=100, random_state=42)
        rf_model.fit(X_train, y_train)

        # Get feature importances
        feature_importance = pd.DataFrame({
            'feature': feature_names,
            'importance': rf_model.feature_importances_
        }).sort_values('importance', ascending=False)

        # Log top features
        top_features = feature_importance.head(10)
        logger.info("  📊 Top 10 most important features:")
        for _, row in top_features.iterrows():
            logger.info(f"    {row['feature']}: {row['importance']:.4f}")

        # Track feature importance
        tracker.log_metrics({
            'top_feature_importance': top_features.iloc[0]['importance'],
            'mean_feature_importance': feature_importance['importance'].mean()
        })

    except Exception as e:
        logger.error(f"❌ Feature analysis failed: {e}")

    # =====================================================
    # STEP 7: Generate Intervention Recommendations
    # =====================================================
    logger.info("💡 Step 7: Generating intervention recommendations...")

    try:
        # Identify aging-related features that are highly important
        aging_genes = []
        for category, genes in targets.items():
            aging_genes.extend(genes)

        # Find overlap between important features and known aging genes
        important_aging_features = []
        if 'feature_importance' in locals():
            for _, row in feature_importance.head(20).iterrows():
                if row['feature'] in aging_genes:
                    important_aging_features.append({
                        'gene': row['feature'],
                        'importance': row['importance'],
                        'category': [cat for cat, genes in targets.items()
                                  if row['feature'] in genes][0]
                    })

        if important_aging_features:
            logger.info("  🎯 High-priority intervention targets identified:")
            for target in important_aging_features[:5]:
                logger.info(f"    {target['gene']} ({target['category']}) - Importance: {target['importance']:.4f}")

            tracker.log_params({
                'priority_targets': [t['gene'] for t in important_aging_features[:5]]
            })
        else:
            logger.info("  ℹ️ No overlap found between important features and known aging genes")

    except Exception as e:
        logger.error(f"❌ Intervention analysis failed: {e}")

    # =====================================================
    # STEP 8: Summary and Results
    # =====================================================
    logger.info("📋 Step 8: Generating summary report...")

    try:
        # Calculate overall pipeline success metrics
        pipeline_success = {
            'data_loaded': X.shape[0] > 0,
            'models_trained': len(ensemble_results) > 0,
            'best_performance': max([r['r2_score'] for r in ensemble_results.values()]) if ensemble_results else 0
        }

        # Log pipeline summary
        tracker.log_metrics({
            'pipeline_success_rate': sum(pipeline_success.values()) / len(pipeline_success),
            'total_models_evaluated': len(ensemble_results)
        })

        # Generate summary
        logger.info("🎉 Phase 1 Integration Pipeline Complete!")
        logger.info("="*60)
        logger.info(f"📈 Dataset: {X.shape[0]:,} samples, {X.shape[1]:,} features")
        logger.info(f"🔬 Experiment tracked in MLflow run: {run_id}")
        logger.info(f"🤖 Ensemble models created: {len(ensemble_results)}")

        if ensemble_results:
            best_ensemble = max(ensemble_results.keys(),
                              key=lambda k: ensemble_results[k]['r2_score'])
            best_score = ensemble_results[best_ensemble]['r2_score']
            logger.info(f"🏆 Best performing ensemble: {best_ensemble} (R² = {best_score:.4f})")

        logger.info("="*60)

        # End experiment run
        tracker.end_run()
        logger.info("✅ Experiment tracking completed")

    except Exception as e:
        logger.error(f"❌ Summary generation failed: {e}")

    logger.info("🎯 Phase 1 integration demo completed successfully!")


if __name__ == "__main__":
    main()

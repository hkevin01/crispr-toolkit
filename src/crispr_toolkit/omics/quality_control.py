"""
Multi-Omics Quality Control Pipeline for CRISPR Toolkit Phase 3
===============================================================

Comprehensive quality control and data validation pipeline for
multi-omics aging intervention research, including batch effect
detection, outlier identification, and data normalization.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import levene, normaltest
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest
from sklearn.preprocessing import RobustScaler, StandardScaler

logger = logging.getLogger(__name__)


@dataclass
class QualityMetrics:
    """Quality metrics for omics data."""
    completeness: float
    distribution_normality: float
    variance_homogeneity: float
    outlier_fraction: float
    batch_effect_strength: float
    technical_noise: float
    overall_quality_score: float


@dataclass
class BatchEffectResult:
    """Results from batch effect analysis."""
    batch_variance_explained: float
    significant_batch_features: List[str]
    batch_correction_recommended: bool
    correction_method: str
    correction_strength: float


@dataclass
class OutlierDetectionResult:
    """Results from outlier detection analysis."""
    outlier_samples: List[str]
    outlier_features: List[str]
    outlier_scores: Dict[str, float]
    removal_recommended: bool
    outlier_categories: Dict[str, List[str]]


@dataclass
class NormalizationResult:
    """Results from data normalization."""
    normalization_method: str
    pre_normalization_stats: Dict[str, float]
    post_normalization_stats: Dict[str, float]
    normalization_success: bool
    quality_improvement: float


class MultiOmicsQualityControl:
    """
    Comprehensive quality control pipeline for multi-omics data.

    This class implements advanced quality control procedures including
    batch effect detection, outlier identification, normalization,
    and overall data quality assessment for aging research.
    """

    def __init__(self, qc_stringency: str = "standard"):
        """
        Initialize multi-omics quality control pipeline.

        Args:
            qc_stringency: "lenient", "standard", "strict"
        """
        self.qc_stringency = qc_stringency
        self.omics_data = {}
        self.qc_results = {}
        self.quality_thresholds = self._set_quality_thresholds(qc_stringency)

        logger.info(f"Initialized MultiOmicsQualityControl with {qc_stringency} stringency")

    def _set_quality_thresholds(self, stringency: str) -> Dict[str, float]:
        """Set quality control thresholds based on stringency level."""

        if stringency == "lenient":
            thresholds = {
                'min_completeness': 0.7,
                'max_outlier_fraction': 0.15,
                'min_normality_pvalue': 0.001,
                'max_batch_variance': 0.3,
                'min_quality_score': 0.5
            }
        elif stringency == "strict":
            thresholds = {
                'min_completeness': 0.95,
                'max_outlier_fraction': 0.05,
                'min_normality_pvalue': 0.05,
                'max_batch_variance': 0.1,
                'min_quality_score': 0.8
            }
        else:  # standard
            thresholds = {
                'min_completeness': 0.85,
                'max_outlier_fraction': 0.1,
                'min_normality_pvalue': 0.01,
                'max_batch_variance': 0.2,
                'min_quality_score': 0.7
            }

        return thresholds

    def add_omics_data(self,
                      omics_type: str,
                      data: pd.DataFrame,
                      metadata: Optional[pd.DataFrame] = None,
                      batch_info: Optional[pd.Series] = None):
        """Add omics data for quality control analysis."""

        self.omics_data[omics_type] = {
            'data': data,
            'metadata': metadata if metadata is not None else pd.DataFrame(),
            'batch_info': batch_info
        }

        logger.info(f"Added {omics_type} data: {data.shape}")

    def run_comprehensive_qc(self) -> Dict[str, Dict[str, Any]]:
        """
        Run comprehensive quality control analysis on all omics data.

        Returns:
            Dictionary of QC results for each omics type
        """
        logger.info("Running comprehensive quality control analysis...")

        qc_results = {}

        for omics_type, omics_info in self.omics_data.items():
            logger.info(f"Analyzing quality for {omics_type}...")

            data = omics_info['data']
            metadata = omics_info['metadata']
            batch_info = omics_info['batch_info']

            # Run individual QC components
            quality_metrics = self._calculate_quality_metrics(data)
            batch_results = self._analyze_batch_effects(data, batch_info)
            outlier_results = self._detect_outliers(data)
            normalization_results = self._evaluate_normalization(data)

            # Compile results
            qc_results[omics_type] = {
                'quality_metrics': quality_metrics,
                'batch_effects': batch_results,
                'outliers': outlier_results,
                'normalization': normalization_results,
                'recommendations': self._generate_qc_recommendations(
                    quality_metrics, batch_results, outlier_results, normalization_results
                )
            }

        self.qc_results = qc_results

        logger.info(f"Quality control completed for {len(qc_results)} omics types")
        return qc_results

    def _calculate_quality_metrics(self, data: pd.DataFrame) -> QualityMetrics:
        """Calculate comprehensive quality metrics for omics data."""

        # Data completeness
        total_values = data.shape[0] * data.shape[1]
        missing_values = data.isna().sum().sum()
        completeness = 1.0 - (missing_values / total_values)

        # Distribution normality (sample mean of feature normality tests)
        normality_pvalues = []
        for feature in data.index:
            feature_data = data.loc[feature].dropna()
            if len(feature_data) > 8:  # Need sufficient data for normality test
                try:
                    _, p_value = normaltest(feature_data)
                    normality_pvalues.append(p_value)
                except:
                    continue

        distribution_normality = np.mean(normality_pvalues) if normality_pvalues else 0.0

        # Variance homogeneity across samples
        sample_variances = data.var(axis=0, skipna=True)
        try:
            # Levene's test for equal variances (using a subset if too many samples)
            if len(sample_variances) > 50:
                subset_indices = np.random.choice(len(sample_variances), 50, replace=False)
                subset_data = [data.iloc[:, i].dropna() for i in subset_indices]
            else:
                subset_data = [data.iloc[:, i].dropna() for i in range(len(sample_variances))]

            _, levene_p = levene(*subset_data)
            variance_homogeneity = levene_p
        except:
            variance_homogeneity = 0.0

        # Outlier detection (using isolation forest)
        outlier_fraction = self._calculate_outlier_fraction(data)

        # Batch effect strength (placeholder - would require batch information)
        batch_effect_strength = 0.1  # Default low batch effect

        # Technical noise (coefficient of variation)
        feature_cvs = (data.std(axis=1) / data.mean(axis=1)).fillna(0)
        technical_noise = feature_cvs.mean()

        # Overall quality score (weighted combination)
        quality_components = [
            completeness * 0.3,
            min(distribution_normality * 10, 1.0) * 0.2,  # Scale p-value to 0-1
            min(variance_homogeneity * 10, 1.0) * 0.2,
            (1.0 - outlier_fraction) * 0.2,
            (1.0 - batch_effect_strength) * 0.1
        ]

        overall_quality_score = sum(quality_components)

        return QualityMetrics(
            completeness=completeness,
            distribution_normality=distribution_normality,
            variance_homogeneity=variance_homogeneity,
            outlier_fraction=outlier_fraction,
            batch_effect_strength=batch_effect_strength,
            technical_noise=technical_noise,
            overall_quality_score=overall_quality_score
        )

    def _calculate_outlier_fraction(self, data: pd.DataFrame) -> float:
        """Calculate fraction of outlier samples using isolation forest."""

        try:
            # Transpose for sample-wise outlier detection
            sample_data = data.T.fillna(data.T.mean())

            # Use subset of features if too many
            if sample_data.shape[1] > 500:
                feature_variance = sample_data.var()
                top_features = feature_variance.nlargest(500).index
                sample_data = sample_data[top_features]

            # Isolation forest
            iso_forest = IsolationForest(contamination=0.1, random_state=42)
            outlier_labels = iso_forest.fit_predict(sample_data)

            outlier_fraction = (outlier_labels == -1).sum() / len(outlier_labels)

        except Exception as e:
            logger.warning(f"Outlier detection failed: {e}")
            outlier_fraction = 0.05  # Default small outlier fraction

        return outlier_fraction

    def _analyze_batch_effects(self,
                             data: pd.DataFrame,
                             batch_info: Optional[pd.Series]) -> BatchEffectResult:
        """Analyze batch effects in omics data."""

        if batch_info is None:
            # No batch information available
            return BatchEffectResult(
                batch_variance_explained=0.0,
                significant_batch_features=[],
                batch_correction_recommended=False,
                correction_method="none",
                correction_strength=0.0
            )

        # Find common samples between data and batch info
        common_samples = list(set(data.columns) & set(batch_info.index))

        if len(common_samples) < 10:
            logger.warning("Insufficient samples with batch information")
            return BatchEffectResult(
                batch_variance_explained=0.0,
                significant_batch_features=[],
                batch_correction_recommended=False,
                correction_method="insufficient_data",
                correction_strength=0.0
            )

        # Subset data to common samples
        subset_data = data[common_samples]
        subset_batch = batch_info[common_samples]

        # Calculate batch variance for each feature
        significant_batch_features = []
        batch_variances = []

        for feature in subset_data.index:
            feature_data = subset_data.loc[feature]

            # Group by batch and calculate variance components
            batch_groups = feature_data.groupby(subset_batch)

            if len(batch_groups) > 1:
                # ANOVA-like analysis
                try:
                    # Between-batch variance vs within-batch variance
                    overall_mean = feature_data.mean()

                    between_batch_var = 0
                    within_batch_var = 0
                    total_samples = 0

                    for batch, group_data in batch_groups:
                        group_mean = group_data.mean()
                        group_size = len(group_data)

                        # Between-batch variance
                        between_batch_var += group_size * (group_mean - overall_mean) ** 2

                        # Within-batch variance
                        within_batch_var += ((group_data - group_mean) ** 2).sum()

                        total_samples += group_size

                    # Variance explained by batch
                    total_var = between_batch_var + within_batch_var
                    if total_var > 0:
                        batch_variance_explained = between_batch_var / total_var
                        batch_variances.append(batch_variance_explained)

                        # Consider significant if batch explains >20% of variance
                        if batch_variance_explained > 0.2:
                            significant_batch_features.append(feature)

                except:
                    continue

        # Overall batch effect strength
        mean_batch_variance = np.mean(batch_variances) if batch_variances else 0.0

        # Recommend correction if batch effects are strong
        batch_correction_recommended = mean_batch_variance > self.quality_thresholds['max_batch_variance']

        # Suggest correction method
        if mean_batch_variance > 0.3:
            correction_method = "ComBat"
            correction_strength = mean_batch_variance
        elif mean_batch_variance > 0.15:
            correction_method = "linear_model"
            correction_strength = mean_batch_variance
        else:
            correction_method = "none"
            correction_strength = 0.0

        return BatchEffectResult(
            batch_variance_explained=mean_batch_variance,
            significant_batch_features=significant_batch_features,
            batch_correction_recommended=batch_correction_recommended,
            correction_method=correction_method,
            correction_strength=correction_strength
        )

    def _detect_outliers(self, data: pd.DataFrame) -> OutlierDetectionResult:
        """Detect outliers in omics data using multiple methods."""

        outlier_samples = []
        outlier_features = []
        outlier_scores = {}
        outlier_categories = {
            'sample_outliers': [],
            'feature_outliers': [],
            'extreme_values': []
        }

        # Sample-wise outlier detection
        try:
            sample_data = data.T.fillna(data.T.mean())  # Samples x features

            # Use PCA for dimensionality reduction if needed
            if sample_data.shape[1] > 1000:
                pca = PCA(n_components=min(50, sample_data.shape[0] - 1))
                sample_data_reduced = pca.fit_transform(sample_data)
            else:
                sample_data_reduced = sample_data.values

            # Isolation forest for sample outliers
            iso_forest = IsolationForest(contamination=0.1, random_state=42)
            sample_outlier_labels = iso_forest.fit_predict(sample_data_reduced)

            # Get outlier samples
            sample_outliers = [sample for i, sample in enumerate(sample_data.index)
                             if sample_outlier_labels[i] == -1]
            outlier_samples.extend(sample_outliers)
            outlier_categories['sample_outliers'] = sample_outliers

            # Calculate outlier scores
            outlier_scores_array = iso_forest.decision_function(sample_data_reduced)
            for i, sample in enumerate(sample_data.index):
                outlier_scores[f"sample_{sample}"] = outlier_scores_array[i]

        except Exception as e:
            logger.warning(f"Sample outlier detection failed: {e}")

        # Feature-wise outlier detection
        try:
            for feature in data.index:
                feature_data = data.loc[feature].dropna()

                if len(feature_data) > 10:
                    # Z-score based outlier detection
                    z_scores = np.abs(stats.zscore(feature_data))
                    extreme_samples = feature_data[z_scores > 3].index.tolist()

                    if len(extreme_samples) > 0:
                        outlier_categories['extreme_values'].extend(extreme_samples)

                    # Feature-level quality metrics
                    feature_cv = feature_data.std() / feature_data.mean() if feature_data.mean() != 0 else 0
                    feature_skewness = abs(stats.skew(feature_data))

                    # Consider feature as outlier if extremely variable or skewed
                    if feature_cv > 2.0 or feature_skewness > 2.0:
                        outlier_features.append(feature)
                        outlier_scores[f"feature_{feature}"] = max(feature_cv, feature_skewness)

        except Exception as e:
            logger.warning(f"Feature outlier detection failed: {e}")

        outlier_categories['feature_outliers'] = outlier_features

        # Remove duplicates
        outlier_samples = list(set(outlier_samples))
        outlier_features = list(set(outlier_features))

        # Determine if removal is recommended
        sample_outlier_fraction = len(outlier_samples) / data.shape[1]
        feature_outlier_fraction = len(outlier_features) / data.shape[0]

        removal_recommended = (sample_outlier_fraction > self.quality_thresholds['max_outlier_fraction'] or
                             feature_outlier_fraction > self.quality_thresholds['max_outlier_fraction'])

        return OutlierDetectionResult(
            outlier_samples=outlier_samples,
            outlier_features=outlier_features,
            outlier_scores=outlier_scores,
            removal_recommended=removal_recommended,
            outlier_categories=outlier_categories
        )

    def _evaluate_normalization(self, data: pd.DataFrame) -> NormalizationResult:
        """Evaluate different normalization methods for the data."""

        # Calculate pre-normalization statistics
        pre_stats = {
            'mean_cv': (data.std(axis=1) / data.mean(axis=1)).fillna(0).mean(),
            'mean_skewness': np.mean([abs(stats.skew(data.loc[feat].dropna()))
                                    for feat in data.index]),
            'mean_kurtosis': np.mean([abs(stats.kurtosis(data.loc[feat].dropna()))
                                    for feat in data.index])
        }

        # Test different normalization methods
        normalization_methods = ['log2', 'zscore', 'robust', 'quantile']
        best_method = 'none'
        best_improvement = 0.0
        post_stats = pre_stats.copy()

        for method in normalization_methods:
            try:
                if method == 'log2':
                    # Log2 transformation (add 1 to avoid log(0))
                    normalized_data = np.log2(data + 1)
                elif method == 'zscore':
                    # Z-score normalization
                    scaler = StandardScaler()
                    normalized_data = pd.DataFrame(
                        scaler.fit_transform(data.T).T,
                        index=data.index, columns=data.columns
                    )
                elif method == 'robust':
                    # Robust scaling
                    scaler = RobustScaler()
                    normalized_data = pd.DataFrame(
                        scaler.fit_transform(data.T).T,
                        index=data.index, columns=data.columns
                    )
                elif method == 'quantile':
                    # Quantile normalization (simplified)
                    normalized_data = data.rank(axis=1) / data.shape[1]
                else:
                    continue

                # Calculate post-normalization statistics
                method_stats = {
                    'mean_cv': (normalized_data.std(axis=1) / normalized_data.mean(axis=1)).fillna(0).mean(),
                    'mean_skewness': np.mean([abs(stats.skew(normalized_data.loc[feat].dropna()))
                                            for feat in normalized_data.index]),
                    'mean_kurtosis': np.mean([abs(stats.kurtosis(normalized_data.loc[feat].dropna()))
                                            for feat in normalized_data.index])
                }

                # Calculate improvement score
                cv_improvement = max(0, pre_stats['mean_cv'] - method_stats['mean_cv'])
                skew_improvement = max(0, pre_stats['mean_skewness'] - method_stats['mean_skewness'])
                kurt_improvement = max(0, pre_stats['mean_kurtosis'] - method_stats['mean_kurtosis'])

                improvement = cv_improvement + skew_improvement + kurt_improvement

                if improvement > best_improvement:
                    best_improvement = improvement
                    best_method = method
                    post_stats = method_stats

            except Exception as e:
                logger.warning(f"Normalization method {method} failed: {e}")
                continue

        # Determine if normalization was successful
        normalization_success = best_improvement > 0.1

        return NormalizationResult(
            normalization_method=best_method,
            pre_normalization_stats=pre_stats,
            post_normalization_stats=post_stats,
            normalization_success=normalization_success,
            quality_improvement=best_improvement
        )

    def _generate_qc_recommendations(self,
                                   quality_metrics: QualityMetrics,
                                   batch_results: BatchEffectResult,
                                   outlier_results: OutlierDetectionResult,
                                   normalization_results: NormalizationResult) -> List[str]:
        """Generate quality control recommendations."""

        recommendations = []

        # Data completeness recommendations
        if quality_metrics.completeness < self.quality_thresholds['min_completeness']:
            recommendations.append(
                f"Low data completeness ({quality_metrics.completeness:.2%}) - "
                "consider imputation or additional data collection"
            )

        # Distribution recommendations
        if quality_metrics.distribution_normality < self.quality_thresholds['min_normality_pvalue']:
            recommendations.append("Non-normal distributions detected - consider data transformation")

        # Outlier recommendations
        if outlier_results.removal_recommended:
            recommendations.append(
                f"Significant outliers detected ({len(outlier_results.outlier_samples)} samples, "
                f"{len(outlier_results.outlier_features)} features) - review and consider removal"
            )

        # Batch effect recommendations
        if batch_results.batch_correction_recommended:
            recommendations.append(
                f"Significant batch effects detected - recommend {batch_results.correction_method} correction"
            )

        # Normalization recommendations
        if normalization_results.normalization_success:
            recommendations.append(
                f"Apply {normalization_results.normalization_method} normalization "
                f"for improved data quality"
            )

        # Overall quality recommendations
        if quality_metrics.overall_quality_score < self.quality_thresholds['min_quality_score']:
            recommendations.append(
                f"Overall quality score low ({quality_metrics.overall_quality_score:.2f}) - "
                "comprehensive data cleaning recommended"
            )
        else:
            recommendations.append("Data quality is acceptable for downstream analysis")

        return recommendations

    def apply_quality_corrections(self,
                                omics_type: str,
                                apply_normalization: bool = True,
                                remove_outliers: bool = True,
                                correct_batch_effects: bool = True) -> pd.DataFrame:
        """
        Apply recommended quality corrections to omics data.

        Args:
            omics_type: Type of omics data to correct
            apply_normalization: Whether to apply normalization
            remove_outliers: Whether to remove outlier samples/features
            correct_batch_effects: Whether to correct batch effects

        Returns:
            Corrected data matrix
        """
        logger.info(f"Applying quality corrections to {omics_type} data...")

        if omics_type not in self.qc_results:
            raise ValueError(f"No QC results available for {omics_type}")

        data = self.omics_data[omics_type]['data'].copy()
        qc_result = self.qc_results[omics_type]

        # Remove outlier features
        if remove_outliers and qc_result['outliers'].removal_recommended:
            outlier_features = qc_result['outliers'].outlier_features
            data = data.drop(index=outlier_features, errors='ignore')
            logger.info(f"Removed {len(outlier_features)} outlier features")

        # Remove outlier samples
        if remove_outliers and qc_result['outliers'].removal_recommended:
            outlier_samples = qc_result['outliers'].outlier_samples
            data = data.drop(columns=outlier_samples, errors='ignore')
            logger.info(f"Removed {len(outlier_samples)} outlier samples")

        # Apply normalization
        if apply_normalization and qc_result['normalization'].normalization_success:
            method = qc_result['normalization'].normalization_method

            if method == 'log2':
                data = np.log2(data + 1)
            elif method == 'zscore':
                scaler = StandardScaler()
                data = pd.DataFrame(
                    scaler.fit_transform(data.T).T,
                    index=data.index, columns=data.columns
                )
            elif method == 'robust':
                scaler = RobustScaler()
                data = pd.DataFrame(
                    scaler.fit_transform(data.T).T,
                    index=data.index, columns=data.columns
                )
            elif method == 'quantile':
                data = data.rank(axis=1) / data.shape[1]

            logger.info(f"Applied {method} normalization")

        # Batch effect correction (placeholder - would implement specific methods)
        if correct_batch_effects and qc_result['batch_effects'].batch_correction_recommended:
            logger.info("Batch effect correction recommended but not implemented in this demo")

        logger.info(f"Quality correction completed. Final data shape: {data.shape}")
        return data

    def generate_qc_report(self) -> Dict[str, Any]:
        """Generate comprehensive quality control report."""

        if not self.qc_results:
            raise ValueError("No QC results available. Run run_comprehensive_qc() first.")

        report = {
            'summary': {},
            'omics_quality': {},
            'overall_recommendations': [],
            'quality_scores': {}
        }

        # Overall summary
        total_omics = len(self.qc_results)
        passing_omics = sum([
            1 for result in self.qc_results.values()
            if result['quality_metrics'].overall_quality_score >= self.quality_thresholds['min_quality_score']
        ])

        report['summary'] = {
            'total_omics_analyzed': total_omics,
            'omics_passing_qc': passing_omics,
            'overall_pass_rate': passing_omics / total_omics if total_omics > 0 else 0,
            'qc_stringency': self.qc_stringency
        }

        # Individual omics quality
        for omics_type, qc_result in self.qc_results.items():
            quality_metrics = qc_result['quality_metrics']

            report['omics_quality'][omics_type] = {
                'overall_quality_score': quality_metrics.overall_quality_score,
                'completeness': quality_metrics.completeness,
                'outlier_fraction': quality_metrics.outlier_fraction,
                'batch_effect_strength': quality_metrics.batch_effect_strength,
                'recommendations': qc_result['recommendations']
            }

            report['quality_scores'][omics_type] = quality_metrics.overall_quality_score

        # Overall recommendations
        recommendations = []

        if report['summary']['overall_pass_rate'] < 0.7:
            recommendations.append("Multiple omics datasets have quality issues - comprehensive review needed")

        batch_issues = sum([
            1 for result in self.qc_results.values()
            if result['batch_effects'].batch_correction_recommended
        ])

        if batch_issues > 0:
            recommendations.append(f"Batch effects detected in {batch_issues} omics types - correction recommended")

        outlier_issues = sum([
            1 for result in self.qc_results.values()
            if result['outliers'].removal_recommended
        ])

        if outlier_issues > 0:
            recommendations.append(f"Outliers detected in {outlier_issues} omics types - review recommended")

        if len(recommendations) == 0:
            recommendations.append("Overall data quality is good - proceed with analysis")

        report['overall_recommendations'] = recommendations

        return report


def create_quality_control_demo():
    """Create demonstration of multi-omics quality control."""

    # Create synthetic multi-omics data with quality issues
    n_samples = 100
    n_features = 500
    sample_names = [f"Sample_{i}" for i in range(n_samples)]

    np.random.seed(42)

    # Transcriptomics data with some issues
    transcriptomics_data = pd.DataFrame(
        np.random.lognormal(5, 1, (n_features, n_samples)),
        index=[f"Gene_{i}" for i in range(n_features)],
        columns=sample_names
    )

    # Add missing values
    missing_mask = np.random.random((n_features, n_samples)) < 0.1
    transcriptomics_data[missing_mask] = np.nan

    # Add outlier samples
    outlier_samples = sample_names[:5]
    transcriptomics_data[outlier_samples] *= 3  # Make them outliers

    # Proteomics data
    proteomics_data = pd.DataFrame(
        np.random.lognormal(10, 0.8, (200, n_samples)),
        index=[f"Protein_{i}" for i in range(200)],
        columns=sample_names
    )

    # Create batch information
    batch_info = pd.Series(
        np.random.choice(['Batch_A', 'Batch_B', 'Batch_C'], n_samples),
        index=sample_names
    )

    # Add batch effects to some features
    for i in range(50):
        batch_effect = batch_info.map({'Batch_A': 0, 'Batch_B': 2, 'Batch_C': -1})
        transcriptomics_data.iloc[i] += batch_effect

    # Initialize QC pipeline
    qc_pipeline = MultiOmicsQualityControl(qc_stringency="standard")

    # Add omics data
    qc_pipeline.add_omics_data("transcriptomics", transcriptomics_data, batch_info=batch_info)
    qc_pipeline.add_omics_data("proteomics", proteomics_data)

    # Run comprehensive QC
    qc_results = qc_pipeline.run_comprehensive_qc()

    # Generate report
    qc_report = qc_pipeline.generate_qc_report()

    # Apply corrections
    corrected_transcriptomics = qc_pipeline.apply_quality_corrections(
        "transcriptomics",
        apply_normalization=True,
        remove_outliers=True,
        correct_batch_effects=False
    )

    return qc_pipeline, qc_results, qc_report, corrected_transcriptomics


if __name__ == "__main__":
    # Run demonstration
    qc_pipeline, qc_results, qc_report, corrected_data = create_quality_control_demo()

    print("Multi-omics quality control completed!")
    print(f"Omics types analyzed: {qc_report['summary']['total_omics_analyzed']}")
    print(f"Overall pass rate: {qc_report['summary']['overall_pass_rate']:.2%}")

    for omics_type, quality_info in qc_report['omics_quality'].items():
        print(f"{omics_type} quality score: {quality_info['overall_quality_score']:.3f}")
        print(f"  Completeness: {quality_info['completeness']:.2%}")
        print(f"  Outlier fraction: {quality_info['outlier_fraction']:.2%}")

    print("Overall recommendations:")
    for rec in qc_report['overall_recommendations']:
        print(f"  - {rec}")

    print(f"Corrected data shape: {corrected_data.shape}")

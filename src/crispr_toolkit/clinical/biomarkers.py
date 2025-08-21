"""
Real-time Biomarker Discovery and Validation for CRISPR Toolkit Phase 3
========================================================================

Advanced biomarker discovery system for aging intervention clinical trials
with real-time validation, multi-omics integration, and automated
biomarker qualification for regulatory submission.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, pearsonr, spearmanr
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.feature_selection import SelectKBest, f_regression, mutual_info_regression
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler

logger = logging.getLogger(__name__)


@dataclass
class BiomarkerCandidate:
    """Individual biomarker candidate with validation metrics."""
    biomarker_id: str
    biomarker_name: str
    biomarker_type: str  # 'protein', 'metabolite', 'gene', 'composite'
    omics_source: str
    discovery_method: str
    statistical_significance: float
    effect_size: float
    validation_score: float
    clinical_relevance: float
    regulatory_qualification: str  # 'exploratory', 'probable', 'known'
    reproducibility: float


@dataclass
class BiomarkerPanel:
    """Panel of validated biomarkers for specific indication."""
    panel_id: str
    indication: str
    biomarkers: List[BiomarkerCandidate]
    combined_performance: Dict[str, float]
    validation_studies: List[str]
    regulatory_status: str
    clinical_utility: Dict[str, float]


@dataclass
class BiomarkerValidationResult:
    """Results from biomarker validation analysis."""
    discovered_biomarkers: List[BiomarkerCandidate]
    validated_panels: List[BiomarkerPanel]
    discovery_statistics: Dict[str, Any]
    validation_metrics: Dict[str, float]
    regulatory_readiness: Dict[str, str]
    longitudinal_stability: Dict[str, float]


class RealTimeBiomarkerEngine:
    """
    Real-time biomarker discovery and validation engine.

    This class implements advanced biomarker discovery workflows
    for aging interventions with continuous validation, regulatory
    qualification, and clinical utility assessment.
    """

    def __init__(self, discovery_mode: str = "comprehensive"):
        """
        Initialize real-time biomarker discovery engine.

        Args:
            discovery_mode: "rapid", "comprehensive", "regulatory"
        """
        self.discovery_mode = discovery_mode
        self.omics_data = {}
        self.clinical_data = None
        self.outcomes_data = None
        self.longitudinal_data = {}
        self.biomarker_candidates = []
        self.validated_panels = []

        logger.info(f"Initialized RealTimeBiomarkerEngine in {discovery_mode} mode")

    def load_study_data(self,
                       clinical_data: pd.DataFrame,
                       outcomes_data: pd.DataFrame,
                       omics_data: Dict[str, pd.DataFrame],
                       longitudinal_data: Optional[Dict[str, pd.DataFrame]] = None):
        """Load study data for biomarker discovery."""

        self.clinical_data = clinical_data
        self.outcomes_data = outcomes_data
        self.omics_data = omics_data

        if longitudinal_data:
            self.longitudinal_data = longitudinal_data

        logger.info("Loaded study data:")
        logger.info(f"  Clinical: {clinical_data.shape}")
        logger.info(f"  Outcomes: {outcomes_data.shape}")
        logger.info(f"  Omics layers: {list(omics_data.keys())}")
        if longitudinal_data:
            logger.info(f"  Longitudinal timepoints: {len(longitudinal_data)}")

    def discover_biomarkers(self,
                          target_outcome: str,
                          discovery_criteria: Optional[Dict[str, Any]] = None) -> BiomarkerValidationResult:
        """
        Perform comprehensive biomarker discovery and validation.

        Args:
            target_outcome: Primary outcome variable for biomarker discovery
            discovery_criteria: Criteria for biomarker selection and validation

        Returns:
            BiomarkerValidationResult with discovered and validated biomarkers
        """
        logger.info(f"Starting biomarker discovery for outcome: {target_outcome}")

        if discovery_criteria is None:
            discovery_criteria = self._get_default_discovery_criteria()

        # Phase 1: Discovery
        discovered_biomarkers = self._perform_biomarker_discovery(target_outcome, discovery_criteria)

        # Phase 2: Validation
        validated_biomarkers = self._validate_biomarker_candidates(discovered_biomarkers, target_outcome)

        # Phase 3: Panel Construction
        biomarker_panels = self._construct_biomarker_panels(validated_biomarkers, target_outcome)

        # Phase 4: Regulatory Assessment
        regulatory_assessment = self._assess_regulatory_readiness(validated_biomarkers, biomarker_panels)

        # Phase 5: Longitudinal Stability
        longitudinal_stability = self._assess_longitudinal_stability(validated_biomarkers)

        # Compile results
        discovery_statistics = self._compile_discovery_statistics(discovered_biomarkers)
        validation_metrics = self._compile_validation_metrics(validated_biomarkers)

        result = BiomarkerValidationResult(
            discovered_biomarkers=validated_biomarkers,
            validated_panels=biomarker_panels,
            discovery_statistics=discovery_statistics,
            validation_metrics=validation_metrics,
            regulatory_readiness=regulatory_assessment,
            longitudinal_stability=longitudinal_stability
        )

        logger.info("Biomarker discovery completed:")
        logger.info(f"  Discovered: {len(discovered_biomarkers)} candidates")
        logger.info(f"  Validated: {len(validated_biomarkers)} biomarkers")
        logger.info(f"  Panels: {len(biomarker_panels)} constructed")

        return result

    def _perform_biomarker_discovery(self,
                                   target_outcome: str,
                                   criteria: Dict[str, Any]) -> List[BiomarkerCandidate]:
        """Perform initial biomarker discovery across all omics layers."""

        discovered_candidates = []

        # Get target outcome data
        if target_outcome not in self.outcomes_data.columns:
            raise ValueError(f"Target outcome {target_outcome} not found in outcomes data")

        target_values = self.outcomes_data[target_outcome]

        # Discovery across each omics layer
        for omics_type, omics_df in self.omics_data.items():
            logger.info(f"Discovering biomarkers in {omics_type} data...")

            omics_candidates = self._discover_omics_biomarkers(
                omics_df, target_values, omics_type, criteria
            )
            discovered_candidates.extend(omics_candidates)

        # Clinical biomarker discovery
        if self.clinical_data is not None:
            clinical_candidates = self._discover_clinical_biomarkers(
                self.clinical_data, target_values, criteria
            )
            discovered_candidates.extend(clinical_candidates)

        # Composite biomarker discovery
        composite_candidates = self._discover_composite_biomarkers(
            discovered_candidates, target_values, criteria
        )
        discovered_candidates.extend(composite_candidates)

        # Filter by discovery criteria
        filtered_candidates = self._filter_discovery_candidates(discovered_candidates, criteria)

        logger.info(f"Discovery phase completed: {len(filtered_candidates)} candidates")
        return filtered_candidates

    def _discover_omics_biomarkers(self,
                                 omics_data: pd.DataFrame,
                                 target_values: pd.Series,
                                 omics_type: str,
                                 criteria: Dict[str, Any]) -> List[BiomarkerCandidate]:
        """Discover biomarkers from single omics layer."""

        candidates = []

        # Align data
        common_samples = list(set(omics_data.columns) & set(target_values.index))
        if len(common_samples) < 10:
            logger.warning(f"Insufficient overlapping samples for {omics_type}: {len(common_samples)}")
            return candidates

        omics_aligned = omics_data[common_samples]
        target_aligned = target_values[common_samples]

        # Statistical discovery methods
        discovery_methods = {
            'correlation': self._correlation_discovery,
            'differential': self._differential_discovery,
            'feature_selection': self._feature_selection_discovery,
            'machine_learning': self._ml_feature_discovery
        }

        for method_name, method_func in discovery_methods.items():
            try:
                method_candidates = method_func(omics_aligned, target_aligned, omics_type, method_name)
                candidates.extend(method_candidates)
            except Exception as e:
                logger.warning(f"Discovery method {method_name} failed for {omics_type}: {e}")
                continue

        return candidates

    def _correlation_discovery(self,
                             omics_data: pd.DataFrame,
                             target_values: pd.Series,
                             omics_type: str,
                             method_name: str) -> List[BiomarkerCandidate]:
        """Discover biomarkers using correlation analysis."""

        candidates = []

        for feature in omics_data.index:
            feature_values = omics_data.loc[feature]

            # Pearson correlation
            pearson_r, pearson_p = pearsonr(feature_values, target_values)

            # Spearman correlation
            spearman_r, spearman_p = spearmanr(feature_values, target_values)

            # Use stronger correlation
            if abs(pearson_r) > abs(spearman_r):
                correlation, p_value = pearson_r, pearson_p
            else:
                correlation, p_value = spearman_r, spearman_p

            # Create candidate if significant
            if p_value < 0.05 and abs(correlation) > 0.3:
                candidate = BiomarkerCandidate(
                    biomarker_id=f"{omics_type}_{feature}",
                    biomarker_name=str(feature),
                    biomarker_type=omics_type,
                    omics_source=omics_type,
                    discovery_method=f"{method_name}_correlation",
                    statistical_significance=p_value,
                    effect_size=abs(correlation),
                    validation_score=0.0,  # To be filled during validation
                    clinical_relevance=0.0,  # To be assessed
                    regulatory_qualification='exploratory',
                    reproducibility=0.0  # To be assessed
                )
                candidates.append(candidate)

        return candidates

    def _differential_discovery(self,
                              omics_data: pd.DataFrame,
                              target_values: pd.Series,
                              omics_type: str,
                              method_name: str) -> List[BiomarkerCandidate]:
        """Discover biomarkers using differential analysis."""

        candidates = []

        # Convert continuous target to binary if needed
        if target_values.dtype in ['float64', 'float32']:
            # Use median split
            median_val = target_values.median()
            binary_target = target_values > median_val
        else:
            binary_target = target_values.astype(bool)

        group1_samples = target_values[binary_target].index
        group2_samples = target_values[~binary_target].index

        for feature in omics_data.index:
            try:
                group1_values = omics_data.loc[feature, group1_samples]
                group2_values = omics_data.loc[feature, group2_samples]

                # Mann-Whitney U test (non-parametric)
                statistic, p_value = mannwhitneyu(group1_values, group2_values, alternative='two-sided')

                # Effect size (Cohen's d approximation)
                pooled_std = np.sqrt((group1_values.var() + group2_values.var()) / 2)
                cohens_d = abs(group1_values.mean() - group2_values.mean()) / pooled_std if pooled_std > 0 else 0

                # Create candidate if significant
                if p_value < 0.05 and cohens_d > 0.5:
                    candidate = BiomarkerCandidate(
                        biomarker_id=f"{omics_type}_{feature}",
                        biomarker_name=str(feature),
                        biomarker_type=omics_type,
                        omics_source=omics_type,
                        discovery_method=f"{method_name}_differential",
                        statistical_significance=p_value,
                        effect_size=cohens_d,
                        validation_score=0.0,
                        clinical_relevance=0.0,
                        regulatory_qualification='exploratory',
                        reproducibility=0.0
                    )
                    candidates.append(candidate)

            except Exception:
                continue

        return candidates

    def _feature_selection_discovery(self,
                                   omics_data: pd.DataFrame,
                                   target_values: pd.Series,
                                   omics_type: str,
                                   method_name: str) -> List[BiomarkerCandidate]:
        """Discover biomarkers using feature selection methods."""

        candidates = []

        try:
            # Prepare data
            X = omics_data.T  # Samples x features
            y = target_values

            # Feature selection methods
            feature_selectors = {
                'f_regression': SelectKBest(score_func=f_regression, k=min(50, X.shape[1])),
                'mutual_info': SelectKBest(score_func=mutual_info_regression, k=min(50, X.shape[1]))
            }

            for selector_name, selector in feature_selectors.items():
                try:
                    # Fit selector
                    X_selected = selector.fit_transform(X, y)
                    selected_features = X.columns[selector.get_support()]
                    feature_scores = selector.scores_[selector.get_support()]

                    # Create candidates for selected features
                    for feature, score in zip(selected_features, feature_scores):
                        # Normalize score to 0-1 range
                        normalized_score = min(score / np.max(feature_scores), 1.0) if np.max(feature_scores) > 0 else 0

                        candidate = BiomarkerCandidate(
                            biomarker_id=f"{omics_type}_{feature}",
                            biomarker_name=str(feature),
                            biomarker_type=omics_type,
                            omics_source=omics_type,
                            discovery_method=f"{method_name}_{selector_name}",
                            statistical_significance=0.05,  # Placeholder
                            effect_size=normalized_score,
                            validation_score=0.0,
                            clinical_relevance=0.0,
                            regulatory_qualification='exploratory',
                            reproducibility=0.0
                        )
                        candidates.append(candidate)

                except Exception as e:
                    logger.warning(f"Feature selector {selector_name} failed: {e}")
                    continue

        except Exception as e:
            logger.warning(f"Feature selection discovery failed: {e}")

        return candidates

    def _ml_feature_discovery(self,
                            omics_data: pd.DataFrame,
                            target_values: pd.Series,
                            omics_type: str,
                            method_name: str) -> List[BiomarkerCandidate]:
        """Discover biomarkers using machine learning feature importance."""

        candidates = []

        try:
            # Prepare data
            X = omics_data.T  # Samples x features
            y = target_values

            # Scale features
            scaler = StandardScaler()
            X_scaled = scaler.fit_transform(X)
            X_scaled = pd.DataFrame(X_scaled, index=X.index, columns=X.columns)

            # ML models for feature importance
            ml_models = {
                'random_forest': RandomForestRegressor(n_estimators=100, random_state=42),
                'gradient_boosting': GradientBoostingRegressor(random_state=42),
                'elastic_net': ElasticNet(random_state=42)
            }

            for model_name, model in ml_models.items():
                try:
                    # Fit model
                    model.fit(X_scaled, y)

                    # Get feature importance
                    if hasattr(model, 'feature_importances_'):
                        importance_scores = model.feature_importances_
                    elif hasattr(model, 'coef_'):
                        importance_scores = np.abs(model.coef_)
                    else:
                        continue

                    # Create candidates for top features
                    top_indices = np.argsort(importance_scores)[-20:]  # Top 20 features

                    for idx in top_indices:
                        feature = X.columns[idx]
                        importance = importance_scores[idx]

                        # Normalize importance
                        normalized_importance = importance / np.max(importance_scores) if np.max(importance_scores) > 0 else 0

                        if normalized_importance > 0.1:  # Threshold for inclusion
                            candidate = BiomarkerCandidate(
                                biomarker_id=f"{omics_type}_{feature}",
                                biomarker_name=str(feature),
                                biomarker_type=omics_type,
                                omics_source=omics_type,
                                discovery_method=f"{method_name}_{model_name}",
                                statistical_significance=0.05,  # Placeholder
                                effect_size=normalized_importance,
                                validation_score=0.0,
                                clinical_relevance=0.0,
                                regulatory_qualification='exploratory',
                                reproducibility=0.0
                            )
                            candidates.append(candidate)

                except Exception as e:
                    logger.warning(f"ML model {model_name} failed: {e}")
                    continue

        except Exception as e:
            logger.warning(f"ML feature discovery failed: {e}")

        return candidates

    def _discover_clinical_biomarkers(self,
                                    clinical_data: pd.DataFrame,
                                    target_values: pd.Series,
                                    criteria: Dict[str, Any]) -> List[BiomarkerCandidate]:
        """Discover clinical biomarkers from routine clinical measurements."""

        candidates = []

        # Find common samples
        common_samples = list(set(clinical_data.index) & set(target_values.index))
        if len(common_samples) < 10:
            return candidates

        clinical_aligned = clinical_data.loc[common_samples]
        target_aligned = target_values[common_samples]

        # Test each clinical variable
        for clinical_var in clinical_aligned.columns:
            if clinical_aligned[clinical_var].dtype in ['float64', 'float32', 'int64', 'int32']:
                try:
                    # Correlation analysis
                    correlation, p_value = pearsonr(clinical_aligned[clinical_var], target_aligned)

                    if p_value < 0.05 and abs(correlation) > 0.2:
                        candidate = BiomarkerCandidate(
                            biomarker_id=f"clinical_{clinical_var}",
                            biomarker_name=clinical_var,
                            biomarker_type='clinical',
                            omics_source='clinical',
                            discovery_method='clinical_correlation',
                            statistical_significance=p_value,
                            effect_size=abs(correlation),
                            validation_score=0.0,
                            clinical_relevance=0.8,  # High for clinical measures
                            regulatory_qualification='probable',  # Clinical measures often validated
                            reproducibility=0.0
                        )
                        candidates.append(candidate)

                except Exception:
                    continue

        return candidates

    def _discover_composite_biomarkers(self,
                                     individual_candidates: List[BiomarkerCandidate],
                                     target_values: pd.Series,
                                     criteria: Dict[str, Any]) -> List[BiomarkerCandidate]:
        """Discover composite biomarkers from combinations of individual biomarkers."""

        composite_candidates = []

        # Select top candidates for combination
        top_candidates = sorted(individual_candidates, key=lambda x: x.effect_size, reverse=True)[:10]

        if len(top_candidates) < 2:
            return composite_candidates

        try:
            # Prepare data matrix
            feature_data = []
            feature_names = []

            for candidate in top_candidates:
                # Extract candidate data based on type
                if candidate.omics_source in self.omics_data:
                    omics_df = self.omics_data[candidate.omics_source]
                    if candidate.biomarker_name in omics_df.index:
                        common_samples = list(set(omics_df.columns) & set(target_values.index))
                        candidate_values = omics_df.loc[candidate.biomarker_name, common_samples]
                        feature_data.append(candidate_values.values)
                        feature_names.append(candidate.biomarker_id)
                elif candidate.omics_source == 'clinical' and self.clinical_data is not None:
                    if candidate.biomarker_name in self.clinical_data.columns:
                        common_samples = list(set(self.clinical_data.index) & set(target_values.index))
                        candidate_values = self.clinical_data.loc[common_samples, candidate.biomarker_name]
                        feature_data.append(candidate_values.values)
                        feature_names.append(candidate.biomarker_id)

            if len(feature_data) >= 2:
                # Create feature matrix
                X = np.column_stack(feature_data)
                y = target_values[common_samples].values

                # Train composite biomarker model
                composite_model = ElasticNet(alpha=0.1, random_state=42)
                composite_model.fit(X, y)

                # Calculate composite score
                composite_scores = composite_model.predict(X)
                composite_correlation, composite_p = pearsonr(composite_scores, y)

                if composite_p < 0.05 and abs(composite_correlation) > 0.4:
                    # Create composite biomarker
                    composite_candidate = BiomarkerCandidate(
                        biomarker_id=f"composite_{'_'.join(feature_names[:3])}",
                        biomarker_name=f"Composite_{len(feature_names)}_biomarkers",
                        biomarker_type='composite',
                        omics_source='multi_omics',
                        discovery_method='elastic_net_combination',
                        statistical_significance=composite_p,
                        effect_size=abs(composite_correlation),
                        validation_score=0.0,
                        clinical_relevance=0.6,
                        regulatory_qualification='exploratory',
                        reproducibility=0.0
                    )
                    composite_candidates.append(composite_candidate)

        except Exception as e:
            logger.warning(f"Composite biomarker discovery failed: {e}")

        return composite_candidates

    def _filter_discovery_candidates(self,
                                   candidates: List[BiomarkerCandidate],
                                   criteria: Dict[str, Any]) -> List[BiomarkerCandidate]:
        """Filter discovery candidates based on criteria."""

        filtered = []

        # Remove duplicates (same biomarker from different methods)
        unique_biomarkers = {}
        for candidate in candidates:
            if candidate.biomarker_id not in unique_biomarkers:
                unique_biomarkers[candidate.biomarker_id] = candidate
            else:
                # Keep the one with better effect size
                if candidate.effect_size > unique_biomarkers[candidate.biomarker_id].effect_size:
                    unique_biomarkers[candidate.biomarker_id] = candidate

        # Apply filtering criteria
        min_effect_size = criteria.get('min_effect_size', 0.2)
        max_p_value = criteria.get('max_p_value', 0.05)

        for candidate in unique_biomarkers.values():
            if (candidate.effect_size >= min_effect_size and
                candidate.statistical_significance <= max_p_value):
                filtered.append(candidate)

        # Sort by effect size
        filtered.sort(key=lambda x: x.effect_size, reverse=True)

        # Limit number of candidates
        max_candidates = criteria.get('max_candidates', 100)
        return filtered[:max_candidates]

    def _validate_biomarker_candidates(self,
                                     candidates: List[BiomarkerCandidate],
                                     target_outcome: str) -> List[BiomarkerCandidate]:
        """Validate biomarker candidates using cross-validation and resampling."""

        validated_candidates = []

        for candidate in candidates:
            try:
                validation_score = self._perform_biomarker_validation(candidate, target_outcome)

                # Update candidate with validation score
                candidate.validation_score = validation_score

                # Assess clinical relevance
                clinical_relevance = self._assess_clinical_relevance(candidate)
                candidate.clinical_relevance = clinical_relevance

                # Assess reproducibility
                reproducibility = self._assess_reproducibility(candidate, target_outcome)
                candidate.reproducibility = reproducibility

                # Update regulatory qualification
                candidate.regulatory_qualification = self._determine_regulatory_qualification(candidate)

                # Include if passes validation threshold
                if validation_score > 0.6:
                    validated_candidates.append(candidate)

            except Exception as e:
                logger.warning(f"Validation failed for {candidate.biomarker_id}: {e}")
                continue

        logger.info(f"Validated {len(validated_candidates)} of {len(candidates)} candidates")
        return validated_candidates

    def _perform_biomarker_validation(self,
                                    candidate: BiomarkerCandidate,
                                    target_outcome: str) -> float:
        """Perform cross-validation for individual biomarker."""

        try:
            # Get biomarker data
            if candidate.omics_source in self.omics_data:
                omics_df = self.omics_data[candidate.omics_source]
                if candidate.biomarker_name in omics_df.index:
                    biomarker_data = omics_df.loc[candidate.biomarker_name]
                else:
                    return 0.0
            elif candidate.omics_source == 'clinical' and self.clinical_data is not None:
                if candidate.biomarker_name in self.clinical_data.columns:
                    biomarker_data = self.clinical_data[candidate.biomarker_name]
                else:
                    return 0.0
            else:
                return 0.0

            # Get target data
            target_data = self.outcomes_data[target_outcome]

            # Find common samples
            common_samples = list(set(biomarker_data.index) & set(target_data.index))
            if len(common_samples) < 10:
                return 0.0

            # Align data
            X = biomarker_data[common_samples].values.reshape(-1, 1)
            y = target_data[common_samples].values

            # Cross-validation
            cv_scores = cross_val_score(
                RandomForestRegressor(n_estimators=50, random_state=42),
                X, y, cv=5, scoring='r2'
            )

            # Return mean CV score (clipped to 0-1)
            return max(0.0, np.mean(cv_scores))

        except Exception as e:
            logger.warning(f"Validation failed for {candidate.biomarker_id}: {e}")
            return 0.0

    def _assess_clinical_relevance(self, candidate: BiomarkerCandidate) -> float:
        """Assess clinical relevance of biomarker."""

        # Scoring factors
        score = 0.0

        # Type-based scoring
        if candidate.biomarker_type == 'clinical':
            score += 0.4  # Clinical measures are inherently relevant
        elif candidate.biomarker_type == 'protein':
            score += 0.3  # Proteins are clinically actionable
        elif candidate.biomarker_type == 'composite':
            score += 0.2  # Composite may be more robust
        else:
            score += 0.1  # Other omics types

        # Effect size contribution
        score += min(candidate.effect_size, 0.3)

        # Statistical significance contribution
        score += max(0, 0.3 - candidate.statistical_significance * 10)

        return min(score, 1.0)

    def _assess_reproducibility(self,
                              candidate: BiomarkerCandidate,
                              target_outcome: str) -> float:
        """Assess biomarker reproducibility using bootstrap resampling."""

        try:
            # Get biomarker and target data
            if candidate.omics_source in self.omics_data:
                omics_df = self.omics_data[candidate.omics_source]
                biomarker_data = omics_df.loc[candidate.biomarker_name]
            elif candidate.omics_source == 'clinical':
                biomarker_data = self.clinical_data[candidate.biomarker_name]
            else:
                return 0.0

            target_data = self.outcomes_data[target_outcome]
            common_samples = list(set(biomarker_data.index) & set(target_data.index))

            if len(common_samples) < 20:
                return 0.0

            # Bootstrap resampling
            n_bootstrap = 100
            correlations = []

            for _ in range(n_bootstrap):
                # Sample with replacement
                bootstrap_indices = np.random.choice(common_samples, len(common_samples), replace=True)

                bootstrap_biomarker = biomarker_data[bootstrap_indices]
                bootstrap_target = target_data[bootstrap_indices]

                # Calculate correlation
                corr, p_val = pearsonr(bootstrap_biomarker, bootstrap_target)
                if not np.isnan(corr):
                    correlations.append(corr)

            if len(correlations) > 0:
                # Reproducibility based on stability of correlation
                correlation_std = np.std(correlations)
                reproducibility = max(0, 1.0 - correlation_std)
                return min(reproducibility, 1.0)
            else:
                return 0.0

        except Exception:
            return 0.0

    def _determine_regulatory_qualification(self, candidate: BiomarkerCandidate) -> str:
        """Determine regulatory qualification level for biomarker."""

        # Qualification criteria
        if (candidate.validation_score > 0.8 and
            candidate.clinical_relevance > 0.7 and
            candidate.reproducibility > 0.7):
            return 'known'
        elif (candidate.validation_score > 0.6 and
              candidate.clinical_relevance > 0.5):
            return 'probable'
        else:
            return 'exploratory'

    def _construct_biomarker_panels(self,
                                  validated_biomarkers: List[BiomarkerCandidate],
                                  target_outcome: str) -> List[BiomarkerPanel]:
        """Construct optimized biomarker panels."""

        panels = []

        if len(validated_biomarkers) < 2:
            return panels

        try:
            # Group biomarkers by type
            biomarker_groups = {}
            for biomarker in validated_biomarkers:
                group = biomarker.biomarker_type
                if group not in biomarker_groups:
                    biomarker_groups[group] = []
                biomarker_groups[group].append(biomarker)

            # Create panels
            panel_configs = [
                {'name': 'top_performers', 'biomarkers': validated_biomarkers[:5]},
                {'name': 'multi_omics', 'biomarkers': self._select_multi_omics_panel(validated_biomarkers)},
                {'name': 'clinical_focused', 'biomarkers': [b for b in validated_biomarkers if b.biomarker_type == 'clinical'][:3]}
            ]

            for config in panel_configs:
                if len(config['biomarkers']) >= 2:
                    panel_performance = self._evaluate_panel_performance(config['biomarkers'], target_outcome)

                    panel = BiomarkerPanel(
                        panel_id=f"{target_outcome}_{config['name']}_panel",
                        indication=target_outcome,
                        biomarkers=config['biomarkers'],
                        combined_performance=panel_performance,
                        validation_studies=[],  # Would be populated with actual studies
                        regulatory_status='development',
                        clinical_utility=self._assess_panel_clinical_utility(config['biomarkers'])
                    )
                    panels.append(panel)

        except Exception as e:
            logger.warning(f"Panel construction failed: {e}")

        return panels

    def _select_multi_omics_panel(self, biomarkers: List[BiomarkerCandidate]) -> List[BiomarkerCandidate]:
        """Select diverse biomarkers for multi-omics panel."""

        # Group by omics type
        omics_groups = {}
        for biomarker in biomarkers:
            omics_type = biomarker.omics_source
            if omics_type not in omics_groups:
                omics_groups[omics_type] = []
            omics_groups[omics_type].append(biomarker)

        # Select top biomarker from each omics type
        multi_omics_panel = []
        for omics_type, group_biomarkers in omics_groups.items():
            # Sort by validation score
            sorted_biomarkers = sorted(group_biomarkers, key=lambda x: x.validation_score, reverse=True)
            multi_omics_panel.append(sorted_biomarkers[0])

        return multi_omics_panel[:5]  # Limit to 5 biomarkers

    def _evaluate_panel_performance(self,
                                   panel_biomarkers: List[BiomarkerCandidate],
                                   target_outcome: str) -> Dict[str, float]:
        """Evaluate combined performance of biomarker panel."""

        try:
            # Collect biomarker data
            panel_data = []
            feature_names = []

            for biomarker in panel_biomarkers:
                if biomarker.omics_source in self.omics_data:
                    omics_df = self.omics_data[biomarker.omics_source]
                    if biomarker.biomarker_name in omics_df.index:
                        panel_data.append(omics_df.loc[biomarker.biomarker_name])
                        feature_names.append(biomarker.biomarker_id)
                elif biomarker.omics_source == 'clinical':
                    if biomarker.biomarker_name in self.clinical_data.columns:
                        panel_data.append(self.clinical_data[biomarker.biomarker_name])
                        feature_names.append(biomarker.biomarker_id)

            if len(panel_data) == 0:
                return {'combined_r2': 0.0, 'individual_mean': 0.0}

            # Create feature matrix
            common_samples = set(panel_data[0].index)
            for data in panel_data[1:]:
                common_samples &= set(data.index)
            common_samples &= set(self.outcomes_data.index)
            common_samples = list(common_samples)

            if len(common_samples) < 10:
                return {'combined_r2': 0.0, 'individual_mean': 0.0}

            X = np.column_stack([data[common_samples].values for data in panel_data])
            y = self.outcomes_data.loc[common_samples, target_outcome].values

            # Evaluate combined performance
            model = RandomForestRegressor(n_estimators=100, random_state=42)
            cv_scores = cross_val_score(model, X, y, cv=5, scoring='r2')
            combined_r2 = max(0.0, np.mean(cv_scores))

            # Individual biomarker mean performance
            individual_scores = [biomarker.validation_score for biomarker in panel_biomarkers]
            individual_mean = np.mean(individual_scores)

            return {
                'combined_r2': combined_r2,
                'individual_mean': individual_mean,
                'improvement': max(0, combined_r2 - individual_mean),
                'n_biomarkers': len(panel_biomarkers)
            }

        except Exception as e:
            logger.warning(f"Panel performance evaluation failed: {e}")
            return {'combined_r2': 0.0, 'individual_mean': 0.0}

    def _assess_panel_clinical_utility(self, panel_biomarkers: List[BiomarkerCandidate]) -> Dict[str, float]:
        """Assess clinical utility of biomarker panel."""

        return {
            'ease_of_measurement': np.mean([b.clinical_relevance for b in panel_biomarkers]),
            'cost_effectiveness': 0.7,  # Placeholder
            'turnaround_time': 0.8,  # Placeholder
            'regulatory_readiness': np.mean([1.0 if b.regulatory_qualification == 'known' else
                                           0.7 if b.regulatory_qualification == 'probable' else 0.3
                                           for b in panel_biomarkers])
        }

    def _assess_regulatory_readiness(self,
                                   biomarkers: List[BiomarkerCandidate],
                                   panels: List[BiomarkerPanel]) -> Dict[str, str]:
        """Assess regulatory readiness of biomarkers and panels."""

        # Count biomarkers by qualification level
        known_count = sum(1 for b in biomarkers if b.regulatory_qualification == 'known')
        probable_count = sum(1 for b in biomarkers if b.regulatory_qualification == 'probable')
        exploratory_count = sum(1 for b in biomarkers if b.regulatory_qualification == 'exploratory')

        # Overall readiness assessment
        if known_count >= 3:
            overall_readiness = 'submission_ready'
        elif known_count + probable_count >= 5:
            overall_readiness = 'validation_ready'
        elif exploratory_count >= 10:
            overall_readiness = 'discovery_complete'
        else:
            overall_readiness = 'early_discovery'

        return {
            'overall_readiness': overall_readiness,
            'known_biomarkers': known_count,
            'probable_biomarkers': probable_count,
            'exploratory_biomarkers': exploratory_count,
            'recommended_path': self._recommend_regulatory_path(overall_readiness)
        }

    def _recommend_regulatory_path(self, readiness_level: str) -> str:
        """Recommend regulatory development path."""

        paths = {
            'submission_ready': 'FDA biomarker qualification submission',
            'validation_ready': 'Conduct validation studies for qualification',
            'discovery_complete': 'Initiate analytical validation studies',
            'early_discovery': 'Continue discovery and replication studies'
        }

        return paths.get(readiness_level, 'Continue research and development')

    def _assess_longitudinal_stability(self, biomarkers: List[BiomarkerCandidate]) -> Dict[str, float]:
        """Assess longitudinal stability of biomarkers."""

        stability_scores = {}

        if not self.longitudinal_data:
            return stability_scores

        for biomarker in biomarkers:
            try:
                # Collect longitudinal data for biomarker
                timepoint_correlations = []

                for timepoint, timepoint_data in self.longitudinal_data.items():
                    if biomarker.omics_source in timepoint_data:
                        omics_df = timepoint_data[biomarker.omics_source]
                        if biomarker.biomarker_name in omics_df.index:
                            # Compare with baseline if available
                            if 'baseline' in self.longitudinal_data:
                                baseline_data = self.longitudinal_data['baseline'][biomarker.omics_source]
                                if biomarker.biomarker_name in baseline_data.index:
                                    common_samples = list(set(omics_df.columns) & set(baseline_data.columns))
                                    if len(common_samples) > 5:
                                        corr, _ = pearsonr(
                                            omics_df.loc[biomarker.biomarker_name, common_samples],
                                            baseline_data.loc[biomarker.biomarker_name, common_samples]
                                        )
                                        if not np.isnan(corr):
                                            timepoint_correlations.append(abs(corr))

                # Calculate stability score
                if timepoint_correlations:
                    stability_score = np.mean(timepoint_correlations)
                    stability_scores[biomarker.biomarker_id] = stability_score

            except Exception:
                continue

        return stability_scores

    def _compile_discovery_statistics(self, candidates: List[BiomarkerCandidate]) -> Dict[str, Any]:
        """Compile discovery phase statistics."""

        # Group by omics type
        omics_counts = {}
        for candidate in candidates:
            omics_type = candidate.omics_source
            omics_counts[omics_type] = omics_counts.get(omics_type, 0) + 1

        # Group by discovery method
        method_counts = {}
        for candidate in candidates:
            method = candidate.discovery_method
            method_counts[method] = method_counts.get(method, 0) + 1

        return {
            'total_candidates': len(candidates),
            'omics_distribution': omics_counts,
            'method_distribution': method_counts,
            'effect_size_stats': {
                'mean': np.mean([c.effect_size for c in candidates]),
                'median': np.median([c.effect_size for c in candidates]),
                'max': np.max([c.effect_size for c in candidates]) if candidates else 0
            }
        }

    def _compile_validation_metrics(self, validated_biomarkers: List[BiomarkerCandidate]) -> Dict[str, float]:
        """Compile validation phase metrics."""

        if not validated_biomarkers:
            return {}

        return {
            'mean_validation_score': np.mean([b.validation_score for b in validated_biomarkers]),
            'mean_clinical_relevance': np.mean([b.clinical_relevance for b in validated_biomarkers]),
            'mean_reproducibility': np.mean([b.reproducibility for b in validated_biomarkers]),
            'high_confidence_count': sum(1 for b in validated_biomarkers
                                       if b.validation_score > 0.8 and b.reproducibility > 0.7)
        }

    def _get_default_discovery_criteria(self) -> Dict[str, Any]:
        """Get default discovery criteria."""

        return {
            'min_effect_size': 0.2,
            'max_p_value': 0.05,
            'max_candidates': 100,
            'min_samples': 10,
            'cross_validation_folds': 5
        }


def create_biomarker_discovery_demo():
    """Create demonstration of biomarker discovery system."""

    # Generate synthetic study data
    n_patients = 150
    patient_ids = [f"patient_{i:03d}" for i in range(n_patients)]

    np.random.seed(42)

    # Clinical data
    clinical_data = pd.DataFrame({
        'age': np.random.normal(65, 12, n_patients),
        'bmi': np.random.normal(25, 4, n_patients),
        'systolic_bp': np.random.normal(130, 20, n_patients),
        'cholesterol': np.random.normal(200, 40, n_patients),
        'inflammation_marker': np.random.lognormal(1, 0.8, n_patients)
    }, index=patient_ids)

    # Outcomes data
    # Create synthetic aging outcome with realistic correlations
    aging_score = (
        0.3 * (clinical_data['age'] - 65) / 12 +
        0.2 * (clinical_data['bmi'] - 25) / 4 +
        0.2 * (clinical_data['systolic_bp'] - 130) / 20 +
        0.3 * np.log(clinical_data['inflammation_marker']) +
        np.random.normal(0, 0.5, n_patients)
    )

    outcomes_data = pd.DataFrame({
        'aging_score': aging_score,
        'intervention_response': np.random.choice([0, 1], n_patients, p=[0.6, 0.4])
    }, index=patient_ids)

    # Multi-omics data
    omics_data = {
        'transcriptomics': pd.DataFrame(
            np.random.normal(0, 1, (500, n_patients)),
            columns=patient_ids,
            index=[f"gene_{i}" for i in range(500)]
        ),
        'proteomics': pd.DataFrame(
            np.random.lognormal(5, 1, (200, n_patients)),
            columns=patient_ids,
            index=[f"protein_{i}" for i in range(200)]
        ),
        'metabolomics': pd.DataFrame(
            np.random.lognormal(3, 0.8, (100, n_patients)),
            columns=patient_ids,
            index=[f"metabolite_{i}" for i in range(100)]
        )
    }

    # Add some realistic correlations to aging score
    for i in range(20):  # Top 20 features correlated with aging
        correlation_strength = np.random.uniform(0.3, 0.7)
        noise = np.random.normal(0, 0.3, n_patients)
        omics_data['transcriptomics'].iloc[i] = correlation_strength * aging_score + noise

    # Initialize biomarker discovery engine
    biomarker_engine = RealTimeBiomarkerEngine(discovery_mode="comprehensive")

    # Load study data
    biomarker_engine.load_study_data(
        clinical_data=clinical_data,
        outcomes_data=outcomes_data,
        omics_data=omics_data
    )

    # Perform biomarker discovery
    discovery_result = biomarker_engine.discover_biomarkers(
        target_outcome='aging_score'
    )

    # Display results
    print("BIOMARKER DISCOVERY RESULTS")
    print("=" * 50)
    print(f"Total biomarkers discovered: {len(discovery_result.discovered_biomarkers)}")
    print(f"Biomarker panels created: {len(discovery_result.validated_panels)}")
    print(f"Regulatory readiness: {discovery_result.regulatory_readiness['overall_readiness']}")

    print("\nTOP VALIDATED BIOMARKERS:")
    top_biomarkers = sorted(discovery_result.discovered_biomarkers,
                          key=lambda x: x.validation_score, reverse=True)[:10]

    for i, biomarker in enumerate(top_biomarkers, 1):
        print(f"{i:2d}. {biomarker.biomarker_name} ({biomarker.biomarker_type})")
        print(f"    Validation score: {biomarker.validation_score:.3f}")
        print(f"    Effect size: {biomarker.effect_size:.3f}")
        print(f"    Clinical relevance: {biomarker.clinical_relevance:.3f}")
        print(f"    Regulatory status: {biomarker.regulatory_qualification}")

    print("\nBIOMARKER PANELS:")
    for panel in discovery_result.validated_panels:
        print(f"Panel: {panel.panel_id}")
        print(f"  Biomarkers: {len(panel.biomarkers)}")
        print(f"  Combined R²: {panel.combined_performance.get('combined_r2', 0):.3f}")
        print(f"  Regulatory readiness: {panel.clinical_utility.get('regulatory_readiness', 0):.3f}")

    print("\nDISCOVERY STATISTICS:")
    stats = discovery_result.discovery_statistics
    print(f"  Total candidates: {stats['total_candidates']}")
    print(f"  Omics distribution: {stats['omics_distribution']}")
    print(f"  Mean effect size: {stats['effect_size_stats']['mean']:.3f}")

    print("\nVALIDATION METRICS:")
    metrics = discovery_result.validation_metrics
    for metric, value in metrics.items():
        print(f"  {metric}: {value:.3f}")

    return biomarker_engine, discovery_result


if __name__ == "__main__":
    # Run demonstration
    engine, results = create_biomarker_discovery_demo()
    print("\nReal-time biomarker discovery completed!")
    print("\nReal-time biomarker discovery completed!")

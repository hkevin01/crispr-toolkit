"""
Automated Adverse Event Detection for CRISPR Toolkit Phase 3
===========================================================

Advanced automated adverse event (AE) detection and safety monitoring
system for aging intervention clinical trials with real-time signal
detection, causality assessment, and regulatory reporting.
"""

import logging
from dataclasses import dataclass
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import IsolationForest

logger = logging.getLogger(__name__)


@dataclass
class AdverseEvent:
    """Individual adverse event record."""
    event_id: str
    patient_id: str
    event_term: str
    onset_date: datetime
    severity: str  # "Mild", "Moderate", "Severe", "Life-threatening"
    seriousness: str  # "Serious", "Non-serious"
    causality: str  # "Unrelated", "Unlikely", "Possible", "Probable", "Definite"
    outcome: str  # "Recovered", "Recovering", "Not recovered", "Fatal"
    action_taken: str  # "None", "Dose reduction", "Drug discontinued", etc.
    narrative: str
    reporter: str
    report_date: datetime


@dataclass
class SafetySignal:
    """Detected safety signal."""
    signal_id: str
    event_term: str
    detection_method: str
    signal_strength: float
    statistical_threshold: float
    affected_patients: int
    total_exposure: int
    confidence_level: float
    detection_date: datetime
    alert_level: str  # "Low", "Medium", "High", "Critical"
    recommended_actions: List[str]


@dataclass
class SafetyAssessment:
    """Overall safety assessment result."""
    assessment_id: str
    assessment_date: datetime
    total_patients: int
    total_aes: int
    serious_aes: int
    unexpected_aes: int
    drug_related_aes: int
    detected_signals: List[SafetySignal]
    risk_benefit_ratio: float
    overall_safety_status: str
    recommendations: List[str]


class AdverseEventDetector:
    """
    Automated adverse event detection and safety monitoring system.

    This class implements sophisticated algorithms for real-time
    detection of safety signals in clinical trial data.
    """

    def __init__(self, detection_methods: Optional[List[str]] = None):
        """
        Initialize adverse event detector.

        Args:
            detection_methods: List of detection methods to use
        """
        if detection_methods is None:
            detection_methods = [
                "disproportionality",
                "temporal_clustering",
                "dose_response",
                "anomaly_detection",
                "bayesian_monitoring"
            ]

        self.detection_methods = detection_methods
        self.ae_data = None
        self.exposure_data = None
        self.historical_data = None
        self.signal_history = []

        # Detection thresholds
        self.thresholds = {
            "disproportionality_ror": 2.0,
            "temporal_clustering_pvalue": 0.01,
            "dose_response_correlation": 0.7,
            "anomaly_score": 0.1,
            "bayesian_posterior": 0.95
        }

        logger.info(f"Initialized AdverseEventDetector with methods: {detection_methods}")

    def load_safety_data(self,
                        ae_data: pd.DataFrame,
                        exposure_data: Optional[pd.DataFrame] = None,
                        historical_data: Optional[pd.DataFrame] = None):
        """Load adverse event and exposure data."""

        self.ae_data = ae_data
        if exposure_data is not None:
            self.exposure_data = exposure_data
        if historical_data is not None:
            self.historical_data = historical_data

        logger.info("Loaded safety data:")
        logger.info(f"  Adverse events: {ae_data.shape}")
        if exposure_data is not None:
            logger.info(f"  Exposure data: {exposure_data.shape}")
        if historical_data is not None:
            logger.info(f"  Historical data: {historical_data.shape}")

    def detect_safety_signals(self,
                             analysis_date: Optional[datetime] = None) -> List[SafetySignal]:
        """
        Detect safety signals using multiple methods.

        Args:
            analysis_date: Date of analysis (default: current date)

        Returns:
            List of detected safety signals
        """
        if analysis_date is None:
            analysis_date = datetime.now()

        logger.info(f"Starting safety signal detection at {analysis_date}")

        if self.ae_data is None:
            logger.warning("No AE data loaded")
            return []

        detected_signals = []

        # Apply each detection method
        for method in self.detection_methods:
            try:
                method_signals = self._apply_detection_method(method, analysis_date)
                detected_signals.extend(method_signals)
                logger.info(f"{method}: detected {len(method_signals)} signals")
            except Exception as e:
                logger.error(f"Error in {method} detection: {str(e)}")

        # Remove duplicates and rank signals
        unique_signals = self._consolidate_signals(detected_signals)
        ranked_signals = self._rank_signals(unique_signals)

        # Update signal history
        self.signal_history.extend(ranked_signals)

        logger.info(f"Total unique signals detected: {len(ranked_signals)}")
        return ranked_signals

    def assess_causality(self,
                        ae: AdverseEvent,
                        patient_data: Optional[Dict[str, Any]] = None) -> str:
        """
        Assess causality relationship between intervention and AE.

        Args:
            ae: Adverse event to assess
            patient_data: Additional patient data for assessment

        Returns:
            Causality assessment
        """
        logger.debug(f"Assessing causality for AE: {ae.event_id}")

        # Naranjo algorithm implementation
        score = 0

        # Temporal relationship
        if hasattr(ae, 'onset_date') and hasattr(ae, 'exposure_start'):
            days_to_onset = (ae.onset_date - ae.exposure_start).days
            if 0 <= days_to_onset <= 30:
                score += 2
            elif 31 <= days_to_onset <= 90:
                score += 1

        # Dose-response relationship
        if patient_data and 'dose_changes' in patient_data:
            dose_changes = patient_data['dose_changes']
            if self._check_dose_response_relationship(ae, dose_changes):
                score += 1

        # Dechallenge/rechallenge
        if patient_data and 'drug_discontinuation' in patient_data:
            if patient_data['drug_discontinuation']:
                score += 1

        # Known reaction pattern
        if self._check_known_reaction(ae.event_term):
            score += 1

        # Alternative causes
        if patient_data and 'concomitant_meds' in patient_data:
            if not patient_data['concomitant_meds']:
                score += 1

        # Convert score to causality assessment
        if score >= 9:
            causality = "Definite"
        elif score >= 5:
            causality = "Probable"
        elif score >= 1:
            causality = "Possible"
        else:
            causality = "Unlikely"

        return causality

    def generate_safety_report(self,
                             analysis_date: Optional[datetime] = None) -> SafetyAssessment:
        """
        Generate comprehensive safety assessment report.

        Args:
            analysis_date: Date of assessment

        Returns:
            SafetyAssessment with comprehensive safety data
        """
        if analysis_date is None:
            analysis_date = datetime.now()

        logger.info(f"Generating safety assessment for {analysis_date}")

        # Detect current safety signals
        detected_signals = self.detect_safety_signals(analysis_date)

        # Calculate safety metrics
        safety_metrics = self._calculate_safety_metrics()

        # Assess overall safety status
        safety_status = self._assess_overall_safety_status(safety_metrics, detected_signals)

        # Generate recommendations
        recommendations = self._generate_safety_recommendations(safety_metrics, detected_signals)

        assessment = SafetyAssessment(
            assessment_id=f"safety_assessment_{analysis_date.strftime('%Y%m%d_%H%M')}",
            assessment_date=analysis_date,
            total_patients=safety_metrics.get('total_patients', 0),
            total_aes=safety_metrics.get('total_aes', 0),
            serious_aes=safety_metrics.get('serious_aes', 0),
            unexpected_aes=safety_metrics.get('unexpected_aes', 0),
            drug_related_aes=safety_metrics.get('drug_related_aes', 0),
            detected_signals=detected_signals,
            risk_benefit_ratio=safety_metrics.get('risk_benefit_ratio', 1.0),
            overall_safety_status=safety_status,
            recommendations=recommendations
        )

        logger.info(f"Safety assessment completed: {safety_status}")
        return assessment

    def _apply_detection_method(self,
                              method: str,
                              analysis_date: datetime) -> List[SafetySignal]:
        """Apply specific detection method."""

        if method == "disproportionality":
            return self._disproportionality_analysis(analysis_date)
        elif method == "temporal_clustering":
            return self._temporal_clustering_analysis(analysis_date)
        elif method == "dose_response":
            return self._dose_response_analysis(analysis_date)
        elif method == "anomaly_detection":
            return self._anomaly_detection_analysis(analysis_date)
        elif method == "bayesian_monitoring":
            return self._bayesian_monitoring_analysis(analysis_date)
        else:
            logger.warning(f"Unknown detection method: {method}")
            return []

    def _disproportionality_analysis(self, analysis_date: datetime) -> List[SafetySignal]:
        """Perform disproportionality analysis using reporting odds ratio."""

        signals = []

        if self.ae_data is None:
            return signals

        # Group AEs by event term
        ae_counts = self.ae_data['event_term'].value_counts()

        # Calculate total exposure
        total_patients = len(self.ae_data['patient_id'].unique())

        for event_term, count in ae_counts.items():
            if count < 3:  # Minimum count threshold
                continue

            # Calculate reporting odds ratio (ROR)
            # Compare to historical data or expected background rate
            if self.historical_data is not None:
                background_rate = self._get_background_rate(event_term)
            else:
                background_rate = 0.01  # Default 1% background rate

            observed_rate = count / total_patients
            expected_count = background_rate * total_patients

            if expected_count > 0:
                ror = observed_rate / background_rate

                # Calculate confidence interval
                se_log_ror = np.sqrt(1/count + 1/expected_count)
                ci_lower = np.exp(np.log(ror) - 1.96 * se_log_ror)
                ci_upper = np.exp(np.log(ror) + 1.96 * se_log_ror)

                # Check if signal threshold is met
                if ror >= self.thresholds["disproportionality_ror"] and ci_lower > 1:
                    alert_level = self._determine_alert_level(ror, count)

                    signal = SafetySignal(
                        signal_id=f"disp_{event_term}_{analysis_date.strftime('%Y%m%d')}",
                        event_term=event_term,
                        detection_method="disproportionality",
                        signal_strength=ror,
                        statistical_threshold=self.thresholds["disproportionality_ror"],
                        affected_patients=count,
                        total_exposure=total_patients,
                        confidence_level=0.95,
                        detection_date=analysis_date,
                        alert_level=alert_level,
                        recommended_actions=self._get_recommended_actions(alert_level, event_term)
                    )
                    signals.append(signal)

        return signals

    def _temporal_clustering_analysis(self, analysis_date: datetime) -> List[SafetySignal]:
        """Detect temporal clustering of adverse events."""

        signals = []

        if self.ae_data is None or 'onset_date' not in self.ae_data.columns:
            return signals

        # Convert onset dates to datetime
        ae_data = self.ae_data.copy()
        ae_data['onset_date'] = pd.to_datetime(ae_data['onset_date'])

        # Group by event term
        for event_term in ae_data['event_term'].unique():
            event_data = ae_data[ae_data['event_term'] == event_term]

            if len(event_data) < 3:
                continue

            # Calculate time intervals between events
            onset_times = sorted(event_data['onset_date'])
            intervals = [(onset_times[i+1] - onset_times[i]).days
                        for i in range(len(onset_times)-1)]

            if len(intervals) > 0:
                # Test for clustering using exponential distribution
                mean_interval = np.mean(intervals)

                # Kolmogorov-Smirnov test for exponential distribution
                # If p-value < threshold, suggests clustering
                statistic, p_value = stats.kstest(intervals,
                                                lambda x: stats.expon.cdf(x, scale=mean_interval))

                if p_value < self.thresholds["temporal_clustering_pvalue"]:
                    signal_strength = 1 - p_value
                    alert_level = self._determine_alert_level(signal_strength, len(event_data))

                    signal = SafetySignal(
                        signal_id=f"temp_{event_term}_{analysis_date.strftime('%Y%m%d')}",
                        event_term=event_term,
                        detection_method="temporal_clustering",
                        signal_strength=signal_strength,
                        statistical_threshold=self.thresholds["temporal_clustering_pvalue"],
                        affected_patients=len(event_data),
                        total_exposure=len(self.ae_data['patient_id'].unique()),
                        confidence_level=0.99,
                        detection_date=analysis_date,
                        alert_level=alert_level,
                        recommended_actions=self._get_recommended_actions(alert_level, event_term)
                    )
                    signals.append(signal)

        return signals

    def _dose_response_analysis(self, analysis_date: datetime) -> List[SafetySignal]:
        """Analyze dose-response relationships for safety signals."""

        signals = []

        if (self.ae_data is None or
            'dose' not in self.ae_data.columns or
            'severity' not in self.ae_data.columns):
            return signals

        # Encode severity as numeric
        severity_mapping = {'Mild': 1, 'Moderate': 2, 'Severe': 3, 'Life-threatening': 4}
        ae_data = self.ae_data.copy()
        ae_data['severity_numeric'] = ae_data['severity'].map(severity_mapping)

        # Group by event term
        for event_term in ae_data['event_term'].unique():
            event_data = ae_data[ae_data['event_term'] == event_term]

            if len(event_data) < 5:  # Need minimum sample size
                continue

            # Calculate correlation between dose and severity
            valid_data = event_data.dropna(subset=['dose', 'severity_numeric'])

            if len(valid_data) >= 5:
                correlation, p_value = stats.pearsonr(valid_data['dose'],
                                                    valid_data['severity_numeric'])

                if (correlation >= self.thresholds["dose_response_correlation"] and
                    p_value < 0.05):

                    alert_level = self._determine_alert_level(correlation, len(valid_data))

                    signal = SafetySignal(
                        signal_id=f"dose_{event_term}_{analysis_date.strftime('%Y%m%d')}",
                        event_term=event_term,
                        detection_method="dose_response",
                        signal_strength=correlation,
                        statistical_threshold=self.thresholds["dose_response_correlation"],
                        affected_patients=len(valid_data),
                        total_exposure=len(self.ae_data['patient_id'].unique()),
                        confidence_level=0.95,
                        detection_date=analysis_date,
                        alert_level=alert_level,
                        recommended_actions=self._get_recommended_actions(alert_level, event_term)
                    )
                    signals.append(signal)

        return signals

    def _anomaly_detection_analysis(self, analysis_date: datetime) -> List[SafetySignal]:
        """Use machine learning for anomaly detection in AE patterns."""

        signals = []

        if self.ae_data is None:
            return signals

        try:
            # Prepare features for anomaly detection
            features = self._prepare_anomaly_features()

            if len(features) < 10:  # Minimum data requirement
                return signals

            # Isolation Forest for anomaly detection
            isolation_forest = IsolationForest(contamination=0.1, random_state=42)
            anomaly_scores = isolation_forest.fit_predict(features)

            # Identify anomalous patterns
            anomaly_indices = np.where(anomaly_scores == -1)[0]

            if len(anomaly_indices) > 0:
                # Analyze anomalous cases
                anomalous_data = self.ae_data.iloc[anomaly_indices]

                for event_term in anomalous_data['event_term'].unique():
                    term_anomalies = anomalous_data[anomalous_data['event_term'] == event_term]

                    if len(term_anomalies) >= 2:
                        anomaly_score = len(term_anomalies) / len(anomalous_data)

                        if anomaly_score >= self.thresholds["anomaly_score"]:
                            alert_level = self._determine_alert_level(anomaly_score, len(term_anomalies))

                            signal = SafetySignal(
                                signal_id=f"anom_{event_term}_{analysis_date.strftime('%Y%m%d')}",
                                event_term=event_term,
                                detection_method="anomaly_detection",
                                signal_strength=anomaly_score,
                                statistical_threshold=self.thresholds["anomaly_score"],
                                affected_patients=len(term_anomalies),
                                total_exposure=len(self.ae_data['patient_id'].unique()),
                                confidence_level=0.90,
                                detection_date=analysis_date,
                                alert_level=alert_level,
                                recommended_actions=self._get_recommended_actions(alert_level, event_term)
                            )
                            signals.append(signal)

        except Exception as e:
            logger.error(f"Error in anomaly detection: {str(e)}")

        return signals

    def _bayesian_monitoring_analysis(self, analysis_date: datetime) -> List[SafetySignal]:
        """Bayesian sequential monitoring for safety signals."""

        signals = []

        if self.ae_data is None:
            return signals

        # Group by event term
        ae_counts = self.ae_data['event_term'].value_counts()
        total_patients = len(self.ae_data['patient_id'].unique())

        for event_term, observed_count in ae_counts.items():
            if observed_count < 2:
                continue

            # Bayesian analysis with Beta-Binomial model
            # Prior: Beta(1, 1) - uninformative prior
            alpha_prior = 1
            beta_prior = 1

            # Posterior: Beta(alpha_prior + observed, beta_prior + n - observed)
            alpha_posterior = alpha_prior + observed_count
            beta_posterior = beta_prior + total_patients - observed_count

            # Posterior probability that rate > threshold
            threshold_rate = 0.05  # 5% threshold
            posterior_prob = 1 - stats.beta.cdf(threshold_rate, alpha_posterior, beta_posterior)

            if posterior_prob >= self.thresholds["bayesian_posterior"]:
                alert_level = self._determine_alert_level(posterior_prob, observed_count)

                signal = SafetySignal(
                    signal_id=f"bayes_{event_term}_{analysis_date.strftime('%Y%m%d')}",
                    event_term=event_term,
                    detection_method="bayesian_monitoring",
                    signal_strength=posterior_prob,
                    statistical_threshold=self.thresholds["bayesian_posterior"],
                    affected_patients=observed_count,
                    total_exposure=total_patients,
                    confidence_level=0.95,
                    detection_date=analysis_date,
                    alert_level=alert_level,
                    recommended_actions=self._get_recommended_actions(alert_level, event_term)
                )
                signals.append(signal)

        return signals

    def _prepare_anomaly_features(self) -> np.ndarray:
        """Prepare features for anomaly detection."""

        # Create feature matrix from AE data
        features = []

        # Group by patient
        for patient_id in self.ae_data['patient_id'].unique():
            patient_aes = self.ae_data[self.ae_data['patient_id'] == patient_id]

            # Feature vector for this patient
            feature_vector = [
                len(patient_aes),  # Number of AEs
                len(patient_aes['event_term'].unique()),  # Number of unique AE terms
                len(patient_aes[patient_aes['severity'].isin(['Severe', 'Life-threatening'])]),  # Severe AEs
                len(patient_aes[patient_aes['seriousness'] == 'Serious']),  # Serious AEs
            ]

            # Add onset timing features if available
            if 'onset_date' in patient_aes.columns:
                onset_dates = pd.to_datetime(patient_aes['onset_date'])
                if len(onset_dates) > 1:
                    time_span = (onset_dates.max() - onset_dates.min()).days
                    feature_vector.append(time_span)
                else:
                    feature_vector.append(0)

            features.append(feature_vector)

        return np.array(features)

    def _consolidate_signals(self, signals: List[SafetySignal]) -> List[SafetySignal]:
        """Remove duplicate signals and consolidate similar ones."""

        # Group signals by event term
        signal_groups = {}
        for signal in signals:
            if signal.event_term not in signal_groups:
                signal_groups[signal.event_term] = []
            signal_groups[signal.event_term].append(signal)

        consolidated = []

        for event_term, term_signals in signal_groups.items():
            if len(term_signals) == 1:
                consolidated.append(term_signals[0])
            else:
                # Take the strongest signal for this event term
                strongest_signal = max(term_signals, key=lambda s: s.signal_strength)
                consolidated.append(strongest_signal)

        return consolidated

    def _rank_signals(self, signals: List[SafetySignal]) -> List[SafetySignal]:
        """Rank signals by importance and alert level."""

        # Define ranking criteria
        def signal_priority(signal):
            priority_score = 0

            # Alert level weight
            alert_weights = {"Critical": 4, "High": 3, "Medium": 2, "Low": 1}
            priority_score += alert_weights.get(signal.alert_level, 0) * 10

            # Signal strength weight
            priority_score += signal.signal_strength * 5

            # Number of affected patients
            priority_score += min(signal.affected_patients, 10)

            return priority_score

        return sorted(signals, key=signal_priority, reverse=True)

    def _calculate_safety_metrics(self) -> Dict[str, Any]:
        """Calculate overall safety metrics."""

        if self.ae_data is None:
            return {}

        total_patients = len(self.ae_data['patient_id'].unique())
        total_aes = len(self.ae_data)

        # Count serious AEs
        serious_aes = len(self.ae_data[self.ae_data['seriousness'] == 'Serious'])

        # Count unexpected AEs (simplified - would need reference safety info)
        unexpected_aes = 0  # Placeholder

        # Count drug-related AEs
        drug_related = self.ae_data['causality'].isin(['Possible', 'Probable', 'Definite'])
        drug_related_aes = len(self.ae_data[drug_related])

        # Calculate rates
        ae_rate = total_aes / total_patients if total_patients > 0 else 0
        serious_ae_rate = serious_aes / total_patients if total_patients > 0 else 0

        # Risk-benefit ratio (simplified)
        risk_benefit_ratio = 1.0  # Would require efficacy data

        return {
            'total_patients': total_patients,
            'total_aes': total_aes,
            'serious_aes': serious_aes,
            'unexpected_aes': unexpected_aes,
            'drug_related_aes': drug_related_aes,
            'ae_rate': ae_rate,
            'serious_ae_rate': serious_ae_rate,
            'risk_benefit_ratio': risk_benefit_ratio
        }

    def _assess_overall_safety_status(self,
                                    metrics: Dict[str, Any],
                                    signals: List[SafetySignal]) -> str:
        """Assess overall safety status."""

        # Count critical and high-priority signals
        critical_signals = len([s for s in signals if s.alert_level == "Critical"])
        high_signals = len([s for s in signals if s.alert_level == "High"])

        # Assess based on signals and metrics
        if critical_signals > 0:
            return "Critical safety concern"
        elif high_signals > 2:
            return "High safety concern"
        elif metrics.get('serious_ae_rate', 0) > 0.1:  # >10% serious AE rate
            return "Moderate safety concern"
        elif len(signals) > 0:
            return "Low safety concern"
        else:
            return "No significant safety concerns"

    def _generate_safety_recommendations(self,
                                       metrics: Dict[str, Any],
                                       signals: List[SafetySignal]) -> List[str]:
        """Generate safety recommendations."""

        recommendations = []

        # Based on detected signals
        critical_signals = [s for s in signals if s.alert_level == "Critical"]
        high_signals = [s for s in signals if s.alert_level == "High"]

        if critical_signals:
            recommendations.append("Consider immediate trial suspension pending safety review")
            recommendations.append("Notify regulatory authorities within 24 hours")

        if high_signals:
            recommendations.append("Conduct urgent safety committee review")
            recommendations.append("Consider protocol amendments or enhanced monitoring")

        # Based on AE rates
        if metrics.get('serious_ae_rate', 0) > 0.15:
            recommendations.append("Review inclusion/exclusion criteria")
            recommendations.append("Enhanced patient monitoring recommended")

        # General recommendations
        if len(signals) > 0:
            recommendations.append("Update safety documentation")
            recommendations.append("Consider additional safety endpoints")

        if not recommendations:
            recommendations.append("Continue current safety monitoring")

        return recommendations

    def _get_background_rate(self, event_term: str) -> float:
        """Get background rate for event term from historical data."""

        if self.historical_data is None:
            return 0.01  # Default 1% background rate

        # Look up in historical data
        if 'event_term' in self.historical_data.columns:
            historical_events = self.historical_data[
                self.historical_data['event_term'] == event_term
            ]
            if len(historical_events) > 0:
                return historical_events['rate'].mean()

        return 0.01

    def _determine_alert_level(self, signal_strength: float, affected_patients: int) -> str:
        """Determine alert level based on signal strength and impact."""

        if signal_strength > 0.9 and affected_patients >= 5:
            return "Critical"
        elif signal_strength > 0.7 and affected_patients >= 3:
            return "High"
        elif signal_strength > 0.5 or affected_patients >= 2:
            return "Medium"
        else:
            return "Low"

    def _get_recommended_actions(self, alert_level: str, event_term: str) -> List[str]:
        """Get recommended actions based on alert level."""

        actions = []

        if alert_level == "Critical":
            actions.extend([
                f"Immediate review of all {event_term} cases",
                "Consider trial suspension",
                "Emergency safety committee meeting",
                "Regulatory notification required"
            ])
        elif alert_level == "High":
            actions.extend([
                f"Urgent review of {event_term} cases",
                "Safety committee consultation",
                "Consider enhanced monitoring",
                "Prepare regulatory update"
            ])
        elif alert_level == "Medium":
            actions.extend([
                f"Detailed review of {event_term} cases",
                "Document in safety update",
                "Monitor closely"
            ])
        else:  # Low
            actions.extend([
                f"Continue monitoring {event_term}",
                "Include in routine safety report"
            ])

        return actions

    def _check_dose_response_relationship(self,
                                        ae: AdverseEvent,
                                        dose_changes: List[Dict]) -> bool:
        """Check for dose-response relationship."""

        # Simplified check - would need more sophisticated analysis
        return len(dose_changes) > 0

    def _check_known_reaction(self, event_term: str) -> bool:
        """Check if event is a known reaction."""

        # Simplified - would check against reference safety information
        known_reactions = [
            "injection site reaction",
            "fatigue",
            "headache",
            "nausea"
        ]

        return any(known in event_term.lower() for known in known_reactions)


def create_ae_detection_demo():
    """Create demonstration of automated AE detection."""

    # Generate synthetic AE data
    np.random.seed(42)

    n_patients = 200
    n_aes = 150

    patient_ids = [f"patient_{i:03d}" for i in range(n_patients)]

    # Event terms with different frequencies
    event_terms = [
        "Headache", "Fatigue", "Nausea", "Dizziness", "Rash",
        "Injection site reaction", "Fever", "Muscle pain", "Insomnia",
        "Rare serious event"  # Low frequency event
    ]

    event_probabilities = [0.2, 0.15, 0.12, 0.1, 0.08, 0.15, 0.05, 0.07, 0.04, 0.04]

    # Generate AE data
    ae_data = []
    for i in range(n_aes):
        event_term = np.random.choice(event_terms, p=event_probabilities)
        patient_id = np.random.choice(patient_ids)

        # Severity based on event type
        if event_term == "Rare serious event":
            severity = np.random.choice(["Severe", "Life-threatening"], p=[0.7, 0.3])
            seriousness = "Serious"
        else:
            severity = np.random.choice(["Mild", "Moderate", "Severe"], p=[0.6, 0.3, 0.1])
            seriousness = np.random.choice(["Non-serious", "Serious"], p=[0.9, 0.1])

        # Causality assessment
        causality = np.random.choice(
            ["Unrelated", "Unlikely", "Possible", "Probable", "Definite"],
            p=[0.2, 0.2, 0.3, 0.2, 0.1]
        )

        ae_record = {
            'ae_id': f"ae_{i:03d}",
            'patient_id': patient_id,
            'event_term': event_term,
            'onset_date': datetime.now() - timedelta(days=np.random.randint(1, 90)),
            'severity': severity,
            'seriousness': seriousness,
            'causality': causality,
            'outcome': np.random.choice(["Recovered", "Recovering", "Not recovered"], p=[0.7, 0.2, 0.1]),
            'dose': np.random.uniform(10, 100)  # Dose in mg
        }
        ae_data.append(ae_record)

    ae_df = pd.DataFrame(ae_data)

    # Historical data for background rates
    historical_data = pd.DataFrame({
        'event_term': event_terms,
        'rate': [0.05, 0.03, 0.02, 0.02, 0.01, 0.08, 0.01, 0.02, 0.005, 0.001]
    })

    # Initialize AE detector
    ae_detector = AdverseEventDetector()

    # Load data
    ae_detector.load_safety_data(
        ae_data=ae_df,
        historical_data=historical_data
    )

    # Detect safety signals
    detected_signals = ae_detector.detect_safety_signals()

    # Generate safety assessment
    safety_assessment = ae_detector.generate_safety_report()

    # Display results
    print("AUTOMATED ADVERSE EVENT DETECTION RESULTS")
    print("=" * 50)
    print(f"Total patients: {safety_assessment.total_patients}")
    print(f"Total AEs: {safety_assessment.total_aes}")
    print(f"Serious AEs: {safety_assessment.serious_aes}")
    print(f"Drug-related AEs: {safety_assessment.drug_related_aes}")
    print(f"Overall safety status: {safety_assessment.overall_safety_status}")

    print(f"\nDETECTED SAFETY SIGNALS ({len(detected_signals)}):")
    for signal in detected_signals:
        print(f"\n{signal.event_term}:")
        print(f"  Detection method: {signal.detection_method}")
        print(f"  Signal strength: {signal.signal_strength:.3f}")
        print(f"  Alert level: {signal.alert_level}")
        print(f"  Affected patients: {signal.affected_patients}")
        print(f"  Recommended actions: {signal.recommended_actions[:2]}")

    print("\nSAFETY RECOMMENDATIONS:")
    for rec in safety_assessment.recommendations:
        print(f"  • {rec}")

    # Test causality assessment
    if len(ae_data) > 0:
        test_ae = AdverseEvent(
            event_id="test_ae_001",
            patient_id="patient_001",
            event_term="Headache",
            onset_date=datetime.now(),
            severity="Moderate",
            seriousness="Non-serious",
            causality="Unknown",
            outcome="Recovered",
            action_taken="None",
            narrative="Patient reported headache 2 days after treatment",
            reporter="Investigator",
            report_date=datetime.now()
        )

        causality = ae_detector.assess_causality(test_ae)
        print("\nCAUSALITY ASSESSMENT EXAMPLE:")
        print(f"Event: {test_ae.event_term}")
        print(f"Assessed causality: {causality}")

    return ae_detector, detected_signals, safety_assessment


if __name__ == "__main__":
    # Run demonstration
    detector, signals, assessment = create_ae_detection_demo()
    print("\nAutomated adverse event detection completed!")
    print("\nAutomated adverse event detection completed!")

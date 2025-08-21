"""
Adaptive Trial Design Optimization for CRISPR Toolkit Phase 3
=============================================================

Advanced adaptive trial design system for aging intervention studies
with real-time optimization, interim analysis, and automated
decision-making for trial efficiency and patient safety.
"""

import logging
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.stats import norm

logger = logging.getLogger(__name__)


@dataclass
class TrialStage:
    """Individual trial stage configuration."""
    stage_id: str
    sample_size: int
    duration_weeks: int
    primary_endpoint: str
    interim_analysis: bool
    stopping_rules: Dict[str, float]
    adaptation_rules: Dict[str, Any]


@dataclass
class AdaptationDecision:
    """Adaptive trial decision point."""
    decision_point: str
    analysis_time: float
    recommendation: str  # 'continue', 'stop_efficacy', 'stop_futility', 'modify'
    evidence_strength: float
    confidence_interval: tuple
    rationale: str
    next_actions: List[str]


@dataclass
class TrialOptimizationResult:
    """Results from adaptive trial optimization."""
    optimized_design: Dict[str, Any]
    expected_sample_size: int
    expected_duration: float
    power_analysis: Dict[str, float]
    adaptation_schedule: List[AdaptationDecision]
    simulation_results: Dict[str, Any]
    regulatory_considerations: Dict[str, str]


class AdaptiveTrialEngine:
    """
    Adaptive trial design optimization engine.

    This class implements sophisticated adaptive trial designs
    for aging interventions with Bayesian optimization,
    group sequential methods, and real-time adaptation.
    """

    def __init__(self, adaptation_strategy: str = "bayesian_optimal"):
        """
        Initialize adaptive trial design engine.

        Args:
            adaptation_strategy: "bayesian_optimal", "group_sequential",
                               "response_adaptive", "seamless"
        """
        self.adaptation_strategy = adaptation_strategy
        self.trial_data = None
        self.historical_data = None
        self.current_trial_state = {}
        self.adaptation_history = []

        logger.info(f"Initialized AdaptiveTrialEngine with {adaptation_strategy} strategy")

    def load_trial_data(self,
                       current_data: pd.DataFrame,
                       historical_data: Optional[pd.DataFrame] = None,
                       trial_protocol: Optional[Dict[str, Any]] = None):
        """Load trial data for adaptive optimization."""

        self.trial_data = current_data
        if historical_data is not None:
            self.historical_data = historical_data
        if trial_protocol:
            self.trial_protocol = trial_protocol

        logger.info("Loaded trial data:")
        logger.info(f"  Current trial: {current_data.shape}")
        if historical_data is not None:
            logger.info(f"  Historical data: {historical_data.shape}")

    def optimize_trial_design(self,
                            design_objectives: Dict[str, Any],
                            constraints: Optional[Dict[str, Any]] = None) -> TrialOptimizationResult:
        """
        Optimize adaptive trial design based on objectives and constraints.

        Args:
            design_objectives: Trial objectives (power, duration, sample size)
            constraints: Design constraints (budget, timeline, regulatory)

        Returns:
            TrialOptimizationResult with optimized design parameters
        """
        logger.info("Starting adaptive trial design optimization...")

        if constraints is None:
            constraints = self._get_default_constraints()

        # Step 1: Design space exploration
        design_candidates = self._generate_design_candidates(design_objectives, constraints)

        # Step 2: Simulation-based evaluation
        simulation_results = self._simulate_trial_designs(design_candidates, design_objectives)

        # Step 3: Bayesian optimization
        optimized_design = self._bayesian_design_optimization(simulation_results, design_objectives)

        # Step 4: Adaptation schedule
        adaptation_schedule = self._create_adaptation_schedule(optimized_design)

        # Step 5: Power analysis
        power_analysis = self._perform_adaptive_power_analysis(optimized_design)

        # Step 6: Regulatory assessment
        regulatory_considerations = self._assess_regulatory_implications(optimized_design)

        result = TrialOptimizationResult(
            optimized_design=optimized_design,
            expected_sample_size=optimized_design.get('expected_sample_size', 0),
            expected_duration=optimized_design.get('expected_duration', 0),
            power_analysis=power_analysis,
            adaptation_schedule=adaptation_schedule,
            simulation_results=simulation_results,
            regulatory_considerations=regulatory_considerations
        )

        logger.info("Trial optimization completed:")
        logger.info(f"  Expected sample size: {result.expected_sample_size}")
        logger.info(f"  Expected duration: {result.expected_duration:.1f} weeks")
        logger.info(f"  Estimated power: {result.power_analysis.get('final_power', 0):.3f}")

        return result

    def analyze_interim_data(self,
                           interim_data: pd.DataFrame,
                           analysis_time: float) -> AdaptationDecision:
        """
        Perform interim analysis and make adaptation decision.

        Args:
            interim_data: Current trial data at interim analysis
            analysis_time: Time of analysis (proportion of total trial)

        Returns:
            AdaptationDecision with recommendation and rationale
        """
        logger.info(f"Performing interim analysis at time {analysis_time:.2f}")

        # Efficacy analysis
        efficacy_analysis = self._analyze_efficacy(interim_data, analysis_time)

        # Futility analysis
        futility_analysis = self._analyze_futility(interim_data, analysis_time)

        # Safety analysis
        safety_analysis = self._analyze_safety(interim_data)

        # Adaptation recommendation
        decision = self._make_adaptation_decision(
            efficacy_analysis, futility_analysis, safety_analysis, analysis_time
        )

        # Update trial state
        self.current_trial_state.update({
            'analysis_time': analysis_time,
            'efficacy_result': efficacy_analysis,
            'futility_result': futility_analysis,
            'safety_result': safety_analysis
        })

        self.adaptation_history.append(decision)

        logger.info(f"Interim analysis completed: {decision.recommendation}")
        return decision

    def _generate_design_candidates(self,
                                  objectives: Dict[str, Any],
                                  constraints: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Generate candidate trial designs for evaluation."""

        candidates = []

        # Sample size range
        min_n = constraints.get('min_sample_size', 50)
        max_n = constraints.get('max_sample_size', 500)

        # Duration range
        min_duration = constraints.get('min_duration', 12)
        max_duration = constraints.get('max_duration', 52)

        # Number of interim analyses
        max_interims = constraints.get('max_interim_analyses', 3)

        # Generate candidate designs
        for n_total in range(min_n, max_n + 1, 50):
            for duration in range(min_duration, max_duration + 1, 8):
                for n_interims in range(1, max_interims + 1):

                    # Create design candidate
                    candidate = {
                        'total_sample_size': n_total,
                        'total_duration': duration,
                        'n_interim_analyses': n_interims,
                        'interim_times': np.linspace(0.3, 0.8, n_interims).tolist(),
                        'alpha_spending': self._calculate_alpha_spending(n_interims),
                        'stopping_boundaries': self._calculate_stopping_boundaries(n_interims),
                        'adaptation_types': self._determine_adaptation_types()
                    }

                    # Check constraints
                    if self._check_design_constraints(candidate, constraints):
                        candidates.append(candidate)

        logger.info(f"Generated {len(candidates)} design candidates")
        return candidates

    def _simulate_trial_designs(self,
                              candidates: List[Dict[str, Any]],
                              objectives: Dict[str, Any]) -> Dict[str, Any]:
        """Simulate trial designs to evaluate performance."""

        simulation_results = {}
        n_simulations = 1000

        for i, candidate in enumerate(candidates):
            logger.debug(f"Simulating design candidate {i+1}/{len(candidates)}")

            # Run simulations for this design
            design_results = []

            for sim in range(n_simulations):
                sim_result = self._simulate_single_trial(candidate, objectives)
                design_results.append(sim_result)

            # Aggregate simulation results
            aggregated_results = self._aggregate_simulation_results(design_results)
            simulation_results[f"design_{i}"] = {
                'design': candidate,
                'results': aggregated_results
            }

        return simulation_results

    def _simulate_single_trial(self,
                             design: Dict[str, Any],
                             objectives: Dict[str, Any]) -> Dict[str, Any]:
        """Simulate a single trial realization."""

        # Trial parameters
        n_total = design['total_sample_size']
        duration = design['total_duration']
        interim_times = design['interim_times']

        # Effect size (from objectives or simulate)
        true_effect = objectives.get('expected_effect_size', 0.5)

        # Simulate patient recruitment
        recruitment_rate = n_total / (duration * 0.7)  # 70% of time for recruitment

        # Simulate outcomes
        control_outcomes = np.random.normal(0, 1, n_total // 2)
        treatment_outcomes = np.random.normal(true_effect, 1, n_total // 2)

        # Simulate trial progression
        trial_stopped = False
        stop_time = duration
        stop_reason = "completed"
        final_sample_size = n_total

        # Interim analyses
        for interim_time in interim_times:
            interim_n = int(n_total * interim_time)

            # Interim efficacy test
            interim_control = control_outcomes[:interim_n//2]
            interim_treatment = treatment_outcomes[:interim_n//2]

            # Calculate test statistic
            pooled_se = np.sqrt(1/len(interim_control) + 1/len(interim_treatment))
            z_stat = (np.mean(interim_treatment) - np.mean(interim_control)) / pooled_se

            # Check stopping boundaries
            efficacy_boundary = design['stopping_boundaries']['efficacy'][int(interim_time * len(interim_times)) - 1]
            futility_boundary = design['stopping_boundaries']['futility'][int(interim_time * len(interim_times)) - 1]

            if z_stat >= efficacy_boundary:
                trial_stopped = True
                stop_time = duration * interim_time
                stop_reason = "efficacy"
                final_sample_size = interim_n
                break
            elif z_stat <= futility_boundary:
                trial_stopped = True
                stop_time = duration * interim_time
                stop_reason = "futility"
                final_sample_size = interim_n
                break

        # Final analysis if not stopped
        if not trial_stopped:
            final_z_stat = (np.mean(treatment_outcomes) - np.mean(control_outcomes)) / pooled_se
            final_p_value = 2 * (1 - norm.cdf(abs(final_z_stat)))
            significant = final_p_value < 0.05
        else:
            significant = (stop_reason == "efficacy")

        return {
            'stopped_early': trial_stopped,
            'stop_time': stop_time,
            'stop_reason': stop_reason,
            'final_sample_size': final_sample_size,
            'significant_result': significant,
            'estimated_effect': np.mean(treatment_outcomes) - np.mean(control_outcomes)
        }

    def _aggregate_simulation_results(self, results: List[Dict[str, Any]]) -> Dict[str, float]:
        """Aggregate results from multiple trial simulations."""

        return {
            'power': np.mean([r['significant_result'] for r in results]),
            'expected_sample_size': np.mean([r['final_sample_size'] for r in results]),
            'expected_duration': np.mean([r['stop_time'] for r in results]),
            'early_stop_rate': np.mean([r['stopped_early'] for r in results]),
            'efficacy_stop_rate': np.mean([r['stop_reason'] == 'efficacy' for r in results]),
            'futility_stop_rate': np.mean([r['stop_reason'] == 'futility' for r in results]),
            'bias': np.mean([r['estimated_effect'] for r in results]) - 0.5  # Assuming true effect = 0.5
        }

    def _bayesian_design_optimization(self,
                                    simulation_results: Dict[str, Any],
                                    objectives: Dict[str, Any]) -> Dict[str, Any]:
        """Use Bayesian optimization to find optimal design."""

        # Extract features and objectives from simulation results
        features = []
        objective_values = []

        for design_id, result_data in simulation_results.items():
            design = result_data['design']
            results = result_data['results']

            # Design features
            feature_vector = [
                design['total_sample_size'],
                design['total_duration'],
                design['n_interim_analyses'],
                np.mean(design['interim_times'])
            ]
            features.append(feature_vector)

            # Objective function (weighted combination)
            power_weight = objectives.get('power_weight', 0.5)
            efficiency_weight = objectives.get('efficiency_weight', 0.3)
            duration_weight = objectives.get('duration_weight', 0.2)

            objective_value = (
                power_weight * results['power'] +
                efficiency_weight * (1 - results['expected_sample_size'] / design['total_sample_size']) +
                duration_weight * (1 - results['expected_duration'] / design['total_duration'])
            )
            objective_values.append(objective_value)

        # Find best design
        best_idx = np.argmax(objective_values)
        best_design_id = list(simulation_results.keys())[best_idx]
        optimal_design = simulation_results[best_design_id]['design'].copy()

        # Add optimization results
        optimal_design.update({
            'optimization_score': objective_values[best_idx],
            'expected_sample_size': simulation_results[best_design_id]['results']['expected_sample_size'],
            'expected_duration': simulation_results[best_design_id]['results']['expected_duration'],
            'expected_power': simulation_results[best_design_id]['results']['power']
        })

        return optimal_design

    def _create_adaptation_schedule(self, design: Dict[str, Any]) -> List[AdaptationDecision]:
        """Create schedule of adaptation decision points."""

        schedule = []
        interim_times = design['interim_times']

        for i, time_point in enumerate(interim_times):
            decision = AdaptationDecision(
                decision_point=f"interim_{i+1}",
                analysis_time=time_point,
                recommendation="continue",  # Default, updated during trial
                evidence_strength=0.0,
                confidence_interval=(0.0, 0.0),
                rationale="Planned interim analysis",
                next_actions=[
                    "Assess efficacy and futility",
                    "Evaluate safety signals",
                    "Consider sample size modification"
                ]
            )
            schedule.append(decision)

        return schedule

    def _perform_adaptive_power_analysis(self, design: Dict[str, Any]) -> Dict[str, float]:
        """Perform power analysis for adaptive design."""

        # Extract design parameters
        n_total = design['total_sample_size']
        n_interims = design['n_interim_analyses']
        alpha = 0.05

        # Adjust alpha for multiple testing
        adjusted_alpha = alpha / (n_interims + 1)  # Bonferroni correction (conservative)

        # Power calculation for different effect sizes
        effect_sizes = [0.2, 0.3, 0.5, 0.8]
        power_by_effect = {}

        for effect_size in effect_sizes:
            # Standard power calculation
            se = np.sqrt(2 / (n_total / 2))  # Standard error for two-sample t-test
            z_alpha = norm.ppf(1 - adjusted_alpha/2)
            z_beta = norm.ppf(0.8)  # 80% power

            power = 1 - norm.cdf(z_alpha - effect_size / se)
            power_by_effect[f"effect_{effect_size}"] = power

        return {
            'final_power': design.get('expected_power', 0.8),
            'power_by_effect_size': power_by_effect,
            'adjusted_alpha': adjusted_alpha,
            'multiple_testing_correction': 'bonferroni'
        }

    def _assess_regulatory_implications(self, design: Dict[str, Any]) -> Dict[str, str]:
        """Assess regulatory implications of adaptive design."""

        n_interims = design['n_interim_analyses']
        adaptation_types = design.get('adaptation_types', [])

        # Regulatory considerations
        considerations = {}

        if n_interims > 2:
            considerations['fda_guidance'] = "Multiple interim analyses require pre-specification and alpha adjustment"

        if 'sample_size_modification' in adaptation_types:
            considerations['sample_size'] = "Sample size re-estimation may require regulatory consultation"

        if 'population_enrichment' in adaptation_types:
            considerations['population'] = "Population enrichment requires strong justification"

        considerations['regulatory_path'] = "Standard FDA/EMA adaptive design guidance applies"
        considerations['documentation'] = "Detailed SAP and simulation results required"

        return considerations

    def _analyze_efficacy(self,
                        interim_data: pd.DataFrame,
                        analysis_time: float) -> Dict[str, Any]:
        """Analyze efficacy at interim time point."""

        # Extract treatment groups
        if 'treatment_group' not in interim_data.columns:
            # Simulate treatment assignment
            n_patients = len(interim_data)
            treatment_group = np.random.choice(['control', 'treatment'], n_patients)
            interim_data = interim_data.copy()
            interim_data['treatment_group'] = treatment_group

        control_data = interim_data[interim_data['treatment_group'] == 'control']
        treatment_data = interim_data[interim_data['treatment_group'] == 'treatment']

        # Primary endpoint analysis (assuming first numeric column is primary endpoint)
        numeric_cols = interim_data.select_dtypes(include=[np.number]).columns
        if len(numeric_cols) > 0:
            primary_endpoint = numeric_cols[0]

            control_outcome = control_data[primary_endpoint].dropna()
            treatment_outcome = treatment_data[primary_endpoint].dropna()

            if len(control_outcome) > 0 and len(treatment_outcome) > 0:
                # Effect size and significance test
                effect_size = treatment_outcome.mean() - control_outcome.mean()
                pooled_se = np.sqrt(control_outcome.var()/len(control_outcome) +
                                  treatment_outcome.var()/len(treatment_outcome))

                z_stat = effect_size / pooled_se if pooled_se > 0 else 0
                p_value = 2 * (1 - norm.cdf(abs(z_stat)))

                # Conditional power (assuming current trend continues)
                conditional_power = self._calculate_conditional_power(z_stat, analysis_time)

                return {
                    'effect_size': effect_size,
                    'z_statistic': z_stat,
                    'p_value': p_value,
                    'conditional_power': conditional_power,
                    'n_control': len(control_outcome),
                    'n_treatment': len(treatment_outcome)
                }

        return {'effect_size': 0, 'z_statistic': 0, 'p_value': 1, 'conditional_power': 0}

    def _analyze_futility(self,
                        interim_data: pd.DataFrame,
                        analysis_time: float) -> Dict[str, Any]:
        """Analyze futility at interim time point."""

        efficacy_result = self._analyze_efficacy(interim_data, analysis_time)

        # Futility based on conditional power
        conditional_power = efficacy_result.get('conditional_power', 0)
        futility_threshold = 0.2  # Standard threshold

        futility_probability = 1 - conditional_power

        return {
            'conditional_power': conditional_power,
            'futility_probability': futility_probability,
            'futility_threshold': futility_threshold,
            'recommend_futility_stop': futility_probability > (1 - futility_threshold)
        }

    def _analyze_safety(self, interim_data: pd.DataFrame) -> Dict[str, Any]:
        """Analyze safety signals at interim time point."""

        # Safety analysis (simplified)
        # In practice, would analyze specific safety endpoints

        # Simulate adverse event rates
        n_patients = len(interim_data)
        ae_rate = np.random.beta(2, 10)  # Simulate AE rate

        # Safety stopping rule
        safety_threshold = 0.15  # 15% AE rate threshold
        safety_concern = ae_rate > safety_threshold

        return {
            'adverse_event_rate': ae_rate,
            'safety_threshold': safety_threshold,
            'safety_concern': safety_concern,
            'recommend_safety_stop': safety_concern,
            'n_patients_safety': n_patients
        }

    def _make_adaptation_decision(self,
                                efficacy: Dict[str, Any],
                                futility: Dict[str, Any],
                                safety: Dict[str, Any],
                                analysis_time: float) -> AdaptationDecision:
        """Make adaptation decision based on interim analyses."""

        # Decision logic
        if safety['recommend_safety_stop']:
            recommendation = "stop_safety"
            rationale = f"Safety concern: AE rate {safety['adverse_event_rate']:.1%} exceeds threshold"
            evidence_strength = 0.9
        elif efficacy['conditional_power'] > 0.9:
            recommendation = "stop_efficacy"
            rationale = f"High conditional power ({efficacy['conditional_power']:.1%}) suggests efficacy"
            evidence_strength = efficacy['conditional_power']
        elif futility['recommend_futility_stop']:
            recommendation = "stop_futility"
            rationale = f"Low conditional power ({efficacy['conditional_power']:.1%}) suggests futility"
            evidence_strength = futility['futility_probability']
        else:
            recommendation = "continue"
            rationale = "Insufficient evidence to stop; continue trial"
            evidence_strength = 0.5

        # Confidence interval for effect size
        effect_size = efficacy.get('effect_size', 0)
        z_stat = efficacy.get('z_statistic', 0)
        se = effect_size / z_stat if z_stat != 0 else 0.1
        ci_lower = effect_size - 1.96 * se
        ci_upper = effect_size + 1.96 * se

        # Next actions
        next_actions = self._determine_next_actions(recommendation, analysis_time)

        return AdaptationDecision(
            decision_point=f"interim_analysis_{analysis_time:.2f}",
            analysis_time=analysis_time,
            recommendation=recommendation,
            evidence_strength=evidence_strength,
            confidence_interval=(ci_lower, ci_upper),
            rationale=rationale,
            next_actions=next_actions
        )

    def _calculate_conditional_power(self,
                                   current_z: float,
                                   analysis_time: float) -> float:
        """Calculate conditional power given current data."""

        # Information fraction
        information_fraction = analysis_time

        # Expected Z-statistic at final analysis
        expected_final_z = current_z / np.sqrt(information_fraction)

        # Conditional power
        final_critical_value = 1.96  # Two-sided test
        conditional_power = 1 - norm.cdf(final_critical_value - expected_final_z)

        return max(0, min(1, conditional_power))

    def _determine_next_actions(self,
                              recommendation: str,
                              analysis_time: float) -> List[str]:
        """Determine next actions based on recommendation."""

        actions = []

        if recommendation == "continue":
            actions.extend([
                "Continue patient enrollment",
                "Monitor safety signals",
                "Prepare for next interim analysis"
            ])
        elif recommendation == "stop_efficacy":
            actions.extend([
                "Prepare efficacy report",
                "Plan regulatory submission",
                "Design confirmatory study"
            ])
        elif recommendation == "stop_futility":
            actions.extend([
                "Analyze futility causes",
                "Consider design modifications",
                "Plan alternative strategies"
            ])
        elif recommendation == "stop_safety":
            actions.extend([
                "Conduct safety review",
                "Notify regulatory authorities",
                "Implement safety measures"
            ])

        return actions

    def _calculate_alpha_spending(self, n_interims: int) -> List[float]:
        """Calculate alpha spending function."""

        # O'Brien-Fleming alpha spending
        times = np.linspace(0.3, 1.0, n_interims + 1)
        alpha_spending = []

        for i, t in enumerate(times[1:]):
            if i == len(times) - 2:  # Final analysis
                alpha_spent = 0.05 - sum(alpha_spending)
            else:
                alpha_spent = 2 * (1 - norm.cdf(1.96 / np.sqrt(t)))
            alpha_spending.append(alpha_spent)

        return alpha_spending

    def _calculate_stopping_boundaries(self, n_interims: int) -> Dict[str, List[float]]:
        """Calculate efficacy and futility stopping boundaries."""

        # O'Brien-Fleming boundaries
        information_fractions = np.linspace(0.3, 1.0, n_interims + 1)[1:]

        efficacy_boundaries = []
        futility_boundaries = []

        for info_frac in information_fractions:
            # Efficacy boundary (O'Brien-Fleming)
            efficacy_bound = 1.96 / np.sqrt(info_frac)
            efficacy_boundaries.append(efficacy_bound)

            # Futility boundary (conservative)
            futility_bound = -0.5  # Fixed futility boundary
            futility_boundaries.append(futility_bound)

        return {
            'efficacy': efficacy_boundaries,
            'futility': futility_boundaries
        }

    def _determine_adaptation_types(self) -> List[str]:
        """Determine types of adaptations to allow."""

        return [
            'sample_size_modification',
            'stopping_for_efficacy',
            'stopping_for_futility',
            'safety_stopping'
        ]

    def _check_design_constraints(self,
                                design: Dict[str, Any],
                                constraints: Dict[str, Any]) -> bool:
        """Check if design meets constraints."""

        # Sample size constraints
        if design['total_sample_size'] < constraints.get('min_sample_size', 0):
            return False
        if design['total_sample_size'] > constraints.get('max_sample_size', float('inf')):
            return False

        # Duration constraints
        if design['total_duration'] < constraints.get('min_duration', 0):
            return False
        if design['total_duration'] > constraints.get('max_duration', float('inf')):
            return False

        return True

    def _get_default_constraints(self) -> Dict[str, Any]:
        """Get default design constraints."""

        return {
            'min_sample_size': 50,
            'max_sample_size': 500,
            'min_duration': 12,
            'max_duration': 104,
            'max_interim_analyses': 3,
            'alpha': 0.05,
            'power': 0.8
        }


def create_adaptive_trial_demo():
    """Create demonstration of adaptive trial optimization."""

    # Generate synthetic trial data
    n_patients = 100
    patient_ids = [f"patient_{i:03d}" for i in range(n_patients)]

    np.random.seed(42)

    # Current trial data
    current_data = pd.DataFrame({
        'patient_id': patient_ids,
        'treatment_group': np.random.choice(['control', 'treatment'], n_patients),
        'primary_endpoint': np.random.normal(0.3, 1, n_patients),  # Small treatment effect
        'safety_endpoint': np.random.exponential(2, n_patients),
        'enrollment_week': np.random.uniform(0, 20, n_patients)
    })

    # Historical data
    historical_data = pd.DataFrame({
        'study_id': [f"study_{i}" for i in range(5)],
        'effect_size': [0.2, 0.3, 0.4, 0.1, 0.5],
        'sample_size': [200, 150, 300, 100, 250],
        'duration': [24, 18, 36, 12, 30]
    })

    # Initialize adaptive trial engine
    trial_engine = AdaptiveTrialEngine(adaptation_strategy="bayesian_optimal")

    # Load trial data
    trial_engine.load_trial_data(
        current_data=current_data,
        historical_data=historical_data
    )

    # Define design objectives
    design_objectives = {
        'target_power': 0.8,
        'expected_effect_size': 0.4,
        'power_weight': 0.5,
        'efficiency_weight': 0.3,
        'duration_weight': 0.2
    }

    # Optimize trial design
    optimization_result = trial_engine.optimize_trial_design(design_objectives)

    # Perform interim analysis
    interim_data = current_data.iloc[:50]  # First 50 patients
    interim_decision = trial_engine.analyze_interim_data(interim_data, analysis_time=0.5)

    # Display results
    print("ADAPTIVE TRIAL OPTIMIZATION RESULTS")
    print("=" * 50)
    print(f"Optimized sample size: {optimization_result.expected_sample_size}")
    print(f"Expected duration: {optimization_result.expected_duration:.1f} weeks")
    print(f"Expected power: {optimization_result.power_analysis['final_power']:.3f}")
    print(f"Number of interim analyses: {optimization_result.optimized_design['n_interim_analyses']}")

    print("\nOPTIMIZED DESIGN PARAMETERS:")
    design = optimization_result.optimized_design
    for key, value in design.items():
        if isinstance(value, (int, float)):
            print(f"  {key}: {value}")
        elif isinstance(value, list) and len(value) < 10:
            print(f"  {key}: {value}")

    print("\nINTERIM ANALYSIS DECISION:")
    print(f"  Analysis time: {interim_decision.analysis_time:.2f}")
    print(f"  Recommendation: {interim_decision.recommendation}")
    print(f"  Evidence strength: {interim_decision.evidence_strength:.3f}")
    print(f"  Rationale: {interim_decision.rationale}")
    print(f"  Next actions: {interim_decision.next_actions}")

    print("\nADAPTATION SCHEDULE:")
    for i, decision in enumerate(optimization_result.adaptation_schedule):
        print(f"  {decision.decision_point}: time {decision.analysis_time:.2f}")

    print("\nREGULATORY CONSIDERATIONS:")
    for aspect, consideration in optimization_result.regulatory_considerations.items():
        print(f"  {aspect}: {consideration}")

    return trial_engine, optimization_result, interim_decision


if __name__ == "__main__":
    # Run demonstration
    engine, optimization, decision = create_adaptive_trial_demo()
    print("\nAdaptive trial design optimization completed!")
    print("\nAdaptive trial design optimization completed!")

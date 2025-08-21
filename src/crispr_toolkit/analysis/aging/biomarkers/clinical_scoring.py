"""
Clinical Aging Scorer Module

Provides clinical interpretation and scoring of aging biomarkers with specific
focus on ReHMGB1/RAGE signaling pathways and senescence research applications.

Clinical Features:
- Risk stratification for aging-related diseases
- Mortality prediction scoring
- Intervention recommendation algorithms
- Biomarker interpretation for clinical translation
- ReHMGB1 pathway clinical significance scoring

Author: CRISPR Toolkit Development Team
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import pandas as pd


@dataclass
class ClinicalRiskScore:
    """Container for clinical risk assessment results."""
    score: float
    category: str
    percentile: Optional[float]
    interpretation: str
    recommendations: List[str]
    confidence: float


@dataclass
class AgingProfile:
    """Container for comprehensive aging profile."""
    chronological_age: float
    biological_age: float
    age_acceleration: float
    mortality_risk: Optional[float]
    senescence_burden: float
    inflammation_score: float
    pathway_scores: Dict[str, float]


class ClinicalAgingScorer:
    """
    Clinical aging assessment and scoring system.

    Provides clinically interpretable scores and recommendations based on
    aging biomarker analysis with focus on ReHMGB1/senescence pathways.
    """

    def __init__(self, verbose: bool = True):
        """
        Initialize clinical aging scorer.

        Args:
            verbose: Enable verbose logging
        """
        self.verbose = verbose

        # Load clinical reference data
        self.clinical_thresholds = self._load_clinical_thresholds()
        self.risk_models = self._load_risk_models()
        self.intervention_guidelines = self._load_intervention_guidelines()

        if verbose:
            print("🏥 Clinical aging scorer initialized")

    def _load_clinical_thresholds(self) -> Dict[str, Dict]:
        """Load clinical thresholds for aging biomarkers."""
        return {
            'age_acceleration': {
                'very_low': -10,    # >10 years younger
                'low': -5,          # 5-10 years younger
                'normal': 5,        # ±5 years
                'high': 10,         # 5-10 years older
                'very_high': float('inf')  # >10 years older
            },
            'mortality_risk': {
                'very_low': 0.01,   # <1% 10-year mortality
                'low': 0.05,        # 1-5% 10-year mortality
                'moderate': 0.15,   # 5-15% 10-year mortality
                'high': 0.30,       # 15-30% 10-year mortality
                'very_high': float('inf')  # >30% 10-year mortality
            },
            'senescence_burden': {
                'minimal': 0.2,     # Minimal senescence markers
                'low': 0.4,         # Low senescence burden
                'moderate': 0.6,    # Moderate senescence burden
                'high': 0.8,        # High senescence burden
                'severe': float('inf')  # Severe senescence burden
            },
            'inflammation_score': {
                'low': 0.3,         # Low inflammatory aging
                'moderate': 0.6,    # Moderate inflammatory aging
                'high': 0.8,        # High inflammatory aging
                'severe': float('inf')  # Severe inflammatory aging
            }
        }

    def _load_risk_models(self) -> Dict[str, Dict]:
        """Load clinical risk prediction models."""
        return {
            'cardiovascular_risk': {
                'description': 'Cardiovascular disease risk from aging biomarkers',
                'coefficients': {
                    'age_acceleration': 0.15,
                    'inflammation_score': 0.25,
                    'senescence_burden': 0.20
                },
                'baseline_risk': 0.05
            },
            'all_cause_mortality': {
                'description': 'All-cause mortality risk prediction',
                'coefficients': {
                    'age_acceleration': 0.12,
                    'mortality_risk_biomarker': 0.40,
                    'senescence_burden': 0.18,
                    'inflammation_score': 0.15
                },
                'baseline_risk': 0.02
            },
            'cognitive_decline': {
                'description': 'Cognitive decline risk assessment',
                'coefficients': {
                    'age_acceleration': 0.10,
                    'inflammation_score': 0.30,
                    'senescence_burden': 0.25
                },
                'baseline_risk': 0.03
            },
            'rehmgb1_related_aging': {
                'description': 'ReHMGB1/RAGE pathway aging risk',
                'coefficients': {
                    'rage_signaling_score': 0.35,
                    'inflammation_score': 0.30,
                    'senescence_burden': 0.20,
                    'nfkb_activity': 0.15
                },
                'baseline_risk': 0.04
            }
        }

    def _load_intervention_guidelines(self) -> Dict[str, Dict]:
        """Load intervention recommendation guidelines."""
        return {
            'lifestyle_interventions': {
                'high_age_acceleration': [
                    'Caloric restriction or intermittent fasting',
                    'Regular aerobic and resistance exercise',
                    'Stress reduction techniques (meditation, yoga)',
                    'Optimize sleep quality (7-9 hours)',
                    'Mediterranean or anti-inflammatory diet'
                ],
                'high_inflammation': [
                    'Anti-inflammatory diet (omega-3 rich)',
                    'Regular exercise to reduce inflammation',
                    'Stress management and adequate sleep',
                    'Consider curcumin or other anti-inflammatory supplements',
                    'Avoid pro-inflammatory foods (processed foods, sugar)'
                ],
                'high_senescence_burden': [
                    'Senolytic interventions (under medical supervision)',
                    'Exercise to promote cellular renewal',
                    'Antioxidant-rich diet',
                    'Consider fasting-mimicking diets',
                    'Regular monitoring of senescence markers'
                ]
            },
            'pharmaceutical_interventions': {
                'high_mortality_risk': [
                    'Consider metformin (anti-aging effects)',
                    'Rapamycin (under medical supervision)',
                    'NAD+ precursors (NMN, NR)',
                    'Senolytics for high senescence burden',
                    'Regular comprehensive health monitoring'
                ],
                'rehmgb1_pathway_dysregulation': [
                    'Anti-RAGE therapeutics (if available)',
                    'NF-κB pathway modulators',
                    'JAK/STAT pathway inhibitors',
                    'Anti-inflammatory medications',
                    'Antioxidant supplementation'
                ]
            },
            'monitoring_recommendations': {
                'all_patients': [
                    'Annual aging biomarker assessment',
                    'Regular cardiovascular risk evaluation',
                    'Cognitive function monitoring',
                    'Inflammatory marker tracking',
                    'Quality of life assessments'
                ],
                'high_risk_patients': [
                    'Quarterly aging biomarker monitoring',
                    'Enhanced cardiovascular surveillance',
                    'Monthly inflammatory marker tracking',
                    'Consider continuous glucose monitoring',
                    'Regular senescence marker evaluation'
                ]
            }
        }

    def calculate_clinical_risk_score(
        self,
        aging_data: Dict[str, Any],
        chronological_age: float,
        risk_type: str = 'all_cause_mortality'
    ) -> ClinicalRiskScore:
        """
        Calculate clinical risk score based on aging biomarkers.

        Args:
            aging_data: Dictionary of aging biomarker data
            chronological_age: Patient's chronological age
            risk_type: Type of risk to calculate

        Returns:
            Clinical risk score with interpretation
        """
        if risk_type not in self.risk_models:
            available_types = list(self.risk_models.keys())
            raise ValueError(f"Risk type '{risk_type}' not available. "
                           f"Choose from: {available_types}")

        model = self.risk_models[risk_type]

        # Calculate risk score
        risk_score = model['baseline_risk']

        for factor, coefficient in model['coefficients'].items():
            if factor in aging_data:
                risk_score += coefficient * aging_data[factor]

        # Apply age adjustment
        age_factor = max(0, (chronological_age - 50) / 50)  # Increase risk with age
        risk_score *= (1 + age_factor * 0.5)

        # Ensure score is between 0 and 1
        risk_score = max(0, min(1, risk_score))

        # Categorize risk
        category = self._categorize_risk(risk_score, risk_type)

        # Generate interpretation
        interpretation = self._interpret_risk_score(risk_score, risk_type, category)

        # Generate recommendations
        recommendations = self._generate_risk_recommendations(
            category, risk_type, aging_data
        )

        # Calculate confidence (based on data completeness)
        data_completeness = len(aging_data) / len(model['coefficients'])
        confidence = min(1.0, data_completeness * 1.2)

        return ClinicalRiskScore(
            score=risk_score,
            category=category,
            percentile=None,  # Would need population data
            interpretation=interpretation,
            recommendations=recommendations,
            confidence=confidence
        )

    def _categorize_risk(self, score: float, risk_type: str) -> str:
        """Categorize risk score into clinical categories."""
        if risk_type == 'all_cause_mortality':
            if score < 0.05:
                return 'low'
            elif score < 0.15:
                return 'moderate'
            elif score < 0.30:
                return 'high'
            else:
                return 'very_high'
        else:
            # Generic categorization
            if score < 0.1:
                return 'low'
            elif score < 0.3:
                return 'moderate'
            elif score < 0.6:
                return 'high'
            else:
                return 'very_high'

    def _interpret_risk_score(
        self,
        score: float,
        risk_type: str,
        category: str
    ) -> str:
        """Generate clinical interpretation of risk score."""
        interpretations = {
            'all_cause_mortality': {
                'low': f"Low mortality risk ({score:.1%} estimated 10-year risk). "
                       "Continue current health practices.",
                'moderate': f"Moderate mortality risk ({score:.1%} estimated 10-year risk). "
                           "Consider lifestyle interventions.",
                'high': f"High mortality risk ({score:.1%} estimated 10-year risk). "
                        "Recommend comprehensive intervention strategy.",
                'very_high': f"Very high mortality risk ({score:.1%} estimated 10-year risk). "
                            "Urgent intervention required."
            },
            'cardiovascular_risk': {
                'low': f"Low cardiovascular risk ({score:.1%}). "
                       "Maintain heart-healthy lifestyle.",
                'moderate': f"Moderate cardiovascular risk ({score:.1%}). "
                           "Focus on cardiovascular health optimization.",
                'high': f"High cardiovascular risk ({score:.1%}). "
                        "Consider cardiology consultation.",
                'very_high': f"Very high cardiovascular risk ({score:.1%}). "
                            "Immediate cardiovascular intervention needed."
            },
            'cognitive_decline': {
                'low': f"Low cognitive decline risk ({score:.1%}). "
                       "Continue cognitive health practices.",
                'moderate': f"Moderate cognitive decline risk ({score:.1%}). "
                           "Implement cognitive protection strategies.",
                'high': f"High cognitive decline risk ({score:.1%}). "
                        "Consider neurology consultation.",
                'very_high': f"Very high cognitive decline risk ({score:.1%}). "
                            "Urgent cognitive health intervention needed."
            },
            'rehmgb1_related_aging': {
                'low': f"Low ReHMGB1 pathway aging burden ({score:.1%}). "
                       "RAGE signaling appears well-controlled.",
                'moderate': f"Moderate ReHMGB1 pathway aging ({score:.1%}). "
                           "Consider anti-inflammatory interventions.",
                'high': f"High ReHMGB1 pathway aging burden ({score:.1%}). "
                        "Recommend targeted RAGE pathway interventions.",
                'very_high': f"Very high ReHMGB1 pathway dysfunction ({score:.1%}). "
                            "Urgent pathway-specific intervention required."
            }
        }

        return interpretations.get(risk_type, {}).get(
            category,
            f"{category.title()} risk ({score:.1%}) for {risk_type}"
        )

    def _generate_risk_recommendations(
        self,
        category: str,
        risk_type: str,
        aging_data: Dict[str, Any]
    ) -> List[str]:
        """Generate personalized recommendations based on risk assessment."""
        recommendations = []

        # Base recommendations by risk category
        if category in ['high', 'very_high']:
            if 'age_acceleration' in aging_data and aging_data['age_acceleration'] > 5:
                recommendations.extend(
                    self.intervention_guidelines['lifestyle_interventions']['high_age_acceleration']
                )

            if 'inflammation_score' in aging_data and aging_data['inflammation_score'] > 0.6:
                recommendations.extend(
                    self.intervention_guidelines['lifestyle_interventions']['high_inflammation']
                )

            if 'senescence_burden' in aging_data and aging_data['senescence_burden'] > 0.6:
                recommendations.extend(
                    self.intervention_guidelines['lifestyle_interventions']['high_senescence_burden']
                )

        # Risk-type specific recommendations
        if risk_type == 'all_cause_mortality' and category in ['high', 'very_high']:
            recommendations.extend(
                self.intervention_guidelines['pharmaceutical_interventions']['high_mortality_risk']
            )

        if risk_type == 'rehmgb1_related_aging' and category in ['high', 'very_high']:
            recommendations.extend(
                self.intervention_guidelines['pharmaceutical_interventions']['rehmgb1_pathway_dysregulation']
            )

        # Monitoring recommendations
        if category in ['high', 'very_high']:
            recommendations.extend(
                self.intervention_guidelines['monitoring_recommendations']['high_risk_patients']
            )
        else:
            recommendations.extend(
                self.intervention_guidelines['monitoring_recommendations']['all_patients']
            )

        # Remove duplicates while preserving order
        seen = set()
        unique_recommendations = []
        for rec in recommendations:
            if rec not in seen:
                seen.add(rec)
                unique_recommendations.append(rec)

        return unique_recommendations

    def create_aging_profile(
        self,
        aging_predictions: pd.DataFrame,
        chronological_age: float,
        expression_data: Optional[pd.DataFrame] = None
    ) -> AgingProfile:
        """
        Create comprehensive aging profile for an individual.

        Args:
            aging_predictions: Aging clock predictions
            chronological_age: Chronological age
            expression_data: Gene expression data (optional)

        Returns:
            Comprehensive aging profile
        """
        # Calculate biological age (mean of available clocks)
        if 'predicted_age' in aging_predictions.columns:
            biological_age = aging_predictions['predicted_age'].mean()
        else:
            biological_age = aging_predictions['predicted_value'].mean()

        # Calculate age acceleration
        age_acceleration = biological_age - chronological_age

        # Extract mortality risk if available
        mortality_risk = None
        mortality_clocks = ['GrimAge', 'grimage', 'PhenoAge', 'phenoage']
        mortality_preds = aging_predictions[
            aging_predictions['clock'].isin(mortality_clocks)
        ]
        if not mortality_preds.empty:
            if 'predicted_age' in mortality_preds.columns:
                mortality_risk = mortality_preds['predicted_age'].mean()
            else:
                mortality_risk = mortality_preds['predicted_value'].mean()

        # Calculate senescence burden
        senescence_burden = self._calculate_senescence_burden(
            aging_predictions, expression_data
        )

        # Calculate inflammation score
        inflammation_score = self._calculate_inflammation_score(
            aging_predictions, expression_data
        )

        # Calculate pathway-specific scores
        pathway_scores = self._calculate_pathway_scores(
            aging_predictions, expression_data
        )

        return AgingProfile(
            chronological_age=chronological_age,
            biological_age=biological_age,
            age_acceleration=age_acceleration,
            mortality_risk=mortality_risk,
            senescence_burden=senescence_burden,
            inflammation_score=inflammation_score,
            pathway_scores=pathway_scores
        )

    def _calculate_senescence_burden(
        self,
        aging_predictions: pd.DataFrame,
        expression_data: Optional[pd.DataFrame]
    ) -> float:
        """Calculate senescence burden score."""
        # Use senescence-related clocks
        senescence_clocks = ['PhenoAge', 'DunedinPACE', 'phenoage', 'dunedinpace']
        senescence_preds = aging_predictions[
            aging_predictions['clock'].isin(senescence_clocks)
        ]

        if senescence_preds.empty:
            return 0.5  # Default moderate burden

        # Normalize predictions to 0-1 scale
        if 'predicted_age' in senescence_preds.columns:
            pred_values = senescence_preds['predicted_age']
        else:
            pred_values = senescence_preds['predicted_value']

        # Convert to burden score (higher values = higher burden)
        burden_score = (pred_values.mean() - 30) / 50  # Assume age range 30-80
        burden_score = max(0, min(1, burden_score))

        return burden_score

    def _calculate_inflammation_score(
        self,
        aging_predictions: pd.DataFrame,
        expression_data: Optional[pd.DataFrame]
    ) -> float:
        """Calculate inflammation aging score."""
        # Use inflammation-related clocks
        inflam_clocks = ['GrimAge', 'grimage', 'inflammage']
        inflam_preds = aging_predictions[
            aging_predictions['clock'].isin(inflam_clocks)
        ]

        if inflam_preds.empty:
            return 0.5  # Default moderate inflammation

        # Normalize to 0-1 scale
        if 'predicted_age' in inflam_preds.columns:
            pred_values = inflam_preds['predicted_age']
        else:
            pred_values = inflam_preds['predicted_value']

        inflammation_score = (pred_values.mean() - 30) / 50
        inflammation_score = max(0, min(1, inflammation_score))

        return inflammation_score

    def _calculate_pathway_scores(
        self,
        aging_predictions: pd.DataFrame,
        expression_data: Optional[pd.DataFrame]
    ) -> Dict[str, float]:
        """Calculate pathway-specific aging scores."""
        pathway_scores = {}

        # Clock-based pathway scores
        clock_pathways = {
            'dna_methylation': ['Horvath2013', 'Hannum2013', 'horvath2013', 'hannum2013'],
            'mortality_prediction': ['GrimAge', 'PhenoAge', 'grimage', 'phenoage'],
            'aging_pace': ['DunedinPACE', 'dunedinpace']
        }

        for pathway, clocks in clock_pathways.items():
            pathway_preds = aging_predictions[
                aging_predictions['clock'].isin(clocks)
            ]

            if not pathway_preds.empty:
                if 'predicted_age' in pathway_preds.columns:
                    pathway_score = pathway_preds['predicted_age'].mean()
                else:
                    pathway_score = pathway_preds['predicted_value'].mean()

                # Normalize to 0-1 scale
                pathway_scores[pathway] = max(0, min(1, (pathway_score - 30) / 50))

        return pathway_scores

    def generate_clinical_report(
        self,
        aging_profile: AgingProfile,
        risk_scores: Dict[str, ClinicalRiskScore],
        patient_id: str = "Patient"
    ) -> str:
        """
        Generate comprehensive clinical report.

        Args:
            aging_profile: Patient's aging profile
            risk_scores: Dictionary of risk assessments
            patient_id: Patient identifier

        Returns:
            Formatted clinical report
        """
        report = f"""
CLINICAL AGING ASSESSMENT REPORT
{'=' * 50}

Patient ID: {patient_id}
Report Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}

AGING PROFILE SUMMARY
{'-' * 30}
Chronological Age: {aging_profile.chronological_age:.1f} years
Biological Age: {aging_profile.biological_age:.1f} years
Age Acceleration: {aging_profile.age_acceleration:+.1f} years

Senescence Burden: {aging_profile.senescence_burden:.1%} ({self._categorize_burden(aging_profile.senescence_burden)})
Inflammation Score: {aging_profile.inflammation_score:.1%} ({self._categorize_inflammation(aging_profile.inflammation_score)})
"""

        if aging_profile.mortality_risk:
            report += f"Mortality Risk Biomarker: {aging_profile.mortality_risk:.1f} years\n"

        report += f"\nPATHWAY-SPECIFIC SCORES\n{'-' * 30}\n"
        for pathway, score in aging_profile.pathway_scores.items():
            report += f"{pathway.replace('_', ' ').title()}: {score:.1%}\n"

        report += f"\nRISK ASSESSMENTS\n{'-' * 30}\n"
        for risk_type, risk_score in risk_scores.items():
            report += f"\n{risk_type.replace('_', ' ').title()}:\n"
            report += f"  Risk Score: {risk_score.score:.1%}\n"
            report += f"  Category: {risk_score.category.upper()}\n"
            report += f"  Confidence: {risk_score.confidence:.1%}\n"
            report += f"  Interpretation: {risk_score.interpretation}\n"

        # Compile all recommendations
        all_recommendations = set()
        for risk_score in risk_scores.values():
            all_recommendations.update(risk_score.recommendations)

        report += f"\nCLINICAL RECOMMENDATIONS\n{'-' * 30}\n"
        for i, rec in enumerate(sorted(all_recommendations), 1):
            report += f"{i}. {rec}\n"

        report += f"\nREHMGB1 PATHWAY CONSIDERATIONS\n{'-' * 30}\n"
        report += self._generate_rehmgb1_commentary(aging_profile, risk_scores)

        report += f"\nDISCLAIMER\n{'-' * 30}\n"
        report += """This report is for research purposes only and should not be used
for clinical decision-making without validation by qualified healthcare
professionals. Aging biomarkers are still under development and their
clinical utility is being established."""

        return report

    def _categorize_burden(self, burden: float) -> str:
        """Categorize senescence burden."""
        thresholds = self.clinical_thresholds['senescence_burden']
        if burden < thresholds['minimal']:
            return 'Minimal'
        elif burden < thresholds['low']:
            return 'Low'
        elif burden < thresholds['moderate']:
            return 'Moderate'
        elif burden < thresholds['high']:
            return 'High'
        else:
            return 'Severe'

    def _categorize_inflammation(self, inflammation: float) -> str:
        """Categorize inflammation score."""
        thresholds = self.clinical_thresholds['inflammation_score']
        if inflammation < thresholds['low']:
            return 'Low'
        elif inflammation < thresholds['moderate']:
            return 'Moderate'
        elif inflammation < thresholds['high']:
            return 'High'
        else:
            return 'Severe'

    def _generate_rehmgb1_commentary(
        self,
        aging_profile: AgingProfile,
        risk_scores: Dict[str, ClinicalRiskScore]
    ) -> str:
        """Generate ReHMGB1-specific clinical commentary."""
        commentary = ""

        # High inflammation commentary
        if aging_profile.inflammation_score > 0.6:
            commentary += """High inflammation score suggests potential RAGE pathway activation.
ReHMGB1 release may be contributing to systemic inflammatory aging.
Consider anti-inflammatory interventions and RAGE pathway modulators.\n\n"""

        # Age acceleration commentary
        if aging_profile.age_acceleration > 5:
            commentary += """Significant age acceleration may involve ReHMGB1-mediated senescence.
RAGE signaling could be accelerating cellular aging processes.
Monitor senescence markers and consider senolytic interventions.\n\n"""

        # ReHMGB1-specific risk assessment
        if 'rehmgb1_related_aging' in risk_scores:
            rehmgb1_risk = risk_scores['rehmgb1_related_aging']
            if rehmgb1_risk.category in ['high', 'very_high']:
                commentary += f"""ReHMGB1 pathway shows {rehmgb1_risk.category} dysfunction.
This indicates significant RAGE-mediated aging processes.
Prioritize pathway-specific interventions and regular monitoring.\n\n"""

        if not commentary:
            commentary = """ReHMGB1 pathway markers appear within normal ranges.
Continue monitoring as part of comprehensive aging assessment.\n"""

        return commentary.strip()


def create_demo_clinical_scoring() -> str:
    """Create demonstration script for clinical aging scoring."""

    demo_script = '''#!/usr/bin/env python3
"""
Clinical Aging Scoring Demo

Demonstrates clinical interpretation and scoring of aging biomarkers
with focus on ReHMGB1/RAGE pathway clinical significance.
"""

import pandas as pd
import numpy as np
from crispr_toolkit.analysis.aging.biomarkers import ClinicalAgingScorer

def main():
    print("🏥 Clinical Aging Scoring Demo")
    print("=" * 50)

    # Initialize clinical scorer
    print("\\n1. Initializing clinical aging scorer...")
    scorer = ClinicalAgingScorer(verbose=True)

    # Create sample aging predictions
    print("\\n2. Generating sample aging biomarker data...")

    # Simulate aging predictions for a 55-year-old patient
    aging_predictions = pd.DataFrame({
        'sample_id': ['Patient_001'] * 4,
        'clock': ['Horvath2013', 'PhenoAge', 'GrimAge', 'DunedinPACE'],
        'predicted_age': [62.5, 65.2, 68.1, 1.15],  # Last one is pace, not age
        'predicted_value': [62.5, 65.2, 68.1, 1.15]
    })

    chronological_age = 55.0

    print(f"   📊 Sample patient: {chronological_age} years old")
    print("   🔬 Aging clock predictions:")
    for _, row in aging_predictions.iterrows():
        value = row['predicted_age'] if row['clock'] != 'DunedinPACE' else f"{row['predicted_age']:.2f}x pace"
        print(f"      {row['clock']}: {value}")

    # Create aging profile
    print("\\n3. Creating comprehensive aging profile...")
    aging_profile = scorer.create_aging_profile(
        aging_predictions=aging_predictions,
        chronological_age=chronological_age
    )

    print(f"   📈 Biological age: {aging_profile.biological_age:.1f} years")
    print(f"   ⚡ Age acceleration: {aging_profile.age_acceleration:+.1f} years")
    print(f"   🔥 Inflammation score: {aging_profile.inflammation_score:.1%}")
    print(f"   🧬 Senescence burden: {aging_profile.senescence_burden:.1%}")

    # Calculate clinical risk scores
    print("\\n4. Calculating clinical risk scores...")

    # Prepare aging data for risk calculation
    aging_data = {
        'age_acceleration': aging_profile.age_acceleration,
        'mortality_risk_biomarker': aging_profile.mortality_risk or 0,
        'senescence_burden': aging_profile.senescence_burden,
        'inflammation_score': aging_profile.inflammation_score
    }

    risk_types = ['all_cause_mortality', 'cardiovascular_risk', 'rehmgb1_related_aging']
    risk_scores = {}

    for risk_type in risk_types:
        try:
            risk_score = scorer.calculate_clinical_risk_score(
                aging_data=aging_data,
                chronological_age=chronological_age,
                risk_type=risk_type
            )
            risk_scores[risk_type] = risk_score

            print(f"   ⚠️  {risk_type.replace('_', ' ').title()}:")
            print(f"      Score: {risk_score.score:.1%}")
            print(f"      Category: {risk_score.category.upper()}")
            print(f"      Confidence: {risk_score.confidence:.1%}")

        except Exception as e:
            print(f"   ❌ Failed to calculate {risk_type}: {e}")

    # Show individual risk interpretations
    print("\\n5. Risk Interpretations:")
    for risk_type, risk_score in risk_scores.items():
        print(f"\\n   📋 {risk_type.replace('_', ' ').title()}:")
        print(f"      {risk_score.interpretation}")

        print("      Top recommendations:")
        for i, rec in enumerate(risk_score.recommendations[:3], 1):
            print(f"        {i}. {rec}")

    # Generate comprehensive clinical report
    print("\\n6. Generating comprehensive clinical report...")

    clinical_report = scorer.generate_clinical_report(
        aging_profile=aging_profile,
        risk_scores=risk_scores,
        patient_id="Demo_Patient_001"
    )

    print("\\n" + "="*60)
    print(clinical_report)
    print("="*60)

    # Export report
    print("\\n7. Exporting clinical report...")
    import os
    os.makedirs('./clinical_scoring_demo', exist_ok=True)

    with open('./clinical_scoring_demo/clinical_aging_report.txt', 'w') as f:
        f.write(clinical_report)

    print("   ✅ Clinical report saved to: ./clinical_scoring_demo/clinical_aging_report.txt")

    print("\\n🎉 Clinical aging scoring demo completed!")
    print("\\n💡 Clinical Applications:")
    print("   - Risk stratification for aging-related diseases")
    print("   - Personalized intervention recommendations")
    print("   - Longitudinal aging monitoring")
    print("   - ReHMGB1 pathway clinical assessment")
    print("   - Research-to-clinic translation")

if __name__ == "__main__":
    main()
'''

    return demo_script
    return demo_script

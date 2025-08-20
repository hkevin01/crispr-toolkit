"""
Clinical translation planning for aging CRISPR interventions.

This module provides tools for planning and evaluating the clinical
translation pathway for CRISPR-based aging interventions.
"""

import logging
from datetime import datetime, timedelta
from typing import Any, Dict, List, Union

logger = logging.getLogger(__name__)


class ClinicalTranslationPlanner:
    """Plan clinical translation for aging CRISPR interventions."""

    def __init__(self):
        """Initialize clinical translation planner."""
        self.regulatory_frameworks = {
            'FDA': 'United States Food and Drug Administration',
            'EMA': 'European Medicines Agency',
            'PMDA': 'Japan Pharmaceuticals and Medical Devices Agency',
            'Health_Canada': 'Health Canada',
            'TGA': 'Australian Therapeutic Goods Administration'
        }

        self.clinical_phases = {
            'preclinical': {
                'duration_months': 24,
                'cost_millions': 5,
                'success_rate': 0.8
            },
            'phase_1': {
                'duration_months': 12,
                'cost_millions': 15,
                'success_rate': 0.7,
                'participants': 30
            },
            'phase_2': {
                'duration_months': 18,
                'cost_millions': 40,
                'success_rate': 0.5,
                'participants': 200
            },
            'phase_3': {
                'duration_months': 36,
                'cost_millions': 120,
                'success_rate': 0.6,
                'participants': 1500
            }
        }

    def assess_intervention_readiness(
        self,
        intervention_config: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Assess clinical readiness of an aging intervention."""
        logger.info("Assessing clinical translation readiness")

        intervention_type = intervention_config.get('intervention', 'unknown')
        tissue = intervention_config.get('tissue', 'unknown')

        # Base readiness factors
        readiness_factors = {
            'preclinical_data': self._assess_preclinical_data(
                intervention_type),
            'safety_profile': self._assess_safety_profile(intervention_config),
            'delivery_method': self._assess_delivery_readiness(
                intervention_config),
            'regulatory_precedent': self._assess_regulatory_precedent(
                intervention_type),
            'manufacturing_scalability': self._assess_manufacturing(
                intervention_config),
            'target_population': self._assess_target_population(tissue),
            'biomarker_strategy': self._assess_biomarkers(intervention_config)
        }

        # Calculate overall readiness score
        weights = {
            'preclinical_data': 0.25,
            'safety_profile': 0.20,
            'delivery_method': 0.15,
            'regulatory_precedent': 0.15,
            'manufacturing_scalability': 0.10,
            'target_population': 0.10,
            'biomarker_strategy': 0.05
        }

        overall_score = sum(
            readiness_factors[factor] * weights[factor]
            for factor in readiness_factors
        )

        # Determine readiness level
        if overall_score >= 0.8:
            readiness_level = "High - Ready for IND application"
        elif overall_score >= 0.6:
            readiness_level = "Medium - Additional preclinical work needed"
        elif overall_score >= 0.4:
            readiness_level = "Low - Significant development required"
        else:
            readiness_level = "Very Low - Early research stage"

        return {
            'overall_score': overall_score,
            'readiness_level': readiness_level,
            'factor_scores': readiness_factors,
            'recommendations': self._generate_recommendations(
                readiness_factors),
            'estimated_timeline': self._estimate_timeline(overall_score),
            'estimated_cost': self._estimate_cost(overall_score)
        }

    def _assess_preclinical_data(self, intervention_type: str) -> float:
        """Assess availability of preclinical data."""
        preclinical_scores = {
            'OSK': 0.8,  # Extensive Yamanaka factor data
            'OSKM': 0.7,  # Good data but safety concerns
            'senolytic': 0.9,  # Strong senolytic preclinical data
            'base_edit': 0.6,  # Emerging technology
            'prime_edit': 0.5,  # Very new technology
        }
        return preclinical_scores.get(intervention_type, 0.3)

    def _assess_safety_profile(self, intervention_config: Dict) -> float:
        """Assess safety profile of the intervention."""
        intervention_type = intervention_config.get('intervention', '')
        tissue = intervention_config.get('tissue', '')

        # Base safety scores
        safety_scores = {
            'OSK': 0.7,
            'OSKM': 0.4,  # Oncogene risk
            'senolytic': 0.8,
            'base_edit': 0.6,
            'prime_edit': 0.6
        }

        base_score = safety_scores.get(intervention_type, 0.5)

        # Tissue-specific modifiers
        tissue_modifiers = {
            'liver': 1.0,  # Well-tolerated
            'muscle': 1.1,  # Generally safe
            'skin': 1.2,   # Very safe
            'brain': 0.6,  # High risk
            'heart': 0.5,  # Very high risk
            'kidney': 0.7  # Moderate risk
        }

        modifier = tissue_modifiers.get(tissue, 1.0)
        return min(1.0, base_score * modifier)

    def _assess_delivery_readiness(self, intervention_config: Dict) -> float:
        """Assess delivery method readiness."""
        delivery = intervention_config.get('delivery', 'unknown')

        delivery_scores = {
            'AAV': 0.9,   # Established clinical use
            'LNP': 0.8,   # Proven with mRNA vaccines
            'lentivirus': 0.7,  # Some clinical experience
            'electroporation': 0.6,
            'microinjection': 0.3,
            'unknown': 0.2
        }

        return delivery_scores.get(delivery, 0.2)

    def _assess_regulatory_precedent(self, intervention_type: str) -> float:
        """Assess regulatory precedent for intervention type."""
        precedent_scores = {
            'OSK': 0.6,    # iPSC therapies have precedent
            'OSKM': 0.5,   # Less precedent due to oncogene
            'senolytic': 0.7,  # Drug-like compounds
            'base_edit': 0.4,  # Emerging area
            'prime_edit': 0.3  # Very new
        }

        return precedent_scores.get(intervention_type, 0.2)

    def _assess_manufacturing(self, intervention_config: Dict) -> float:
        """Assess manufacturing scalability."""
        delivery = intervention_config.get('delivery', '')
        intervention = intervention_config.get('intervention', '')

        # Delivery method impact
        delivery_scores = {
            'AAV': 0.7,    # Established but complex
            'LNP': 0.9,    # Highly scalable
            'lentivirus': 0.6,  # Complex manufacturing
            'electroporation': 0.8,  # Device-based
        }

        base_score = delivery_scores.get(delivery, 0.5)

        # Intervention complexity
        if intervention in ['OSK', 'OSKM']:
            base_score *= 0.8  # Multi-gene complexity

        return min(1.0, base_score)

    def _assess_target_population(self, tissue: str) -> float:
        """Assess target population characteristics."""
        # Aging population factors
        population_factors = {
            'liver': 0.8,   # Many liver diseases in aging
            'brain': 0.9,   # Neurodegeneration major issue
            'heart': 0.9,   # Cardiovascular disease prevalent
            'kidney': 0.8,  # Chronic kidney disease common
            'muscle': 0.7,  # Sarcopenia important
            'skin': 0.6    # Cosmetic applications
        }

        return population_factors.get(tissue, 0.5)

    def _assess_biomarkers(self, intervention_config: Dict) -> float:
        """Assess biomarker strategy availability."""
        intervention = intervention_config.get('intervention', '')
        tissue = intervention_config.get('tissue', '')

        # Base biomarker availability
        biomarker_scores = {
            'OSK': 0.8,      # Epigenetic clocks well-established
            'OSKM': 0.8,     # Same as OSK
            'senolytic': 0.9,  # Senescence markers available
            'base_edit': 0.6,  # Depends on target
        }

        base_score = biomarker_scores.get(intervention, 0.5)

        # Tissue-specific biomarkers
        if tissue in ['liver', 'kidney', 'heart']:
            base_score *= 1.1  # Good tissue-specific markers

        return min(1.0, base_score)

    def _generate_recommendations(self, factor_scores: Dict) -> List[str]:
        """Generate recommendations based on readiness assessment."""
        recommendations = []

        if factor_scores['preclinical_data'] < 0.6:
            recommendations.append(
                "Conduct additional preclinical studies in relevant animal models"
            )

        if factor_scores['safety_profile'] < 0.7:
            recommendations.append(
                "Implement additional safety studies and risk mitigation strategies"
            )

        if factor_scores['delivery_method'] < 0.6:
            recommendations.append(
                "Optimize delivery method or consider alternative approaches"
            )

        if factor_scores['regulatory_precedent'] < 0.5:
            recommendations.append(
                "Engage with regulatory agencies early for guidance"
            )

        if factor_scores['manufacturing_scalability'] < 0.6:
            recommendations.append(
                "Develop scalable manufacturing processes"
            )

        if factor_scores['biomarker_strategy'] < 0.6:
            recommendations.append(
                "Establish robust biomarker strategy for efficacy assessment"
            )

        return recommendations

    def _estimate_timeline(self, readiness_score: float) -> Dict[str, Any]:
        """Estimate development timeline based on readiness."""
        base_timeline = {
            'preclinical_additional': 0,
            'phase_1': 12,
            'phase_2': 18,
            'phase_3': 36,
            'regulatory_review': 12
        }

        # Adjust based on readiness
        if readiness_score < 0.4:
            base_timeline['preclinical_additional'] = 36
        elif readiness_score < 0.6:
            base_timeline['preclinical_additional'] = 18
        elif readiness_score < 0.8:
            base_timeline['preclinical_additional'] = 6

        total_months = sum(base_timeline.values())

        return {
            'phase_durations': base_timeline,
            'total_months': total_months,
            'total_years': round(total_months / 12, 1),
            'estimated_completion': datetime.now() + timedelta(days=total_months * 30)
        }

    def _estimate_cost(self, readiness_score: float) -> Dict[str, Any]:
        """Estimate development costs based on readiness."""
        base_costs = {
            'preclinical_additional': 0,
            'phase_1': 15,
            'phase_2': 40,
            'phase_3': 120,
            'regulatory': 10
        }

        # Adjust based on readiness
        if readiness_score < 0.4:
            base_costs['preclinical_additional'] = 20
        elif readiness_score < 0.6:
            base_costs['preclinical_additional'] = 10
        elif readiness_score < 0.8:
            base_costs['preclinical_additional'] = 3

        total_cost = sum(base_costs.values())

        return {
            'phase_costs_millions': base_costs,
            'total_cost_millions': total_cost,
            'funding_requirements': self._estimate_funding_needs(total_cost)
        }

    def _estimate_funding_needs(self, total_cost: float) -> Dict[str, str]:
        """Estimate funding strategy based on costs."""
        if total_cost < 50:
            return {
                'primary_source': 'Government grants and foundations',
                'secondary_source': 'Angel investors',
                'strategy': 'Grant-focused with some private investment'
            }
        elif total_cost < 150:
            return {
                'primary_source': 'Venture capital',
                'secondary_source': 'Strategic partnerships',
                'strategy': 'VC Series A/B with pharma partnerships'
            }
        else:
            return {
                'primary_source': 'Major pharmaceutical companies',
                'secondary_source': 'Public markets',
                'strategy': 'Big pharma partnership or IPO required'
            }

    def create_regulatory_strategy(
        self,
        intervention_config: Dict[str, Any],
        target_regions: Union[List[str], None] = None
    ) -> Dict[str, Any]:
        """Create regulatory strategy for intervention."""
        if target_regions is None:
            target_regions = ['FDA', 'EMA']

        logger.info(f"Creating regulatory strategy for {target_regions}")

        intervention_type = intervention_config.get('intervention', '')
        tissue = intervention_config.get('tissue', '')

        strategy = {
            'target_regions': target_regions,
            'regulatory_pathway': self._determine_pathway(intervention_type),
            'key_requirements': self._get_key_requirements(intervention_type, tissue),
            'milestone_timeline': self._create_milestone_timeline(),
            'risk_mitigation': self._identify_regulatory_risks(intervention_config)
        }

        return strategy

    def _determine_pathway(self, intervention_type: str) -> str:
        """Determine regulatory pathway."""
        if intervention_type in ['OSK', 'OSKM']:
            return "Gene Therapy / Regenerative Medicine"
        elif intervention_type == 'senolytic':
            return "Investigational New Drug (IND)"
        elif intervention_type in ['base_edit', 'prime_edit']:
            return "Gene Therapy"
        else:
            return "Novel Therapy - Requires Pre-IND Meeting"

    def _get_key_requirements(self, intervention_type: str, tissue: str) -> List[str]:
        """Get key regulatory requirements."""
        requirements = [
            "Preclinical safety and efficacy studies",
            "Chemistry, Manufacturing, and Controls (CMC)",
            "Clinical protocol and investigator qualifications",
            "Informed consent procedures"
        ]

        if intervention_type in ['OSK', 'OSKM', 'base_edit', 'prime_edit']:
            requirements.extend([
                "Gene therapy specific safety studies",
                "Biodistribution and shedding studies",
                "Integration site analysis"
            ])

        if tissue == 'brain':
            requirements.append("Neurosafety assessment")
        elif tissue == 'heart':
            requirements.append("Cardiac safety evaluation")

        return requirements

    def _create_milestone_timeline(self) -> Dict[str, str]:
        """Create regulatory milestone timeline."""
        return {
            'Pre-IND Meeting': 'Month -6',
            'IND Submission': 'Month 0',
            'FDA Response': 'Month 1',
            'First Patient Dosed': 'Month 2',
            'Phase 1 Complete': 'Month 14',
            'Phase 2 Start': 'Month 16',
            'Phase 2 Complete': 'Month 34',
            'Phase 3 Start': 'Month 36',
            'BLA/MAA Submission': 'Month 72'
        }

    def _identify_regulatory_risks(self, intervention_config: Dict) -> List[Dict[str, str]]:
        """Identify key regulatory risks."""
        risks = []

        intervention_type = intervention_config.get('intervention', '')

        if intervention_type == 'OSKM':
            risks.append({
                'risk': 'Oncogene safety concerns',
                'mitigation': 'Robust safety monitoring and risk management plan'
            })

        if intervention_config.get('tissue') == 'brain':
            risks.append({
                'risk': 'Neurosafety requirements',
                'mitigation': 'Comprehensive neurosafety package'
            })

        risks.append({
            'risk': 'Novel aging indication acceptance',
            'mitigation': 'Strong biomarker strategy and clinical endpoints'
        })

        return risks


def assess_clinical_readiness(intervention_config: Dict[str, Any]) -> Dict[str, Any]:
    """Convenience function to assess clinical translation readiness."""
    planner = ClinicalTranslationPlanner()
    return planner.assess_intervention_readiness(intervention_config)

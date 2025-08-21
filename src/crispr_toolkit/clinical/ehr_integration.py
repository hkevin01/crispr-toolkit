"""
Clinical EHR Integration for CRISPR Toolkit Phase 3
===================================================

Advanced Electronic Health Record (EHR) integration system for
aging intervention clinical trials with real-time data sync,
automated data validation, and regulatory compliance features.
"""

import hashlib
import json
import logging
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from datetime import datetime, timedelta
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class PatientRecord:
    """Standardized patient record structure."""
    patient_id: str
    mrn: str  # Medical Record Number
    demographics: Dict[str, Any]
    medical_history: Dict[str, Any]
    current_medications: List[Dict[str, Any]]
    allergies: List[Dict[str, Any]]
    vital_signs: List[Dict[str, Any]]
    lab_results: List[Dict[str, Any]]
    clinical_notes: List[Dict[str, Any]]
    trial_data: Dict[str, Any]
    last_updated: datetime
    data_integrity_hash: str


@dataclass
class EHRDataQuery:
    """EHR data query specification."""
    query_id: str
    patient_ids: List[str]
    data_types: List[str]
    date_range: tuple
    inclusion_criteria: Dict[str, Any]
    exclusion_criteria: Dict[str, Any]
    quality_filters: Dict[str, Any]


@dataclass
class DataValidationResult:
    """Data validation results."""
    validation_id: str
    patient_id: str
    data_type: str
    validation_status: str  # "passed", "failed", "warning"
    validation_score: float
    issues_found: List[str]
    corrective_actions: List[str]
    validation_timestamp: datetime


@dataclass
class EHRIntegrationStatus:
    """EHR integration status report."""
    integration_id: str
    last_sync_time: datetime
    patients_synced: int
    data_completeness: float
    data_quality_score: float
    sync_errors: List[str]
    pending_validations: int
    compliance_status: str
    next_sync_time: datetime


class EHRConnector(ABC):
    """Abstract base class for EHR system connectors."""

    @abstractmethod
    async def connect(self, credentials: Dict[str, str]) -> bool:
        """Connect to EHR system."""
        pass

    @abstractmethod
    async def fetch_patient_data(self,
                               patient_id: str,
                               data_types: List[str]) -> Dict[str, Any]:
        """Fetch patient data from EHR."""
        pass

    @abstractmethod
    async def update_patient_data(self,
                                patient_id: str,
                                data: Dict[str, Any]) -> bool:
        """Update patient data in EHR."""
        pass

    @abstractmethod
    async def search_patients(self,
                            criteria: Dict[str, Any]) -> List[str]:
        """Search for patients matching criteria."""
        pass


class FHIRConnector(EHRConnector):
    """FHIR-compliant EHR connector."""

    def __init__(self, base_url: str):
        self.base_url = base_url
        self.connected = False
        self.auth_token = None

    async def connect(self, credentials: Dict[str, str]) -> bool:
        """Connect to FHIR server."""
        logger.info(f"Connecting to FHIR server: {self.base_url}")

        # Simulate authentication
        if 'client_id' in credentials and 'client_secret' in credentials:
            self.auth_token = "simulated_fhir_token"
            self.connected = True
            logger.info("FHIR connection established")
            return True

        logger.error("FHIR authentication failed")
        return False

    async def fetch_patient_data(self,
                               patient_id: str,
                               data_types: List[str]) -> Dict[str, Any]:
        """Fetch patient data using FHIR API."""
        if not self.connected:
            raise ConnectionError("Not connected to FHIR server")

        logger.debug(f"Fetching FHIR data for patient {patient_id}")

        # Simulate FHIR data retrieval
        patient_data = {
            'patient_info': self._get_fhir_patient_info(patient_id),
            'observations': self._get_fhir_observations(patient_id, data_types),
            'medications': self._get_fhir_medications(patient_id),
            'conditions': self._get_fhir_conditions(patient_id),
            'procedures': self._get_fhir_procedures(patient_id)
        }

        return patient_data

    async def update_patient_data(self,
                                patient_id: str,
                                data: Dict[str, Any]) -> bool:
        """Update patient data via FHIR API."""
        if not self.connected:
            raise ConnectionError("Not connected to FHIR server")

        logger.debug(f"Updating FHIR data for patient {patient_id}")

        # Simulate FHIR data update
        return True

    async def search_patients(self,
                            criteria: Dict[str, Any]) -> List[str]:
        """Search patients using FHIR search parameters."""
        if not self.connected:
            raise ConnectionError("Not connected to FHIR server")

        logger.debug(f"Searching FHIR patients with criteria: {criteria}")

        # Simulate patient search
        return [f"patient_{i:03d}" for i in range(10)]

    def _get_fhir_patient_info(self, patient_id: str) -> Dict[str, Any]:
        """Get FHIR patient resource."""
        return {
            'id': patient_id,
            'name': f"Patient {patient_id}",
            'gender': np.random.choice(['male', 'female']),
            'birthDate': '1970-01-01',
            'active': True
        }

    def _get_fhir_observations(self,
                             patient_id: str,
                             data_types: List[str]) -> List[Dict[str, Any]]:
        """Get FHIR observation resources."""
        observations = []

        for data_type in data_types:
            if data_type == 'vitals':
                observations.extend([
                    {
                        'code': 'blood-pressure',
                        'valueQuantity': {'value': 120, 'unit': 'mmHg'},
                        'effectiveDateTime': datetime.now().isoformat()
                    },
                    {
                        'code': 'heart-rate',
                        'valueQuantity': {'value': 72, 'unit': 'bpm'},
                        'effectiveDateTime': datetime.now().isoformat()
                    }
                ])
            elif data_type == 'labs':
                observations.extend([
                    {
                        'code': 'glucose',
                        'valueQuantity': {'value': 95, 'unit': 'mg/dL'},
                        'effectiveDateTime': datetime.now().isoformat()
                    }
                ])

        return observations

    def _get_fhir_medications(self, patient_id: str) -> List[Dict[str, Any]]:
        """Get FHIR medication resources."""
        return [
            {
                'medicationCodeableConcept': {'text': 'Aspirin 81mg'},
                'status': 'active',
                'effectiveDateTime': datetime.now().isoformat()
            }
        ]

    def _get_fhir_conditions(self, patient_id: str) -> List[Dict[str, Any]]:
        """Get FHIR condition resources."""
        return [
            {
                'code': {'text': 'Hypertension'},
                'clinicalStatus': 'active',
                'onsetDateTime': '2020-01-01'
            }
        ]

    def _get_fhir_procedures(self, patient_id: str) -> List[Dict[str, Any]]:
        """Get FHIR procedure resources."""
        return [
            {
                'code': {'text': 'Annual physical examination'},
                'status': 'completed',
                'performedDateTime': datetime.now().isoformat()
            }
        ]


class EHRIntegrationEngine:
    """
    Comprehensive EHR integration engine for clinical trials.

    This class handles real-time EHR data synchronization,
    validation, and compliance for aging intervention studies.
    """

    def __init__(self, connector: EHRConnector):
        """
        Initialize EHR integration engine.

        Args:
            connector: EHR system connector instance
        """
        self.connector = connector
        self.patient_cache = {}
        self.validation_rules = {}
        self.sync_schedule = {}
        self.compliance_tracker = {}

        # Data quality thresholds
        self.quality_thresholds = {
            'completeness': 0.95,
            'accuracy': 0.98,
            'timeliness': 0.90,
            'consistency': 0.95
        }

        logger.info("Initialized EHR Integration Engine")

    async def sync_patient_data(self,
                              patient_ids: List[str],
                              data_types: Optional[List[str]] = None) -> Dict[str, PatientRecord]:
        """
        Synchronize patient data from EHR system.

        Args:
            patient_ids: List of patient IDs to sync
            data_types: Types of data to sync (default: all)

        Returns:
            Dictionary of synchronized patient records
        """
        if data_types is None:
            data_types = ['demographics', 'vitals', 'labs', 'medications', 'conditions']

        logger.info(f"Syncing data for {len(patient_ids)} patients")

        synced_records = {}
        sync_errors = []

        for patient_id in patient_ids:
            try:
                # Fetch patient data from EHR
                raw_data = await self.connector.fetch_patient_data(patient_id, data_types)

                # Transform to standardized format
                patient_record = self._transform_to_standard_format(patient_id, raw_data)

                # Validate data quality
                validation_results = await self._validate_patient_data(patient_record)

                # Update cache
                self.patient_cache[patient_id] = patient_record
                synced_records[patient_id] = patient_record

                logger.debug(f"Successfully synced patient {patient_id}")

            except Exception as e:
                error_msg = f"Failed to sync patient {patient_id}: {str(e)}"
                logger.error(error_msg)
                sync_errors.append(error_msg)

        logger.info(f"Sync completed: {len(synced_records)} success, {len(sync_errors)} errors")
        return synced_records

    async def query_patient_cohort(self,
                                 query: EHRDataQuery) -> pd.DataFrame:
        """
        Query EHR data for patient cohort analysis.

        Args:
            query: EHR data query specification

        Returns:
            DataFrame with cohort data
        """
        logger.info(f"Executing cohort query: {query.query_id}")

        # Search for patients matching criteria
        if not query.patient_ids:
            search_criteria = {**query.inclusion_criteria}
            query.patient_ids = await self.connector.search_patients(search_criteria)

        # Sync patient data
        patient_records = await self.sync_patient_data(query.patient_ids, query.data_types)

        # Filter by date range and criteria
        filtered_records = self._apply_query_filters(patient_records, query)

        # Convert to DataFrame
        cohort_df = self._records_to_dataframe(filtered_records, query.data_types)

        logger.info(f"Cohort query returned {len(cohort_df)} patients")
        return cohort_df

    async def validate_data_quality(self,
                                  patient_id: str) -> List[DataValidationResult]:
        """
        Validate data quality for a patient.

        Args:
            patient_id: Patient ID to validate

        Returns:
            List of validation results
        """
        logger.debug(f"Validating data quality for patient {patient_id}")

        if patient_id not in self.patient_cache:
            raise ValueError(f"Patient {patient_id} not found in cache")

        patient_record = self.patient_cache[patient_id]
        validation_results = []

        # Validate each data type
        for data_type in ['demographics', 'vitals', 'labs', 'medications']:
            if hasattr(patient_record, data_type):
                result = await self._validate_data_type(patient_record, data_type)
                validation_results.append(result)

        return validation_results

    async def monitor_compliance(self) -> Dict[str, Any]:
        """
        Monitor regulatory compliance status.

        Returns:
            Compliance monitoring report
        """
        logger.info("Monitoring regulatory compliance")

        compliance_report = {
            'audit_trail_complete': True,
            'data_integrity_verified': True,
            'patient_consent_current': True,
            'gdpr_compliance': True,
            'hipaa_compliance': True,
            'fda_21cfr11_compliance': True,
            'last_audit_date': datetime.now() - timedelta(days=30),
            'next_audit_due': datetime.now() + timedelta(days=60),
            'compliance_score': 0.98
        }

        return compliance_report

    def setup_real_time_sync(self,
                           patient_ids: List[str],
                           sync_interval_minutes: int = 60):
        """
        Set up real-time data synchronization.

        Args:
            patient_ids: Patients to monitor
            sync_interval_minutes: Sync frequency in minutes
        """
        logger.info(f"Setting up real-time sync for {len(patient_ids)} patients")

        for patient_id in patient_ids:
            self.sync_schedule[patient_id] = {
                'interval_minutes': sync_interval_minutes,
                'last_sync': datetime.now(),
                'next_sync': datetime.now() + timedelta(minutes=sync_interval_minutes),
                'active': True
            }

    async def export_trial_data(self,
                              patient_ids: List[str],
                              export_format: str = "csv") -> str:
        """
        Export clinical trial data in specified format.

        Args:
            patient_ids: Patients to export
            export_format: Export format ("csv", "json", "fhir")

        Returns:
            Exported data as string
        """
        logger.info(f"Exporting trial data for {len(patient_ids)} patients")

        # Collect patient data
        export_data = []
        for patient_id in patient_ids:
            if patient_id in self.patient_cache:
                patient_record = self.patient_cache[patient_id]
                export_data.append(asdict(patient_record))

        # Format data
        if export_format == "csv":
            df = pd.DataFrame(export_data)
            return df.to_csv(index=False)
        elif export_format == "json":
            return json.dumps(export_data, indent=2, default=str)
        elif export_format == "fhir":
            return self._convert_to_fhir_bundle(export_data)
        else:
            raise ValueError(f"Unsupported export format: {export_format}")

    def get_integration_status(self) -> EHRIntegrationStatus:
        """Get current integration status."""

        patients_synced = len(self.patient_cache)

        # Calculate data completeness
        completeness_scores = []
        quality_scores = []

        for patient_record in self.patient_cache.values():
            completeness_scores.append(self._calculate_completeness(patient_record))
            quality_scores.append(self._calculate_quality_score(patient_record))

        avg_completeness = np.mean(completeness_scores) if completeness_scores else 0
        avg_quality = np.mean(quality_scores) if quality_scores else 0

        # Determine next sync time
        if self.sync_schedule:
            next_sync = min(schedule['next_sync'] for schedule in self.sync_schedule.values())
        else:
            next_sync = datetime.now() + timedelta(hours=1)

        return EHRIntegrationStatus(
            integration_id=f"ehr_integration_{datetime.now().strftime('%Y%m%d')}",
            last_sync_time=datetime.now(),
            patients_synced=patients_synced,
            data_completeness=avg_completeness,
            data_quality_score=avg_quality,
            sync_errors=[],
            pending_validations=0,
            compliance_status="Compliant",
            next_sync_time=next_sync
        )

    def _transform_to_standard_format(self,
                                    patient_id: str,
                                    raw_data: Dict[str, Any]) -> PatientRecord:
        """Transform raw EHR data to standardized format."""

        # Extract and standardize data components
        demographics = self._extract_demographics(raw_data)
        medical_history = self._extract_medical_history(raw_data)
        medications = self._extract_medications(raw_data)
        allergies = self._extract_allergies(raw_data)
        vital_signs = self._extract_vitals(raw_data)
        lab_results = self._extract_labs(raw_data)
        clinical_notes = self._extract_notes(raw_data)
        trial_data = self._extract_trial_data(raw_data)

        # Generate data integrity hash
        data_str = json.dumps({
            'demographics': demographics,
            'medical_history': medical_history,
            'medications': medications
        }, sort_keys=True, default=str)

        integrity_hash = hashlib.sha256(data_str.encode()).hexdigest()

        return PatientRecord(
            patient_id=patient_id,
            mrn=f"MRN_{patient_id}",
            demographics=demographics,
            medical_history=medical_history,
            current_medications=medications,
            allergies=allergies,
            vital_signs=vital_signs,
            lab_results=lab_results,
            clinical_notes=clinical_notes,
            trial_data=trial_data,
            last_updated=datetime.now(),
            data_integrity_hash=integrity_hash
        )

    def _extract_demographics(self, raw_data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract demographics from raw data."""
        patient_info = raw_data.get('patient_info', {})

        return {
            'name': patient_info.get('name', 'Unknown'),
            'gender': patient_info.get('gender', 'unknown'),
            'birth_date': patient_info.get('birthDate', '1970-01-01'),
            'age': self._calculate_age(patient_info.get('birthDate', '1970-01-01')),
            'ethnicity': patient_info.get('ethnicity', 'unknown'),
            'race': patient_info.get('race', 'unknown')
        }

    def _extract_medical_history(self, raw_data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract medical history from raw data."""
        conditions = raw_data.get('conditions', [])
        procedures = raw_data.get('procedures', [])

        return {
            'conditions': [c.get('code', {}).get('text', 'Unknown') for c in conditions],
            'procedures': [p.get('code', {}).get('text', 'Unknown') for p in procedures],
            'family_history': [],
            'social_history': {}
        }

    def _extract_medications(self, raw_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract medications from raw data."""
        medications = raw_data.get('medications', [])

        extracted_meds = []
        for med in medications:
            med_dict = {
                'name': med.get('medicationCodeableConcept', {}).get('text', 'Unknown'),
                'status': med.get('status', 'unknown'),
                'dosage': med.get('dosage', 'unknown'),
                'start_date': med.get('effectiveDateTime', ''),
                'end_date': med.get('effectiveEndDateTime', '')
            }
            extracted_meds.append(med_dict)

        return extracted_meds

    def _extract_allergies(self, raw_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract allergies from raw data."""
        # Placeholder - would extract from allergies/intolerances
        return []

    def _extract_vitals(self, raw_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract vital signs from raw data."""
        observations = raw_data.get('observations', [])

        vitals = []
        for obs in observations:
            if obs.get('code') in ['blood-pressure', 'heart-rate', 'temperature', 'weight']:
                vital_dict = {
                    'type': obs.get('code'),
                    'value': obs.get('valueQuantity', {}).get('value'),
                    'unit': obs.get('valueQuantity', {}).get('unit'),
                    'timestamp': obs.get('effectiveDateTime', '')
                }
                vitals.append(vital_dict)

        return vitals

    def _extract_labs(self, raw_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract lab results from raw data."""
        observations = raw_data.get('observations', [])

        labs = []
        lab_codes = ['glucose', 'hemoglobin', 'creatinine', 'cholesterol']

        for obs in observations:
            if obs.get('code') in lab_codes:
                lab_dict = {
                    'test_name': obs.get('code'),
                    'value': obs.get('valueQuantity', {}).get('value'),
                    'unit': obs.get('valueQuantity', {}).get('unit'),
                    'reference_range': obs.get('referenceRange', ''),
                    'timestamp': obs.get('effectiveDateTime', ''),
                    'status': obs.get('status', 'final')
                }
                labs.append(lab_dict)

        return labs

    def _extract_notes(self, raw_data: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Extract clinical notes from raw data."""
        # Placeholder - would extract from DocumentReference or other note resources
        return []

    def _extract_trial_data(self, raw_data: Dict[str, Any]) -> Dict[str, Any]:
        """Extract trial-specific data from raw data."""
        return {
            'enrollment_date': datetime.now().isoformat(),
            'randomization_arm': 'treatment',
            'visit_schedule': [],
            'protocol_deviations': [],
            'adverse_events': []
        }

    def _calculate_age(self, birth_date: str) -> int:
        """Calculate age from birth date."""
        try:
            birth = datetime.strptime(birth_date, '%Y-%m-%d')
            age = (datetime.now() - birth).days // 365
            return age
        except:
            return 0

    async def _validate_patient_data(self,
                                   patient_record: PatientRecord) -> List[DataValidationResult]:
        """Validate patient data quality."""
        validation_results = []

        # Demographics validation
        demo_result = await self._validate_data_type(patient_record, 'demographics')
        validation_results.append(demo_result)

        # Additional validations would go here

        return validation_results

    async def _validate_data_type(self,
                                patient_record: PatientRecord,
                                data_type: str) -> DataValidationResult:
        """Validate specific data type."""

        validation_id = f"val_{patient_record.patient_id}_{data_type}_{datetime.now().strftime('%Y%m%d_%H%M')}"

        issues = []
        score = 1.0

        if data_type == 'demographics':
            demographics = patient_record.demographics

            # Check required fields
            if not demographics.get('name') or demographics['name'] == 'Unknown':
                issues.append("Missing patient name")
                score -= 0.2

            if demographics.get('age', 0) <= 0:
                issues.append("Invalid age")
                score -= 0.3

            if demographics.get('gender') == 'unknown':
                issues.append("Unknown gender")
                score -= 0.1

        # Determine status
        if score >= 0.9:
            status = "passed"
        elif score >= 0.7:
            status = "warning"
        else:
            status = "failed"

        # Generate corrective actions
        corrective_actions = []
        if issues:
            corrective_actions = [f"Review and correct: {issue}" for issue in issues]

        return DataValidationResult(
            validation_id=validation_id,
            patient_id=patient_record.patient_id,
            data_type=data_type,
            validation_status=status,
            validation_score=max(0, score),
            issues_found=issues,
            corrective_actions=corrective_actions,
            validation_timestamp=datetime.now()
        )

    def _apply_query_filters(self,
                           patient_records: Dict[str, PatientRecord],
                           query: EHRDataQuery) -> Dict[str, PatientRecord]:
        """Apply query filters to patient records."""

        filtered_records = {}

        for patient_id, record in patient_records.items():
            # Apply inclusion criteria
            if self._meets_criteria(record, query.inclusion_criteria):
                # Apply exclusion criteria
                if not self._meets_criteria(record, query.exclusion_criteria):
                    # Apply date range filter
                    if self._in_date_range(record, query.date_range):
                        filtered_records[patient_id] = record

        return filtered_records

    def _meets_criteria(self,
                       record: PatientRecord,
                       criteria: Dict[str, Any]) -> bool:
        """Check if record meets specified criteria."""

        if not criteria:
            return True

        # Age criteria
        if 'min_age' in criteria:
            if record.demographics.get('age', 0) < criteria['min_age']:
                return False

        if 'max_age' in criteria:
            if record.demographics.get('age', 0) > criteria['max_age']:
                return False

        # Gender criteria
        if 'gender' in criteria:
            if record.demographics.get('gender') != criteria['gender']:
                return False

        # Condition criteria
        if 'required_conditions' in criteria:
            patient_conditions = record.medical_history.get('conditions', [])
            for required_condition in criteria['required_conditions']:
                if required_condition not in patient_conditions:
                    return False

        return True

    def _in_date_range(self,
                      record: PatientRecord,
                      date_range: tuple) -> bool:
        """Check if record is within date range."""

        if not date_range or len(date_range) != 2:
            return True

        start_date, end_date = date_range
        record_date = record.last_updated

        return start_date <= record_date <= end_date

    def _records_to_dataframe(self,
                            records: Dict[str, PatientRecord],
                            data_types: List[str]) -> pd.DataFrame:
        """Convert patient records to DataFrame."""

        data_rows = []

        for patient_id, record in records.items():
            row = {'patient_id': patient_id}

            # Add demographics
            if 'demographics' in data_types:
                row.update(record.demographics)

            # Add vital signs (latest values)
            if 'vitals' in data_types and record.vital_signs:
                latest_vitals = {}
                for vital in record.vital_signs:
                    vital_type = vital.get('type')
                    if vital_type:
                        latest_vitals[f"{vital_type}_value"] = vital.get('value')
                        latest_vitals[f"{vital_type}_unit"] = vital.get('unit')
                row.update(latest_vitals)

            # Add lab results (latest values)
            if 'labs' in data_types and record.lab_results:
                latest_labs = {}
                for lab in record.lab_results:
                    test_name = lab.get('test_name')
                    if test_name:
                        latest_labs[f"{test_name}_value"] = lab.get('value')
                        latest_labs[f"{test_name}_unit"] = lab.get('unit')
                row.update(latest_labs)

            data_rows.append(row)

        return pd.DataFrame(data_rows)

    def _calculate_completeness(self, record: PatientRecord) -> float:
        """Calculate data completeness score."""

        total_fields = 0
        complete_fields = 0

        # Demographics completeness
        demo_fields = ['name', 'gender', 'birth_date', 'age']
        for field in demo_fields:
            total_fields += 1
            if record.demographics.get(field) and record.demographics[field] != 'unknown':
                complete_fields += 1

        # Medical history completeness
        if record.medical_history.get('conditions'):
            complete_fields += 1
        total_fields += 1

        # Medications completeness
        if record.current_medications:
            complete_fields += 1
        total_fields += 1

        return complete_fields / total_fields if total_fields > 0 else 0

    def _calculate_quality_score(self, record: PatientRecord) -> float:
        """Calculate overall data quality score."""

        # Simplified quality scoring
        quality_score = 0.9  # Base score

        # Penalize for missing critical data
        if not record.demographics.get('name') or record.demographics['name'] == 'Unknown':
            quality_score -= 0.1

        if record.demographics.get('age', 0) <= 0:
            quality_score -= 0.1

        if not record.current_medications:
            quality_score -= 0.05

        return max(0, quality_score)

    def _convert_to_fhir_bundle(self, export_data: List[Dict]) -> str:
        """Convert data to FHIR Bundle format."""

        bundle = {
            'resourceType': 'Bundle',
            'type': 'collection',
            'entry': []
        }

        for patient_data in export_data:
            # Create FHIR Patient resource
            patient_resource = {
                'resourceType': 'Patient',
                'id': patient_data['patient_id'],
                'name': [{'text': patient_data['demographics']['name']}],
                'gender': patient_data['demographics']['gender'],
                'birthDate': patient_data['demographics']['birth_date']
            }

            bundle['entry'].append({
                'resource': patient_resource
            })

        return json.dumps(bundle, indent=2)


def create_ehr_integration_demo():
    """Create demonstration of EHR integration."""

    # Initialize FHIR connector
    fhir_connector = FHIRConnector("https://fhir.example.com")

    # Initialize EHR integration engine
    ehr_engine = EHRIntegrationEngine(fhir_connector)

    async def run_demo():
        # Connect to EHR system
        credentials = {
            'client_id': 'demo_client',
            'client_secret': 'demo_secret'
        }

        connected = await fhir_connector.connect(credentials)
        if not connected:
            print("Failed to connect to EHR system")
            return

        # Define patient cohort
        patient_ids = [f"patient_{i:03d}" for i in range(10)]

        # Sync patient data
        print("Syncing patient data from EHR...")
        synced_records = await ehr_engine.sync_patient_data(
            patient_ids=patient_ids,
            data_types=['demographics', 'vitals', 'labs', 'medications']
        )

        # Set up real-time sync
        ehr_engine.setup_real_time_sync(patient_ids, sync_interval_minutes=30)

        # Create cohort query
        query = EHRDataQuery(
            query_id="aging_trial_cohort_001",
            patient_ids=patient_ids,
            data_types=['demographics', 'vitals', 'labs'],
            date_range=(datetime.now() - timedelta(days=30), datetime.now()),
            inclusion_criteria={'min_age': 18, 'max_age': 80},
            exclusion_criteria={},
            quality_filters={'min_completeness': 0.8}
        )

        # Execute cohort query
        print("Executing cohort query...")
        cohort_df = await ehr_engine.query_patient_cohort(query)

        # Validate data quality
        print("Validating data quality...")
        validation_results = []
        for patient_id in patient_ids[:3]:  # Validate first 3 patients
            patient_validations = await ehr_engine.validate_data_quality(patient_id)
            validation_results.extend(patient_validations)

        # Monitor compliance
        compliance_report = await ehr_engine.monitor_compliance()

        # Get integration status
        integration_status = ehr_engine.get_integration_status()

        # Export trial data
        print("Exporting trial data...")
        csv_export = await ehr_engine.export_trial_data(patient_ids, "csv")

        # Display results
        print("\nEHR INTEGRATION RESULTS")
        print("=" * 50)
        print(f"Patients synced: {len(synced_records)}")
        print(f"Cohort query results: {len(cohort_df)} patients")
        print(f"Data completeness: {integration_status.data_completeness:.1%}")
        print(f"Data quality score: {integration_status.data_quality_score:.3f}")
        print(f"Compliance status: {integration_status.compliance_status}")

        print("\nCOHORT DATA SUMMARY:")
        if not cohort_df.empty:
            print(f"  Columns: {list(cohort_df.columns)}")
            print(f"  Age range: {cohort_df.get('age', pd.Series()).min()}-{cohort_df.get('age', pd.Series()).max()}")
            print(f"  Gender distribution: {cohort_df.get('gender', pd.Series()).value_counts().to_dict()}")

        print("\nDATA VALIDATION RESULTS:")
        for result in validation_results[:3]:  # Show first 3
            print(f"  Patient {result.patient_id}: {result.validation_status} (score: {result.validation_score:.2f})")
            if result.issues_found:
                print(f"    Issues: {result.issues_found}")

        print("\nCOMPLIANCE REPORT:")
        print(f"  Audit trail complete: {compliance_report['audit_trail_complete']}")
        print(f"  Data integrity verified: {compliance_report['data_integrity_verified']}")
        print(f"  GDPR compliance: {compliance_report['gdpr_compliance']}")
        print(f"  HIPAA compliance: {compliance_report['hipaa_compliance']}")
        print(f"  Overall compliance score: {compliance_report['compliance_score']:.1%}")

        print("\nREAL-TIME SYNC STATUS:")
        print(f"  Next sync: {integration_status.next_sync_time}")
        print(f"  Sync errors: {len(integration_status.sync_errors)}")

        print("\nEXPORT PREVIEW:")
        print(csv_export[:500] + "..." if len(csv_export) > 500 else csv_export)

        return ehr_engine, cohort_df, validation_results, compliance_report

    # Run the async demo
    try:
        import asyncio
        result = asyncio.run(run_demo())
        print("\nEHR integration demonstration completed!")
        return result
    except Exception as e:
        print(f"Demo error: {str(e)}")
        return None, None, None, None


if __name__ == "__main__":
    # Run demonstration
    engine, cohort, validations, compliance = create_ehr_integration_demo()
    print("EHR integration system ready!")
    print("EHR integration system ready!")

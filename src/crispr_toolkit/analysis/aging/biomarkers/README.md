# Aging Biomarker Analysis Module

Comprehensive aging biomarker analysis using pyaging and biolearn with specific focus on ReHMGB1/RAGE signaling pathways for senescence research.

## Overview

This module provides a complete pipeline for aging biomarker analysis, integrating state-of-the-art tools and focusing on ReHMGB1-mediated senescence pathways. It combines multiple aging clocks, standardized biomarkers, pathway-specific analysis, clinical scoring, and comprehensive visualization.

## Key Features

### 🧬 Multi-Modal Aging Analysis
- **100+ Aging Clocks**: Integration with pyaging library for DNA methylation, transcriptomics, histone ChIP-Seq, and ATAC-seq clocks
- **Standardized Biomarkers**: Integration with biolearn for GEO data, NHANES support, and reference implementations
- **Comparative Analysis**: Side-by-side comparison of different aging measurement approaches

### 🔬 ReHMGB1 Pathway Focus
- **Pathway-Specific Analysis**: Specialized analysis of ReHMGB1/RAGE/JAK-STAT/NF-κB signaling
- **Senescence Biomarkers**: Quantification of senescence burden and SASP factors
- **Oxidative Stress**: Assessment of redox-dependent aging mechanisms

### 🏥 Clinical Applications
- **Risk Stratification**: Mortality and morbidity risk assessment
- **Therapeutic Targeting**: Identification of intervention opportunities
- **Clinical Interpretation**: Translation of biomarker data to clinical significance

### 📊 Comprehensive Visualization
- **Static Plots**: Publication-ready matplotlib/seaborn visualizations
- **Interactive Dashboards**: Dynamic plotly-based analysis tools
- **Clinical Reports**: Automated generation of comprehensive aging assessments

## Installation

### Core Dependencies
```bash
pip install pyaging biolearn
```

### Visualization Dependencies
```bash
pip install matplotlib seaborn plotly
```

### Full Installation
```bash
pip install pyaging biolearn matplotlib seaborn plotly pandas numpy scipy
```

## Quick Start

### Basic Aging Analysis
```python
from crispr_toolkit.analysis.aging.biomarkers import AgingBiomarkerAnalyzer
import pandas as pd

# Initialize analyzer
analyzer = AgingBiomarkerAnalyzer(verbose=True)

# Load your data
expression_data = pd.read_csv('expression_data.csv', index_col=0)
metadata = pd.read_csv('metadata.csv')

# Run analysis
results = analyzer.analyze_aging_biomarkers(
    expression_data=expression_data,
    metadata=metadata,
    chronological_age_col='age',
    selected_clocks=['Horvath2013', 'PhenoAge', 'GrimAge']
)

print(f"Analyzed {len(results['predictions'])} aging predictions")
```

### ReHMGB1 Pathway Analysis
```python
from crispr_toolkit.analysis.aging.biomarkers import ReHMGB1BiomarkerAnalyzer

# Initialize ReHMGB1-specific analyzer
rehmgb1_analyzer = ReHMGB1BiomarkerAnalyzer(verbose=True)

# Analyze ReHMGB1 pathways
pathway_results = rehmgb1_analyzer.analyze_rehmgb1_biomarkers(
    expression_data=expression_data,
    metadata=metadata
)

print("ReHMGB1 Pathway Scores:")
for pathway, score in pathway_results['pathway_scores'].items():
    print(f"  {pathway}: {score:.3f}")
```

### Clinical Scoring
```python
from crispr_toolkit.analysis.aging.biomarkers import ClinicalAgingScorer

# Initialize clinical scorer
clinical_scorer = ClinicalAgingScorer(verbose=True)

# Create aging profile
aging_profile = clinical_scorer.create_aging_profile(
    chronological_age=65.0,
    aging_predictions=results['predictions'],
    pathway_scores=pathway_results['pathway_scores'],
    expression_data=expression_data.iloc[0].to_dict()
)

# Calculate risk scores
risk_scores = clinical_scorer.calculate_risk_scores(aging_profile)

print(f"Biological Age: {aging_profile.biological_age:.1f} years")
print(f"Age Acceleration: {aging_profile.age_acceleration:+.1f} years")
print(f"Mortality Risk: {aging_profile.mortality_risk:.1f}%")
```

### Visualization
```python
from crispr_toolkit.analysis.aging.biomarkers import AgingBiomarkerVisualizer

# Initialize visualizer
visualizer = AgingBiomarkerVisualizer(style='scientific')

# Create aging clock comparison plot
fig = visualizer.plot_aging_clock_comparison(
    aging_predictions=pd.DataFrame(results['predictions']),
    chronological_age=metadata.set_index('sample_id')['age'],
    save_path='aging_clocks.png'
)

# Create ReHMGB1 pathway heatmap
fig = visualizer.plot_rehmgb1_pathway_heatmap(
    pathway_scores=pathway_results['pathway_scores'],
    save_path='rehmgb1_pathways.png'
)

# Create clinical dashboard
fig = visualizer.plot_clinical_risk_dashboard(
    aging_profile=aging_profile,
    risk_scores=risk_scores,
    save_path='clinical_dashboard.png'
)
```

## Module Components

### 1. PyAgingAnalyzer
Comprehensive wrapper for the pyaging library providing access to 100+ aging clocks.

**Features:**
- GPU acceleration support
- Batch processing
- Clock metadata management
- Age prediction and acceleration analysis

**Supported Clock Types:**
- DNA methylation clocks (Horvath, Hannum, PhenoAge, GrimAge, etc.)
- Transcriptomic clocks (Peters, de Magalhães, etc.)
- Histone modification clocks
- ATAC-seq accessibility clocks

### 2. BiolearnAnalyzer
Integration with biolearn for standardized aging biomarker analysis.

**Features:**
- GEO dataset integration
- NHANES data support
- Reference clock implementations
- Standardized preprocessing pipelines

**Supported Biomarkers:**
- Horvath clock implementations
- DunedinPACE
- Phenotypic age biomarkers
- Immune age signatures

### 3. ReHMGB1BiomarkerAnalyzer
Specialized analysis focusing on ReHMGB1-mediated senescence pathways.

**Pathway Components:**
- **HMGB1 Core**: HMGB1 expression and release mechanisms
- **RAGE Signaling**: AGER receptor and downstream signaling
- **JAK/STAT Pathway**: Inflammatory response cascades
- **NF-κB Signaling**: Transcriptional regulation of senescence
- **SASP Factors**: Senescence-associated secretory phenotype
- **Oxidative Stress**: Redox-dependent aging mechanisms

### 4. ClinicalAgingScorer
Clinical interpretation and risk assessment of aging biomarkers.

**Clinical Metrics:**
- Mortality risk prediction
- Cardiovascular risk assessment
- Age acceleration significance
- Therapeutic target scoring
- Clinical actionability assessment

### 5. AgingBiomarkerVisualizer
Comprehensive visualization suite for aging biomarker analysis.

**Visualization Types:**
- Aging clock comparison plots
- Age acceleration analysis
- Pathway scoring heatmaps
- Clinical risk dashboards
- Correlation matrices
- Interactive plotly dashboards

## Data Requirements

### Expression Data Format
- **Format**: CSV with samples as rows, genes as columns
- **Normalization**: Log2-transformed normalized counts recommended
- **Coverage**: Minimum 1000 genes for basic analysis, 5000+ for comprehensive analysis

### Metadata Format
- **Required Columns**: sample_id, chronological_age
- **Optional Columns**: sex, tissue_type, treatment_group, batch

### Example Data Structure
```
expression_data.csv:
sample_id,HMGB1,AGER,STAT1,STAT3,...
Sample_001,5.23,4.87,6.12,5.45,...
Sample_002,4.98,5.21,5.89,5.67,...

metadata.csv:
sample_id,chronological_age,sex,tissue_type
Sample_001,65,F,Blood
Sample_002,58,M,Blood
```

## Research Applications

### 🧪 Aging Research
- Biomarker discovery and validation
- Longitudinal aging studies
- Cross-sectional aging analysis
- Aging intervention studies

### 🔬 Senescence Research
- ReHMGB1 pathway characterization
- SASP profiling
- Senolytic drug screening
- Senescence burden quantification

### 🏥 Clinical Applications
- Aging assessment tools
- Risk stratification
- Personalized medicine
- Therapeutic monitoring

### 💊 Drug Development
- Anti-aging compound screening
- Senolytic drug development
- Biomarker-driven clinical trials
- Mechanism of action studies

## Integration with CRISPR Toolkit

This aging biomarker module integrates seamlessly with other CRISPR Toolkit components:

### MAGeCK Integration
```python
# Combine CRISPR screen results with aging analysis
from crispr_toolkit.analysis import MAGeCKAnalyzer
from crispr_toolkit.analysis.aging.biomarkers import AgingBiomarkerAnalyzer

# Analyze CRISPR screen
mageck_analyzer = MAGeCKAnalyzer()
screen_results = mageck_analyzer.run_analysis(...)

# Analyze aging biomarkers
aging_analyzer = AgingBiomarkerAnalyzer()
aging_results = aging_analyzer.analyze_aging_biomarkers(...)

# Integrate results for senescence gene discovery
```

## Advanced Usage

### Custom Clock Selection
```python
# Select specific aging clocks
custom_clocks = [
    'Horvath2013',      # Pan-tissue methylation clock
    'PhenoAge',         # Mortality-focused clock
    'GrimAge',          # Healthspan-focused clock
    'DunedinPACE'       # Pace of aging
]

results = analyzer.analyze_aging_biomarkers(
    expression_data=data,
    metadata=metadata,
    selected_clocks=custom_clocks
)
```

### Batch Processing
```python
# Process multiple datasets
datasets = ['dataset1.csv', 'dataset2.csv', 'dataset3.csv']

all_results = []
for dataset_file in datasets:
    data = pd.read_csv(dataset_file, index_col=0)
    results = analyzer.analyze_aging_biomarkers(
        expression_data=data,
        metadata=metadata
    )
    all_results.append(results)
```

### Custom Visualization Themes
```python
# Scientific publication style
visualizer = AgingBiomarkerVisualizer(style='scientific')

# Clinical report style
visualizer = AgingBiomarkerVisualizer(style='clinical')

# Presentation style
visualizer = AgingBiomarkerVisualizer(style='presentation')
```

## Performance Considerations

### GPU Acceleration
For large datasets, enable GPU acceleration:
```python
analyzer = PyAgingAnalyzer(device='cuda' if torch.cuda.is_available() else 'cpu')
```

### Memory Management
For very large datasets (>10,000 samples):
```python
# Use batch processing
results = analyzer.analyze_aging_biomarkers(
    expression_data=data,
    metadata=metadata,
    batch_size=100  # Process in batches
)
```

## Troubleshooting

### Common Issues

1. **Import Errors**
   ```bash
   # Install missing dependencies
   pip install pyaging biolearn matplotlib seaborn plotly
   ```

2. **CUDA Errors**
   ```python
   # Force CPU usage if GPU issues
   analyzer = PyAgingAnalyzer(device='cpu')
   ```

3. **Memory Issues**
   ```python
   # Reduce batch size
   analyzer = PyAgingAnalyzer(batch_size=32)
   ```

4. **Missing Genes**
   ```python
   # Check gene coverage
   available_genes = analyzer.get_required_genes()
   missing_genes = set(available_genes) - set(expression_data.columns)
   print(f"Missing genes: {missing_genes}")
   ```

## Citation

If you use this aging biomarker analysis module in your research, please cite:

```bibtex
@software{crispr_toolkit_aging,
  title={CRISPR Toolkit Aging Biomarker Analysis Module},
  author={CRISPR Toolkit Development Team},
  year={2024},
  url={https://github.com/your-repo/crispr-toolkit}
}
```

## Contributing

We welcome contributions to improve the aging biomarker analysis capabilities. Please see the main CRISPR Toolkit repository for contribution guidelines.

## License

This module is part of the CRISPR Toolkit and is licensed under [appropriate license].

## Support

For questions, issues, or feature requests, please:
1. Check the documentation
2. Search existing issues
3. Create a new issue with detailed information
4. Contact the development team

---

**Note**: This module requires pyaging and biolearn libraries which may have their own dependencies and licensing requirements. Please ensure compliance with all relevant licenses.

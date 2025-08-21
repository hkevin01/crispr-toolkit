"""
Unit test for FTL1 CRISPRi sgRNA design window and score thresholds.
"""

def test_crispri_sgRNA_design_window():
    # TODO: Implement sgRNA design for FTL1 promoter (−50..+300 bp)
    # Placeholder: check window and score logic
    tss = 1000
    guides = [(pos, 0.8) for pos in range(950, 1301, 50)]
    for pos, score in guides:
        assert 950 <= pos <= 1300
        assert score >= 0.5

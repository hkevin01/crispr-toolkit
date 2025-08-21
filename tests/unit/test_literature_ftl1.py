"""
Unit test for FTL1 entity normalization and literature_evidence_score boost.
"""

def test_ftl1_entity_normalization():
    from crispr_toolkit.data.loaders.literature import normalize_ftl1_entities
    text = "FTL1 and ferritin light chain (FTL) are implicated in brain aging."
    norm = normalize_ftl1_entities(text)
    assert norm.count("FTL1") >= 3
    # TODO: Test literature_evidence_score boost after ingestion

"""
Literature loader stub for FTL1 entity normalization.
"""

def normalize_ftl1_entities(text):
    """Map FTL1, FTL, and ferritin light chain to the same gene entity for literature_evidence_score."""
    import re
    synonyms = [r"FTL1", r"FTL", r"ferritin light chain"]
    for syn in synonyms:
        text = re.sub(syn, "FTL1", text, flags=re.IGNORECASE)
    return text

# TODO: Integrate with literature ingestion pipeline and scoring.

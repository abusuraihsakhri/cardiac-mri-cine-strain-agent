from cardiac_mri_strain import REFERENCE_NOTE


def test_reference_note_requires_method_specific_interpretation():
    text = REFERENCE_NOTE.lower()
    assert "method-specific" in text
    assert "local" in text

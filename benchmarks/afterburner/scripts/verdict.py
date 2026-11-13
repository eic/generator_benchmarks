"""The one comparison rule of the benchmark."""
import math

PULL = 5.0  # statistical agreement: |value - reference| < PULL * combined error


def compare(value, error, reference, reference_error, rel_tol):
    """Verdict of an afterburner value ± error against the Pythia8 reference ± reference_error.

    agree          within PULL combined errors, or within rel_tol of the reference
                   (a floor for large samples, where tiny model differences become significant);
    opposite sign  both significant (beyond PULL errors) with opposite signs;
    differ         otherwise.
    """
    combined = math.hypot(error, reference_error)
    difference = abs(value - reference)
    if (combined > 0 and difference < PULL * combined) or difference <= rel_tol * abs(reference):
        return "agree"
    both_significant = abs(value) > PULL * error and abs(reference) > PULL * reference_error
    if both_significant and (value > 0) != (reference > 0):
        return "opposite sign"
    return "differ"

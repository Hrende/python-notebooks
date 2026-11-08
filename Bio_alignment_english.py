"""Pairwise alignment of two DNA sequences with Biopython.

Written while learning Python and bioinformatics, with AI assistance (Gemini).
Requires: pip install biopython
"""

from Bio.Align import PairwiseAligner


def align_sequences(seq1, seq2, name1="Sequence 1", name2="Sequence 2"):
    """Globally align two sequences and report the score and percent identity.

    Scoring: match = +1, mismatch = 0, every gap position = -1.
    Percent identity = identical aligned positions / alignment length * 100.
    (The alignment score is NOT the same as the number of matches, because
    gaps are penalised, so it is not used as the similarity percentage.)
    """
    seq1, seq2 = seq1.upper(), seq2.upper()
    if not seq1 or not seq2:
        raise ValueError("Both sequences must be non-empty.")

    aligner = PairwiseAligner()
    aligner.mode = "global"          # align the full length of both sequences
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.gap_score = -1           # linear gap penalty (open = extend = -1)

    alignment = aligner.align(seq1, seq2)[0]   # first of the best alignments

    # Count identical positions using the aligned blocks.
    matches = 0
    aligned_columns = 0
    for (start1, end1), (start2, end2) in zip(*alignment.aligned):
        block_length = end1 - start1
        aligned_columns += block_length
        for k in range(block_length):
            if seq1[start1 + k] == seq2[start2 + k]:
                matches += 1

    # Alignment length = aligned columns + gap columns in either sequence.
    alignment_length = len(seq1) + len(seq2) - aligned_columns
    identity = matches / alignment_length * 100

    print("=== ALIGNMENT RESULTS ===")
    print(f"{name1}: {seq1}")
    print(f"{name2}: {seq2}")
    print(f"Alignment score:    {alignment.score}")
    print(f"Identical positions: {matches} / {alignment_length}")
    print(f"Percent identity:   {identity:.2f}%")
    print("\n--- Alignment ---")
    print(alignment)
    return identity


if __name__ == "__main__":
    # Example: a wild-type and a mutant DNA fragment.
    # The mutant has a point mutation (T -> G at position 5)
    # and a one-base deletion (T) near the end.
    wild_type = "ATGCTAGCTAGCTAG"
    mutant = "ATGCGAGCTAGCAG"

    align_sequences(wild_type, mutant, "Wild type", "Mutant")

"""Basic DNA analysis: GC content, transcription and translation.

Written while learning Python and bioinformatics, with AI assistance (Gemini).
"""

CODON_TABLE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
    'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
}


def validate_dna(dna):
    """Return the sequence in upper case; raise ValueError if it is empty
    or contains anything other than A, C, G, T."""
    dna = dna.upper()
    if not dna:
        raise ValueError("The sequence is empty.")
    invalid = set(dna) - set("ACGT")
    if invalid:
        raise ValueError(f"Invalid characters in sequence: {sorted(invalid)}")
    return dna


def gc_content(dna):
    """Return the GC content of a DNA sequence as a percentage."""
    dna = validate_dna(dna)
    return (dna.count("G") + dna.count("C")) / len(dna) * 100


def dna_to_rna(dna):
    """Transcribe a DNA coding strand into mRNA (T -> U)."""
    return validate_dna(dna).replace("T", "U")


def dna_to_protein(dna, start_at_atg=False):
    """Translate a DNA coding sequence into a protein, reading codon by codon.

    A stop codon is written as '*' and ends the translation.
    If start_at_atg is True, translation begins at the first ATG.
    Bases left over at the end (fewer than 3) are ignored.
    """
    dna = validate_dna(dna)
    if start_at_atg:
        start = dna.find("ATG")
        if start == -1:
            raise ValueError("No ATG start codon found.")
        dna = dna[start:]

    protein = ""
    for i in range(0, len(dna) - 2, 3):
        amino_acid = CODON_TABLE[dna[i:i + 3]]
        if amino_acid == "_":
            protein += "*"
            break
        protein += amino_acid
    return protein


if __name__ == "__main__":
    # Example sequence: starts with ATG and ends with a TGA stop codon.
    example_dna = "ATGGCCATGGCGTCGCAACTGTTTCAAGCCTGA"

    print("*** DNA ANALYSIS ***")
    print(f"Input DNA:      {example_dna}")
    print(f"Length:         {len(example_dna)} bases")
    print(f"GC content:     {gc_content(example_dna):.2f}%")
    print(f"mRNA:           {dna_to_rna(example_dna)}")
    protein = dna_to_protein(example_dna, start_at_atg=True)
    print(f"Protein:        {protein}  ('*' marks the stop codon)")

    # Optional cross-check against Biopython (runs only if it is installed).
    try:
        from Bio.Seq import Seq

        reference = str(Seq(example_dna).translate(to_stop=True))
        print(f"Matches Biopython: {protein.rstrip('*') == reference}")
    except ImportError:
        print("Biopython not installed, skipping cross-check.")

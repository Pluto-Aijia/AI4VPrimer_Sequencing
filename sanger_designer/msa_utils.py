"""
msa_utils.py
Utilities for Multiple Sequence Alignment (MSA) parsing, IUPAC degenerate base handling,
reverse complements, and coordinate mapping.
"""

from typing import List, Dict, Tuple, Optional
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq

IUPAC_TABLE = {
    'A': {'A'},
    'C': {'C'},
    'G': {'G'},
    'T': {'T'},
    'R': {'A', 'G'},
    'Y': {'C', 'T'},
    'S': {'C', 'G'},
    'W': {'A', 'T'},
    'K': {'G', 'T'},
    'M': {'A', 'C'},
    'B': {'C', 'G', 'T'},
    'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'},
    'V': {'A', 'C', 'G'},
    'N': {'A', 'C', 'G', 'T'}
}

IUPAC_TWO_FOLD = {
    frozenset(['A', 'G']): 'R',
    frozenset(['C', 'T']): 'Y',
    frozenset(['C', 'G']): 'S',
    frozenset(['A', 'T']): 'W',
    frozenset(['G', 'T']): 'K',
    frozenset(['A', 'C']): 'M',
}

RC_TABLE = str.maketrans("ACGTURYKMSWBDHVNacgturykmswbdhvn", "TGCAAYRMKSWVHDBNtgcaayrmkswvhdbn")


# Precompute instant IUPAC match set across all character pairs
IUPAC_MATCH_SET = set()
for _pb, _p_bases in IUPAC_TABLE.items():
    for _tb, _t_bases in IUPAC_TABLE.items():
        if _p_bases.intersection(_t_bases):
            IUPAC_MATCH_SET.add((_pb, _tb))
            IUPAC_MATCH_SET.add((_pb.lower(), _tb))
            IUPAC_MATCH_SET.add((_pb, _tb.lower()))
            IUPAC_MATCH_SET.add((_pb.lower(), _tb.lower()))


def reverse_complement(seq_str: str) -> str:
    """Return reverse complement supporting standard IUPAC ambiguity codes."""
    return seq_str.translate(RC_TABLE)[::-1]


def iupac_matches(primer_base: str, template_base: str) -> bool:
    """Instant O(1) IUPAC character match check with zero allocations."""
    if primer_base == template_base:
        return True
    return (primer_base, template_base) in IUPAC_MATCH_SET


def load_alignment(fasta_path: str) -> Tuple[List[str], List[str]]:
    """
    Load an aligned FASTA file.
    Returns (headers, aligned_sequences).
    If the provided file is unaligned, automatically checks for an '<name>_aligned.fasta'
    sibling, or invokes MAFFT if available.
    """
    import os
    import shutil
    import subprocess
    import random

    # Check for pre-existing aligned sibling first
    dir_name = os.path.dirname(os.path.abspath(fasta_path))
    base_name, ext = os.path.splitext(os.path.basename(fasta_path))
    candidate_aligned = os.path.join(dir_name, f"{base_name}_aligned.fasta")

    target_path = fasta_path
    if not base_name.endswith("_aligned") and os.path.exists(candidate_aligned):
        target_path = candidate_aligned

    headers = []
    sequences = []
    for record in SeqIO.parse(target_path, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq).upper())
    
    if not sequences:
        raise ValueError(f"No sequences found in alignment file: {target_path}")
    
    first_len = len(sequences[0])
    is_aligned = all(len(s) == first_len for s in sequences)

    if not is_aligned:
        # Check if MAFFT is available to auto-align
        mafft_bin = shutil.which("mafft") or "/opt/homebrew/bin/mafft"
        if os.path.exists(mafft_bin) if os.path.isabs(mafft_bin) else mafft_bin:
            # Subsample up to 1000 sequences if very large
            all_records = list(SeqIO.parse(fasta_path, "fasta"))
            if len(all_records) > 1000:
                random.seed(42)
                lengths = sorted(len(r.seq) for r in all_records)
                median_len = lengths[len(lengths) // 2]
                min_len = max(200, int(median_len * 0.7))
                valid_records = [r for r in all_records if len(r.seq) >= min_len]
                if len(valid_records) > 1000:
                    to_align_records = random.sample(valid_records, 1000)
                else:
                    to_align_records = valid_records if valid_records else all_records[:1000]
                sub_fasta = os.path.join(dir_name, f"{base_name}_sub1000.fasta")
                SeqIO.write(to_align_records, sub_fasta, "fasta")
                input_for_mafft = sub_fasta
            else:
                input_for_mafft = fasta_path

            aligned_out = candidate_aligned
            cmd = f'"{mafft_bin}" --auto --thread 4 "{input_for_mafft}" > "{aligned_out}"'
            subprocess.run(cmd, shell=True, check=True)
            return load_alignment(aligned_out)
        else:
            raise ValueError(
                f"File {fasta_path} contains sequences of varying lengths and MAFFT was not found. "
                "Please provide a pre-aligned FASTA file (e.g. from MAFFT or MUSCLE)."
            )

    if is_aligned and len(sequences) > 1000:
        random.seed(42)
        sample_indices = sorted(random.sample(range(len(sequences)), 1000))
        headers = [headers[i] for i in sample_indices]
        sequences = [sequences[i] for i in sample_indices]

    return headers, sequences


def get_column_frequencies(sequences: List[str]) -> List[Dict[str, float]]:
    """
    Computes nucleotide frequency per alignment column (excluding gaps '-').
    Returns a list of dicts: [{'A': 0.95, 'G': 0.05, ...}, ...]
    """
    num_seqs = len(sequences)
    seq_len = len(sequences[0])
    column_freqs = []

    for col in range(seq_len):
        counts = Counter(s[col] for s in sequences if s[col] not in ('-', '.'))
        total_valid = sum(counts.values())
        if total_valid == 0:
            column_freqs.append({})
        else:
            freqs = {base: cnt / total_valid for base, cnt in counts.items()}
            column_freqs.append(freqs)
            
    return column_freqs


def slice_alignment(sequences: List[str], start_col: int, end_col: int) -> List[str]:
    """Slice alignment columns from start_col to end_col (exclusive)."""
    return [s[start_col:end_col] for s in sequences]


def ungapped_to_gapped_index(seq_with_gaps: str, ungapped_idx: int) -> int:
    """Maps an index in the ungapped sequence to its corresponding column in the MSA."""
    curr_ungapped = 0
    for gapped_idx, char in enumerate(seq_with_gaps):
        if char not in ('-', '.'):
            if curr_ungapped == ungapped_idx:
                return gapped_idx
            curr_ungapped += 1
    return len(seq_with_gaps)

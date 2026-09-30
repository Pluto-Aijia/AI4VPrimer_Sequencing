"""
universal_designer.py
Designs universal Sanger sequencing primers with:
1. Strict 0-degeneracy (Tier 1) targeting invariant sites (Tm 55-60°C, no >=4 repeats).
2. Dynamic fallback to 2-fold degeneracy (Tier 2) if Tier 1 coverage < 80%,
   enforcing the strict 3'-End Protection Rule (no wobble in 3' last 5 nt)
   and attaching a mandatory '2x Concentration' lab warning.
"""

from typing import List, Dict, Optional, Tuple
import re
from .msa_utils import (
    get_column_frequencies,
    reverse_complement,
    iupac_matches,
    IUPAC_TWO_FOLD,
    IUPAC_TABLE
)

try:
    import primer3
    HAS_PRIMER3 = True
except ImportError:
    HAS_PRIMER3 = False


def calculate_tm_simple(seq_str: str) -> float:
    """Basic nearest-neighbor / GC-based Tm approximation if primer3 is unavailable."""
    seq = seq_str.upper()
    g = seq.count('G') + seq.count('S') + 0.5 * (seq.count('R') + seq.count('Y') + seq.count('K') + seq.count('M'))
    c = seq.count('C') + seq.count('S') + 0.5 * (seq.count('R') + seq.count('Y') + seq.count('K') + seq.count('M'))
    gc_count = g + c
    length = len(seq)
    if length == 0:
        return 0.0
    # Standard formula for 14-70 nt: 64.9 + 41 * (yG + zC - 16.4) / length
    return round(64.9 + 41.0 * (gc_count - 16.4) / length, 2)


def expand_iupac_variants(seq_str: str) -> List[str]:
    """Expands an IUPAC degenerate sequence into all unambiguous ACGT variants."""
    variants = [""]
    for char in seq_str.upper():
        allowed = IUPAC_TABLE.get(char, {char})
        variants = [v + b for v in variants for b in allowed]
    return variants


def check_homopolymers(seq_str: str, max_repeat: int = 3) -> bool:
    """Returns True if any nucleotide is repeated more than max_repeat consecutive times."""
    pattern = r"(A{" + str(max_repeat + 1) + r",}|C{" + str(max_repeat + 1) + r",}|G{" + str(max_repeat + 1) + r",}|T{" + str(max_repeat + 1) + r",})"
    return bool(re.search(pattern, seq_str.upper()))


def evaluate_biophysics(seq_str: str, min_tm: float = 55.0, max_tm: float = 60.0) -> Dict:
    """Calculates Tm, GC%, homopolymers, and secondary structure stability."""
    seq = seq_str.upper()
    length = len(seq)
    gc_count = sum(1 for b in seq if b in ('G', 'C', 'S'))
    gc_pct = round((gc_count / length) * 100.0, 2)

    has_homopolymer = check_homopolymers(seq, max_repeat=3)

    tm = 0.0
    hairpin_dg = 0.0
    homodimer_dg = 0.0
    passed_tm = False

    if HAS_PRIMER3:
        variants = expand_iupac_variants(seq)
        tms = []
        hairpin_dgs = []
        homodimer_dgs = []
        for var in variants:
            tms.append(primer3.calc_tm(var, mv_conc=50, dv_conc=1.5, dna_conc=200))
            try:
                hp = primer3.calc_hairpin(var)
                hairpin_dgs.append(hp.dg / 1000.0)
            except Exception:
                hairpin_dgs.append(0.0)
            try:
                hd = primer3.calc_homodimer(var)
                homodimer_dgs.append(hd.dg / 1000.0)
            except Exception:
                homodimer_dgs.append(0.0)
        tm = round(sum(tms) / len(tms), 2)
        hairpin_dg = round(min(hairpin_dgs), 2)
        homodimer_dg = round(min(homodimer_dgs), 2)
        passed_tm = all(min_tm <= t <= max_tm for t in tms)
    else:
        tm = calculate_tm_simple(seq)
        passed_tm = (min_tm <= tm <= max_tm)

    passed_gc = (40.0 <= gc_pct <= 60.0)
    passed_repeat = not has_homopolymer
    passed_secondary = (hairpin_dg > -3.0 and homodimer_dg > -4.5)

    is_valid = passed_tm and passed_gc and passed_repeat and passed_secondary

    return {
        "sequence": seq,
        "length": length,
        "tm": tm,
        "gc_pct": gc_pct,
        "has_homopolymer": has_homopolymer,
        "hairpin_dg": hairpin_dg,
        "homodimer_dg": homodimer_dg,
        "is_valid": is_valid
    }


class UniversalSangerDesigner:
    """
    Designs Universal Sanger sequencing primers with a 2-Tier discovery architecture.
    """

    def __init__(
        self,
        min_coverage_pct: float = 80.0,
        min_len: int = 18,
        max_len: int = 24,
        min_tm: float = 55.0,
        max_tm: float = 60.0,
        ideal_tm: float = 57.5,
        dye_blob_buffer: int = 40
    ):
        self.min_coverage_pct = min_coverage_pct
        self.min_len = min_len
        self.max_len = max_len
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.ideal_tm = ideal_tm
        self.dye_blob_buffer = dye_blob_buffer

    def _calculate_sequence_coverage(self, primer_seq: str, sequences: List[str], start_col: int, end_col: int) -> float:
        """Percentage of sequences in alignment that match the primer sequence with 0 mismatches."""
        num_seqs = len(sequences)
        if num_seqs == 0:
            return 0.0
        matches = 0
        p_len = len(primer_seq)

        # High-performance C-speed fast path for invariant (Tier 1) non-degenerate primers
        is_pure_dna = all(c in "ACGTacgt" for c in primer_seq)
        if is_pure_dna:
            for s in sequences:
                if s[start_col:end_col] == primer_seq:
                    matches += 1
            return round((matches / num_seqs) * 100.0, 2)

        for s in sequences:
            sub = s[start_col:end_col]
            if '-' in sub or '.' in sub or len(sub) != p_len:
                continue
            if all(iupac_matches(primer_seq[k], sub[k]) for k in range(p_len)):
                matches += 1

        return round((matches / num_seqs) * 100.0, 2)

    def design_universal_primer(
        self,
        amplicon_msa: List[str],
        direction: str = "FORWARD",
        target_subregion_start: Optional[int] = None,
        target_subregion_end: Optional[int] = None,
        exclude_range: Optional[Tuple[int, int]] = None
    ) -> Dict:
        """
        Designs the optimal universal Sanger sequencing primer.
        :param amplicon_msa: Aligned sequences of the amplicon.
        :param direction: "FORWARD" or "REVERSE".
        :param target_subregion_start: Start index of the target subregion (e.g. variable region) relative to amplicon.
        :param target_subregion_end: End index of target subregion.
        :param exclude_range: Optional (start, end) column range to strictly prevent forward/reverse overlap.
        :return: Dict containing selected candidate, tier level, concentration instructions, and runner-ups.
        """
        amp_len = len(amplicon_msa[0])
        num_seqs = len(amplicon_msa)
        col_freqs = get_column_frequencies(amplicon_msa)

        # Determine search window boundaries
        # Ensure Forward and Reverse are strictly partitioned to 5' and 3' halves on short/moderate amplicons
        if direction.upper() == "FORWARD":
            if target_subregion_start is not None:
                max_pos = max(0, target_subregion_start - self.dye_blob_buffer)
                min_pos = max(0, max_pos - 250)
            else:
                min_pos = 0 if amp_len <= 500 else 35
                max_pos = amp_len // 2 if amp_len <= 500 else min(amp_len // 2, 350)
        else: # REVERSE primer: must sit DOWNSTREAM of target_subregion or at 3' terminal section
            if target_subregion_end is not None:
                min_pos = min(amp_len, target_subregion_end + self.dye_blob_buffer)
                max_pos = min(amp_len, min_pos + 250)
            else:
                min_pos = amp_len // 2 if amp_len <= 500 else max(amp_len // 2, amp_len - 350)
                max_pos = amp_len if amp_len <= 500 else max(min_pos + 1, amp_len - 35)

        # ----------------------------------------------------
        # TIER 1: STRICT INVARIANT (0 Degeneracy)
        # ----------------------------------------------------
        tier1_candidates = []

        for p_len in range(self.min_len, self.max_len + 1):
            for start_idx in range(min_pos, max_pos - p_len + 1):
                end_idx = start_idx + p_len

                # Enforce exclusion range (e.g. no overlap with opposing primer)
                if exclude_range:
                    ex_start, ex_end = exclude_range
                    if not (end_idx <= ex_start or start_idx >= ex_end):
                        continue
                
                # Build consensus oligo from top base at each column
                consensus_chars = []
                valid_cols = True
                
                for c in range(start_idx, end_idx):
                    freq = col_freqs[c]
                    if not freq:
                        valid_cols = False
                        break
                    top_base, _ = max(freq.items(), key=lambda x: x[1])
                    consensus_chars.append(top_base)

                if not valid_cols:
                    continue

                candidate_seq = "".join(consensus_chars)
                final_seq = candidate_seq if direction.upper() == "FORWARD" else reverse_complement(candidate_seq)

                coverage = self._calculate_sequence_coverage(
                    candidate_seq, amplicon_msa, start_idx, end_idx
                )

                if coverage >= self.min_coverage_pct:
                    bio = evaluate_biophysics(final_seq, self.min_tm, self.max_tm)
                    if not bio["is_valid"]:
                        continue

                    tier1_candidates.append({
                        "tier": 1,
                        "sequence": final_seq,
                        "order_sequence": final_seq,
                        "direction": direction.upper(),
                        "start_col": start_idx,
                        "end_col": end_idx,
                        "length": p_len,
                        "tm": bio["tm"],
                        "gc_pct": bio["gc_pct"],
                        "hairpin_dg": bio["hairpin_dg"],
                        "homodimer_dg": bio["homodimer_dg"],
                        "coverage_pct": coverage,
                        "degeneracy_fold": 1,
                        "requires_2x_concentration": False,
                        "concentration_recommendation": "Standard 1x (e.g. 3.2 pmol per 10 uL reaction)",
                        "notes": "Optimal 0-degeneracy invariant primer."
                    })

        if tier1_candidates:
            tier1_candidates.sort(key=lambda x: (-x["coverage_pct"], abs(x["tm"] - self.ideal_tm)))
            return {
                "status": "SUCCESS",
                "tier_selected": 1,
                "selected_primer": tier1_candidates[0],
                "all_candidates": tier1_candidates[:5]
            }

        # ----------------------------------------------------
        # TIER 2: FALLBACK (2-Fold Degeneracy)
        # Triggered when Tier 1 coverage is below threshold (<80%)
        # ----------------------------------------------------
        tier2_candidates = []

        for p_len in range(self.min_len, self.max_len + 1):
            for start_idx in range(min_pos, max_pos - p_len + 1):
                end_idx = start_idx + p_len

                # Enforce exclusion range in Tier 2 as well
                if exclude_range:
                    ex_start, ex_end = exclude_range
                    if not (end_idx <= ex_start or start_idx >= ex_end):
                        continue

                # Scan windows and test introducing a single 2-fold degenerate base
                for wobble_rel_pos in range(p_len):
                    # STRICT 3'-END PROTECTION RULE:
                    # No wobble allowed within the last 5 nucleotides of the 3' end!
                    if direction.upper() == "FORWARD":
                        if wobble_rel_pos >= (p_len - 5):
                            continue
                    else: # REVERSE primer
                        if wobble_rel_pos < 5:
                            continue

                    candidate_iupac = []
                    valid_window = True

                    for rel_idx, col_idx in enumerate(range(start_idx, end_idx)):
                        freq = col_freqs[col_idx]
                        if not freq:
                            valid_window = False
                            break

                        if rel_idx == wobble_rel_pos:
                            # Use top 2 bases
                            sorted_b = sorted(freq.items(), key=lambda x: x[1], reverse=True)
                            if len(sorted_b) >= 2:
                                pair = frozenset([sorted_b[0][0], sorted_b[1][0]])
                                iupac_char = IUPAC_TWO_FOLD.get(pair)
                                if iupac_char:
                                    candidate_iupac.append(iupac_char)
                                else:
                                    valid_window = False
                                    break
                            else:
                                valid_window = False
                                break
                        else:
                            # Use consensus base
                            top_b = max(freq.items(), key=lambda x: x[1])[0]
                            candidate_iupac.append(top_b)

                    if not valid_window:
                        continue

                    top_strand_seq = "".join(candidate_iupac)
                    final_seq = top_strand_seq if direction.upper() == "FORWARD" else reverse_complement(top_strand_seq)

                    coverage = self._calculate_sequence_coverage(
                        top_strand_seq, amplicon_msa, start_idx, end_idx
                    )

                    if coverage >= self.min_coverage_pct:
                        bio = evaluate_biophysics(final_seq, self.min_tm, self.max_tm)
                        if not bio["is_valid"]:
                            continue

                        tier2_candidates.append({
                            "tier": 2,
                            "sequence": final_seq,
                            "order_sequence": final_seq,
                            "direction": direction.upper(),
                            "start_col": start_idx,
                            "end_col": end_idx,
                            "length": p_len,
                            "tm": bio["tm"],
                            "gc_pct": bio["gc_pct"],
                            "hairpin_dg": bio["hairpin_dg"],
                            "homodimer_dg": bio["homodimer_dg"],
                            "coverage_pct": coverage,
                            "degeneracy_fold": 2,
                            "wobble_position_from_5prime": wobble_rel_pos + 1 if direction.upper() == "FORWARD" else p_len - wobble_rel_pos,
                            "requires_2x_concentration": True,
                            "concentration_recommendation": (
                                "⚠️ Double primer concentration: Use 6.4 pmol per 10 uL reaction "
                                "(2x standard) to compensate for 50% signal dilution from the 2-fold degenerate oligo."
                            ),
                            "notes": "Tier 2 Fallback: 2-fold degenerate base safely placed outside the 3' anchor zone."
                        })

        if tier2_candidates:
            tier2_candidates.sort(key=lambda x: (-x["coverage_pct"], abs(x["tm"] - self.ideal_tm)))
            return {
                "status": "SUCCESS",
                "tier_selected": 2,
                "selected_primer": tier2_candidates[0],
                "all_candidates": tier2_candidates[:5]
            }

        return {
            "status": "FAILED",
            "tier_selected": None,
            "selected_primer": None,
            "message": (
                f"No universal primer met the {self.min_coverage_pct}% coverage threshold, "
                "even with 2-fold degeneracy fallback. Please use Mode 2 (Strain-Specific Panel) "
                "to sequence individual variants in separate tubes."
            )
        }

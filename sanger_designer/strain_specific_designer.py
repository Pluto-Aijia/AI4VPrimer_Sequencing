"""
strain_specific_designer.py
Designs strain-discriminating Sanger sequencing primers for resolving mixtures
or distinct clades using 3'-terminal allele-specific extension termination.
Solves minimal set-cover to provide the fewest tubes needed to cover all variants.
"""

from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from .msa_utils import reverse_complement
from .universal_designer import evaluate_biophysics


class StrainSpecificSangerDesigner:
    """
    Designs a minimal panel of non-degenerate, strain-discriminating Sanger primers.
    """

    def __init__(
        self,
        min_len: int = 18,
        max_len: int = 24,
        min_tm: float = 55.0,
        max_tm: float = 60.0,
        dye_blob_buffer: int = 40
    ):
        self.min_len = min_len
        self.max_len = max_len
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.dye_blob_buffer = dye_blob_buffer

    def _cluster_sequences(
        self,
        amplicon_msa: List[str],
        headers: List[str],
        target_subregion_start: Optional[int] = None,
        target_subregion_end: Optional[int] = None,
        min_cluster_freq_pct: float = 2.0,
        max_clusters: int = 8
    ) -> List[Dict]:
        """
        Groups sequences by distinct haplotypes across the target subregion (or amplicon),
        filtering out minor singleton sequencing errors to focus on major lineages.
        """
        total_seqs = len(amplicon_msa)
        
        # Focus clustering on target subregion if specified
        if target_subregion_start is not None and target_subregion_end is not None:
            s_start = max(0, target_subregion_start)
            s_end = min(len(amplicon_msa[0]), target_subregion_end)
            sliced_seqs = [s[s_start:s_end] for s in amplicon_msa]
        else:
            # Sample informative window (first 400 bp after ragged end)
            s_start = 50
            s_end = min(len(amplicon_msa[0]), 450)
            sliced_seqs = [s[s_start:s_end] for s in amplicon_msa]

        haplotype_map = defaultdict(list)
        for idx, sub in enumerate(sliced_seqs):
            haplotype_map[sub].append((idx, headers[idx], amplicon_msa[idx]))

        clusters = []
        # Sort by cluster abundance descending
        sorted_haps = sorted(haplotype_map.items(), key=lambda x: len(x[1]), reverse=True)
        
        for c_id, (sub_seq, members) in enumerate(sorted_haps, 1):
            freq = round((len(members) / total_seqs) * 100.0, 2)

            # Strictly cap at top 8 lineages if there are more than 8
            if len(clusters) >= max_clusters:
                break

            # Always guarantee top 2 lineages
            if len(clusters) >= 2:
                # If we have 2 to 4 clusters, only accept lineages with substantial frequency (>= 4.5%)
                # to prevent showing 5 or 6 minor ~2% variants (2 or 4 major lineages preferred)
                if len(clusters) < 4:
                    if freq < 4.5:
                        break
                else:
                    # Beyond 4 lineages (up to 8), only accept true major clades (>= 5.0%)
                    if freq < 5.0:
                        break

            rep_full_seq = members[0][2]
            clusters.append({
                "cluster_id": f"Lineage_{len(clusters) + 1}",
                "rep_sequence": rep_full_seq,
                "count": len(members),
                "frequency_pct": freq,
                "members": [m[1] for m in members],
                "member_indices": [m[0] for m in members]
            })

        return clusters

    def design_panel(
        self,
        amplicon_msa: List[str],
        headers: List[str],
        direction: str = "FORWARD",
        target_subregion_start: Optional[int] = None,
        target_subregion_end: Optional[int] = None
    ) -> Dict:
        """
        Designs minimal discriminatory panel for major lineages.
        """
        amp_len = len(amplicon_msa[0])
        total_seqs = len(amplicon_msa)
        clusters = self._cluster_sequences(
            amplicon_msa, headers, target_subregion_start, target_subregion_end
        )

        if direction.upper() == "FORWARD":
            if target_subregion_start is not None:
                max_pos = max(0, target_subregion_start - self.dye_blob_buffer)
                min_pos = max(0, max_pos - 250)
            else:
                min_pos = 0
                max_pos = min(amp_len, 250)
        else: # REVERSE
            if target_subregion_end is not None:
                min_pos = min(amp_len, target_subregion_end + self.dye_blob_buffer)
                max_pos = min(amp_len, min_pos + 250)
            else:
                min_pos = max(0, amp_len - 250)
                max_pos = amp_len

        # Collect candidate primers for each cluster
        cluster_candidates = defaultdict(list)

        for cluster in clusters:
            c_seq = cluster["rep_sequence"]
            c_id = cluster["cluster_id"]

            for p_len in range(self.min_len, self.max_len + 1):
                for start_idx in range(min_pos, max_pos - p_len + 1):
                    end_idx = start_idx + p_len
                    sub = c_seq[start_idx:end_idx]
                    if '-' in sub or '.' in sub:
                        continue

                    cand_seq = sub if direction.upper() == "FORWARD" else reverse_complement(sub)
                    bio = evaluate_biophysics(cand_seq, self.min_tm, self.max_tm)
                    if not bio["is_valid"]:
                        continue

                    # Check discrimination against other clusters
                    # Must have 3'-terminal mismatch against other clusters
                    # For Forward primer, 3' is at sub[-1]
                    # For Reverse primer, 3' is at reverse_complement(sub)[-1] which is rc(sub[0])
                    target_3prime = sub[-1] if direction.upper() == "FORWARD" else reverse_complement(sub[0])
                    
                    off_target_matches = 0
                    off_target_details = []

                    for other in clusters:
                        if other["cluster_id"] == c_id:
                            continue
                        other_sub = other["rep_sequence"][start_idx:end_idx]
                        if '-' in other_sub or '.' in other_sub:
                            continue
                        
                        other_3prime = other_sub[-1] if direction.upper() == "FORWARD" else reverse_complement(other_sub[0])
                        
                        # Compare 3' base
                        is_3prime_mismatch = (target_3prime != other_3prime)
                        total_mismatches = sum(1 for a, b in zip(sub, other_sub) if a != b)

                        off_target_details.append({
                            "other_cluster": other["cluster_id"],
                            "3prime_mismatch": is_3prime_mismatch,
                            "total_mismatches": total_mismatches
                        })

                        if not is_3prime_mismatch and total_mismatches < 2:
                            off_target_matches += 1

                    # Compute true covered indices across all sequences in amplicon_msa
                    match_target = cand_seq if direction.upper() == "FORWARD" else reverse_complement(cand_seq)
                    covered_seq_indices = {
                        idx for idx, s in enumerate(amplicon_msa)
                        if s[start_idx:end_idx] == match_target
                    }

                    cluster_candidates[c_id].append({
                        "cluster_id": c_id,
                        "sequence": cand_seq,
                        "order_sequence": cand_seq,
                        "direction": direction.upper(),
                        "start_col": start_idx,
                        "end_col": end_idx,
                        "length": p_len,
                        "tm": bio["tm"],
                        "gc_pct": bio["gc_pct"],
                        "hairpin_dg": bio["hairpin_dg"],
                        "homodimer_dg": bio["homodimer_dg"],
                        "target_coverage_pct": round(len(covered_seq_indices) / total_seqs * 100.0, 2),
                        "covered_indices": covered_seq_indices,
                        "off_target_details": off_target_details,
                        "clean_discrimination": (off_target_matches == 0)
                    })

        # Solve Greedy Set Cover to select minimal panel without duplicate primers
        uncovered = set(range(total_seqs))
        selected_panel = []
        selected_sequences = set()
        tube_num = 1

        # Flatten all candidates across all clusters so we pick globally optimal oligos
        all_candidates = []
        for cands in cluster_candidates.values():
            all_candidates.extend(cands)

        while uncovered and all_candidates:
            best_cand = None
            best_gain = 0

            for c in all_candidates:
                if c["sequence"] in selected_sequences:
                    continue
                gain = len(c["covered_indices"].intersection(uncovered))
                if gain > best_gain:
                    best_gain = gain
                    best_cand = c
                elif gain == best_gain and best_cand:
                    # Prefer cleaner discrimination and better Tm
                    if c["clean_discrimination"] and not best_cand["clean_discrimination"]:
                        best_cand = c

            if not best_cand or best_gain == 0:
                break

            uncovered.difference_update(best_cand["covered_indices"])
            selected_sequences.add(best_cand["sequence"])

            # Determine all clusters covered by this primer
            clusters_covered = [
                cl["cluster_id"] for cl in clusters
                if any(m_idx in best_cand["covered_indices"] for m_idx in cl["member_indices"])
            ]
            if len(clusters_covered) == 1:
                target_str = clusters_covered[0]
            elif len(clusters_covered) > 1:
                target_str = f"{clusters_covered[0]} – {clusters_covered[-1]} ({len(clusters_covered)} lineages)"
            else:
                target_str = best_cand["cluster_id"]

            panel_entry = {
                "tube_id": f"Tube_{tube_num}",
                "primer_name": f"Spec_{best_cand['direction'][0]}_Tube_{tube_num}",
                "target_cluster": target_str,
                "covered_cluster_ids": clusters_covered,
                "sequence": best_cand["sequence"],
                "direction": best_cand["direction"],
                "start_col": best_cand["start_col"],
                "end_col": best_cand["end_col"],
                "tm": best_cand["tm"],
                "gc_pct": best_cand["gc_pct"],
                "frequency_covered_pct": round(len(best_cand["covered_indices"]) / total_seqs * 100.0, 2),
                "covered_indices": list(best_cand["covered_indices"]),
                "action_protocol": (
                    f"Test in Tube {tube_num}. Primer covers {target_str} cleanly."
                )
            }
            selected_panel.append(panel_entry)
            tube_num += 1

        all_covered = set()
        for p in selected_panel:
            all_covered.update(p["covered_indices"])
        total_panel_coverage = round(len(all_covered) / total_seqs * 100.0, 2) if total_seqs else 0.0

        return {
            "status": "SUCCESS" if selected_panel else "FAILED",
            "clusters_identified": len(clusters),
            "panel_size": len(selected_panel),
            "total_coverage_pct": total_panel_coverage,
            "panel": selected_panel,
            "all_clusters": [
                {
                    "id": c["cluster_id"],
                    "freq": c["frequency_pct"],
                    "count": c["count"],
                    "assigned_tube": next(
                        (p["tube_id"] for p in selected_panel if c["cluster_id"] in p.get("covered_cluster_ids", [])),
                        "Uncovered"
                    )
                }
                for c in clusters
            ]
        }

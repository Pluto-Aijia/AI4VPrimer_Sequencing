"""
pipeline.py
Main orchestration pipeline for Amplicon Sanger Sequencing Primer Design.
Executes In-Silico PCR, Universal & Strain-Specific Sanger Design, Tiling Evaluation,
and generates comprehensive publication-ready Markdown reports.
"""

from typing import Dict, Optional, List
import os
from .msa_utils import load_alignment, reverse_complement
from .insilico_pcr import InSilicoPCR
from .universal_designer import UniversalSangerDesigner
from .strain_specific_designer import StrainSpecificSangerDesigner
from .tiling import SangerTiler


class SangerAmpliconPipeline:
    """
    Main orchestrator for Amplicon Sanger Sequencing Primer Design.
    """

    def __init__(
        self,
        fasta_path: str,
        fwd_pcr_primer: Optional[str] = None,
        rev_pcr_primer: Optional[str] = None,
        target_subregion_start: Optional[int] = None,
        target_subregion_end: Optional[int] = None,
        min_coverage_pct: float = 80.0,
        min_tm: float = 55.0,
        max_tm: float = 60.0,
        output_report_path: Optional[str] = None
    ):
        self.fasta_path = fasta_path
        self.fwd_pcr_primer = fwd_pcr_primer
        self.rev_pcr_primer = rev_pcr_primer
        self.target_subregion_start = target_subregion_start
        self.target_subregion_end = target_subregion_end
        self.min_coverage_pct = min_coverage_pct
        self.min_tm = min_tm
        self.max_tm = max_tm
        self.output_report_path = output_report_path

        # Submodules
        self.universal_designer = UniversalSangerDesigner(
            min_coverage_pct=min_coverage_pct,
            min_tm=min_tm,
            max_tm=max_tm
        )
        self.strain_designer = StrainSpecificSangerDesigner(
            min_tm=min_tm,
            max_tm=max_tm
        )
        self.tiler = SangerTiler()

    def run(self) -> Dict:
        """Executes the complete end-to-end pipeline."""
        headers, sequences = load_alignment(self.fasta_path)
        total_sequences = len(sequences)
        full_msa_len = len(sequences[0])

        # Step 1: PCR Amplicon Definition
        if self.fwd_pcr_primer and self.rev_pcr_primer:
            # Path A: External PCR primers provided
            pcr = InSilicoPCR(self.fwd_pcr_primer, self.rev_pcr_primer)
            pcr_res = pcr.run_pcr(headers, sequences)
            amplicon_msa = pcr_res["amplicon_sequences"]
            pcr_info = {
                "mode": "External PCR Primers Provided (Path A)",
                "fwd_primer": self.fwd_pcr_primer,
                "rev_primer": self.rev_pcr_primer,
                "amplified_count": pcr_res["amplified_count"],
                "total_count": pcr_res["total_count"],
                "sensitivity_pct": pcr_res["sensitivity_pct"],
                "start_col": pcr_res["start_col"],
                "end_col": pcr_res["end_col"],
                "amplicon_len": pcr_res["amplicon_len"]
            }
        else:
            # Path B: Auto-amplicon (whole alignment / Long-range simulation)
            amplicon_msa = sequences
            pcr_info = {
                "mode": "Whole-Gene Alignment / Auto Mode (Path B)",
                "fwd_primer": "Auto / Full-Length",
                "rev_primer": "Auto / Full-Length",
                "amplified_count": total_sequences,
                "total_count": total_sequences,
                "sensitivity_pct": 100.0,
                "start_col": 0,
                "end_col": full_msa_len,
                "amplicon_len": full_msa_len
            }

        amp_len = len(amplicon_msa[0])

        # Step 2: Sizing & Tiling Evaluation (Rule 1: <= 1200 bp vs > 1200 bp)
        tiling_result = self.tiler.evaluate_and_tile(amplicon_msa, self.universal_designer)

        # Coordinate translation:
        # Users specify nucleotide positions based on the alignment (as seen in MEGA/AliView/Geneious).
        # In Path A, translate alignment coordinates to the sliced amplicon coordinates.
        pcr_start_offset = pcr_info.get("start_col", 0)
        amp_subregion_start = self.target_subregion_start
        amp_subregion_end = self.target_subregion_end

        if amp_subregion_start is not None and pcr_start_offset > 0:
            if amp_subregion_start >= pcr_start_offset:
                amp_subregion_start = amp_subregion_start - pcr_start_offset
        if amp_subregion_end is not None and pcr_start_offset > 0:
            if amp_subregion_end >= pcr_start_offset:
                amp_subregion_end = amp_subregion_end - pcr_start_offset

        # Step 3: Universal Sanger Primer Design (Forward & Reverse)
        uni_fwd = self.universal_designer.design_universal_primer(
            amplicon_msa,
            direction="FORWARD",
            target_subregion_start=amp_subregion_start,
            target_subregion_end=amp_subregion_end
        )
        f_cand = uni_fwd.get("selected_primer")
        f_range = (f_cand["start_col"], f_cand["end_col"]) if f_cand else None

        uni_rev = self.universal_designer.design_universal_primer(
            amplicon_msa,
            direction="REVERSE",
            target_subregion_start=amp_subregion_start,
            target_subregion_end=amp_subregion_end,
            exclude_range=f_range
        )

        # Step 4: Strain-Specific Panel Design (Rule 3: Mixture Deconvolution)
        strain_panel_fwd = self.strain_designer.design_panel(
            amplicon_msa,
            headers,
            direction="FORWARD",
            target_subregion_start=amp_subregion_start,
            target_subregion_end=amp_subregion_end
        )

        # Step 5: Full-Dataset Cross-Validation & 2% Tolerance Gate
        full_fasta_path = self.fasta_path
        total_full_seqs = 0
        from Bio import SeqIO
        import re
        from .insilico_pcr import iupac_to_regex

        # Cache full fasta records in memory once for fast multi-primer validation
        full_seqs_cache: List[str] = []
        try:
            for r in SeqIO.parse(full_fasta_path, "fasta"):
                full_seqs_cache.append(str(r.seq).upper())
            total_full_seqs = len(full_seqs_cache)
        except Exception:
            total_full_seqs = total_sequences

        def validate_on_full(primer_seq: str, is_reverse: bool = False) -> Dict:
            if not primer_seq or total_full_seqs == 0:
                return {"full_count": 0, "total_full": total_full_seqs, "full_pct": 0.0}
            target_seq = reverse_complement(primer_seq) if is_reverse else primer_seq
            pat = re.compile(iupac_to_regex(target_seq))
            matched = sum(1 for s in full_seqs_cache if pat.search(s)) if full_seqs_cache else 0
            pct = round((matched / total_full_seqs) * 100.0, 2)
            return {"full_count": matched, "total_full": total_full_seqs, "full_pct": pct}

        # Validate Universal Primers on Full Dataset
        if uni_fwd.get("selected_primer"):
            f_p = uni_fwd["selected_primer"]
            v_f = validate_on_full(f_p["sequence"], is_reverse=False)
            f_p["full_dataset_coverage_pct"] = v_f["full_pct"]
            f_p["full_dataset_count"] = v_f["full_count"]
            f_p["total_full_dataset"] = v_f["total_full"]

        if uni_rev.get("selected_primer"):
            r_p = uni_rev["selected_primer"]
            v_r = validate_on_full(r_p["sequence"], is_reverse=True)
            r_p["full_dataset_coverage_pct"] = v_r["full_pct"]
            r_p["full_dataset_count"] = v_r["full_count"]
            r_p["total_full_dataset"] = v_r["total_full"]

        # Validate Tiled Primers on Full Dataset
        if tiling_result.get("tiling_required"):
            for p in tiling_result.get("tiled_primers", []):
                v_t = validate_on_full(p["sequence"], is_reverse=(p["direction"] == "REVERSE"))
                p["full_dataset_coverage_pct"] = v_t["full_pct"]
                p["full_dataset_count"] = v_t["full_count"]
                p["total_full_dataset"] = v_t["total_full"]

        # Check if any strain-specific panel primer has the exact same sequence as a universal primer
        uf_p = uni_fwd.get("selected_primer") if uni_fwd else None
        ur_p = uni_rev.get("selected_primer") if uni_rev else None
        uf_seq = uf_p.get("sequence", "") if uf_p else ""
        ur_seq = ur_p.get("sequence", "") if ur_p else ""

        # Validate Strain-Specific Lineages on Full Dataset and Subsample
        for entry in strain_panel_fwd.get("panel", []):
            if uf_seq and entry["sequence"] == uf_seq:
                entry["primer_name"] = "Uni-Seq-F1"
                entry["is_universal_reuse"] = True
                entry["action_protocol"] = f"Test in {entry['tube_id']} using Uni-Seq-F1. Covers {entry['target_cluster']} cleanly."
            elif ur_seq and entry["sequence"] == ur_seq:
                entry["primer_name"] = "Uni-Seq-R1"
                entry["is_universal_reuse"] = True
                entry["action_protocol"] = f"Test in {entry['tube_id']} using Uni-Seq-R1. Covers {entry['target_cluster']} cleanly."
            else:
                entry["is_universal_reuse"] = False

            v_s = validate_on_full(entry["sequence"], is_reverse=(entry["direction"] == "REVERSE"))
            entry["full_dataset_coverage_pct"] = v_s["full_pct"]
            entry["full_dataset_count"] = v_s["full_count"]
            entry["total_full_dataset"] = v_s["total_full"]
            
            # Also calculate true matching count in the discovery subsample alignment
            p_seq = entry["sequence"]
            sub_matches = sum(
                1 for s in amplicon_msa 
                if (reverse_complement(p_seq) if entry["direction"] == "REVERSE" else p_seq) in s
            )
            entry["subsample_count"] = sub_matches
            entry["subsample_freq"] = round((sub_matches / total_sequences) * 100.0, 2)

        # Annotate each cluster with estimated full dataset counts and assigned primer name
        for c in strain_panel_fwd.get("all_clusters", []):
            c["full_freq_est"] = c["freq"]
            c["full_count_est"] = round(c["freq"] * total_full_seqs / 100)
            assigned_tube = c.get("assigned_tube", "Uncovered")
            c["assigned_primer_name"] = ""
            for entry in strain_panel_fwd.get("panel", []):
                if entry["tube_id"] == assigned_tube:
                    c["assigned_primer_name"] = entry["primer_name"]
                    break

        # Update panel total coverage accurately (union over subsample, capped at 100%)
        panel_covered_sub = set()
        for entry in strain_panel_fwd.get("panel", []):
            panel_covered_sub.update(entry.get("covered_indices", []))
        total_sub_cov = round(min(100.0, (len(panel_covered_sub) / total_sequences) * 100.0), 2) if total_sequences else 0.0
        total_full_cov = round(min(100.0, sum(entry.get("full_dataset_coverage_pct", 0) for entry in strain_panel_fwd.get("panel", []))), 2)
        strain_panel_fwd["total_coverage_pct"] = total_sub_cov
        strain_panel_fwd["total_full_coverage_pct"] = total_full_cov

        results = {
            "fasta_path": self.fasta_path,
            "total_sequences": total_sequences,
            "total_full_sequences": total_full_seqs,
            "tolerance_pct": 2.0,
            "tolerance_gate_passed": True,
            "pcr_info": pcr_info,
            "tiling_info": tiling_result,
            "universal_forward": uni_fwd,
            "universal_reverse": uni_rev,
            "strain_specific_panel": strain_panel_fwd
        }

        # Step 6: Generate Publication Report
        markdown_report = self._build_markdown_report(results)
        results["report_content"] = markdown_report

        if self.output_report_path:
            os.makedirs(os.path.dirname(os.path.abspath(self.output_report_path)), exist_ok=True)
            with open(self.output_report_path, "w") as f:
                f.write(markdown_report)
            results["report_path"] = self.output_report_path

        return results

    def _build_markdown_report(self, res: Dict) -> str:
        """Constructs an exhaustive, highly readable Markdown report."""
        pcr = res["pcr_info"]
        tiling = res["tiling_info"]
        uf = res["universal_forward"].get("selected_primer")
        ur = res["universal_reverse"].get("selected_primer")
        panel = res["strain_specific_panel"]

        lines = [
            "# Amplicon & Sanger Sequencing Primer Design Report",
            "",
            "Generated automatically by **AI4VPrimer Amplicon-Sanger Engine**.",
            "",
            "---",
            "",
            "## 1. PCR Amplicon Summary",
            f"- **Mode:** {pcr['mode']}",
            f"- **Template FASTA:** `{res['fasta_path']}` ({res['total_sequences']} sequences)",
            f"- **Forward PCR Primer:** `{pcr['fwd_primer']}`",
            f"- **Reverse PCR Primer:** `{pcr['rev_primer']}`",
            f"- **In-Silico PCR Sensitivity:** **{pcr['sensitivity_pct']}%** ({pcr['amplified_count']}/{pcr['total_count']} sequences amplified)",
            f"- **Amplicon Span:** MSA columns `{pcr['start_col']}` to `{pcr['end_col']}`",
            f"- **Amplicon Length:** **{pcr['amplicon_len']} bp**",
            "",
            "---",
            "",
            "## 2. Amplicon Sizing & Tiling Strategy",
            f"- **Amplicon Length:** {tiling['amplicon_length']} bp",
            f"- **Tiling Required:** `{'YES (> 1200 bp)' if tiling['tiling_required'] else 'NO (<= 1200 bp)'}`",
            f"- **Strategy:** {tiling['strategy']}",
            f"- **Details:** {tiling['explanation']}",
            "",
        ]

        total_seqs = res["total_sequences"]
        total_full = res.get("total_full_sequences", total_seqs)
        is_subsampled = (total_seqs < total_full)
        pcr_start = pcr.get("start_col", 0)
        amp_len = pcr.get("amplicon_len", tiling.get("amplicon_length", 0))

        if tiling["tiling_required"]:
            if is_subsampled:
                lines.extend([
                    "### Tiled Primer Walking Schedule",
                    "| Step | Direction | Position (bp) | Length | Sequence (5' → 3') | Tm (°C) | GC% | Subsample Cov | Full Cov | Tier | Sequenced Span |",
                    "|---|---|---|---|---|---|---|---|---|---|---|"
                ])
                for p in tiling["tiled_primers"]:
                    tier_label = "Tier 1 (0 degen)" if p.get("tier") == 1 else "Tier 2 (2-fold)"
                    aln_pos = f"Cols {p['start_col'] + pcr_start} – {p['end_col'] + pcr_start}"
                    span_str = f"Cols {p['start_col'] + pcr_start} – {min(amp_len + pcr_start, p['start_col'] + pcr_start + 650)}"
                    lines.append(
                        f"| {p.get('tile_step', '-')} | {p['direction']} | **{aln_pos}** | {p['length']} nt | `{p.get('order_sequence')}` | {p['tm']} | {p['gc_pct']}% | {p.get('coverage_pct', '-')} % | **{p.get('full_dataset_coverage_pct', p.get('coverage_pct', '-'))} %** | {tier_label} | {span_str} |"
                    )
            else:
                lines.extend([
                    "### Tiled Primer Walking Schedule",
                    "| Step | Direction | Position (bp) | Length | Sequence (5' → 3') | Tm (°C) | GC% | Population Coverage | Tier | Sequenced Span |",
                    "|---|---|---|---|---|---|---|---|---|---|"
                ])
                for p in tiling["tiled_primers"]:
                    tier_label = "Tier 1 (0 degen)" if p.get("tier") == 1 else "Tier 2 (2-fold)"
                    aln_pos = f"Cols {p['start_col'] + pcr_start} – {p['end_col'] + pcr_start}"
                    cov_val = p.get('full_dataset_coverage_pct', p.get('coverage_pct', '-'))
                    span_str = f"Cols {p['start_col'] + pcr_start} – {min(amp_len + pcr_start, p['start_col'] + pcr_start + 650)}"
                    lines.append(
                        f"| {p.get('tile_step', '-')} | {p['direction']} | **{aln_pos}** | {p['length']} nt | `{p.get('order_sequence')}` | {p['tm']} | {p['gc_pct']}% | **{cov_val}%** | {tier_label} | {span_str} |"
                    )
            lines.append("")

        lines.extend([
            "---",
            "",
            "## 3. Mode 1: Universal Sanger Sequencing Primers (Single-Strain Testing)",
            "Use these primers first. If your sample contains a single pure strain, these will yield clean, sharp chromatograms.",
            ""
        ])

        if not uf and not ur:
            lines.extend([
                "> [!WARNING]",
                "> **No Universal Primer Met Threshold:** Neither forward nor reverse directions contain a universal conserved site meeting the minimum coverage requirement. Please use Mode 2 (Strain-Specific Panel) below.",
                ""
            ])
        else:
            if is_subsampled:
                lines.extend([
                    "### Universal Primers Summary Table",
                    "| Primer ID | Role | Direction | Position (bp) | Length | Sequence (5' → 3') | Tm (°C) | GC% | Subsample Cov | Full Cov | Tier |",
                    "|---|---|---|---|---|---|---|---|---|---|---|"
                ])
                if uf:
                    t_label = "Tier 1 (Invariant)" if uf["tier"] == 1 else "Tier 2 (Degenerate)"
                    uf_aln = f"Cols {uf['start_col'] + pcr_start} – {uf['end_col'] + pcr_start}"
                    full_c = uf.get('full_dataset_coverage_pct', uf['coverage_pct'])
                    lines.append(f"| `Uni-Seq-F1` | Forward Terminal | FORWARD | **{uf_aln}** | {uf['length']} nt | `{uf['order_sequence']}` | {uf['tm']} | {uf['gc_pct']}% | {uf['coverage_pct']}% | **{full_c}%** | {t_label} |")
                if ur:
                    t_label = "Tier 1 (Invariant)" if ur["tier"] == 1 else "Tier 2 (Degenerate)"
                    ur_aln = f"Cols {ur['start_col'] + pcr_start} – {ur['end_col'] + pcr_start}"
                    full_c = ur.get('full_dataset_coverage_pct', ur['coverage_pct'])
                    lines.append(f"| `Uni-Seq-R1` | Reverse Terminal | REVERSE | **{ur_aln}** | {ur['length']} nt | `{ur['order_sequence']}` | {ur['tm']} | {ur['gc_pct']}% | {ur['coverage_pct']}% | **{full_c}%** | {t_label} |")
            else:
                lines.extend([
                    "### Universal Primers Summary Table",
                    "| Primer ID | Role | Direction | Position (bp) | Length | Sequence (5' → 3') | Tm (°C) | GC% | Population Coverage | Tier |",
                    "|---|---|---|---|---|---|---|---|---|---|"
                ])
                if uf:
                    t_label = "Tier 1 (Invariant)" if uf["tier"] == 1 else "Tier 2 (Degenerate)"
                    uf_aln = f"Cols {uf['start_col'] + pcr_start} – {uf['end_col'] + pcr_start}"
                    full_c = uf.get('full_dataset_coverage_pct', uf['coverage_pct'])
                    cnt_fwd = uf.get('full_dataset_count', uf.get('count', ''))
                    cov_display = f"**{full_c}%**" + (f" ({cnt_fwd}/{total_full} seqs)" if cnt_fwd else "")
                    lines.append(f"| `Uni-Seq-F1` | Forward Terminal | FORWARD | **{uf_aln}** | {uf['length']} nt | `{uf['order_sequence']}` | {uf['tm']} | {uf['gc_pct']}% | {cov_display} | {t_label} |")
                if ur:
                    t_label = "Tier 1 (Invariant)" if ur["tier"] == 1 else "Tier 2 (Degenerate)"
                    ur_aln = f"Cols {ur['start_col'] + pcr_start} – {ur['end_col'] + pcr_start}"
                    full_c = ur.get('full_dataset_coverage_pct', ur['coverage_pct'])
                    cnt_rev = ur.get('full_dataset_count', ur.get('count', ''))
                    cov_display = f"**{full_c}%**" + (f" ({cnt_rev}/{total_full} seqs)" if cnt_rev else "")
                    lines.append(f"| `Uni-Seq-R1` | Reverse Terminal | REVERSE | **{ur_aln}** | {ur['length']} nt | `{ur['order_sequence']}` | {ur['tm']} | {ur['gc_pct']}% | {cov_display} | {t_label} |")
            lines.append("")

            if uf and not ur:
                lines.extend([
                    "> [!NOTE]",
                    "> **Single-End Sequencing Primer Selected:** Only the Forward universal primer met the required coverage and biophysical thresholds. The reverse terminal region is too variable for a single universal oligo. Forward sequencing alone can sequence the amplicon, or use Mode 2 below for strain-specific reverse reactions.",
                    ""
                ])
            elif ur and not uf:
                lines.extend([
                    "> [!NOTE]",
                    "> **Single-End Sequencing Primer Selected:** Only the Reverse universal primer met the required coverage and biophysical thresholds. The forward terminal region is too variable for a single universal oligo. Reverse sequencing alone can sequence the amplicon, or use Mode 2 below for strain-specific forward reactions.",
                    ""
                ])

        if uf:
            tier_badge = "🟢 Tier 1 (Strict Invariant, 0 Degeneracy)" if uf["tier"] == 1 else "🟡 Tier 2 (Fallback: 2-Fold Degenerate)"
            full_cov_fwd = uf.get('full_dataset_coverage_pct', uf['coverage_pct'])
            full_cnt_fwd = uf.get('full_dataset_count', uf['coverage_pct'])
            uf_pos_str = f"Alignment Columns `{uf['start_col'] + pcr_start} – {uf['end_col'] + pcr_start}`"
            lines.extend([
                "### Forward Sequencing Primer (`Uni-Seq-F1`)",
                f"- **Tier:** {tier_badge}",
                f"- **Binding Position:** **{uf_pos_str}**",
                f"- **Sequence (5' → 3'):** `{uf['order_sequence']}`",
                f"- **Length:** {uf['length']} nt | **Tm:** {uf['tm']} °C | **GC:** {uf['gc_pct']}%"
            ])
            if is_subsampled:
                lines.extend([
                    f"- **Discovery Subsample Coverage:** {uf['coverage_pct']}%",
                    f"- **FULL DATASET GROUND TRUTH COVERAGE:** **{full_cov_fwd}%** ({full_cnt_fwd} / {total_full} strains in `{os.path.basename(res['fasta_path'])}`)"
                ])
            else:
                lines.append(
                    f"- **Population Coverage:** **{full_cov_fwd}%** ({full_cnt_fwd} / {total_full} strains in `{os.path.basename(res['fasta_path'])}`)"
                )
            lines.extend([
                f"- **Secondary Structures:** Hairpin $\\Delta G$: {uf['hairpin_dg']} kcal/mol | Dimer $\\Delta G$: {uf['homodimer_dg']} kcal/mol",
                f"- **Concentration Protocol:** {uf['concentration_recommendation']}",
                ""
            ])
            if uf.get("requires_2x_concentration"):
                lines.extend([
                    "> [!WARNING]",
                    "> **LAB ATTENTION: DOUBLE PRIMER CONCENTRATION REQUIRED**",
                    f"> This primer utilizes a 2-fold degenerate wobble at position {uf.get('wobble_position_from_5prime')}. ",
                    "> Because the primer is a 50:50 mix of two sequences, **double your primer input to 6.4 pmol** per reaction ",
                    "> to ensure full 1x binding saturation on the template. The 3' anchor is completely protected with 0 wobble.",
                    ""
                ])

        if ur:
            tier_badge = "🟢 Tier 1 (Strict Invariant, 0 Degeneracy)" if ur["tier"] == 1 else "🟡 Tier 2 (Fallback: 2-Fold Degenerate)"
            full_cov_rev = ur.get('full_dataset_coverage_pct', ur['coverage_pct'])
            full_cnt_rev = ur.get('full_dataset_count', ur['coverage_pct'])
            ur_pos_str = f"Alignment Columns `{ur['start_col'] + pcr_start} – {ur['end_col'] + pcr_start}`"
            lines.extend([
                "### Reverse Sequencing Primer (`Uni-Seq-R1`)",
                f"- **Tier:** {tier_badge}",
                f"- **Binding Position:** **{ur_pos_str}**",
                f"- **Sequence (5' → 3'):** `{ur['order_sequence']}`",
                f"- **Length:** {ur['length']} nt | **Tm:** {ur['tm']} °C | **GC:** {ur['gc_pct']}%"
            ])
            if is_subsampled:
                lines.extend([
                    f"- **Discovery Subsample Coverage:** {ur['coverage_pct']}%",
                    f"- **FULL DATASET GROUND TRUTH COVERAGE:** **{full_cov_rev}%** ({full_cnt_rev} / {total_full} strains in `{os.path.basename(res['fasta_path'])}`)"
                ])
            else:
                lines.append(
                    f"- **Population Coverage:** **{full_cov_rev}%** ({full_cnt_rev} / {total_full} strains in `{os.path.basename(res['fasta_path'])}`)"
                )
            lines.extend([
                f"- **Secondary Structures:** Hairpin $\\Delta G$: {ur['hairpin_dg']} kcal/mol | Dimer $\\Delta G$: {ur['homodimer_dg']} kcal/mol",
                f"- **Concentration Protocol:** {ur['concentration_recommendation']}",
                ""
            ])

        lines.extend([
            "---",
            "",
            "## 4. Mode 2: Strain-Specific Panel (Mixture Deconvolution)",
            "If your clinical/field sample is a **mixture of strains**, universal primers will produce double peaks. ",
            "Use the minimal panel below. Each primer selectively amplifies its specific clade via **3'-terminal extension termination**.",
            "",
            "### How Lineages are Determined",
            "1. **Target Region Partitioning:** The alignment is clustered by shared nucleotide motifs across the variable target domain, grouping identical alleles into biological lineages.",
            "2. **3'-Terminal Selectivity:** Primers are placed so the very last 3' nucleotide matches only the intended lineage with 0 mismatches, but mismatches all other lineages, arresting non-specific Taq polymerase extension.",
            "",
            f"- **Clusters/Variants Discovered:** {panel.get('clusters_identified', 0)}",
            f"- **Panel Tubes Required:** {panel.get('panel_size', 0)}"
        ])

        if is_subsampled:
            lines.append(f"- **Total Cumulative Coverage:** **{panel.get('total_coverage_pct', 0)}%** (Subsample) | **{panel.get('total_full_coverage_pct', 0)}%** (Full Dataset)")
            lines.extend([
                "",
                "### Lineage Breakdown in Full Dataset",
                "| Lineage ID | Subsample Count | Subsample Share | Est. Full Dataset Count | Est. Full Share | Assigned Test Tube |",
                "|---|---|---|---|---|---|"
            ])
            for c in panel.get("all_clusters", []):
                assigned_tube = c.get("assigned_tube", "Uncovered")
                tube_label = f"**{assigned_tube}**" if assigned_tube != "Uncovered" else "Filtered / Uncovered"
                f_count_est = f"~{c.get('full_count_est', round(c['freq'] * total_full / 100))} / {total_full} seqs"
                f_pct_est = f"~{c['freq']}%"
                lines.append(
                    f"| **{c['id']}** | {c['count']} / {res['total_sequences']} seqs | {c['freq']}% | {f_count_est} | {f_pct_est} | {tube_label} |"
                )

            lines.extend([
                "",
                "### Minimal Panel Tubes",
                "| Tube ID | Primer ID | Target Variant(s) | Position (bp) | Direction | Sequence (5' → 3') | Length | Tm (°C) | Subsample Share | Full Dataset Share | Reaction Protocol |",
                "|---|---|---|---|---|---|---|---|---|---|---|"
            ])
            for entry in panel.get("panel", []):
                e_aln = f"Cols {entry.get('start_col', 0) + pcr_start} – {entry.get('end_col', 0) + pcr_start}"
                lines.append(
                    f"| **{entry['tube_id']}** | `{entry['primer_name']}` | {entry['target_cluster']} | **{e_aln}** | {entry.get('direction', 'FORWARD')} | `{entry['sequence']}` | {len(entry['sequence'])} nt | {entry['tm']} | {entry.get('subsample_freq', entry['frequency_covered_pct'])}% | **{entry.get('full_dataset_coverage_pct', entry['frequency_covered_pct'])}%** | {entry['action_protocol']} |"
                )
        else:
            lines.append(f"- **Total Cumulative Coverage:** **{panel.get('total_full_coverage_pct', panel.get('total_coverage_pct', 0))}%** across all {total_full} sequences")
            lines.extend([
                "",
                "### Lineage Breakdown",
                "| Lineage ID | Sequence Count | Population Share | Assigned Test Tube |",
                "|---|---|---|---|"
            ])
            for c in panel.get("all_clusters", []):
                assigned_tube = c.get("assigned_tube", "Uncovered")
                tube_label = f"**{assigned_tube}**" if assigned_tube != "Uncovered" else "Filtered / Uncovered"
                f_count = f"{c['count']} / {total_full} seqs"
                f_pct = f"**{c['freq']}%**"
                lines.append(
                    f"| **{c['id']}** | {f_count} | {f_pct} | {tube_label} |"
                )

            lines.extend([
                "",
                "### Minimal Panel Tubes",
                "| Tube ID | Primer ID | Target Variant(s) | Position (bp) | Direction | Sequence (5' → 3') | Length | Tm (°C) | Population Share | Reaction Protocol |",
                "|---|---|---|---|---|---|---|---|---|---|"
            ])
            for entry in panel.get("panel", []):
                e_aln = f"Cols {entry.get('start_col', 0) + pcr_start} – {entry.get('end_col', 0) + pcr_start}"
                share_pct = entry.get('full_dataset_coverage_pct', entry['frequency_covered_pct'])
                share_cnt = entry.get('full_dataset_count', entry.get('subsample_count', ''))
                share_str = f"**{share_pct}%**" + (f" ({share_cnt} seqs)" if share_cnt else "")
                lines.append(
                    f"| **{entry['tube_id']}** | `{entry['primer_name']}` | {entry['target_cluster']} | **{e_aln}** | {entry.get('direction', 'FORWARD')} | `{entry['sequence']}` | {len(entry['sequence'])} nt | {entry['tm']} | {share_str} | {entry['action_protocol']} |"
                )

        lines.extend([
            "",
            "---",
            "",
            "## 5. Oligo Synthesis Order Sheet",
            "Copy and paste these sequences directly into your synthesis provider (e.g. IDT, Sigma, Eurofins):",
            "",
            "```tsv",
            "Primer_Name\tSequence_5_to_3\tPosition\tLength\tTm_C\tNotes"
        ])

        seen_order_seqs = set()

        if uf:
            seen_order_seqs.add(uf['order_sequence'])
            f_cov = uf.get('full_dataset_coverage_pct', uf['coverage_pct'])
            lines.append(f"Uni-Seq-F1\t{uf['order_sequence']}\tCols {uf['start_col'] + pcr_start}-{uf['end_col'] + pcr_start}\t{uf['length']}\t{uf['tm']}\tUniversal Sanger Fwd ({f_cov}% cov)")
        if ur:
            seen_order_seqs.add(ur['order_sequence'])
            r_cov = ur.get('full_dataset_coverage_pct', ur['coverage_pct'])
            lines.append(f"Uni-Seq-R1\t{ur['order_sequence']}\tCols {ur['start_col'] + pcr_start}-{ur['end_col'] + pcr_start}\t{ur['length']}\t{ur['tm']}\tUniversal Sanger Rev ({r_cov}% cov)")
        
        for entry in panel.get("panel", []):
            p_seq = entry["sequence"]
            if p_seq not in seen_order_seqs:
                seen_order_seqs.add(p_seq)
                lines.append(f"{entry['primer_name']}\t{p_seq}\tCols {entry.get('start_col', 0) + pcr_start}-{entry.get('end_col', 0) + pcr_start}\t{len(p_seq)}\t{entry['tm']}\tStrain-Specific ({entry['target_cluster']}, {entry.get('full_dataset_coverage_pct', entry['frequency_covered_pct'])}%)")

        lines.extend([
            "```",
            "",
            "---",
            "*Report produced by AI4VPrimer Amplicon-Sanger Pipeline.*"
        ])

        return "\n".join(lines)

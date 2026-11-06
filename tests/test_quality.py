#!/usr/bin/env python3
"""Tests of NCBI's quality reports in Eukore's genome ranking (eukore/quality.py): an assembly that is mostly another
organism's (by its RefSeq copy's FCS-GX summary) or fails NCBI's ANI species check never taken; one contaminated by
another eukaryote, or inconclusive, flagged; of two otherwise equal, the one without contaminant taken.

    python3 -m unittest discover -s tests -t tests
"""

import gzip
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import genomes, ncbi as ncbi_mod, quality, sources  # noqa: E402
import fixtures as fx  # noqa: E402

FCS_HEAD = ("assembly_accession", "fcs_gx_software_rev", "fcs_gx_db_rev", "fcs_gx_run_date", "asserted_division",
            "top_inferred_division", "fcs_gx_aggregate_coverage", "total_contaminant_accessions",
            "total_contaminant_length", "fraction_contaminated", "exclude_ct", "fix_ct", "trim_ct", "rare_ct",
            "review_ct", "info_ct", "adaptor_vector_ct", "organelle_ct", "primates", "anml_minus_primates", "fung",
            "plnt", "prst", "prok", "arch", "virs", "synt", "adaptor_vector")
ANI_HEAD = ("genbank-accession", "refseq-accession", "taxid", "species-taxid", "organism-name", "species-name",
            "assembly-name", "assembly-type-category", "excluded-from-refseq", "declared-type-assembly",
            "declared-type-organism-name", "declared-type-category", "declared-type-ANI", "declared-type-qcoverage",
            "declared-type-scoverage", "best-match-type-assembly", "best-match-species-taxid", "best-match-species-name",
            "best-match-type-category", "best-match-type-ANI", "best-match-type-qcoverage",
            "best-match-type-scoverage", "best-match-status", "comment", "taxonomy-check-status")


def fcs_row(acc, fraction=0.0, bp=0, **kinds):
    r = dict.fromkeys(FCS_HEAD, "0")
    r.update(assembly_accession=acc, fraction_contaminated=str(fraction), total_contaminant_length=str(bp),
             **{k: str(v) for k, v in kinds.items()})
    return "\t".join(r[c] for c in FCS_HEAD) + "\n"


def reports(src_dir):
    folder = quality.folder_of(src_dir)
    os.makedirs(folder)
    with gzip.open(os.path.join(folder, "fcs_summary_genbank.txt.gz"), "wt") as fh:
        fh.write("#" + "\t".join(FCS_HEAD) + "\n")
        fh.write(fcs_row("GCA_000146045.2"))                         # the yeast's GenBank copy: clean,
        fh.write(fcs_row("GCA_000002516.1", 0.0001, 500))            # K. lactis' newer: a little contaminant
        fh.write(fcs_row("GCA_000182965.3", 0.001, 200_000, fung=200_000))  # another fungus' contigs
    with gzip.open(os.path.join(folder, "fcs_summary_refseq.txt.gz"), "wt") as fh:
        fh.write("#" + "\t".join(FCS_HEAD) + "\n")
        fh.write(fcs_row("GCF_000146045.2", 0.6, 7_000_000, prok=7_000_000))  # but its RefSeq copy (fetched) not
    with open(os.path.join(folder, "ANI_report_saccharomycotina.txt"), "w") as fh:
        fh.write("# " + "\t".join(ANI_HEAD) + "\n")
        for acc, status in (("GCA_000000013.1", "Failed"), ("GCA_000182965.3", "Inconclusive"),
                            ("GCA_000002515.1", "OK")):
            r = dict.fromkeys(ANI_HEAD, "na")
            r.update({"genbank-accession": acc, "taxonomy-check-status": status})
            fh.write("\t".join(r[c] for c in ANI_HEAD) + "\n")


class QualityTest(unittest.TestCase):
    def test_flags_and_ranking(self):
        with tempfile.TemporaryDirectory() as tmp:
            ncbi = ncbi_mod.Ncbi(fx.taxdump(os.path.join(tmp, "taxdump.tar.gz")))
            gb = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "genbank_fungi.txt"), fx.GENBANK)])
            rs = sources.read_assemblies([fx.assembly_table(os.path.join(tmp, "refseq_fungi.txt"), fx.REFSEQ)])
            before = {sp: best["accession"] for sp, (best, _o) in genomes.representatives(gb, ncbi).items()}
            self.assertEqual((before[28985], before[27291], before[4932]),
                             ("GCA_000002516.1", "GCA_000000013.1", "GCA_000146045.2"))
            src = os.path.join(tmp, "sources")
            reports(src)
            self.assertEqual(len(quality.paths(src)), 3)
            counts = quality.annotate(gb, rs, src)
            self.assertEqual(counts, {quality.FCS_FOREIGN: 1, quality.FCS_CONTAMINATED: 1, quality.ANI_FAILED: 1,
                                      quality.ANI_INCONCLUSIVE: 1})
            row = {r["accession"]: r for r in gb}
            self.assertEqual(row["GCA_000182965.3"]["flags"], [quality.ANI_INCONCLUSIVE, quality.FCS_CONTAMINATED])
            self.assertEqual((row["GCA_000002516.1"]["fcs_bp"], row["GCA_000002515.1"]["fcs_bp"]), (500, None))
            after = {sp: best["accession"] for sp, (best, _o) in genomes.representatives(gb, ncbi).items()}
            self.assertEqual(after[28985], "GCA_000002515.1")   # of two equal, the one without contaminant
            self.assertEqual(after[27291], "GCA_000000012.1")   # the ANI failure never; a demoted one instead
            self.assertEqual(after[4932], "GCA_000000003.1")    # its RefSeq copy mostly another organism's
            self.assertEqual(after[5476], "GCA_000182965.3")    # flagged, but the species' only one
            jobs = genomes.jobs(ncbi, gb, rs)
            self.assertIn(quality.FCS_CONTAMINATED, jobs["GCA_000182965.3"]["flags"])

    def test_fetch(self):
        with tempfile.TemporaryDirectory() as tmp:
            got = []

            def download(url, dest):
                got.append(url)
                open(dest, "w").close()
            quality.fetch(tmp, download, lambda _t: None)
            quality.fetch(tmp, download, lambda _t: None)  # kept when there
            self.assertEqual(got, [f"{quality.REPORTS}/{n}" for n in quality.FILES])
            self.assertEqual(quality.annotate([], [], os.path.join(tmp, "none")), {})


if __name__ == "__main__":
    unittest.main()

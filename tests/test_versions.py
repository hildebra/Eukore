#!/usr/bin/env python3
"""Tests of the steps' versions (versions.py; 2026-11-19): a kept result is made again where its step's MAJOR changed,
never because the code's text did; the version constants the modules key their kept results by are the registry's,
with the values they had before it (no kept result of a build made again by the change).

    python3 -m unittest discover -s tests -t tests
"""

import os
import re
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from eukore import audit, genomes, markers, organelle, qc, refscreen, rrna, versions  # noqa: E402


class TestVersions(unittest.TestCase):
    def test_every_version_is_major_minor(self):
        for step, v in versions.STEPS.items():
            self.assertRegex(v, r"^[1-9]\d*\.\d+$", step)
            self.assertRegex(step, r"^[a-z0-9_]+$")

    def test_the_constants_keep_their_values(self):
        # the values before 2026-11-19, so that every kept result stays good: a change here makes them again
        import eukore_build  # noqa: PLC0415
        self.assertEqual(genomes.GENOME_RECORD, 4)
        self.assertEqual((qc.SEARCH_RECORD, qc.CM_RECORD, qc.CMALIGN_RECORD, qc.MAP_RECORD, qc.UCHIME_RECORD),
                         (1, 1, 1, 1, 2))
        self.assertEqual(audit.AUDIT_RECORD, 1)
        self.assertEqual(rrna.JOIN_RECORD, 1)
        self.assertEqual(refscreen.VERSION, "1")
        self.assertEqual((organelle.CACHE_VERSION, organelle.CANDIDATE_RECORD, organelle.ONE_ORGANELLE,
                          organelle.CHOICE_RECORD), (5, "candidates1", "one_organelle1", 1))
        self.assertEqual(eukore_build.PUBLISHED_REGIONS, "regions.4.tsv")
        self.assertEqual(markers.RECIPE, 1)

    def test_tag_none_at_major_1(self):
        with mock.patch.dict(versions.STEPS, {"x": "1.7"}):
            self.assertEqual(versions.tag("x"), "")
        with mock.patch.dict(versions.STEPS, {"x": "3.0"}):
            self.assertEqual(versions.tag("x"), ".v3")
            self.assertEqual(versions.tag("x", "_"), "_v3")

    def test_digest_follows_majors_alone(self):
        before = versions.digest()
        step = next(iter(versions.STEPS))
        n = versions.major(step)
        with mock.patch.dict(versions.STEPS, {step: f"{n}.99"}):
            self.assertEqual(versions.digest(), before)  # a MINOR raised: the same
        with mock.patch.dict(versions.STEPS, {step: f"{n + 1}.0"}):
            self.assertNotEqual(versions.digest(), before)  # a MAJOR raised: another

    def test_no_module_hashes_the_code(self):
        # since 2026-11-19 no kept result is keyed by the code's text: no module reads a .py file's bytes
        root = os.path.dirname(HERE)
        paths = [os.path.join(root, "eukore_build.py")] + [
            os.path.join(root, "eukore", n) for n in os.listdir(os.path.join(root, "eukore")) if n.endswith(".py")]
        for path in paths:
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
            self.assertIsNone(re.search(r"glob\([^)]*\*\.py", text), path)
            self.assertNotIn("getsource", text, path)
            self.assertNotIn("__code__", text, path)

    def test_the_fingerprint_follows_majors_not_the_code(self):
        # the process stage's nothing_left (since 2026-11-16) held a digest of every module's text until 2026-11-19:
        # now a MINOR raised, or any edit, leaves it as it was; a MAJOR raised changes it
        import eukore_build  # noqa: PLC0415
        import types  # noqa: PLC0415
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "sources")
            os.makedirs(src)
            os.makedirs(os.path.join(tmp, "genomes"))
            opts = types.SimpleNamespace(work=tmp, genomes=os.path.join(tmp, "genomes"))
            before = eukore_build.process_fingerprint(opts, src, ["--stage", "process"])
            with mock.patch.object(eukore_build, "__file__", os.path.join(tmp, "elsewhere.py")):
                self.assertEqual(eukore_build.process_fingerprint(opts, src, ["--stage", "process"]), before)
            with mock.patch.dict(versions.STEPS, {"genome": f"{versions.major('genome')}.99"}):
                self.assertEqual(eukore_build.process_fingerprint(opts, src, ["--stage", "process"]), before)
            with mock.patch.dict(versions.STEPS, {"genome": f"{versions.major('genome') + 1}.0"}):
                self.assertNotEqual(eukore_build.process_fingerprint(opts, src, ["--stage", "process"]), before)

    def test_the_manifest_and_banner_name_them(self):
        self.assertIn("genome 4.", versions.text())
        self.assertEqual(versions.text().count(","), len(versions.STEPS) - 1)


if __name__ == "__main__":
    unittest.main()

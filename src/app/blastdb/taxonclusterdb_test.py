#!/usr/bin/env python3
"""Integration tests; run from the checkout with BLAST binaries on PATH."""
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import unittest


@unittest.skipUnless(shutil.which("makeblastdb") and shutil.which("taxonclusterdb"),
                     "makeblastdb and taxonclusterdb must be on PATH")
class TestTaxonClusterDb(unittest.TestCase):
    def setUp(self):
        self.work = Path.cwd() / (".taxonclusterdb-test-%d" % os.getpid())
        self.work.mkdir()
        self.addCleanup(shutil.rmtree, self.work)
        (self.work / "input.fa").write_text(
            ">a first\x01a2 alternate\nACGTACGTACGTACGN\n>b\nACGTACGTACGTACGN\n"
            ">c\nACGTACGTACGTACGN\n>d\nACGTACGTACGTACGA\n"
            ">e\nACGTACGTACGTACGN\n")
        (self.work / "taxids.txt").write_text("a 11\na2 11\nb 11\nc 12\nd 11\n")

    def make_db(self, version):
        subprocess.run(["makeblastdb", "-in", "input.fa", "-dbtype", "nucl",
                        "-parse_seqids", "-taxid_map", "taxids.txt",
                        "-blastdb_version", str(version), "-out", "db"],
                       cwd=self.work, check=True, capture_output=True, text=True)

    def run_cluster(self, *args, env=None):
        return subprocess.run(["taxonclusterdb", "-db", "db", *args], cwd=self.work,
                              capture_output=True, text=True, env=env)

    def test_exact_taxonomy_and_singletons(self):
        self.make_db(5)
        first = self.run_cluster()
        self.assertEqual(first.returncode, 0, first.stderr)
        self.assertEqual(first.stdout,
                         "representative_oid\tmember_oid\n"
                         "0\t0\n0\t1\n2\t2\n3\t3\n4\t4\n")
        second = self.run_cluster()
        self.assertEqual(second.returncode, 0, second.stderr)
        self.assertEqual(first.stdout, second.stdout)

    def test_rejects_v4(self):
        self.make_db(4)
        result = self.run_cluster()
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(result.stdout, "")

    def test_near_is_opt_in_and_writes_representatives(self):
        rng = random.Random(42)
        sequence = "".join(rng.choice("ACGT") for _ in range(300))
        substitution = sequence[:100] + ("C" if sequence[100] == "A" else "A") + sequence[101:]
        insertion = sequence[:150] + "G" + sequence[150:]
        (self.work / "input.fa").write_text(
            ">a\n" + sequence + "\n>b\n" + substitution + "\n>c\n" + insertion
            + "\n>d\n" + substitution + "\n")
        (self.work / "taxids.txt").write_text("a 11\nb 11\nc 11\nd 12\n")
        self.make_db(5)
        exact = self.run_cluster()
        self.assertEqual(exact.returncode, 0, exact.stderr)
        self.assertEqual(exact.stdout, "representative_oid\tmember_oid\n"
                         "0\t0\n1\t1\n2\t2\n3\t3\n")
        near = self.run_cluster("-near", "-out_dir", "near-output")
        self.assertEqual(near.returncode, 0, near.stderr)
        self.assertEqual((self.work / "near-output/data/clusters.tsv").read_text(),
                         "representative_oid\tmember_oid\n"
                         "0\t0\n0\t1\n0\t2\n3\t3\n")

    def test_configurable_edit_percent(self):
        rng = random.Random(42)
        sequence = "".join(rng.choice("ACGT") for _ in range(300))
        changed = list(sequence)
        for i in (0, 50, 100, 150):
            changed[i] = "C" if changed[i] == "A" else "A"
        (self.work / "input.fa").write_text(">a\n" + sequence + "\n>b\n"
                                            + "".join(changed) + "\n")
        (self.work / "taxids.txt").write_text("a 11\nb 11\n")
        self.make_db(5)
        strict = self.run_cluster("-near")
        self.assertEqual(strict.returncode, 0, strict.stderr)
        self.assertEqual(strict.stdout, "representative_oid\tmember_oid\n0\t0\n1\t1\n")
        relaxed = self.run_cluster("-near", "-max_edit_percent", "2")
        self.assertEqual(relaxed.returncode, 0, relaxed.stderr)
        self.assertEqual(relaxed.stdout, "representative_oid\tmember_oid\n0\t0\n0\t1\n")
        invalid = self.run_cluster("-near", "-max_edit_percent", "100.1")
        self.assertNotEqual(invalid.returncode, 0)
        self.assertEqual(invalid.stdout, "")

    @unittest.skipUnless(shutil.which("blastdbcmd"), "blastdbcmd must be on PATH")
    def test_taxonomy_lca_and_output_taxids(self):
        (self.work / "nodes.dmp").write_text(
            "1 | 1 | no rank |\n10 | 1 | species |\n11 | 10 | no rank |\n"
            "12 | 11 | no rank |\n13 | 12 | no rank |\n14 | 12 | strain |\n"
            "15 | 11 | strain |\n18 | 13 | no rank |\n")
        (self.work / "input.fa").write_text(
            "".join(">" + name + "\n" + base * 300 + "\n"
                    for name, base in (("a", "A"), ("b", "A"), ("c", "A"),
                                       ("d", "A"), ("e", "A"), ("g", "C"),
                                       ("h", "C"))))
        (self.work / "taxids.txt").write_text(
            "a 13\nb 18\nc 14\nd 15\ne 11\ng 13\nh 18\n")
        self.make_db(5)
        default = self.run_cluster()
        self.assertEqual(default.returncode, 0, default.stderr)
        self.assertIn("1\t1\n", default.stdout)
        result = self.run_cluster("-tax_nodes", "nodes.dmp", "-out_dir", "ranked")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.work / "ranked/data/clusters.tsv").read_text(),
                         "representative_oid\tmember_oid\n"
                         "0\t0\n0\t1\n0\t2\n3\t3\n4\t4\n5\t5\n5\t6\n")
        for accession, expected in (("a", "12"), ("g", "13"), ("d", "15"),
                                    ("e", "11")):
            taxid = subprocess.run(
                ["blastdbcmd", "-db", "ranked/data/db", "-entry", accession,
                 "-outfmt", "%T"], cwd=self.work, check=True,
                capture_output=True, text=True)
            self.assertEqual(taxid.stdout.strip(), expected)

    @unittest.skipUnless(shutil.which("blastdbcmd"), "blastdbcmd must be on PATH")
    def test_database_and_multidefline(self):
        self.make_db(5)
        result = self.run_cluster("-out_dir", "output")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout, "")
        data = self.work / "output" / "data"
        self.assertEqual((data / "clusters.tsv").read_text(),
                         "representative_oid\tmember_oid\n"
                         "0\t0\n0\t1\n2\t2\n3\t3\n4\t4\n")
        db = str(data / "db")
        all_records = subprocess.run(
            ["blastdbcmd", "-db", db, "-entry", "all", "-outfmt", "%o %s %T"],
            cwd=self.work, capture_output=True, text=True, check=True)
        self.assertEqual(len(all_records.stdout.splitlines()), 4)
        alternate = subprocess.run(
            ["blastdbcmd", "-db", db, "-entry", "a2", "-outfmt", "%s"],
            cwd=self.work, capture_output=True, text=True, check=True)
        self.assertEqual(alternate.stdout.strip(), "ACGTACGTACGTACGN")
        singleton = subprocess.run(
            ["blastdbcmd", "-db", db, "-entry", "e", "-outfmt", "%s"],
            cwd=self.work, capture_output=True, text=True, check=True)
        self.assertEqual(singleton.stdout.strip(), "ACGTACGTACGTACGN")

    def test_existing_directory_never_overwritten(self):
        self.make_db(5)
        before = {p.name: p.read_bytes() for p in self.work.glob("db.*")}
        result = self.run_cluster("-out_dir", ".")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.work.glob("db.*")})
        self.assertFalse((self.work / "data").exists())
        (self.work / "existing-link").symlink_to(self.work, target_is_directory=True)
        result = self.run_cluster("-out_dir", "existing-link")
        self.assertNotEqual(result.returncode, 0)
        self.assertEqual(before, {p.name: p.read_bytes() for p in self.work.glob("db.*")})

    def test_writer_failure_not_published(self):
        self.make_db(5)
        env = dict(os.environ, BLASTDB_LMDB_MAP_SIZE="1")
        result = self.run_cluster("-out_dir", "failed", env=env)
        self.assertNotEqual(result.returncode, 0, result.stderr)
        self.assertFalse((self.work / "failed" / "data").exists())

    @unittest.skipUnless(shutil.which("dustmasker") and shutil.which("blastdbcmd"),
                         "dustmasker and blastdbcmd must be on PATH")
    def test_mask_and_custom_algorithm_metadata(self):
        (self.work / "input.fa").write_text(">a\n" + "A" * 100 + "\n"
                                             ">b\n" + "A" * 100 + "\n")
        (self.work / "taxids.txt").write_text("a 11\nb 11\n")
        subprocess.run(["dustmasker", "-in", "input.fa", "-parse_seqids",
                        "-outfmt", "maskinfo_asn1_bin", "-out", "masks.asnb"],
                       cwd=self.work, check=True, capture_output=True, text=True)
        subprocess.run(["makeblastdb", "-in", "input.fa", "-dbtype", "nucl",
                        "-parse_seqids", "-taxid_map", "taxids.txt",
                        "-blastdb_version", "5", "-mask_data", "masks.asnb",
                        "-mask_id", "custom", "-mask_desc", "custom details",
                        "-out", "db"],
                       cwd=self.work, check=True, capture_output=True, text=True)
        result = self.run_cluster("-out_dir", "masked")
        self.assertEqual(result.returncode, 0, result.stderr)

        def masked_ranges(db):
            info = subprocess.run(["blastdbcmd", "-db", db, "-info"],
                                  cwd=self.work, check=True, capture_output=True, text=True)
            match = re.search(r"(?m)^\s*(\d+)\s+custom\s+", info.stdout)
            self.assertIsNotNone(match, info.stdout)
            return subprocess.run(
                ["blastdbcmd", "-db", db, "-entry", "a",
                 "-outfmt", "%m" + match.group(1)],
                cwd=self.work, check=True, capture_output=True, text=True).stdout

        original = masked_ranges("db")
        self.assertNotEqual(original.strip(), "none")
        self.assertEqual(original, masked_ranges("masked/data/db"))


if __name__ == "__main__":
    unittest.main()

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

    def run_cluster(self, *args, env=None, preexec_fn=None):
        return subprocess.run(["taxonclusterdb", "-db", "db", *args], cwd=self.work,
                              capture_output=True, text=True, env=env,
                              preexec_fn=preexec_fn)

    def test_help(self):
        result = subprocess.run(["taxonclusterdb", "-help"], cwd=self.work,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("-tax_nodes", result.stdout)

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

    @unittest.skipUnless(shutil.which("blastdb_aliastool") and shutil.which("blastdbcmd"),
                         "blastdb_aliastool and blastdbcmd must be on PATH")
    def test_filtered_alias_retains_only_selected_oids(self):
        self.make_db(5)
        (self.work / "ids").write_text("b\nc\n")
        subprocess.run(["blastdb_aliastool", "-db", "db", "-dbtype", "nucl",
                        "-seqidlist", "ids", "-out", "selected"],
                       cwd=self.work, check=True, capture_output=True, text=True)
        result = subprocess.run(
            ["taxonclusterdb", "-db", "selected", "-out_dir", "selected-output"],
            cwd=self.work, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.work / "selected-output/data/clusters.tsv").read_text(),
                         "representative_oid\tmember_oid\n1\t1\n2\t2\n")
        records = subprocess.run(
            ["blastdbcmd", "-db", "selected-output/data/db", "-entry", "all",
             "-outfmt", "%a %T"], cwd=self.work, check=True, capture_output=True, text=True)
        self.assertEqual(set(records.stdout.splitlines()), {"b 11", "c 12"})

    @unittest.skipUnless(shutil.which("blastdb_aliastool") and shutil.which("blastdbcmd"),
                         "blastdb_aliastool and blastdbcmd must be on PATH")
    def test_multi_volume_alias(self):
        self.make_db(5)
        (self.work / "input2.fa").write_text(
            ">f\nACGTACGTACGTACGN\n>g\nACGTACGTACGTACGC\n")
        (self.work / "taxids2.txt").write_text("f 11\ng 11\n")
        subprocess.run(["makeblastdb", "-in", "input2.fa", "-dbtype", "nucl",
                        "-parse_seqids", "-taxid_map", "taxids2.txt",
                        "-blastdb_version", "5", "-out", "db2"],
                       cwd=self.work, check=True, capture_output=True, text=True)
        subprocess.run(["blastdb_aliastool", "-dblist", "db db2", "-dbtype", "nucl",
                        "-title", "combined", "-out", "combined"], cwd=self.work, check=True,
                       capture_output=True, text=True)
        result = subprocess.run(
            ["taxonclusterdb", "-db", "combined", "-out_dir", "combined-output"],
            cwd=self.work, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.work / "combined-output/data/clusters.tsv").read_text(),
                         "representative_oid\tmember_oid\n"
                         "0\t0\n0\t1\n2\t2\n3\t3\n4\t4\n0\t5\n6\t6\n")
        info = subprocess.run(
            ["blastdbcmd", "-db", "combined-output/data/db", "-info"],
            cwd=self.work, check=True, capture_output=True, text=True)
        self.assertRegex(info.stdout, r"\b5 sequences\b")

    @unittest.skipUnless(shutil.which("blastn"), "blastn must be on PATH")
    def test_megablast_search_of_output(self):
        rng = random.Random(42)
        first = "".join(rng.choice("ACGT") for _ in range(300))
        second = "".join(rng.choice("ACGT") for _ in range(300))
        changed = first[:100] + ("C" if first[100] == "A" else "A") + first[101:]
        (self.work / "input.fa").write_text(
            ">a\n" + first + "\n>b\n" + changed + "\n>c\n" + second + "\n")
        (self.work / "taxids.txt").write_text("a 11\nb 11\nc 12\n")
        (self.work / "queries.fa").write_text(">q1\n" + first + "\n>q2\n" + second + "\n")
        self.make_db(5)
        result = self.run_cluster("-near", "-out_dir", "search-output")
        self.assertEqual(result.returncode, 0, result.stderr)

        def best_taxids(db):
            output = subprocess.run(
                ["blastn", "-task", "megablast", "-db", db, "-query", "queries.fa",
                 "-dust", "no", "-outfmt", "6 qseqid staxids bitscore"],
                cwd=self.work, check=True, capture_output=True, text=True)
            best = {}
            for line in output.stdout.splitlines():
                query, taxids, score = line.split("\t")
                score = float(score)
                if query not in best or score > best[query][0]:
                    best[query] = (score, set(taxids.split(";")))
                elif score == best[query][0]:
                    best[query][1].update(taxids.split(";"))
            return {query: taxa for query, (_, taxa) in best.items()}

        self.assertEqual(best_taxids("db"), {"q1": {"11"}, "q2": {"12"}})
        self.assertEqual(best_taxids("search-output/data/db"), best_taxids("db"))

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
        oids = [int(line.split()[0]) for line in all_records.stdout.splitlines()]
        self.assertEqual(sorted(set(oids)), [0, 1, 2, 3])
        self.assertEqual(oids.count(0), 2)
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

    @unittest.skipUnless(shutil.which("blastdbcmd"), "blastdbcmd must be on PATH")
    def test_output_path_with_spaces(self):
        self.make_db(5)
        result = self.run_cluster("-out_dir", "output with spaces")
        self.assertEqual(result.returncode, 0, result.stderr)
        records = subprocess.run(
            ["blastdbcmd", "-db", '"output with spaces/data/db"', "-entry", "all",
             "-outfmt", "%o"], cwd=self.work, check=True, capture_output=True, text=True)
        self.assertEqual(set(records.stdout.splitlines()), {"0", "1", "2", "3"})

    def test_input_path_with_spaces(self):
        self.make_db(5)
        directory = self.work / "source with spaces"
        directory.mkdir()
        for path in self.work.glob("db.*"):
            path.rename(directory / path.name)
        result = subprocess.run(
            ["taxonclusterdb", "-db", "source with spaces/db", "-out_dir", "output"],
            cwd=self.work, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.work / "output/data/clusters.tsv").read_text(),
                         "representative_oid\tmember_oid\n"
                         "0\t0\n0\t1\n2\t2\n3\t3\n4\t4\n")

    @unittest.skipUnless(os.name == "posix", "requires POSIX file size limits")
    def test_writer_failure_not_published(self):
        import resource

        self.make_db(5)
        def limit_output_file_size():
            resource.setrlimit(resource.RLIMIT_FSIZE, (256, 256))

        result = self.run_cluster("-out_dir", "failed", preexec_fn=limit_output_file_size)
        self.assertNotEqual(result.returncode, 0, result.stderr)
        self.assertIn("Building a new DB", result.stderr)
        self.assertFalse((self.work / "failed" / "data").exists())

    @unittest.skipUnless(shutil.which("dustmasker") and shutil.which("blastdbcmd"),
                         "dustmasker and blastdbcmd must be on PATH")
    def test_standard_dust_metadata_and_masks(self):
        (self.work / "input.fa").write_text(">a\n" + "A" * 100 + "\n"
                                             ">b\n" + "A" * 100 + "\n")
        (self.work / "taxids.txt").write_text("a 11\nb 11\n")
        subprocess.run(["dustmasker", "-in", "input.fa", "-parse_seqids",
                        "-outfmt", "maskinfo_asn1_bin", "-out", "masks.asnb"],
                       cwd=self.work, check=True, capture_output=True, text=True)
        subprocess.run(["makeblastdb", "-in", "input.fa", "-dbtype", "nucl",
                        "-parse_seqids", "-taxid_map", "taxids.txt",
                        "-blastdb_version", "5", "-mask_data", "masks.asnb", "-out", "db"],
                       cwd=self.work, check=True, capture_output=True, text=True)
        result = self.run_cluster("-out_dir", "standard-mask")
        self.assertEqual(result.returncode, 0, result.stderr)
        info = subprocess.run(["blastdbcmd", "-db", "standard-mask/data/db", "-info"],
                              cwd=self.work, check=True, capture_output=True, text=True)
        self.assertRegex(info.stdout, r"(?m)^\s*\d+\s+dust\s+")

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

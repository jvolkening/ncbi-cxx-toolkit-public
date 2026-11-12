/* This is a United States Government Work, freely available for public use. */
#include "taxon_exact_cluster.hpp"
#include "taxon_lineage.hpp"

#include <cassert>
#include <map>
#include <random>
#include <sstream>
#include <stdexcept>

static std::size_t EditDistance(const std::string& a, const std::string& b)
{
    std::vector<std::size_t> row(b.size() + 1);
    for (std::size_t j = 0; j <= b.size(); ++j) row[j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i) {
        std::size_t diagonal = row[0];
        row[0] = i;
        for (std::size_t j = 1; j <= b.size(); ++j) {
            std::size_t old = row[j];
            row[j] = std::min({row[j] + 1, row[j - 1] + 1,
                               diagonal + (a[i - 1] != b[j - 1])});
            diagonal = old;
        }
    }
    return row.back();
}

int main()
{
    using namespace taxon_exact;
    std::map<int, SRecord> records;
    auto add = [&](int oid, const std::string& seq, std::vector<int> tax,
                   TMasks masks, bool complete = true) {
        SRecord r;
        r.sequence = seq;
        r.taxids = std::move(tax);
        r.masks = std::move(masks);
        r.complete_taxonomy = complete;
        records.emplace(oid, std::move(r));
    };
    add(0, "ACGN", {2, 1}, {{{0, 2}}});
    add(1, "ACGN", {1, 2, 1}, {{{0, 2}}});
    add(2, "ACGN", {1}, {{{0, 2}}});
    add(3, "ACGN", {1, 2}, {{{1, 2}}});
    add(4, "ACGA", {1, 2}, {{{0, 2}}});
    add(5, "ACGN", {1, 2}, {{{0, 2}}}, false);
    add(6, "ACGN", {1, 0}, {{{0, 2}}});
    add(7, "", {1, 2}, {{{0, 2}}});
    add(8, "ACGN", {2, 1}, {{{0, 2}}});
    add(9, "ACGA", {2, 1}, {{{0, 2}}});

    // Force every sequence hash to collide, including unequal sequences.
    CClusterer clusters([&](int oid) { return records.at(oid); },
                        [](const std::string&) { return 0; });
    const int expected[] = {0, 0, 2, 3, 4, 5, 6, 7, 0, 4};
    for (int oid = 0; oid != 10; ++oid) {
        assert(clusters.Add(oid, records.at(oid)) == expected[oid]);
    }
    bool rejected = false;
    try {
        clusters.Add(9, records.at(9));
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    std::mt19937 rng(42);
    std::string sequence;
    for (int i = 0; i < 300; ++i) sequence += "ACGT"[rng() % 4];
    std::string substituted = sequence;
    substituted[100] = substituted[100] == 'A' ? 'C' : 'A';
    std::string inserted = sequence;
    inserted.insert(150, "G");
    std::string divergent = sequence;
    for (int i = 0; i < 4; ++i) {
        divergent[i * 50] = divergent[i * 50] == 'A' ? 'C' : 'A';
    }
    assert(NearEqual(sequence, substituted));
    assert(NearEqual(sequence, inserted));
    assert(!NearEqual(sequence, divergent));
    assert(NearEqual(sequence, divergent, ParseEditPercent("2")));
    assert(!NearEqual(sequence, substituted, ParseEditPercent("0")));
    std::string thousand(1000, 'A'), one_edit = thousand;
    one_edit[500] = 'C';
    assert(NearEqual(thousand, one_edit, ParseEditPercent("0.1")));
    assert(!NearEqual(thousand, one_edit, ParseEditPercent("0.099999")));
    assert(ParseEditPercent("100.000000") == 100000000);
    for (const char* invalid : {"", ".", "-1", "100.000001", "0.0000001", "nan"}) {
        bool bad = false;
        try {
            ParseEditPercent(invalid);
        } catch (const std::invalid_argument&) {
            bad = true;
        }
        assert(bad);
    }
    assert(!NearEqual(std::string(99, 'A'), std::string(98, 'A') + "C"));
    assert(!NearEqual(std::string(10001, 'A'), std::string(10001, 'A')));

    for (int trial = 0; trial < 100; ++trial) {
        std::string a = sequence.substr(0, 100 + rng() % 51), b = a;
        for (unsigned j = 0, edits = rng() % 5; j < edits; ++j) {
            std::size_t pos = rng() % b.size();
            switch (rng() % 3) {
            case 0: b[pos] = "ACGT"[rng() % 4]; break;
            case 1: b.insert(pos, 1, "ACGT"[rng() % 4]); break;
            case 2: b.erase(pos, 1); break;
            }
        }
        std::size_t longest = std::max(a.size(), b.size());
        std::size_t distance = EditDistance(a, b);
        for (const char* percent : {"1", "2", "5", "10", "100"}) {
            std::uint32_t threshold = ParseEditPercent(percent);
            std::size_t limit = longest * static_cast<std::uint64_t>(threshold) / 100000000;
            bool expected = std::min(a.size(), b.size()) >= 100 && limit > 0 &&
                            distance <= limit;
            assert(NearEqual(a, b, threshold) == expected);
        }
    }
    for (int trial = 0; trial < 30; ++trial) {
        std::string a, b;
        for (int i = 0; i < 100 + trial; ++i) {
            a += "ACGT"[rng() % 4];
            b += "ACGT"[rng() % 4];
        }
        std::size_t distance = EditDistance(a, b);
        for (const char* percent : {"20", "40", "60", "80", "100"}) {
            std::uint32_t threshold = ParseEditPercent(percent);
            std::size_t limit = a.size() * static_cast<std::uint64_t>(threshold) / 100000000;
            assert(NearEqual(a, b, threshold) == (distance <= limit));
        }
    }

    std::map<int, SRecord> near_records;
    near_records[0] = {sequence, {1}, {}, true};
    near_records[1] = {substituted, {1}, {}, true};
    near_records[2] = {inserted, {1}, {}, true};
    near_records[3] = {substituted, {2}, {}, true};
    near_records[4] = {substituted, {1}, {{{1, 2}}}, true};
    near_records[5] = {divergent, {1}, {}, true};
    near_records[6] = {substituted, {1}, {}, false};
    CClusterer near_clusters([&](int oid) { return near_records.at(oid); },
                             std::hash<std::string>(), true);
    const int near_expected[] = {0, 0, 0, 3, 4, 5, 6};
    for (int oid = 0; oid != 7; ++oid) {
        assert(near_clusters.Add(oid, near_records.at(oid)) == near_expected[oid]);
    }

    std::istringstream nodes(
        "1 | 1 | no rank |\n"
        "10 | 1 | species |\n"
        "11 | 10 | no rank |\n"
        "12 | 11 | no rank |\n"
        "13 | 12 | no rank |\n"
        "14 | 12 | strain |\n"
        "15 | 11 | strain |\n"
        "18 | 13 | no rank |\n");
    CTaxonomy taxonomy(nodes);
    assert(taxonomy.Group({10}) == 0);
    assert(taxonomy.Group({11}) == 0);
    assert(taxonomy.Group({12, 13, 14, 18}) == 12);
    assert(taxonomy.Group({13, 15}) == 0);
    assert(taxonomy.Group({10, 13}) == 0);
    assert(taxonomy.Group({999}) == 0 && taxonomy.MissingCount() == 1);
    assert(taxonomy.Group({13, 999}) == 0);
    assert(taxonomy.Lca({13, 18}) == 13);
    assert(taxonomy.Lca({13, 14, 18}) == 12);
    assert(taxonomy.Lca({13, 15}) == 11);

    std::map<int, SRecord> clade_records;
    const std::vector<std::vector<int>> taxa = {
        {13}, {18}, {14}, {15}, {11}, {13, 15}, {999}, {999}, {14}, {14}};
    for (std::size_t oid = 0; oid < taxa.size(); ++oid) {
        SRecord r;
        r.sequence = oid == 8 ? substituted : sequence;
        r.taxids = taxa[oid];
        r.complete_taxonomy = true;
        r.clade = taxonomy.Group(r.taxids);
        if (oid == 9) r.masks = {{{1, 2}}};
        clade_records.emplace(static_cast<int>(oid), r);
    }
    CClusterer clades([&](int oid) { return clade_records.at(oid); },
                      std::hash<std::string>(), true);
    const int clade_expected[] = {0, 0, 0, 3, 4, 5, 6, 6, 0, 9};
    for (int oid = 0; oid != 10; ++oid) {
        assert(clades.Add(oid, clade_records.at(oid)) == clade_expected[oid]);
    }

    for (const std::string& bad_nodes : {
             std::string("1 | 1 | no rank |\n2 | 3 | no rank |\n3 | 2 | no rank |\n"),
             std::string("1 | 1 | no rank |\n2 | 99 | no rank |\n")}) {
        std::istringstream input(bad_nodes);
        CTaxonomy bad_taxonomy(input);
        bool bad = false;
        try {
            bad_taxonomy.Group({2});
        } catch (const std::runtime_error&) {
            bad = true;
        }
        assert(bad);
    }
    bool malformed = false;
    try {
        std::istringstream input("1 | 1 | no rank |\n1 | 1 | no rank |\n");
        CTaxonomy bad_taxonomy(input);
    } catch (const std::runtime_error&) {
        malformed = true;
    }
    assert(malformed);
}

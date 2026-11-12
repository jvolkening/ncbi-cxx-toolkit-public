/* This is a United States Government Work, freely available for public use. */
#include "taxon_exact_cluster.hpp"

#include <cassert>
#include <map>
#include <random>
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
        assert(NearEqual(a, b) == (EditDistance(a, b) <= longest / 100));
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
}

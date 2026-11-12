#ifndef TAXON_EXACT_CLUSTER_HPP
#define TAXON_EXACT_CLUSTER_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace taxon_exact {

using TMasks = std::vector<std::vector<std::pair<unsigned int, unsigned int>>>;

struct SRecord {
    std::string sequence;
    std::vector<int> taxids;
    TMasks masks;
    bool complete_taxonomy = false;
    int clade = 0;
};

// Percent is stored in millionths so threshold rounding is exact.
inline std::uint32_t ParseEditPercent(const std::string& text)
{
    if (text.empty()) throw std::invalid_argument("Edit percent is empty");
    std::uint32_t whole = 0, fraction = 0, places = 0;
    bool decimal = false, digit = false;
    for (char c : text) {
        if (c == '.' && !decimal) {
            decimal = true;
        } else if (c >= '0' && c <= '9') {
            digit = true;
            if (decimal) {
                if (++places > 6) throw std::invalid_argument("Edit percent supports at most six decimal places");
                fraction = fraction * 10 + (c - '0');
            } else {
                if (whole > 100 || whole * 10 + (c - '0') > 100) {
                    throw std::invalid_argument("Edit percent must be between 0 and 100");
                }
                whole = whole * 10 + (c - '0');
            }
        } else {
            throw std::invalid_argument("Invalid edit percent");
        }
    }
    if (!digit) throw std::invalid_argument("Invalid edit percent");
    while (places++ < 6) fraction *= 10;
    if (whole == 100 && fraction != 0) {
        throw std::invalid_argument("Edit percent must be between 0 and 100");
    }
    return whole * 1000000 + fraction;
}

// Bounded wavefront edit distance. Exhausting the work budget rejects a
// candidate, so approximate search never accepts an unverified distance.
inline bool NearEqual(const std::string& a, const std::string& b,
                      std::uint32_t edit_percent = 1000000)
{
    const std::size_t longest = std::max(a.size(), b.size());
    if (edit_percent > 100000000) throw std::invalid_argument("Edit percent exceeds 100");
    if (longest < 100 || longest > 10000 || std::min(a.size(), b.size()) < 100) return false;
    const std::size_t limit =
        static_cast<std::uint64_t>(longest) * edit_percent / 100000000;
    if (limit == 0 ||
        longest - std::min(a.size(), b.size()) > limit) {
        return false;
    }
    const int n = static_cast<int>(a.size()), m = static_cast<int>(b.size());
    const int max_edits = static_cast<int>(limit), offset = max_edits + 1;
    std::vector<int> previous(2 * max_edits + 3, -1), current(previous.size(), -1);
    std::size_t work = 0;
    const std::size_t work_limit = 500000;
    auto extend = [&](int x, int diagonal) {
        int y = x - diagonal;
        while (x < n && y < m) {
            if (++work > work_limit) return -1;
            if (a[x] != b[y]) break;
            ++x;
            ++y;
        }
        return x;
    };
    previous[offset] = extend(0, 0);
    if (previous[offset] == n && n == m) return true;
    for (int edits = 1; edits <= max_edits; ++edits) {
        for (int diagonal = -edits; diagonal <= edits; ++diagonal) {
            if (++work > work_limit) return false;
            int best = -1;
            int x = previous[offset + diagonal - 1];
            if (x >= 0 && x < n) best = x + 1;
            x = previous[offset + diagonal + 1];
            if (x >= 0 && x - (diagonal + 1) < m) best = std::max(best, x);
            x = previous[offset + diagonal];
            if (x >= 0 && x < n && x - diagonal < m) best = std::max(best, x + 1);
            if (best >= 0) {
                best = extend(best, diagonal);
                if (best < 0) return false;
                if (best == n && best - diagonal == m) return true;
            }
            current[offset + diagonal] = best;
        }
        previous.swap(current);
    }
    return false;
}

// Bottom 16 distinct 15-mer hashes; non-ACGT bases break a k-mer.
inline std::vector<std::uint64_t> Anchors(const std::string& sequence)
{
    std::set<std::uint64_t> hashes;
    std::uint64_t kmer = 0;
    unsigned valid = 0;
    for (char c : sequence) {
        unsigned value;
        switch (c) {
        case 'A': value = 0; break;
        case 'C': value = 1; break;
        case 'G': value = 2; break;
        case 'T': value = 3; break;
        default: valid = 0; kmer = 0; continue;
        }
        kmer = ((kmer << 2) | value) & ((std::uint64_t(1) << 30) - 1);
        if (++valid < 15) continue;
        // SplitMix64 finalizer makes the chosen anchors independent of DNA order.
        std::uint64_t h = kmer + UINT64_C(0x9e3779b97f4a7c15);
        h = (h ^ (h >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        h = (h ^ (h >> 27)) * UINT64_C(0x94d049bb133111eb);
        h ^= h >> 31;
        hashes.insert(h);
        if (hashes.size() > 16) hashes.erase(--hashes.end());
        valid = 15;
    }
    return std::vector<std::uint64_t>(hashes.begin(), hashes.end());
}

// The reader must return the same record for an OID throughout the run.
class CClusterer {
public:
    using TRead = std::function<SRecord(int)>;
    using THash = std::function<std::size_t(const std::string&)>;

    explicit CClusterer(TRead read, THash hash = std::hash<std::string>(),
                        bool near = false, std::uint32_t edit_percent = 1000000)
        : m_Read(std::move(read)), m_Hash(std::move(hash)), m_Near(near),
          m_EditPercent(edit_percent) {}

    int Add(int oid, const SRecord& record)
    {
        if (oid < 0 || oid <= m_LastOid) {
            throw std::invalid_argument("OIDs must be nonnegative and increasing");
        }
        m_LastOid = oid;
        if (!record.complete_taxonomy || record.taxids.empty() ||
            record.sequence.empty()) {
            return oid;
        }
        auto taxids = record.taxids;
        std::sort(taxids.begin(), taxids.end());
        taxids.erase(std::unique(taxids.begin(), taxids.end()), taxids.end());
        if (taxids.front() <= 0) {
            return oid;
        }
        const std::vector<int> taxon_key = record.clade > 0
            ? std::vector<int>{-record.clade} : taxids;
        TKey key(taxon_key, record.sequence.size(), m_Hash(record.sequence));
        auto exact = m_Representatives.find(key);
        if (exact != m_Representatives.end()) {
            for (int candidate : exact->second) {
                SRecord representative = m_Read(candidate);
                if (record.sequence == representative.sequence &&
                    record.masks == representative.masks) {
                    return candidate;
                }
            }
        }
        if (m_Near && m_EditPercent != 0 && record.sequence.size() >= 100 &&
            record.sequence.size() <= 10000) {
            const auto anchors = Anchors(record.sequence);
            std::set<int> near_candidates;
            for (auto anchor : anchors) {
                auto it = m_Anchors.find(std::make_pair(taxon_key, anchor));
                if (it == m_Anchors.end() ||
                    (!it->second.empty() && it->second.front() == -1)) continue;
                for (int candidate : it->second) {
                    near_candidates.insert(candidate);
                    if (near_candidates.size() == 256) break;
                }
                if (near_candidates.size() == 256) break;
            }
            for (int candidate : near_candidates) {
                SRecord representative = m_Read(candidate);
                if (record.masks == representative.masks &&
                    NearEqual(record.sequence, representative.sequence, m_EditPercent)) {
                    return candidate;
                }
            }
            for (auto anchor : anchors) {
                auto& bucket = m_Anchors[std::make_pair(taxon_key, anchor)];
                if (!bucket.empty() && bucket.front() == -1) continue;
                if (bucket.size() == 64) {
                    bucket.assign(1, -1);
                } else {
                    bucket.push_back(oid);
                }
            }
        }
        m_Representatives[key].push_back(oid);
        return oid;
    }

private:
    using TKey = std::tuple<std::vector<int>, std::size_t, std::size_t>;
    using TAnchorKey = std::pair<std::vector<int>, std::uint64_t>;
    TRead m_Read;
    THash m_Hash;
    bool m_Near;
    std::uint32_t m_EditPercent;
    int m_LastOid = -1;
    std::map<TKey, std::vector<int>> m_Representatives;
    std::map<TAnchorKey, std::vector<int>> m_Anchors;
};

}
#endif

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
};

// Global edit distance, counting substitutions and single-base insertions and
// deletions. Ambiguity codes compare literally.
inline bool NearEqual(const std::string& a, const std::string& b)
{
    const std::size_t longest = std::max(a.size(), b.size());
    const std::size_t limit = longest / 100;
    if (longest > 10000 || limit == 0 ||
        longest - std::min(a.size(), b.size()) > limit) {
        return false;
    }
    const std::size_t inf = limit + 1;
    std::vector<std::size_t> prev(b.size() + 1, inf), cur(b.size() + 1, inf);
    for (std::size_t j = 0; j <= std::min(b.size(), limit); ++j) prev[j] = j;
    for (std::size_t i = 1; i <= a.size(); ++i) {
        const std::size_t lo = i > limit ? i - limit : 0;
        const std::size_t hi = std::min(b.size(), i + limit);
        if (lo > hi) return false;
        if (lo > 0) cur[lo - 1] = inf;
        if (lo == 0) cur[0] = i;
        std::size_t row_min = lo == 0 ? cur[0] : inf;
        for (std::size_t j = std::max<std::size_t>(lo, 1); j <= hi; ++j) {
            cur[j] = std::min({inf, prev[j] + 1, cur[j - 1] + 1,
                               prev[j - 1] + (a[i - 1] != b[j - 1])});
            row_min = std::min(row_min, cur[j]);
        }
        if (hi < b.size()) cur[hi + 1] = inf;
        if (row_min > limit) return false;
        prev.swap(cur);
    }
    return prev[b.size()] <= limit;
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
                        bool near = false)
        : m_Read(std::move(read)), m_Hash(std::move(hash)), m_Near(near) {}

    int Add(int oid, SRecord record)
    {
        if (oid < 0 || oid <= m_LastOid) {
            throw std::invalid_argument("OIDs must be nonnegative and increasing");
        }
        m_LastOid = oid;
        if (!record.complete_taxonomy || record.taxids.empty() ||
            record.sequence.empty()) {
            return oid;
        }
        std::sort(record.taxids.begin(), record.taxids.end());
        record.taxids.erase(std::unique(record.taxids.begin(), record.taxids.end()),
                            record.taxids.end());
        if (record.taxids.front() <= 0) {
            return oid;
        }
        TKey key(record.taxids, record.sequence.size(), m_Hash(record.sequence));
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
        if (m_Near && record.sequence.size() >= 100 &&
            record.sequence.size() <= 10000) {
            const auto anchors = Anchors(record.sequence);
            std::set<int> near_candidates;
            for (auto anchor : anchors) {
                auto it = m_Anchors.find(std::make_pair(record.taxids, anchor));
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
                    NearEqual(record.sequence, representative.sequence)) {
                    return candidate;
                }
            }
            for (auto anchor : anchors) {
                auto& bucket = m_Anchors[std::make_pair(record.taxids, anchor)];
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
    int m_LastOid = -1;
    std::map<TKey, std::vector<int>> m_Representatives;
    std::map<TAnchorKey, std::vector<int>> m_Anchors;
};

}
#endif

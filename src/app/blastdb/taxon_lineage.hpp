#ifndef TAXON_LINEAGE_HPP
#define TAXON_LINEAGE_HPP

#include <algorithm>
#include <istream>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace taxon_exact {

class CTaxonomy {
public:
    explicit CTaxonomy(std::istream& input)
    {
        std::string line;
        std::size_t line_number = 0;
        while (std::getline(input, line)) {
            ++line_number;
            const auto first = line.find('|');
            const auto second = first == std::string::npos ? first : line.find('|', first + 1);
            const auto third = second == std::string::npos ? second : line.find('|', second + 1);
            if (third == std::string::npos) {
                throw std::runtime_error("Malformed nodes.dmp line " + std::to_string(line_number));
            }
            int id = ParseId(line.substr(0, first), line_number);
            int parent = ParseId(line.substr(first + 1, second - first - 1), line_number);
            std::string rank = Trim(line.substr(second + 1, third - second - 1));
            if (rank.empty() || !m_Nodes.emplace(id, SNode{parent, rank == "species"}).second) {
                throw std::runtime_error("Invalid or duplicate node on line " +
                                         std::to_string(line_number));
            }
        }
        if (input.bad() || m_Nodes.empty()) {
            throw std::runtime_error("Cannot read a nonempty nodes.dmp");
        }
    }

    // Returns the ancestor exactly two edges below species if all taxids
    // fall in that subtree. Zero means no eligible common subtree.
    int Group(const std::vector<int>& taxids)
    {
        int group = 0;
        bool compatible = true;
        for (int id : taxids) {
            int current = GroupOne(id);
            if (current == 0 || (group != 0 && current != group)) compatible = false;
            if (group == 0) group = current;
        }
        return compatible ? group : 0;
    }

    int Lca(const std::vector<int>& taxids) const
    {
        if (taxids.empty()) throw std::runtime_error("Cannot find LCA of empty taxid set");
        auto common = Trace(taxids.front());
        std::reverse(common.begin(), common.end());
        for (std::size_t i = 1; i < taxids.size(); ++i) {
            auto path = Trace(taxids[i]);
            std::reverse(path.begin(), path.end());
            std::size_t n = 0;
            while (n < common.size() && n < path.size() && common[n] == path[n]) ++n;
            common.resize(n);
            if (common.empty()) throw std::runtime_error("Taxids have no common root");
        }
        return common.back();
    }

    std::size_t MissingCount() const { return m_Missing.size(); }

private:
    struct SNode {
        int parent;
        bool species;
    };

    static std::string Trim(const std::string& value)
    {
        auto begin = value.find_first_not_of(" \t\r\n");
        if (begin == std::string::npos) return {};
        auto end = value.find_last_not_of(" \t\r\n");
        return value.substr(begin, end - begin + 1);
    }

    static int ParseId(const std::string& value, std::size_t line)
    {
        const auto text = Trim(value);
        unsigned long long number = 0;
        if (text.empty()) throw std::runtime_error("Missing taxid on line " + std::to_string(line));
        for (char c : text) {
            if (c < '0' || c > '9' ||
                number > (std::numeric_limits<int>::max() - (c - '0')) / 10ULL) {
                throw std::runtime_error("Invalid taxid on line " + std::to_string(line));
            }
            number = number * 10 + (c - '0');
        }
        if (number == 0) throw std::runtime_error("Invalid taxid on line " + std::to_string(line));
        return static_cast<int>(number);
    }

    std::vector<int> Trace(int id) const
    {
        std::vector<int> path;
        while (true) {
            auto node = m_Nodes.find(id);
            if (node == m_Nodes.end()) {
                throw std::runtime_error("Missing parent taxid " + std::to_string(id));
            }
            if (std::find(path.begin(), path.end(), id) != path.end()) {
                throw std::runtime_error("Cycle in taxonomy at taxid " + std::to_string(id));
            }
            path.push_back(id);
            if (node->second.parent == id) {
                if (id != 1) throw std::runtime_error("Invalid taxonomy root " + std::to_string(id));
                return path;
            }
            id = node->second.parent;
        }
    }

    int GroupOne(int id)
    {
        auto cached = m_Groups.find(id);
        if (cached != m_Groups.end()) return cached->second;
        if (m_Nodes.find(id) == m_Nodes.end()) {
            m_Missing.insert(id);
            return 0;
        }
        auto path = Trace(id);
        int group = 0;
        for (std::size_t i = 0; i < path.size(); ++i) {
            if (m_Nodes.at(path[i]).species) {
                if (i >= 2) group = path[i - 2];
                break;
            }
        }
        m_Groups.emplace(id, group);
        return group;
    }

    std::unordered_map<int, SNode> m_Nodes;
    std::unordered_map<int, int> m_Groups;
    std::set<int> m_Missing;
};

}
#endif

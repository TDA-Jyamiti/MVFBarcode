#pragma once

#include <algorithm>
#include <charconv>
#include <fstream>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace cmvf {

/*
  SIMPLICES AND BLOCK PARTITIONS
  =============================

  Simplex ids remain fixed throughout an input sequence. Each simplex is a sorted list of vertex
  labels and every codimension-one face appears in the complex. Input partitions are multivector
  fields: each part must be locally closed. A part's canonical label is its smallest simplex id,
  so record labels, record order and vertex order inside a simplex leave the computation unchanged.

  The quotient face graph has an edge from a face's block to a coface's block. A topological order
  of this graph is an admissible block order. Its existence also certifies that every block is
  locally closed: any face chain leaving a block and returning to it would produce a directed
  cycle. A field's finest block partition is obtained by merging the strongly connected components
  of its quotient face graph. Reversing every arrow leaves these components unchanged.
*/

struct partitioned_complex {
  std::vector<std::vector<int>> simplices;
  std::vector<int> block_of;
};

inline std::string trim(std::string text) {
  const auto first = text.find_first_not_of(" \t\r\n"), last = text.find_last_not_of(" \t\r\n");
  return first == std::string::npos ? "" : text.substr(first, last - first + 1);
}

inline int parse_integer(const std::string& text) {
  const auto value = trim(text);
  int result = 0;
  const auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
  if (value.empty() || value.front() == '-' || parsed.ec != std::errc() || parsed.ptr != value.data() + value.size()) throw std::runtime_error("expected a nonnegative integer");
  return result;
}

inline std::vector<int> read_order(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(path + ": cannot open simplex order");
  std::vector<int> result;
  std::string line, value;
  while (std::getline(input, line)) {
    if (const auto comment = line.find('#'); comment != std::string::npos) line.resize(comment);
    std::istringstream entries(line);
    while (entries >> value) result.push_back(parse_integer(value));
  }
  if (input.bad() || result.empty()) throw std::runtime_error(path + ": expected a nonempty simplex order");
  return result;
}

inline std::vector<int> parse_list(const std::string& text) {
  std::vector<int> result;
  for (std::size_t begin = 0;;) {
    const auto end = text.find(',', begin);
    result.push_back(parse_integer(text.substr(begin, end == std::string::npos ? end : end - begin)));
    if (end == std::string::npos) break;
    begin = end + 1;
  }
  std::sort(result.begin(), result.end());
  if (std::adjacent_find(result.begin(), result.end()) != result.end()) throw std::runtime_error("duplicate list member");
  return result;
}

inline std::vector<int> canonical_partition(const std::vector<int>& partition) {
  std::map<int, int> labels;
  std::vector<int> result;
  result.reserve(partition.size());
  for (int i = 0; i < static_cast<int>(partition.size()); ++i) {
    if (partition[i] < 0) throw std::runtime_error("every simplex must belong to a block");
    result.push_back(labels.emplace(partition[i], i).first->second);
  }
  return result;
}

inline bool refines(const std::vector<int>& fine, const std::vector<int>& coarse) {
  if (fine.size() != coarse.size()) return false;
  std::vector<int> image(fine.size(), -1);
  for (int s = 0; s < static_cast<int>(fine.size()); ++s) {
    if (image[fine[s]] >= 0 && image[fine[s]] != coarse[s]) return false;
    image[fine[s]] = coarse[s];
  }
  return true;
}

inline partitioned_complex read_smp(const std::string& path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(path + ": cannot open input");
  partitioned_complex result;
  std::vector<std::vector<int>> blocks;
  std::set<int> labels;
  bool in_blocks = false;
  std::string line;
  std::size_t line_number = 0;
  try {
    while (std::getline(input, line)) {
      ++line_number;
      line = trim(line);
      if (line.empty()) continue;
      if (line.front() == '#') {
        auto header = trim(line.substr(1));
        if (!header.empty() && header.back() == '#') header = trim(header.substr(0, header.size() - 1));
        if (header == "multivectors") { if (in_blocks) throw std::runtime_error("repeated multivectors header"); in_blocks = true; }
        continue;
      }
      if (const auto comment = line.find('#'); comment != std::string::npos) line.resize(comment);
      const auto colon = line.find(':');
      if (colon == std::string::npos) throw std::runtime_error("expected id:members");
      const int label = parse_integer(line.substr(0, colon));
      auto members = parse_list(line.substr(colon + 1));
      if (in_blocks) { if (!labels.insert(label).second) throw std::runtime_error("duplicate multivector label"); blocks.push_back(std::move(members)); }
      else { if (label != static_cast<int>(result.simplices.size())) throw std::runtime_error("simplex ids must be dense and ordered from zero"); result.simplices.push_back(std::move(members)); }
    }
    if (input.bad()) throw std::runtime_error("input read failed");
    if (!in_blocks || result.simplices.empty()) throw std::runtime_error("a nonempty simplex section and a multivectors header are required");
    result.block_of.assign(result.simplices.size(), -1);
    for (const auto& block : blocks) for (int s : block) {
      if (s >= static_cast<int>(result.simplices.size())) throw std::runtime_error("block refers to an unknown simplex");
      if (result.block_of[s] != -1) throw std::runtime_error("blocks overlap");
      result.block_of[s] = block.front();
    }
    for (int s = 0; s < static_cast<int>(result.simplices.size()); ++s) if (result.block_of[s] == -1) result.block_of[s] = s;
  } catch (const std::exception& error) { throw std::runtime_error(path + ":" + std::to_string(line_number) + ": " + error.what()); }
  return result;
}

struct simplicial_complex {
  std::vector<std::vector<int>> simplices, faces;

  explicit simplicial_complex(const std::vector<std::vector<int>>& input) : simplices(input), faces(input.size()) {
    if (input.empty() || input.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) throw std::runtime_error("complex size is out of range");
    std::map<std::vector<int>, int> index;
    for (int s = 0; s < size(); ++s) {
      const auto& vertices = simplices[s];
      if (vertices.empty() || !std::is_sorted(vertices.begin(), vertices.end()) || std::adjacent_find(vertices.begin(), vertices.end()) != vertices.end()) throw std::runtime_error("simplex vertices must be distinct and sorted");
      if (!index.emplace(vertices, s).second) throw std::runtime_error("duplicate simplex");
    }
    for (int s = 0; s < size(); ++s) if (dimension(s)) for (int k = 0; k <= dimension(s); ++k) {
      auto face = simplices[s]; face.erase(face.begin() + k);
      const auto found = index.find(face);
      if (found == index.end()) throw std::runtime_error("missing face of simplex " + std::to_string(s));
      faces[s].push_back(found->second);
    }
  }

  int size() const { return static_cast<int>(simplices.size()); }
  int dimension(int s) const { return static_cast<int>(simplices[s].size()) - 1; }

  std::vector<int> finest_blocks(const std::vector<int>& field) const {
    if (field.size() != simplices.size()) throw std::runtime_error("field size differs from the complex");
    const auto owner = canonical_partition(field);
    std::vector<std::vector<int>> members(size()), edges(size()), reverse(size());
    for (int s = 0; s < size(); ++s) {
      members[owner[s]].push_back(s);
      for (int f : faces[s]) if (owner[s] != owner[f]) { edges[owner[s]].push_back(owner[f]); reverse[owner[f]].push_back(owner[s]); }
    }
    // Iterative Kosaraju: avoid a recursion-depth limit on large fields.
    std::vector<bool> seen(size());
    std::vector<int> finish, component(size(), -1), component_size(size());
    std::vector<std::pair<int, std::size_t>> search;
    for (int b = 0; b < size(); ++b) if (!members[b].empty() && !seen[b]) {
      seen[b] = true; search.emplace_back(b, 0);
      while (!search.empty()) {
        auto& [v, next] = search.back();
        if (next == edges[v].size()) { finish.push_back(v); search.pop_back(); }
        else { const int w = edges[v][next++]; if (!seen[w]) { seen[w] = true; search.emplace_back(w, 0); } }
      }
    }
    std::vector<int> stack;
    for (auto it = finish.rbegin(); it != finish.rend(); ++it) if (component[*it] < 0) {
      const int root = *it;
      component[root] = root; stack.push_back(root);
      while (!stack.empty()) {
        const int v = stack.back(); stack.pop_back(); ++component_size[root];
        for (int w : reverse[v]) if (component[w] < 0) { component[w] = root; stack.push_back(w); }
      }
    }
    // A nonconvex multivector has a face path that leaves it and returns. Such a path lies in
    // its SCC. Acyclic parts and singleton parts are already convex; check the remaining parts.
    std::vector<int> visited(size(), -1);
    for (int b = 0; b < size(); ++b) if (members[b].size() > 1 && component_size[component[b]] > 1) {
      auto visit = [&](int s) { if (component[owner[s]] == component[b] && visited[s] != b) { visited[s] = b; stack.push_back(s); } };
      for (int s : members[b]) for (int f : faces[s]) if (owner[f] != b) visit(f);
      while (!stack.empty()) {
        const int s = stack.back(); stack.pop_back();
        if (owner[s] == b) throw std::runtime_error("multivector containing simplex " + std::to_string(b) + " is not locally closed");
        for (int f : faces[s]) visit(f);
      }
    }
    std::vector<int> result(size());
    for (int s = 0; s < size(); ++s) result[s] = component[owner[s]];
    return canonical_partition(result);
  }

  std::vector<int> block_order(const std::vector<int>& block_of, const std::vector<int>& priority = {}) const {
    if (block_of.size() != simplices.size()) throw std::runtime_error("partition size differs from the complex");
    if (!priority.empty() && priority.size() != simplices.size()) throw std::runtime_error("ordering priority size differs from the complex");
    std::vector<std::vector<int>> edges(size());
    std::vector<int> degree(size()), present(size());
    for (int s = 0; s < size(); ++s) { if (block_of[s] < 0 || block_of[s] >= size()) throw std::runtime_error("block id out of range"); present[block_of[s]] = 1; }
    for (int s = 0; s < size(); ++s) for (int f : faces[s]) if (block_of[s] != block_of[f]) edges[block_of[f]].push_back(block_of[s]);
    for (auto& row : edges) { std::sort(row.begin(), row.end()); row.erase(std::unique(row.begin(), row.end()), row.end()); for (int b : row) ++degree[b]; }
    std::vector<int> first(size(), size());
    for (int s = 0; s < size(); ++s) first[block_of[s]] = std::min(first[block_of[s]], priority.empty() ? s : priority[s]);
    using entry = std::pair<int, int>;
    std::priority_queue<entry, std::vector<entry>, std::greater<entry>> ready;
    for (int b = 0; b < size(); ++b) if (present[b] && !degree[b]) ready.emplace(first[b], b);
    std::vector<int> order;
    while (!ready.empty()) { const int b = ready.top().second; ready.pop(); order.push_back(b); for (int next : edges[b]) if (!--degree[next]) ready.emplace(first[next], next); }
    if (order.size() != static_cast<std::size_t>(std::count(present.begin(), present.end(), 1))) throw std::runtime_error("partition is not an acyclic block partition");
    return order;
  }
};

}

#pragma once

#include <cstdint>
#include <tuple>
#include "bit_matrix.hpp"
#include "partitioned_complex.hpp"

namespace cmvf {

/*
  CONLEY-MORSE BARCODES
  ====================

  D is the boundary matrix of the entire simplicial complex over GF(2). The order places blocks
  in an admissible order, and faces before cofaces inside each block. R = D U throughout. A column
  of R is homogeneous when its pivot belongs to its own block. Homogeneous columns have distinct
  pivots; their pivot rows are targeted. Every remaining column carries one live barcode interval.
  U is upper triangular, has unit diagonal, and has off-diagonal entries only in homogeneous rows
  of the same block. These are the block reduction and uniformity conditions of Sections 4-6.

  block[j] is the first position of the block containing position j. A split or merge changes these
  labels while D remains fixed. A transposition changes both axes of all three matrices and repairs
  the decomposition. The feature slots implement the changes of representatives described in
  Section 6, including the fixed-slot convention for two live columns of the same degree.

  Input fields are replaced by their finest block partitions. Refinement-comparable fields give
  refinement-comparable block partitions. Larger changes are factored into binary splits or merges
  by peeling off children in reverse admissible order (or merging prefixes in admissible order).
  For B = W disjoint-union V, with W closed in B, the simplex-wise expansion in Section 7 puts B
  in order [W | V], splits off
  V from right to left, and merges its singletons from left to right. A merge runs this expansion
  backwards. Whole-block permutations preserve internal orders and implement the admissible
  reordering discussed in Section 8. Every birth and death within one supplied transition receives
  the destination input field index. Unchanged block partitions still occupy an input index.
  Intervals with equal endpoints disappear from the output.

  A supplied initial order places faces before cofaces and makes every initial coalesced block an interval.
  Splits keep the current child order whenever admissible. Adjacent merges preserve matrix order;
  other merges use a topological block order with ties resolved by the current simplex positions.
  Thus cuts of a common admissible order require no transpositions or block shuffles.
*/

struct interval {
  int dimension = 0;
  std::int64_t birth = 0, death = -1;
};

struct barcode_engine {
  simplicial_complex complex;
  binary_matrix D, R, U;
  std::vector<int> order, block, feature_of_column, field;
  std::vector<interval> features;
  std::int64_t partition_index = 0;

  explicit barcode_engine(const partitioned_complex& initial, const std::vector<int>& supplied_order = {}) : complex(initial.simplices), D(complex.size()), R(complex.size()), U(complex.size()) {
    field = canonical_partition(initial.block_of);
    const auto partition = complex.finest_blocks(field);
    const auto blocks = complex.block_order(partition);
    if (!supplied_order.empty()) {
      if (supplied_order.size() != initial.simplices.size()) throw std::runtime_error("simplex order must contain every simplex exactly once");
      order = supplied_order;
      std::vector<int> position(size(), -1);
      for (int j = 0; j < size(); ++j) { const int s = order[j]; if (s < 0 || s >= size() || position[s] >= 0) throw std::runtime_error("simplex order must be a permutation of the simplex ids"); position[s] = j; }
      for (int s = 0; s < size(); ++s) for (int f : complex.faces[s]) if (position[f] >= position[s]) throw std::runtime_error("simplex order must place faces before cofaces");
      check_initial_partition(partition);
      for (int j = 0; j < size(); ++j) block.push_back(j && partition[order[j]] == partition[order[j - 1]] ? block.back() : j);
    } else {
      std::vector<std::vector<int>> members(size());
      for (int s = 0; s < size(); ++s) members[partition[s]].push_back(s);
      for (int b : blocks) {
        auto& cells = members[b];
        std::sort(cells.begin(), cells.end(), [&](int a, int c) { return std::pair(complex.dimension(a), a) < std::pair(complex.dimension(c), c); });
        const int begin = static_cast<int>(order.size());
        for (int s : cells) { order.push_back(s); block.push_back(begin); }
      }
    }
    std::vector<int> position(size());
    for (int j = 0; j < size(); ++j) position[order[j]] = j;
    for (int j = 0; j < size(); ++j) { for (int f : complex.faces[order[j]]) D.set(position[f], j); U.set(j, j); }
    R = D;
    std::vector<int> targets(size(), -1);
    for (int j = 0; j < size(); ++j) {
      while (homogeneous(j) && targets[R.low(j)] >= 0) add(targets[R.low(j)], j);
      if (homogeneous(j)) targets[R.low(j)] = j;
    }
    feature_of_column.assign(size(), -1);
    for (int j = 0; j < size(); ++j) if (!homogeneous(j) && targets[j] < 0) birth(j);
  }

  int size() const { return complex.size(); }
  int block_end(int j) const { const int b = block[j]; while (j < size() && block[j] == b) ++j; return j; }
  bool homogeneous(int j) const { const int p = R.low(j); return p >= 0 && block[p] == block[j]; }
  int targeting(int row) const { for (int j = row + 1; j < block_end(row); ++j) if (homogeneous(j) && R.low(j) == row) return j; return -1; }
  void add(int source, int destination) { R.add(source, destination); U.add(source, destination); }

  void birth(int j) {
    if (feature_of_column[j] >= 0) throw std::logic_error("birth on an occupied feature slot");
    feature_of_column[j] = static_cast<int>(features.size());
    features.push_back({complex.dimension(order[j]), partition_index, -1});
  }

  void die(int j) {
    if (feature_of_column[j] < 0) throw std::logic_error("death on an empty feature slot");
    features[feature_of_column[j]].death = partition_index;
    feature_of_column[j] = -1;
  }

  void move_feature(int source, int destination) {
    if (source == destination) return;
    if (feature_of_column[destination] >= 0) throw std::logic_error("feature destination is occupied");
    feature_of_column[destination] = std::exchange(feature_of_column[source], -1);
  }

  void transpose(int j) {
    const int k = j + 1;
    if (j < 0 || k >= size() || block[j] != block[k] || D.get(j, k)) throw std::runtime_error("transposition requires adjacent incomparable cells in one block");
    const bool j_homogeneous = homogeneous(j), k_homogeneous = homogeneous(k), clear = U.get(j, k);
    const int a = targeting(j), b = targeting(k), low_j = R.low(j), low_k = R.low(k);
    const bool b_hits_j = b >= 0 && R.get(j, b), same_degree = complex.dimension(order[j]) == complex.dimension(order[k]);
    if (clear) add(j, k);
    D.swap_adjacent(j); R.swap_adjacent(j); U.swap_adjacent(j); std::swap(order[j], order[k]);
    if (clear && low_j > low_k) add(j, k);
    if (a >= 0 && b >= 0 && b_hits_j) add(std::min(a, b), std::max(a, b));
    /* In case (ii.1), row k becomes nonhomogeneous. Cancel its contributions to later columns.
       Column k of R has its boundary outside this block, so these additions preserve its internal pivots. */
    if (j_homogeneous && !k_homogeneous && clear) for (int column = k + 1; column < block_end(k); ++column) if (U.get(k, column)) add(k, column);
    if (!same_degree) std::swap(feature_of_column[j], feature_of_column[k]);
    else if (j_homogeneous && !k_homogeneous) { if (!clear) move_feature(k, j); }
    else if (!j_homogeneous && k_homogeneous) move_feature(j, k);
    else if (!j_homogeneous) {
      if (a >= 0 && b < 0) move_feature(k, j);
      else if (a < 0 && b >= 0 && !b_hits_j) move_feature(j, k);
    }
  }

  void split_right(int last) {
    if (last < 0 || last >= size() || block[last] == last || block_end(last) != last + 1) throw std::runtime_error("right split requires the last cell of a nonsingleton block");
    const int pivot = homogeneous(last) ? R.low(last) : -1;
    block[last] = last; R.columns[last] = D.columns[last]; U.clear(last); U.set(last, last);
    if (pivot >= 0) { birth(last); birth(pivot); }
  }

  void merge_right(int singleton) {
    if (singleton <= 0 || singleton >= size() || block[singleton] != singleton || block_end(singleton) != singleton + 1) throw std::runtime_error("right merge requires a singleton immediately after a block");
    block[singleton] = block[singleton - 1];
    while (homogeneous(singleton)) {
      int donor = -1;
      for (int j = block[singleton]; j < singleton; ++j) if (homogeneous(j) && R.low(j) == R.low(singleton)) { donor = j; break; }
      if (donor < 0) break;
      add(donor, singleton);
    }
    if (homogeneous(singleton)) { die(singleton); die(R.low(singleton)); }
  }

  std::vector<int> partition() const {
    std::vector<int> result(size());
    for (int j = 0; j < size(); ++j) result[order[j]] = block[j];
    return canonical_partition(result);
  }

private:
  void check_initial_partition(const std::vector<int>& partition) const {
    std::vector<bool> seen(size());
    int previous = -1;
    for (int s : order) if (partition[s] != previous) {
      const int b = partition[s];
      if (seen[b]) throw std::runtime_error("every initial block must be consecutive in the supplied simplex order");
      seen[b] = true; previous = b;
    }
  }

  void reorder_blocks(const std::vector<int>& cells) {
    if (cells == order) return;
    std::vector<int> position(size()), old_indices(size()), new_block(size()), new_features(size());
    for (int j = 0; j < size(); ++j) position[order[j]] = j;
    for (int j = 0; j < size(); ++j) { old_indices[j] = position[cells[j]]; new_features[j] = feature_of_column[old_indices[j]]; new_block[j] = j && block[old_indices[j]] == block[old_indices[j - 1]] ? new_block[j - 1] : j; }
    D.permute(old_indices); R.permute(old_indices); U.permute(old_indices);
    order = cells; block = std::move(new_block); feature_of_column = std::move(new_features);
  }

  void advance_binary(const std::vector<int>& destination) {
    const auto current = partition();
    const auto destination_order = complex.block_order(destination);
    int current_count = 0;
    for (int j = 0; j < size(); ++j) current_count += block[j] == j;
    const int difference = static_cast<int>(destination_order.size()) - current_count;
    if (difference != 1 && difference != -1) throw std::runtime_error("transition must be exactly one binary block split or merge");
    const bool splitting = difference == 1;
    const auto& fine = splitting ? destination : current;
    const auto& coarse = splitting ? current : destination;
    std::vector<int> image(size(), -1);
    std::vector<std::vector<int>> children(size());
    for (int s = 0; s < size(); ++s) { if (image[fine[s]] >= 0 && image[fine[s]] != coarse[s]) throw std::runtime_error("successive partitions must be related by refinement"); image[fine[s]] = coarse[s]; }
    for (int b = 0; b < size(); ++b) if (image[b] >= 0) children[image[b]].push_back(b);
    int parent = -1;
    for (int b = 0; b < size(); ++b) if (children[b].size() > 1) { if (parent >= 0 || children[b].size() != 2) throw std::runtime_error("transition must change exactly one block"); parent = b; }
    if (parent < 0) throw std::runtime_error("transition has no binary refinement");
    auto closed_in_parent = [&](int child) { for (int s = 0; s < size(); ++s) if (fine[s] == child) for (int f : complex.faces[s]) if (coarse[f] == parent && fine[f] != child) return false; return true; };
    int closed = children[parent][0], open = children[parent][1];
    for (int s : order) if (coarse[s] == parent) { if (fine[s] != closed) std::swap(closed, open); break; }
    if (!closed_in_parent(closed)) std::swap(closed, open);
    if (!closed_in_parent(closed)) throw std::runtime_error("split has no child closed in its parent");

    if (splitting) {
      int begin = 0; while (coarse[order[begin]] != parent) ++begin;
      const int end = block_end(begin);
      for (int j = begin + 1; j < end; ++j) if (fine[order[j]] == closed) for (int k = j; k > begin && fine[order[k - 1]] == open; --k) transpose(k - 1);
      int cut = begin; while (cut < end && fine[order[cut]] == closed) ++cut;
      for (int j = end - 1; j >= cut; --j) split_right(j);
      for (int j = cut + 1; j < end; ++j) merge_right(j);
    } else {
      int closed_begin = 0, open_begin = 0;
      while (fine[order[closed_begin]] != closed) ++closed_begin;
      while (fine[order[open_begin]] != open) ++open_begin;
      if (block_end(closed_begin) != open_begin) {
        std::vector<std::vector<int>> members(size());
        std::vector<int> cells, closed_cells, open_cells, priority(size());
        for (int j = 0; j < size(); ++j) priority[order[j]] = j;
        for (int s : order) { if (fine[s] == closed) closed_cells.push_back(s); else if (fine[s] == open) open_cells.push_back(s); else members[coarse[s]].push_back(s); }
        for (int b : complex.block_order(destination, priority)) {
          if (b == parent) { cells.insert(cells.end(), closed_cells.begin(), closed_cells.end()); cells.insert(cells.end(), open_cells.begin(), open_cells.end()); }
          else cells.insert(cells.end(), members[b].begin(), members[b].end());
        }
        reorder_blocks(cells);
      }
      int cut = 0; while (fine[order[cut]] != open) ++cut;
      const int end = block_end(cut);
      for (int j = end - 1; j > cut; --j) split_right(j);
      for (int j = cut; j < end; ++j) merge_right(j);
    }
    if (partition() != destination) throw std::logic_error("elementary expansion ended at the wrong partition");
  }

public:
  void advance(const partitioned_complex& next) {
    if (next.simplices != complex.simplices) throw std::runtime_error("fields must describe the same simplex-id-preserving complex");
    const auto next_field = canonical_partition(next.block_of);
    const bool splitting = refines(next_field, field);
    if (!splitting && !refines(field, next_field)) throw std::runtime_error("successive multivector fields must be related by refinement");
    const auto current = partition(), destination = complex.finest_blocks(next_field);
    const auto& fine = splitting ? destination : current;
    const auto& coarse = splitting ? current : destination;
    if (!refines(fine, coarse)) throw std::logic_error("coalescing failed to preserve refinement");
    ++partition_index;
    if (current != destination) {
      std::vector<int> priority(size());
      for (int j = 0; j < size(); ++j) priority[order[j]] = j;
      const auto fine_order = complex.block_order(fine, priority), coarse_order = complex.block_order(coarse, priority);
      if (fine_order.size() == coarse_order.size() + 1) advance_binary(destination);
      else {
        std::vector<std::vector<int>> children(size()), members(size());
        for (int s = 0; s < size(); ++s) members[fine[s]].push_back(s);
        for (int b : fine_order) children[coarse[b]].push_back(b);
        auto step = current;
        for (int parent : coarse_order) {
          const auto& parts = children[parent];
          // A prefix is closed within its coarse parent. Contracting it therefore leaves an
          // acyclic partition. Split off the last child, or reverse these prefix contractions.
          for (int k = 1; k < static_cast<int>(parts.size()); ++k) {
            const int cut = splitting ? static_cast<int>(parts.size()) - k : k;
            const int prefix = *std::min_element(parts.begin(), parts.begin() + cut + (splitting ? 0 : 1));
            for (int i = 0; i <= cut; ++i) for (int s : members[parts[i]]) step[s] = splitting && i == cut ? parts[i] : prefix;
            advance_binary(step);
          }
        }
      }
    }
    if (partition() != destination) throw std::logic_error("field expansion ended at the wrong partition");
    field = next_field;
  }

  std::vector<interval> barcodes() const {
    std::vector<interval> result;
    for (const auto& bar : features) if (bar.death < 0 || bar.birth < bar.death) result.push_back(bar);
    std::sort(result.begin(), result.end(), [](const interval& a, const interval& b) { return std::tuple(a.dimension, a.birth, a.death < 0 ? INT64_MAX : a.death) < std::tuple(b.dimension, b.birth, b.death < 0 ? INT64_MAX : b.death); });
    return result;
  }
};

}

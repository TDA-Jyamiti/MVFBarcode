#pragma once

#include <algorithm>
#include <bit>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cmvf {

/*
  BINARY MATRICES
  ===============

  A matrix has one fixed row and column for every simplex of the underlying complex. Columns are
  ordinary vectors of 64-bit words; row i occupies bit i % 64 of word i / 64. Column addition is XOR.
  The largest occupied row is the pivot, with -1 assigned to an empty column. This storage gives
  the global matrices D, R and U the same small interface throughout the computation.

  A permutation lists old indices in their new order. Both axes receive that permutation, so it
  acts on the domain and codomain of a boundary operator together. An adjacent transposition uses
  the same convention and changes just two bits in each column after exchanging the columns.
*/

struct binary_matrix {
  using word = std::uint64_t;
  int size = 0;
  std::vector<std::vector<word>> columns;

  explicit binary_matrix(int n = 0) : size(n) {
    if (n < 0) throw std::invalid_argument("matrix size must be nonnegative");
    columns.assign(n, std::vector<word>((static_cast<std::size_t>(n) + 63) / 64));
  }

  bool get(int i, int j) const { return (columns[j][i / 64] >> (i % 64)) & 1; }
  void set(int i, int j, bool value = true) { const word mask = word(1) << (i % 64); if (value) columns[j][i / 64] |= mask; else columns[j][i / 64] &= ~mask; }
  void clear(int j) { std::fill(columns[j].begin(), columns[j].end(), 0); }
  void add(int source, int destination) { for (std::size_t k = 0; k < columns[source].size(); ++k) columns[destination][k] ^= columns[source][k]; }
  int low(int j) const { for (int k = static_cast<int>(columns[j].size()) - 1; k >= 0; --k) if (columns[j][k]) return 64 * k + 63 - std::countl_zero(columns[j][k]); return -1; }

  template<class function> void each_row(int j, function visit) const {
    for (std::size_t k = 0; k < columns[j].size(); ++k) for (word bits = columns[j][k]; bits; bits &= bits - 1) visit(static_cast<int>(64 * k + std::countr_zero(bits)));
  }

  void swap_adjacent(int j) {
    std::swap(columns[j], columns[j + 1]);
    for (int k = 0; k < size; ++k) { const bool a = get(j, k), b = get(j + 1, k); if (a != b) { set(j, k, b); set(j + 1, k, a); } }
  }

  void permute(const std::vector<int>& old_indices) {
    std::vector<int> position(size);
    for (int i = 0; i < size; ++i) position[old_indices[i]] = i;
    binary_matrix result(size);
    for (int j = 0; j < size; ++j) each_row(old_indices[j], [&](int i) { result.set(position[i], j); });
    *this = std::move(result);
  }
};

}

#include "sw.h"
#include <iostream>
#include <xsimd/xsimd.hpp>

SmithWaterman naive_sw_traceback(const std::string& s1, const std::string& s2, int end_i, int end_j,
    int match, int mismatch, int gap_open, int gap_extend) {
    int len1 = s1.size();
    int len2 = s2.size();
    std::vector<std::vector<int>> H(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> E(len1 + 1, std::vector<int>(len2 + 1, 0));
    std::vector<std::vector<int>> F(len1 + 1, std::vector<int>(len2 + 1, 0));

    // Fill in H, E, F
    for (int i = 1; i <= len1; ++i) {
        for (int j = 1; j <= len2; ++j) {
            E[i][j] = std::max(H[i-1][j] - gap_open, E[i-1][j] - gap_extend);
            F[i][j] = std::max(H[i][j-1] - gap_open, F[i][j-1] - gap_extend);
            int score_diag = H[i-1][j-1] + ((s1[i-1] == s2[j-1]) ? match : -mismatch);
            H[i][j] = std::max(0, std::max({score_diag, E[i][j], F[i][j]}));
        }
    }

    // Traceback from (end_i, end_j)
    int i = end_i;
    int j = end_j;
    std::string aligned1, aligned2, match_line;

    int score = H[i][j];

    while (i > 0 && j > 0 && H[i][j] > 0) {
        if (H[i][j] == H[i-1][j-1] + ((s1[i-1] == s2[j-1]) ? match : -mismatch)) {
            aligned1 += s1[i-1];
            aligned2 += s2[j-1];
            match_line += (s1[i-1] == s2[j-1]) ? '|' : '*';
            --i; --j;
        } else if (H[i][j] == E[i][j]) {
            aligned1 += s1[i-1];
            aligned2 += '-';
            match_line += ' ';
            --i;
        } else {
            aligned1 += '-';
            aligned2 += s2[j-1];
            match_line += ' ';
            --j;
        }
    }

    // Reverse strings
    std::reverse(aligned1.begin(), aligned1.end());
    std::reverse(aligned2.begin(), aligned2.end());
    std::reverse(match_line.begin(), match_line.end());

    SmithWaterman result;
    result.score = score;
    result.aligned_seq1 = aligned1;
    result.aligned_seq2 = aligned2;
    result.match_line = match_line;
    result.start1 = i;
    result.end1 = end_i;
    result.start2 = j;
    result.end2 = end_j;

    return result;      
}

inline int16_t base2num(char b) {
  switch (b) {
  case 'A':
    return 0;
  case 'T':
    return 1;
  case 'C':
    return 2;
  case 'G':
    return 3;
  default:
    return 255;
  }
}

inline char num2base(int16_t n) {
  switch (n) {
  case 0:
    return 'A';
  case 1:
    return 'T';
  case 2:
    return 'C';
  case 3:
    return 'G';
  default:
    return 'N';
  }
}

std::vector<int16_t> seq2vec(const std::string &s) {
  std::vector<int16_t> v;
  for (char c : s)
    v.push_back(base2num(c));
  return v;
}

std::vector<xsimd::batch<int16_t>>
reorder_for_simd(const std::vector<int16_t> &input) {
  using batch = xsimd::batch<int16_t>;
  constexpr int B = batch::size;

  int total = input.size();
  int seq_num = (total + B - 1) / B;

  std::vector<batch> output(seq_num);

  for (int i = 0; i < seq_num; ++i) {
    alignas(alignof(batch)) int16_t temp[B] = {0}; // 一個 batch buffer

    for (int j = 0; j < B; ++j) {
      int idx = j * seq_num + i;
      if (idx < total)
        temp[j] = static_cast<int16_t>(input[idx]);
      else
        temp[j] = -1; // padding
    }

    output[i] = batch::load_aligned(temp);
  }

  return output;
}

xsimd::batch<int16_t> shift_left(const xsimd::batch<int16_t> &input) {
  using batch = xsimd::batch<int16_t>;
  constexpr int B = batch::size;

  alignas(alignof(batch)) int16_t temp[B];
  input.store_aligned(temp);

  for (int i = B - 1; i >= 1; --i)
    temp[i] = temp[i - 1];
  temp[0] = 0;

  return batch::load_aligned(temp);
}

Position striped_smith_waterman(const std::vector<int16_t> &target, const std::vector<int16_t> &query,
        int match_score, int mismatch_penalty, int gap_open_penalty, int gap_extend_penalty) {
  using simd_type = xsimd::batch<int16_t>;
  constexpr int VEC_WIDTH = simd_type::size;

  const int query_len = query.size();
  const int target_len = target.size();
  const int blocks = (query_len + VEC_WIDTH - 1) / VEC_WIDTH;

  std::vector<simd_type> reordered_query = reorder_for_simd(query);
  std::vector<simd_type> curr_H(blocks, simd_type(0));
  std::vector<simd_type> prev_H(blocks, simd_type(0));
  std::vector<simd_type> E_col(blocks, simd_type(0));
  std::vector<simd_type> F_row(blocks, simd_type(0));
  std::vector<simd_type> best_column(blocks, simd_type(0));

  simd_type simd_match = simd_type::broadcast(match_score);
  simd_type simd_mismatch = simd_type::broadcast(mismatch_penalty);
  simd_type simd_gap_open = simd_type::broadcast(gap_open_penalty);
  simd_type simd_gap_ext = simd_type::broadcast(gap_extend_penalty);

  simd_type h_diag = simd_type(0);
  simd_type max_tracker = simd_type(0);
  simd_type col_tracker = simd_type(0);

  int best_i = -1, best_j = -1;

  for (int i = 0; i < target_len; ++i) {
    col_tracker = simd_type(0);
    h_diag = shift_left(curr_H.back());
    std::swap(curr_H, prev_H);
    simd_type ref_nt = simd_type::broadcast(target[i]);
    simd_type f = simd_type(0);

    for (int blk = 0; blk < blocks; ++blk) {
      auto match_mask = reordered_query[blk] == ref_nt;
      simd_type score = xsimd::select(match_mask, simd_match, simd_mismatch);

      h_diag = xsimd::max(h_diag + score, simd_type(0));
      h_diag = xsimd::max(h_diag, (xsimd::max(E_col[blk], f)));
      col_tracker = xsimd::max(col_tracker, h_diag);

      curr_H[blk] = h_diag;
      E_col[blk] = xsimd::max(E_col[blk] + simd_gap_ext, h_diag + simd_gap_open);
      f = xsimd::max(f + simd_gap_ext, h_diag + simd_gap_open);

      h_diag = prev_H[blk];
    }

    int idx = 0;
    while (xsimd::any(f > curr_H[idx] + simd_gap_open)) {
      curr_H[idx] = xsimd::max(curr_H[idx], f);
      col_tracker = xsimd::max(col_tracker, curr_H[idx]);
      f += simd_gap_ext;
      if (++idx >= blocks) {
        idx = 0;
        f = shift_left(f);
      }
    }

    if (xsimd::any(col_tracker > max_tracker)) {
      max_tracker = col_tracker;
      best_column = curr_H;
      best_i = i;
    }
  }

  int best_score = 0;
  for (int blk = 0; blk < blocks; ++blk) {
    alignas(alignof(simd_type)) int16_t temp[VEC_WIDTH];
    best_column[blk].store_aligned(temp);

    for (int k = 0; k < VEC_WIDTH; ++k) {
      int q_pos = k * blocks + blk;
      if (q_pos >= query_len) continue;

      if (temp[k] > best_score) {
        best_score = temp[k];
        best_j = q_pos;
      }
    }
  }

  return {best_i + 1, best_j + 1};
}
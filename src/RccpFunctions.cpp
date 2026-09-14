#include <Rcpp.h>
using namespace Rcpp;
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>

// ------------------------------------------------------------------------------
// 1. allele_diff_paired
//' Count or locate SNPs between paired germline and input sequences
//'
//' Compares each germline sequence in \code{germs} against the input sequence at the
//' same position in \code{inputs}, and returns either the positions of the mismatches
//' or how many there are. The two vectors are paired element by element, which is what
//' distinguishes this from \code{\link{allele_diff_indices}}, where every sequence is
//' compared against the first.
//'
//' A position is ignored when \emph{either} sequence carries a non-mismatch character,
//' by default a gap or an ambiguous base (\code{N}, \code{.}, \code{-}). Sequences of
//' unequal length are padded to the longer with \code{N}, so the padded positions never
//' count as mismatches.
//'
//' @param germs A vector of strings representing germline sequences.
//' @param inputs A vector of strings representing input sequences, paired with
//'   \code{germs} element by element.
//' @param X The position from which mismatches are counted, zero-based (default 0).
//' @param return_count Return the number of mismatches per pair rather than their
//'   positions (default \code{FALSE}).
//' @param non_mismatch_chars_nullable Characters ignored on either side when comparing
//'   (default \code{N}, \code{.}, \code{-}).
//' @return A list of integer vectors of one-based positions when
//'   \code{return_count = FALSE}, or an integer vector of counts when \code{TRUE}.
//'
//' @examples
//' germs  <- c("ACGTACGT", "ACGTACGT")
//' inputs <- c("ACGTTCGT", "ACGTTCGA")
//'
//' allele_diff_paired(germs, inputs, return_count = TRUE)   # 1, 2
//' allele_diff_paired(germs, inputs)                        # list(5), c(5, 8)
//'
//' @name allele_diff_paired
//' @export
// [[Rcpp::export]]
Rcpp::RObject allele_diff_paired(std::vector<std::string> germs,
                                 std::vector<std::string> inputs,
                                 int X = 0,
                                 bool return_count = false,
                                 Rcpp::Nullable<Rcpp::CharacterVector> non_mismatch_chars_nullable = R_NilValue) {
  std::unordered_set<char> non_mismatch_chars = {'N', '.', '-'};
  if (non_mismatch_chars_nullable.isNotNull()) {
    Rcpp::CharacterVector char_vec(non_mismatch_chars_nullable);
    non_mismatch_chars.clear();
    for (const auto& c : char_vec) {
      non_mismatch_chars.insert(Rcpp::as<std::string>(c)[0]);
    }
  }

  if (germs.size() != inputs.size()) {
    Rcpp::stop("The size of germs and inputs must be the same.");
  }

  size_t num_sequences = germs.size();
  auto pad_with_ns = [](std::string& seq, size_t target_length) {
    if (seq.size() < target_length) {
      seq.append(target_length - seq.size(), 'N');
    }
  };

  // One loop, shared by both return shapes. Padding to the longer of the pair is what
  // keeps this in bounds when the two sequences differ in length; because the pad
  // character is itself a non-mismatch character, padded positions are never counted.
  Rcpp::List snp_list(return_count ? 0 : num_sequences);
  std::vector<int> mutation_counts(return_count ? num_sequences : 0);

  for (size_t i = 0; i < num_sequences; ++i) {
    std::string germ = germs[i];
    std::string input = inputs[i];
    size_t max_length = std::max(germ.size(), input.size());
    pad_with_ns(germ, max_length);
    pad_with_ns(input, max_length);

    std::vector<int> snp_indices;
    int count = 0;
    for (size_t j = 0; j < max_length; ++j) {
      if (j >= static_cast<size_t>(X) && germ[j] != input[j] &&
          non_mismatch_chars.find(input[j]) == non_mismatch_chars.end() &&
          non_mismatch_chars.find(germ[j]) == non_mismatch_chars.end()) {
        if (return_count) {
          count++;
        } else {
          snp_indices.push_back(j + 1);
        }
      }
    }
    if (return_count) {
      mutation_counts[i] = count;
    } else {
      snp_list[i] = Rcpp::wrap(snp_indices);
    }
  }

  if (return_count) return Rcpp::wrap(mutation_counts);
  return snp_list;
}

// ------------------------------------------------------------------------------
// 2. insert_gaps
//' Insert gaps into an ungapped sequence based on a gapped reference sequence.
//'
//' This function inserts gaps (e.g., \code{.} or \code{-}) into an ungapped sequence
//' (\code{ungapped}) to match the positions of gaps in a reference sequence
//' (\code{gapped}). It ensures that the aligned sequence has the same gap structure as
//' the reference. Vectorised over both arguments, which are paired element by element.
//'
//' @param gapped A vector of strings representing the reference sequences with gaps.
//' @param ungapped A vector of strings representing the sequences without gaps.
//' @return A vector of strings with gaps inserted to match the gapped reference.
//'
//' @examples
//' gapped <- c("caggtc..aact", "caggtc---aact")
//' ungapped <- c("caggtcaact", "caggtcaact")
//'
//' insert_gaps(gapped, ungapped)   # "caggtc..aact", "caggtc---aact"
//'
//' @name insert_gaps
//' @export
// [[Rcpp::export]]
std::vector<std::string> insert_gaps(const std::vector<std::string>& gapped,
                                     const std::vector<std::string>& ungapped) {
  if (gapped.size() != ungapped.size()) {
    Rcpp::stop("The size of gapped and ungapped vectors must be the same.");
  }
  const std::unordered_set<char> gap_chars = {'.', '-'};
  size_t num_sequences = gapped.size();
  std::vector<std::string> results(num_sequences);

  auto process_sequence = [&gap_chars](const std::string& gapped, const std::string& ungapped) -> std::string {
    std::string result;
    size_t ungapped_index = 0;
    for (char gap_char : gapped) {
      if (gap_chars.find(gap_char) != gap_chars.end()) {
        result.push_back(gap_char);
      } else {
        if (ungapped_index < ungapped.size()) {
          result.push_back(ungapped[ungapped_index]);
          ++ungapped_index;
        } else {
          break;
        }
      }
    }
    return result;
  };

  for (size_t i = 0; i < num_sequences; ++i) {
    results[i] = process_sequence(gapped[i], ungapped[i]);
  }

  return results;
}

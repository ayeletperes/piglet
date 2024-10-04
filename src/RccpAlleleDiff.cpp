// ------------------------------------------------------------------------------
// Allele difference
// The Rcpp code for the `allele_diff` function
// the code for the example is written in R.  
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>


//'Calculate differences between characters in columns of germs and return them as a string vector.
//'
//'@param germs A vector of strings representing germ sequences.
//'@param X The threshold index from which to return differences as strings.
//'@return A vector of strings containing differences between characters in columns.
//'
//' @examples
//' {
//' germs = c("ATCG", "ATCC") 
//' X = 3 
//' result = allele_diff_strings(germs, X) 
//' # "A2T", "T3C", "C2G"
//' }
//' 
//' @export
// [[Rcpp::export]]
std::vector<std::string> allele_diff_strings(std::vector<std::string> germs, int X=0) {
  std::vector<std::vector<char>> germs_m;
  for (const std::string& germ : germs) {
    germs_m.push_back(std::vector<char>(germ.begin(), germ.end()));
  }
  
  int max_length = 0;
  for (const auto& germ : germs_m) {
    max_length = std::max(max_length, static_cast<int>(germ.size()));
  }
  
  for (auto& germ : germs_m) {
    germ.resize(max_length, '.'); 
  }
  
  auto setdiff_mat = [](const std::vector<char>& x) -> int {
    std::unordered_set<char> unique_chars(x.begin(), x.end());
    std::unordered_set<char> filter_chars = { '.', 'N', '-', 'X' };
    int diff_count = 0;
    for (const char& c : unique_chars) {
      if (filter_chars.find(c) == filter_chars.end()) {
        diff_count++;
      }
    }
    return diff_count;
  };
  
  std::vector<std::string> idx_strings; // Return as a string vector
  for (int i = 0; i < max_length; i++) {
    std::vector<char> column_chars;
    for (const auto& germ : germs_m) {
      column_chars.push_back(germ[i]);
    }
    int diff_count = setdiff_mat(column_chars);
    if (diff_count > 1 && i >= (X-1)) {
      std::string concatenated_str = std::string(1, column_chars[0]) + std::to_string(i+1) + std::string(1, column_chars[1]);
      idx_strings.push_back(concatenated_str);
    }
  }
  
  return idx_strings;
}

//'Calculate differences between characters in columns of germs and return their indices as an int vector.
//'
//'@param germs A vector of strings representing germ sequences.
//'@param X The threshold index from which to return differences as indices.
//'@return A vector of integers containing indices of differing columns.
//'
//' @examples 
//' 
//' {
//' germs = c("ATCG", "ATCC") 
//' X = 3 
//' result = allele_diff_indices(germs, X)
//' # 1, 2, 3
//' }
//' 
//' @export
// [[Rcpp::export]]
std::vector<int> allele_diff_indices(std::vector<std::string> germs, int X=0) {
  std::vector<std::vector<char>> germs_m;
  for (const std::string& germ : germs) {
    germs_m.push_back(std::vector<char>(germ.begin(), germ.end()));
  }
  
  int max_length = 0;
  for (const auto& germ : germs_m) {
    max_length = std::max(max_length, static_cast<int>(germ.size()));
  }
  
  for (auto& germ : germs_m) {
    germ.resize(max_length, '.'); 
  }
  
  auto setdiff_mat = [](const std::vector<char>& x) -> int {
    std::unordered_set<char> unique_chars(x.begin(), x.end());
    std::unordered_set<char> filter_chars = { '.', 'N', '-', 'X' };
    int diff_count = 0;
    for (const char& c : unique_chars) {
      if (filter_chars.find(c) == filter_chars.end()) {
        diff_count++;
      }
    }
    return diff_count;
  };
  
  std::vector<int> idx_indices; // Return as an int vector with indices
  for (int i = 0; i < max_length; i++) {
    std::vector<char> column_chars;
    for (const auto& germ : germs_m) {
      column_chars.push_back(germ[i]);
    }
    int diff_count = setdiff_mat(column_chars);
    if (diff_count > 1 && i >= (X-1)) {
      idx_indices.push_back(i+1);
    }
  }
  
  return idx_indices;
}

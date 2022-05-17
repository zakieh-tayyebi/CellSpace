/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "data.h"
#include <string>
#include <vector>
#include <fstream>
#include <numeric>
#include <assert.h>

using namespace std;

namespace starspace {

InternDataHandler::InternDataHandler(shared_ptr<Args> args, std::shared_ptr<Dictionary> dict, std::shared_ptr<DataParser> parser) {
  size_ = 0;
  idx_ = -1;
  examples_.clear();
  dict_ = dict;
  parser_ = parser;
  args_= args;
}

void InternDataHandler::errorOnZeroExample(const string& fileName) {
  std::cerr << "ERROR: File '" << fileName
            << "' does not contain any valid example.\n"
            << "Please check: is the file empty? "
            << "Do the examples contain proper feature and label according to the trainMode? "
            << "If your examples are unlabeled, try to set trainMode=5.\n";
  exit(EXIT_FAILURE);
}

inline bool check_header(string header){
    if(header == "%%MatrixMarket matrix coordinate pattern general") return(true);
    else if(header == "%%MatrixMarket matrix coordinate integer general") return(false);
    else {
        cerr << "Unsupported matrix format!" << endl;
        exit(EXIT_FAILURE);
    }
}

void InternDataHandler::loadFromFile() {

  examples_ = vector<ParseResults>(args_->nPeaks_total);
  ngrams_ = vector<vector<Base>>(args_->nPeaks_total);
  // hashes_ = vector<vector<int32_t>>(args_->nPeaks_total);

  unsigned pfn = 0, pnt = 0, nPeaks_cur = 0;
  for(auto peaks: args_->peaks_list){
    auto nPeaks = args_->nPeaks_list[pfn++];
    cout << "Reading " << nPeaks << " peak sequences from \'" << peaks << "\'" << endl;
    ifstream pf(peaks);
    string line, seq;
    int pi = -1;
    while(getline(pf, line))
      if(line[0] == '>'){
        if(pi >= 0) add_kmer_tokens(nPeaks_cur + pi, seq);
        pi++; pnt++;
        seq = "";
      } else seq += line;
    if(pi >= 0 && seq != "") add_kmer_tokens(nPeaks_cur + pi, seq);
    pf.close();

    if(pi != nPeaks - 1){
      cerr << "The number of peak sequences should match the number of columns in the corresponding count matrix!" << endl;
      exit(EXIT_FAILURE);
    }

    nPeaks_cur += nPeaks;
  }
  assert(pnt == args_->nPeaks_total);

  unsigned cfn = 0, nCells_cur = 0; nPeaks_cur = 0;
  for(auto cp_matrix: args_->cp_matrix_list){
    cout << "Reading a " << args_->nCells_list[cfn] << " cell by " << args_->nPeaks_list[cfn]
         << " peak count matrix from \'" << cp_matrix << "\'" << endl;
    string str, header, mat_format = cp_matrix.substr(cp_matrix.find_last_of(".") + 1);
    bool bin;
    unsigned nCells, nPeaks, peak, cell;
    if(mat_format == "gz"){
#ifdef COMPRESS_FILE
        ifstream ifs2(cp_matrix);
        if (!ifs2.good()) exit(EXIT_FAILURE);
        filtering_istream cf;
        cf.push(gzip_decompressor());
        cf.push(ifs2);

        getline(cf, header);
        bin = check_header(header);

        cf >> nCells >> nPeaks >> str;
        assert(nCells == args_->nCells_list[cfn] && nPeaks == args_->nPeaks_list[cfn]);     

        if(bin)
          while(cf >> cell){
            cf >> peak;
            addCellLabel(nPeaks_cur + peak - 1, nCells_cur + cell);
          }
        else
          while(cf >> cell){
            cf >> peak >> str;
            addCellLabel(nPeaks_cur + peak - 1, nCells_cur + cell);
          }

        nCells_cur += nCells;
        nPeaks_cur += nPeaks;

        ifs2.close();
#endif
    } else {
        ifstream cf(cp_matrix);
        if(!cf.good()) exit(EXIT_FAILURE);

        getline(cf, header);
        bin = check_header(header);

        cf >> nCells >> nPeaks >> str;
        assert(nCells == args_->nCells_list[cfn] && nPeaks == args_->nPeaks_list[cfn]);  

        if(bin)
          while(cf >> cell){
            cf >> peak;
            addCellLabel(nPeaks_cur + peak - 1, nCells_cur + cell);
          }
        else
          while(cf >> cell){
            cf >> peak >> str;
            addCellLabel(nPeaks_cur + peak - 1, nCells_cur + cell);
          }

        nCells_cur += nCells;
        nPeaks_cur += nPeaks;

        cf.close();
    }

    for(unsigned pi = nPeaks_cur - args_->nPeaks_list[cfn]; pi < nPeaks_cur; pi++)
      examples_[pi].dataset = cfn;

    // if(args_->batchLabels)
    //   for(unsigned pi = nPeaks_cur - args_->nPeaks_list[cfn]; pi < nPeaks_cur; pi++){
    //     int32_t wid = dict_->getId(args_->label + "B" + std::to_string(cfn + 1));
    //     assert(wid >= 0);
    //     examples_[pi].RHSTokens.push_back(make_pair(wid, 1.0));
    //   }

    cfn ++;
  }

  assert(nCells_cur == args_->nCells_total && nPeaks_cur == args_->nPeaks_total);

  size_ = args_->nPeaks_total * args_->exmpPerPeak;
  cout << "Number of datasets: " << args_->nCells_list.size() << endl
       << "Total number of cells: " << args_->nCells_total << endl
       << "Total number of peaks: " << args_->nPeaks_total << endl
       << "Number of training examples: " << size_ << " (" << args_->exmpPerPeak << " per peak)"
       << "\n-----------------\n" << endl;

// for(auto ex: examples_){
   // for(auto word: ex.RHSTokens) cerr << dict_->getSymbol(word.first) << " ";
   // cerr << ": ";
   // for(auto word: ex.LHSTokens) cerr << dict_->getSymbol(word.first) << " ";
   // cerr << ex.LHSTokens.size();
   // cerr << endl;
// }

}

int first_ngram_idx(int L, int i, int n){
  if(i <= L - n + 1)
    return(i * (n - 1));
  else {
    int first = (L - n + 1) * (n - 1);
    for(int j = L - n + 2; j <= i; j++)
      first += L - j;
    return(first);
  }
}

// Convert an example for training/testing if needed.
// In the case of trainMode=1, a random label from r.h.s will be selected
// as label, and the rest of labels from r.h.s. will be input features
void InternDataHandler::convert(
    const ParseResults& example,
    // const std::vector<int32_t>& hashes,
    const std::vector<Base>& ngrams,
    ParseResults& rslt) const {

#define LHS_SIZE (example.LHSTokens.size())
  int num_tokens = (args_->sampleLen == -1) ? LHS_SIZE : (args_->sampleLen - args_->k + 1);
  if(num_tokens > LHS_SIZE) num_tokens = LHS_SIZE;
  int first = (num_tokens == LHS_SIZE) ? 0 : (rand() % (LHS_SIZE - num_tokens + 1));

  rslt.weight = example.weight;
  rslt.dataset = example.dataset;
  rslt.LHSTokens.clear();
  rslt.RHSTokens.clear();
  auto first_token = example.LHSTokens.begin() + first;
  rslt.LHSTokens.insert(rslt.LHSTokens.end(), first_token, first_token + num_tokens);

  if(args_->ngrams > 1){
    int first_ngram = first_ngram_idx(LHS_SIZE, first, args_->ngrams),
        last_ngram  = first_ngram_idx(LHS_SIZE, first + num_tokens - 2, args_->ngrams);
    if(last_ngram >= first_ngram){
      rslt.LHSTokens.insert(rslt.LHSTokens.end(), ngrams.begin() + first_ngram, ngrams.begin() + last_ngram);
    }
// std::cerr << num_tokens << "\t" << last_ngram - first_ngram + 1 << std::endl;
  }

  if (args_->trainMode == 0) {
    // lhs is the same, pick one random label as rhs
    assert(example.LHSTokens.size() > 0);
    assert(example.RHSTokens.size() > 0);
    auto idx = rand() % example.RHSTokens.size();
    rslt.RHSTokens.push_back(example.RHSTokens[idx]);
  } else {
    assert(example.RHSTokens.size() > 1);
    if (args_->trainMode == 1) {
      // pick one random label as rhs and the rest is lhs
      auto idx = rand() % example.RHSTokens.size();
      for (unsigned int i = 0; i < example.RHSTokens.size(); i++) {
        auto tok = example.RHSTokens[i];
        if (i == idx) {
          rslt.RHSTokens.push_back(tok);
        } else {
          rslt.LHSTokens.push_back(tok);
        }
      }
    } else
    if (args_->trainMode == 2) {
      // pick one random label as lhs and the rest is rhs
      auto idx = rand() % example.RHSTokens.size();
      for (unsigned int i = 0; i < example.RHSTokens.size(); i++) {
        auto tok = example.RHSTokens[i];
        if (i == idx) {
          rslt.LHSTokens.push_back(tok);
        } else {
          rslt.RHSTokens.push_back(tok);
        }
      }
    } else
    if (args_->trainMode == 3) {
      // pick two random labels, one as lhs and the other as rhs
      auto idx = rand() % example.RHSTokens.size();
      unsigned int idx2;
      do {
        idx2 = rand() % example.RHSTokens.size();
      } while (idx2 == idx);
      rslt.LHSTokens.push_back(example.RHSTokens[idx]);
      rslt.RHSTokens.push_back(example.RHSTokens[idx2]);
    } else
    if (args_->trainMode == 4) {
      // the first one as lhs and the second one as rhs
      rslt.LHSTokens.push_back(example.RHSTokens[0]);
      rslt.RHSTokens.push_back(example.RHSTokens[1]);
    }
  }
}

void InternDataHandler::getWordExamples(
    const vector<Base>& doc,
    vector<ParseResults>& rslts) const {

  rslts.clear();
  for (int widx = 0; widx < (int)(doc.size()); widx++) {
    ParseResults rslt;
    rslt.LHSTokens.clear();
    rslt.RHSTokens.clear();
    rslt.RHSTokens.push_back(doc[widx]);
    for (unsigned int i = max(widx - args_->ws, 0);
         i < min(size_t(widx + args_->ws), doc.size()); i++) {
      if ((int)i != widx) {
        rslt.LHSTokens.push_back(doc[i]);
      }
    }
    rslt.weight = args_->wordWeight;
    rslts.emplace_back(rslt);
  }
}

void InternDataHandler::getWordExamples(
    int idx,
    vector<ParseResults>& rslts) const {
  assert(idx >= 0 && idx < size_);
  idx = idx % args_->nPeaks_total;
  const auto& example = examples_[idx];
  getWordExamples(example.LHSTokens, rslts);
}

void InternDataHandler::addExample(const ParseResults& example) {
  examples_.push_back(example);
  size_++;
}

void InternDataHandler::getExampleById(int32_t idx, ParseResults& rslt) const {
  assert(idx >= 0 && idx < size_);
  idx = idx % args_->nPeaks_total;
  convert(examples_[idx], ngrams_[idx], rslt);
}

void InternDataHandler::getNextExample(ParseResults& rslt) {
//  assert(args_->nPeaks_total > 0);
  idx_ = idx_ + 1;
  // go back to the beginning of the examples if we reach the end
  if (idx_ >= args_->nPeaks_total) {
    idx_ = idx_ - args_->nPeaks_total;
  }
  convert(examples_[idx_], ngrams_[idx_], rslt);
}

void InternDataHandler::getRandomExample(ParseResults& rslt) const {
//  assert(args_->nPeaks_total > 0);
  int32_t idx = rand() % args_->nPeaks_total;
  convert(examples_[idx], ngrams_[idx], rslt);
}

void InternDataHandler::getKRandomExamples(int K, vector<ParseResults>& c) {
  auto kSamples = min(K, (int)args_->nPeaks_total);
  for (int i = 0; i < kSamples; i++) {
    ParseResults example;
    getRandomExample(example);
    c.push_back(example);
  }
}

void InternDataHandler::getNextKExamples(int K, vector<ParseResults>& c) {
  auto kSamples = min(K, (int)args_->nPeaks_total);
  for (int i = 0; i < kSamples; i++) {
    idx_ = (idx_ + 1) % args_->nPeaks_total;
    ParseResults example;
    convert(examples_[idx_], ngrams_[idx_], example);
    c.push_back(example);
  }
}

void InternDataHandler::getRandomWord(vector<Base>& result) {
  result.push_back(word_negatives_[word_iter_]);
  word_iter_++;
  if (word_iter_ >= (int)word_negatives_.size()) {
    word_iter_ = 0;
  }
}

void InternDataHandler::initWordNegatives() {
  word_iter_ = 0;
  word_negatives_.clear();
//  assert(args_->nPeaks_total > 0);
  for (int i = 0; i < MAX_WORD_NEGATIVES_SIZE; i++) {
    word_negatives_.emplace_back(genRandomWord());
  }
}

Base InternDataHandler::genRandomWord() const {
//  assert(args_->nPeaks_total > 0);
  auto& ex = examples_[rand() % args_->nPeaks_total];
  int r = rand() % ex.LHSTokens.size();
  return ex.LHSTokens[r];
}

// Randomly sample one example and randomly sample a label from this example
// The result is usually used as negative samples in training
void InternDataHandler::getRandomRHS(vector<Base>& results, unsigned dataset) const {
//  assert(args_->nPeaks_total > 0);
  results.clear();

  unsigned pi = args_->first_peak_idx[dataset] + (rand() % args_->nPeaks_list[dataset]);
  // std::cerr << args_->first_peak_idx[dataset] << " " << pi << " - ";
  auto& ex = examples_[pi];
  unsigned int r = rand() % ex.RHSTokens.size();
  if (args_->trainMode == 2) {
    for (unsigned int i = 0; i < ex.RHSTokens.size(); i++) {
      if (i != r) {
        results.push_back(ex.RHSTokens[i]);
      }
    }
  } else {
    results.push_back(ex.RHSTokens[r]);
  }
}

void InternDataHandler::save(std::ostream& out) {
  out << "data size : " << args_->nPeaks_total << endl;
  for (auto& example : examples_) {
    out << "lhs : ";
    for (auto t : example.LHSTokens) {out << t.first << ':' << t.second << ' ';}
    out << endl;
    out << "rhs : ";
    for (auto t : example.RHSTokens) {out << t.first << ':' << t.second << ' ';}
    out << endl;
  }
}

} // unamespace starspace

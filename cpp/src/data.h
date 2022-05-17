/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#pragma once

#include "dict.h"
#include "parser.h"
#include "utils/utils.h"
#include <string>
#include <vector>
#include <fstream>

namespace starspace {

class InternDataHandler {
public:
  explicit InternDataHandler(std::shared_ptr<Args> args, std::shared_ptr<Dictionary> dict, std::shared_ptr<DataParser> parser);

  virtual void loadFromFile();

  virtual void convert(const ParseResults& example, const std::vector<Base>& ngrams, ParseResults& rslt) const;

  virtual void getRandomRHS(std::vector<Base>& results, unsigned dataset) const;

  virtual void save(std::ostream& out);

  virtual void getWordExamples(int idx, std::vector<ParseResults>& rslt) const;

  void add_kmer_tokens(int pi, std::string seq){
    std::vector<std::string> tokens;
    for(unsigned i = 0; i < seq.length() - args_->k + 1; i++){
      std::string kmer = "", rc_kmer = "";
      bool N = false;
      for(unsigned j = i; (j < i + args_->k) && !N; j++){
        switch(seq[j]){
          case 'A':
          case 'a':
            kmer = kmer + "A";
            rc_kmer = "T" + rc_kmer;
            break;
          case 'C':
          case 'c':
            kmer = kmer + "C";
            rc_kmer = "G" + rc_kmer;
            break;
          case 'T':
          case 't':
            kmer = kmer + "T";
            rc_kmer = "A" + rc_kmer;
            break;
          case 'G':
          case 'g':
            kmer = kmer + "G";
            rc_kmer = "C" + rc_kmer;
            break;
          default:
            N = true;
            i = j;
        }
      }
      if(!N){
         int32_t kmer_wid = dict_->getId(kmer), rc_kmer_wid = dict_->getId(rc_kmer);
         if(!(kmer_wid < 0)){
           tokens.push_back(kmer);
           examples_[pi].LHSTokens.push_back(std::make_pair(kmer_wid, 1.0));
         } else if(!(rc_kmer_wid < 0)){
           tokens.push_back(rc_kmer);
           examples_[pi].LHSTokens.push_back(std::make_pair(rc_kmer_wid, 1.0));
         } else {
           std::cerr << "Invalid DNA k-mer! \'" << kmer << "\'" << std::endl;
           continue;
         }
      }
    }

    assert(examples_[pi].LHSTokens.size() > 0);

    if(args_->ngrams > 1){
      parser_->addNgrams(tokens, ngrams_[pi], args_->ngrams);
      // for(auto token: tokens) hashes_[pi].push_back(dict_->hash(token));
    }
  }

  void getWordExamples(
      const std::vector<Base>& doc,
      std::vector<ParseResults>& rslt) const;

  void addExample(const ParseResults& example);

  void addCellLabel(unsigned peak, unsigned cell){
    int32_t wid = dict_->getId(args_->label + "C" + std::to_string(cell));
    assert(wid >= 0);
    examples_[peak].RHSTokens.push_back(std::make_pair(wid, 1.0));
  }

  void getExampleById(int32_t idx, ParseResults& rslt) const;

  void getNextExample(ParseResults& rslt);

  void getRandomExample(ParseResults& rslt) const;

  void getKRandomExamples(int K, std::vector<ParseResults>& c);

  void getNextKExamples(int K, std::vector<ParseResults>& c);

  size_t getSize() const { return size_; };

  void errorOnZeroExample(const std::string& fileName);

  void initWordNegatives();
  void getRandomWord(std::vector<Base>& result);


protected:
  virtual Base genRandomWord() const;

  static const int32_t MAX_VOCAB_SIZE = 10000000;
  static const int32_t MAX_WORD_NEGATIVES_SIZE = 10000000;

  std::shared_ptr<Args> args_;
  std::shared_ptr<Dictionary> dict_;
  std::shared_ptr<DataParser> parser_;
  std::vector<ParseResults> examples_;
  std::vector<std::vector<Base>> ngrams_;
  // std::vector<std::vector<int32_t>> hashes_;

  int32_t idx_ = -1;
  int32_t size_ = 0;

  int32_t word_iter_;
  std::vector<Base> word_negatives_;
};

}

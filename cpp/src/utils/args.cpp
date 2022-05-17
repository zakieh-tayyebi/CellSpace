/**
 * Copyright (c) Facebook, Inc. and its affiliates.
 *
 * This source code is licensed under the MIT license found in the
 * LICENSE file in the root directory of this source tree.
 */

#include "args.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <cstring>
#include <assert.h>

using namespace std;

namespace starspace {

Args::Args() {
  lr = 0.01;
  termLr = 1e-9;
  norm = 1.0;
  margin = 0.05;
  wordWeight = 0.5;
  initRandSd = 0.001;
  dropoutLHS = 0.0;
  dropoutRHS = 0.0;
  p = 0.5;
  ws = 5;
  maxTrainTime = 60*60*24*100;
  validationPatience = 10;
  thread = 10;
  maxNegSamples = 10;
  negSearchLimit = 50;
  minCount = 1;
  minCountLabel = 1;
  K = 5;
  batchSize = 5;
  verbose = false;
  debug = false;
  adagrad = true;
  normalizeText = false;
  trainMode = 0;
  // fileFormat = "fastText";
  label = "__label__";
  bucket = 2000000;
  isTrain = true;
  shareEmb = true;
  saveEveryEpoch = false;
  saveTempModel = false;
  useWeight = false;
  trainWord = false;
  excludeLHS = false;
  weightSep = ':';
  numGzFile = 1;

  loss = "hinge";
  similarity = "cosine";
  ngrams = 3;
  dim = 30;
  epoch = 20;
  exmpPerPeak = 20;
  k = 8;
  sampleLen = 150;
  // fixedFeatEmb = false;
  // batchLabels = false;
}

bool Args::isTrue(string arg) {
  std::transform(arg.begin(), arg.end(), arg.begin(),
      [&](char c) { return tolower(c); }
  );
  return (arg == "true" || arg == "1");
}

bool valid_file(std::string path){
    if(path.empty()){
      cerr << "\'" << path << "\' cannot be opened!" << endl;
      return false;
    } else {
      std::ifstream fin(path, std::ifstream::in);
      if (!fin.is_open()){
        cerr << "\'" << path << "\' cannot be opened!" << endl;
        return false;
      }
      fin.close();
    }
    return true;
}

inline bool check_header(string header){
    if(header == "%%MatrixMarket matrix coordinate pattern general") return(true);
    else if(header == "%%MatrixMarket matrix coordinate integer general") return(false);
    else {
        cerr << "Unsupported matrix format!" << endl;
        exit(EXIT_FAILURE);
    }
}

void Args::parseArgs(int argc, char** argv) {
  string a0 = argv[0];
  isCellSpace = a0.substr(a0.length() - 9, 9) == "CellSpace";

  if (argc < 4) {
    // cerr << "Usage: need to specify whether it is train or test.\n";
    printHelp();
    exit(EXIT_FAILURE);
  }
  // if (strcmp(argv[1], "train") == 0) {
  //   isTrain = true;
  // } else if (strcmp(argv[1], "test") == 0) {
  //   isTrain = false;
  //   if(isCellSpace){
  //       std::cerr << "CellSpace test has not been implemented yet!" << std::endl;
  //       exit(EXIT_FAILURE);
  //   }
  // } else if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-help") == 0) {
  //   std::cerr << "Here is the help! Usage:" << std::endl;
  //   printHelp();
  //   exit(EXIT_FAILURE);
  // } else {
  //   cerr << "Usage: the first argument should be either train or test.\n";
  //   printHelp();
  //   exit(EXIT_FAILURE);
  // }
  int i = 1;
  while (i < argc) {
    if (argv[i][0] != '-') {
      cout << "Provided argument without a dash! Usage:" << endl;
      printHelp();
      exit(EXIT_FAILURE);
    }

    // handling "--"
    if (strlen(argv[i]) >= 2 && argv[i][1] == '-') {
      argv[i] = argv[i] + 1;
    }

    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "-help") == 0) {
      // std::cerr << "Here is the help! Usage:" << std::endl;
      printHelp();
      exit(EXIT_FAILURE);
    // } else if (strcmp(argv[i], "-trainFile") == 0 && !isCellSpace) {
    //   trainFile = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-validationFile") == 0 && !isCellSpace) {
    //   validationFile = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-testFile") == 0 && !isCellSpace) {
    //   testFile = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-predictionFile") == 0 && !isCellSpace) {
    //   predictionFile = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-basedoc") == 0 && !isCellSpace) {
    //   basedoc = string(argv[i + 1]);
    } else if (strcmp(argv[i], "-output") == 0) {
      model = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-initModel") == 0 && !isCellSpace) {
    //   initModel = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-fileFormat") == 0 && !isCellSpace) {
    //   fileFormat = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-compressFile") == 0 && !isCellSpace) {
    //   compressFile = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-numGzFile") == 0) {
    //   numGzFile = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-label") == 0) {
      label = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-weightSep") == 0 && !isCellSpace) {
    //   weightSep = argv[i + 1][0];
    // } else if (strcmp(argv[i], "-loss") == 0 && !isCellSpace) {
    //   loss = string(argv[i + 1]);
    // } else if (strcmp(argv[i], "-similarity") == 0) {
    //   similarity = string(argv[i + 1]);
    } else if (strcmp(argv[i], "-lr") == 0) {
      lr = atof(argv[i + 1]);
    } else if (strcmp(argv[i], "-p") == 0) {
      p = atof(argv[i + 1]);
    // } else if (strcmp(argv[i], "-termLr") == 0) {
    //   termLr = atof(argv[i + 1]);
    // } else if (strcmp(argv[i], "-norm") == 0) {
    //   norm = atof(argv[i + 1]);
    } else if (strcmp(argv[i], "-margin") == 0) {
      margin = atof(argv[i + 1]);
    } else if (strcmp(argv[i], "-initRandSd") == 0) {
      initRandSd = atof(argv[i + 1]);
    // } else if (strcmp(argv[i], "-dropoutLHS") == 0) {
    //   dropoutLHS = atof(argv[i + 1]);
    // } else if (strcmp(argv[i], "-dropoutRHS") == 0) {
    //   dropoutRHS = atof(argv[i + 1]);
    // } else if (strcmp(argv[i], "-wordWeight") == 0) {
    //   wordWeight = atof(argv[i + 1]);
    } else if (strcmp(argv[i], "-dim") == 0) {
      dim = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-epoch") == 0) {
      epoch = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-ws") == 0) {
    //   ws = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-maxTrainTime") == 0) {
      maxTrainTime = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-validationPatience") == 0) {
    //   validationPatience = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-thread") == 0) {
      thread = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-maxNegSamples") == 0) {
      maxNegSamples = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-negSearchLimit") == 0) {
      negSearchLimit = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-minCount") == 0 && !isCellSpace) {
    //   minCount = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-minCountLabel") == 0 && !isCellSpace) {
    //   minCountLabel = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-bucket") == 0) {
      bucket = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-ngrams") == 0) {
      ngrams = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-K") == 0) {
    //   K = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-batchSize") == 0) {
      batchSize = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-trainMode") == 0 && !isCellSpace) {
    //   trainMode = atoi(argv[i + 1]);
    // } else if (strcmp(argv[i], "-fixedFeatEmb") == 0 && isCellSpace) {
    //   fixedFeatEmb = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-batchLabels") == 0 && isCellSpace) {
    //   batchLabels = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-verbose") == 0) {
    //   verbose = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-debug") == 0) {
    //   debug = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-adagrad") == 0) {
    //   adagrad = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-shareEmb") == 0) {
    //   shareEmb = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-normalizeText") == 0 && !isCellSpace) {
    //   normalizeText = isTrue(string(argv[i + 1]));
    } else if (strcmp(argv[i], "-saveEveryEpoch") == 0) {
      saveEveryEpoch = isTrue(string(argv[i + 1]));
    } else if (strcmp(argv[i], "-saveTempModel") == 0) {
      saveTempModel = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-useWeight") == 0 && !isCellSpace) {
    //   useWeight = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-trainWord") == 0 && !isCellSpace) {
    //   trainWord = isTrue(string(argv[i + 1]));
    // } else if (strcmp(argv[i], "-excludeLHS") == 0) {
    //   excludeLHS = isTrue(string(argv[i + 1]));
    } else if (strcmp(argv[i], "-cpMat") == 0 && isCellSpace) {
      assert(i + 1 < argc);
      while(i + 1 < argc && argv[i + 1][0] != '-')
        cp_matrix_list.push_back(string(argv[(i++) + 1]));
      i--;
    } else if (strcmp(argv[i], "-peaks") == 0 && isCellSpace) {
      assert(i + 1 < argc);
      while(i + 1 < argc && argv[i + 1][0] != '-')
        peaks_list.push_back(string(argv[(i++) + 1]));
      i--;
    // } else if (strcmp(argv[i], "-featEmb") == 0 && isCellSpace) {
    //   feat_emb = string(argv[i + 1]);
    } else if (strcmp(argv[i], "-k") == 0 && isCellSpace) {
      k = atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-sampleLen") == 0 && isCellSpace) {
      sampleLen = (strcmp(argv[i + 1], "given") == 0) ? -1 : atoi(argv[i + 1]);
    } else if (strcmp(argv[i], "-exmpPerPeak") == 0 && isCellSpace) {
      exmpPerPeak = atoi(argv[i + 1]);
    } else {
      cerr << "Unknown argument: " << argv[i] << std::endl;
      printHelp();
      exit(EXIT_FAILURE);
    }
    i += 2;
  }

  if(isCellSpace){
    trainMode = 0;
    minCount = 1;
    minCountLabel = 1;
    trainWord = false;
    useWeight = false;
    normalizeText = false;
    compressFile = "";

    if(peaks_list.size() != cp_matrix_list.size()){
      cerr << "The number of input files for peak sequences should match the number of input files for count matrices!" << endl;
      exit(EXIT_FAILURE);
    }
    nBatches = cp_matrix_list.size();

    nCells_total = 0; nCells_list.clear();
    nPeaks_total = 0; nPeaks_list.clear();
    first_peak_idx.clear();
    for(auto cp_matrix: cp_matrix_list){
      unsigned nCells = 0, nPeaks = 0;
      if(!valid_file(cp_matrix)) exit(EXIT_FAILURE);
      else {
        string header, mat_format = cp_matrix.substr(cp_matrix.find_last_of(".") + 1);
        if(mat_format == "gz"){
#ifdef COMPRESS_FILE
            ifstream ifs2(cp_matrix);
            if(!ifs2.good()) exit(EXIT_FAILURE);
            filtering_istream fin;
            fin.push(gzip_decompressor());
            fin.push(ifs2);
            getline(fin, header);
            bool bin = check_header(header);
            fin >> nCells >> nPeaks;
            ifs2.close();
#endif
        } else {
            std::ifstream fin(cp_matrix, std::ifstream::in);
            getline(fin, header);
            bool bin = check_header(header);
            fin >> nCells >> nPeaks;
            fin.close();
        }
        if(nCells < 2 || nPeaks < 2){
          cerr << "There must be at least two cells and two peaks in a dataset!" << endl;
          exit(EXIT_FAILURE);
        } else {
          first_peak_idx.push_back(nPeaks_total);
          nCells_list.push_back(nCells); nCells_total += nCells;
          nPeaks_list.push_back(nPeaks); nPeaks_total += nPeaks;
          // cerr << nCells_list.size() << ") cells: " << nCells << ", peaks: " << nPeaks << endl;
        }
      }
    }

    for(auto peaks: peaks_list)
      if(!valid_file(peaks)) exit(EXIT_FAILURE);

    // if(feat_emb != "" && !valid_file(feat_emb)) exit(EXIT_FAILURE);

    if(k < 2 || k > 16){
      cerr << "k must be >1 and <17" << endl;
      exit(EXIT_FAILURE);
    }
    if(sampleLen != -1 && sampleLen < k){
      cerr << "sampleLen must be >=k" << endl;
      exit(EXIT_FAILURE);
    }
    if(exmpPerPeak < 1){
      cerr << "exmpPerPeak must be >0" << endl;
      exit(EXIT_FAILURE);
    }
  } else {
    if (isTrain) {
      if (model.empty()) {
        cerr << "Empty output path." << endl;
        printHelp();
        exit(EXIT_FAILURE);
      }
    } else {
      if (testFile.empty() || model.empty()) {
        cerr << "Empty test file or model path." << endl;
        printHelp();
        exit(EXIT_FAILURE);
      }
    }
  }
}

void Args::printHelp() {
    if(isCellSpace){
       cout << "\n"
       << "\"CellSpace ...\"\n"

       << "\nThe following arguments are mandatory:\n"
       << "  -output          output file to write the embedding matrix for cells and k-mers (.tsv)\n"
       << "  -cpMat           sparse cell by peak/tile count matrix (.mtx)\n"
       << "  -peaks           multi-fasta file containing peak/tile DNA sequences with the order they appear in the corresponding count matrix (.fa)\n"
       // << "  -featEmb         optional embedding matrix to initialize k-mer embeddings from a previously trained CellSpace model (.tsv)\n"

       << "\nThe following arguments are optional:\n"
       << "  -dim             size of embedding vectors [" << dim << "]\n"
       << "  -ngrams          max length of k-mer ngram [" << ngrams << "]\n"
       << "  -k               k-mer length [" << k << "]\n"
       << "  -sampleLen       length of the sequences randomly sampled from the peak/tile DNA sequences [" << (sampleLen == -1 ? "given" : to_string(sampleLen)) << "]\n"
       << "  -exmpPerPeak     number of training examples per peak/tile [" << exmpPerPeak << "]\n"
       << "  -epoch           number of epochs [" << epoch << "]\n"
       // << "  -fixedFeatEmb    whether feature (k-mer) embeddings should be kept fixed [" << fixedFeatEmb << "]\n"
       // << "  -batchLabels     whether the batch (i.e. dataset) should be included as a label (i.e. RHS) in training [" << batchLabels << "]\n"
       << "  -margin          margin parameter in hinge loss. [" << margin << "]\n"
       // << "  -similarity      takes value in [cosine, dot]. Whether to use cosine or dot product as similarity function in  hinge loss. [" << similarity << "]\n"
       << "  -bucket          number of buckets [" << bucket << "]\n"
       << "  -label           labels prefix [" << label << "]\n"
       << "  -lr              learning rate [" << lr << "]\n"
       << "  -maxTrainTime    max train time (seconds) [" << maxTrainTime << "]\n"
       << "  -negSearchLimit  number of negative labels sampled per dataset [" << negSearchLimit << "]\n"
       << "  -maxNegSamples   max number of negatives in a batch update [" << maxNegSamples << "]\n"
       << "  -p               the embedding of an entity equals the sum of its M feature embedding vectors devided by M^p. [" << p << "]\n"
       // << "  -adagrad         whether to use adagrad in training [" << adagrad << "]\n"
       // << "  -shareEmb        whether to use the same embedding matrix for LHS and RHS. [" << shareEmb << "]\n"
       // << "  -dropoutLHS      dropout probability for LHS features. [" << dropoutLHS << "]\n"
       // << "  -dropoutRHS      dropout probability for RHS features. [" << dropoutRHS << "]\n"
       << "  -initRandSd      initial values of embeddings are randomly generated from normal distribution with mean=0 and standard deviation=initRandSd. [" << initRandSd << "]\n"
       << "  -batchSize       size of mini batch in training. [" << batchSize << "]\n"
       // << "  -verbose         verbosity level [" << verbose << "]\n"
       // << "  -debug           whether it's in debug mode [" << debug << "]\n"
       // << "  -saveEveryEpoch  save intermediate models after each epoch [" << saveEveryEpoch << "]\n"
       // << "  -saveTempModel   save intermediate models after each epoch with an unique name including epoch number [" << saveTempModel << "]\n"
       << "  -thread          number of threads [" << thread << "]\n\n";
    }
}

void Args::printArgs() {
  cout << "CellSpace Arguments: \n"
       << "dim: " << dim << endl
       << "ngrams: " << ngrams << endl
       << "k: " << k << endl
       << "sampleLen: " << (sampleLen == -1 ? "given" : to_string(sampleLen)) << endl
       << "exmpPerPeak: " << exmpPerPeak << endl
       << "epoch: " << epoch << endl
       << "margin: " << margin << endl
       << "bucket: " << bucket << endl
       << "label: " << label << endl
       << "lr: " << lr << endl
       << "maxTrainTime: " << maxTrainTime << endl
       << "negSearchLimit: " << negSearchLimit << endl
       << "maxNegSamples: " << maxNegSamples << endl
       << "p: " << p << endl
       << "initRandSd: " << initRandSd << endl
       << "batchSize: " << batchSize << endl
       << "thread: " << thread << endl
       << "-----------------" << endl;
}

void Args::save(std::ostream& out) {
  out.write((char*) &(dim), sizeof(int));
  out.write((char*) &(epoch), sizeof(int));
  out.write((char*) &(maxTrainTime), sizeof(int));
  out.write((char*) &(minCount), sizeof(int));
  out.write((char*) &(minCountLabel), sizeof(int));
  out.write((char*) &(maxNegSamples), sizeof(int));
  out.write((char*) &(negSearchLimit), sizeof(int));
  out.write((char*) &(ngrams), sizeof(int));
  out.write((char*) &(bucket), sizeof(int));
  out.write((char*) &(trainMode), sizeof(int));
  out.write((char*) &(shareEmb), sizeof(bool));
  out.write((char*) &(useWeight), sizeof(bool));
  out.write((char*) &(weightSep), sizeof(char));
  size_t size; // = fileFormat.size();
  out.write((char*) &(size), sizeof(size_t));
  // out.write((char*) &(fileFormat[0]), size);
  size = similarity.size();
  out.write((char*) &(size), sizeof(size_t));
  out.write((char*) &(similarity[0]), size);
  out.write((char*) &(batchSize), sizeof(int));
}

void Args::load(std::istream& in) {
  in.read((char*) &(dim), sizeof(int));
  in.read((char*) &(epoch), sizeof(int));
  // in.read((char*) &(maxTrainTime), sizeof(int));
  in.read((char*) &(minCount), sizeof(int));
  in.read((char*) &(minCountLabel), sizeof(int));
  in.read((char*) &(maxNegSamples), sizeof(int));
  in.read((char*) &(negSearchLimit), sizeof(int));
  in.read((char*) &(ngrams), sizeof(int));
  in.read((char*) &(bucket), sizeof(int));
  in.read((char*) &(trainMode), sizeof(int));
  // in.read((char*) &(shareEmb), sizeof(bool));
  in.read((char*) &(useWeight), sizeof(bool));
  in.read((char*) &(weightSep), sizeof(char));
  size_t size;
  in.read((char*) &(size), sizeof(size_t));
  // fileFormat.resize(size);
  // in.read((char*) &(fileFormat[0]), size);
  in.read((char*) &(size), sizeof(size_t));
  similarity.resize(size);
  in.read((char*) &(similarity[0]), size);
  in.read((char*) &(batchSize), sizeof(int));
}

}

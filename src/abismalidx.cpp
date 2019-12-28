/* Copyright (C) 2018 Andrew D. Smith and Meng Zhou
 *
 * Authors: Andrew D. Smith
 *
 * This file is part of ABISMAL.
 *
 * ABISMAL is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ABISMAL is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"
#include "OptionParser.hpp"

#include "AbismalIndex.hpp"
#include <omp.h>

#include <iostream>
#include <fstream>
#include <unordered_set>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>

using std::vector;
using std::unordered_set;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;

void
BuildIndex(const bool VERBOSE, const string &genome_file,
           AbismalIndex &abismal_index) {

  std::ifstream in(genome_file);
  if (!in)
    throw runtime_error("bad genome file: " + genome_file);

  if (VERBOSE)
    cerr << "[loading genome]" << endl;
  load_genome(genome_file, abismal_index.genome, abismal_index.cl);

  if (VERBOSE)
    cerr << "[validating genome]" << endl;

  transform(begin(abismal_index.genome), end(abismal_index.genome),
            begin(abismal_index.genome),
            [](char c) {return to_valid_five_letter[static_cast<size_t>(c)];});

  unordered_set<uint32_t> big_buckets;
  abismal_index.get_bucket_sizes(big_buckets);

  abismal_index.hash_genome(big_buckets);

  abismal_index.sort_buckets();
}

int main(int argc, const char **argv) {
  try {

    bool VERBOSE = false;

    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "build abismal index",
                           "<genome-fasta> <abismal-index-file>", 2);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("too-big", 'B', "ignore buckets bigger than this",
                      false, AbismalIndex::valid_bucket_limit);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string genome_file = leftover_args.front();
    const string outfile = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    omp_set_num_threads(n_threads);

    AbismalIndex::VERBOSE = VERBOSE;

    AbismalIndex abismal_index;
    BuildIndex(VERBOSE, genome_file, abismal_index);

    if (VERBOSE)
      cerr << "[writing abismal index to: " << outfile << "]" << endl;

    abismal_index.write(outfile);

  }
  catch (const std::runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

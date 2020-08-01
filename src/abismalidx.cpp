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
#include "dna_four_bit.hpp"

#include "AbismalIndex.hpp"
#include <omp.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <unordered_set>

using std::vector;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;
using std::unordered_set;

static void
BuildIndex(const bool VERBOSE,
           const uint32_t max_candidates,
           const uint32_t deadzone_freq,
           const string &genome_file,
           AbismalIndex &ai) {

  std::ifstream in(genome_file);
  if (!in)
    throw runtime_error("bad genome file: " + genome_file);

  if (VERBOSE)
    cerr << "[loading genome]" << endl;
  load_genome(genome_file, ai.genome, ai.cl);

  ai.encode_genome();
  ai.get_bucket_sizes();
  ai.hash_genome();

  ai.sort_buckets(four_letter);
  ai.remove_big_buckets(four_letter, deadzone_freq);
  ai.sort_buckets(two_letter);
  if (max_candidates != std::numeric_limits<uint32_t>::max())
    ai.remove_big_buckets(two_letter, max_candidates);

}

int main(int argc, const char **argv) {
  try {
    bool VERBOSE = false;
    size_t n_threads = omp_get_max_threads();
    uint32_t max_candidates = std::numeric_limits<uint32_t>::max();
    uint32_t deadzone_freq = 200;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "build abismal index",
                           "<genome-fasta> <abismal-index-file>", 2);
    opt_parse.set_show_defaults();
    //opt_parse.add_opt("too-big", 'B', "ignore buckets bigger than this",
    //                  false, AbismalIndex::valid_bucket_limit);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("max-candidates", 'c',
                      "maximum candidates used in practice",
                      false, max_candidates);
    opt_parse.add_opt("solid", 's', "number of solid positions", false,
                      seed::n_solid_positions);
    opt_parse.add_opt("deadzone", 'd', "number of times a "
                      + std::to_string(seed::n_solid_positions) + "-mer should "
                      "appear to be excluded from the genome",
                      false, deadzone_freq);
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
    BuildIndex(VERBOSE, max_candidates, deadzone_freq, genome_file,
               abismal_index);

    if (VERBOSE)
      cerr << "[writing abismal index to: " << outfile << "]" << endl;

    abismal_index.write(outfile, seed::n_solid_positions, max_candidates);

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

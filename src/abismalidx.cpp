/* Copyright (C) 2018 Andrew D. Smith and Meng Zhou
 *
 * Authors: Andrew D. Smith and Guilherme de Sena Brandine
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

#include "abismalidx.hpp"
#include "AbismalIndex.hpp"

#include "OptionParser.hpp"
#include "dna_four_bit.hpp"
#include "smithlab_os.hpp"
#include "smithlab_utils.hpp"

#include <config.h>
#include <omp.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using std::cerr;
using std::endl;
using std::runtime_error;
using std::string;
using std::unordered_set;
using std::vector;

int
abismalidx(int argc, const char **argv) {

  try {

    const string description =
      string("build abismal index (v") + VERSION + string(")");

    string target_regions_file;
    bool VERBOSE = false;
    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<genome-fasta> <abismal-index-file>", 2);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("targets", 'A', "target regions", false,
                      target_regions_file);
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
    const double start_time = omp_get_wtime();
    AbismalIndex::VERBOSE = VERBOSE;

    if (!std::filesystem::exists(genome_file))
      throw runtime_error("file not found: " + genome_file);

    if (!std::filesystem::is_regular_file(genome_file))
      throw runtime_error("not regular file: " + genome_file);

    /****************** START BUILDING INDEX *************/
    AbismalIndex abismal_index;
    if (!target_regions_file.empty())
      abismal_index.create_index(target_regions_file, genome_file);
    else
      abismal_index.create_index(genome_file);

    if (VERBOSE)
      cerr << "[writing abismal index to: " << outfile << "]\n";

    abismal_index.write(outfile);
    if (VERBOSE)
      cerr << "[total indexing time: " << omp_get_wtime() - start_time << "]"
           << endl;
    /****************** END BUILDING INDEX *************/
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

/*
  Test strategy:
  reference genomes:
  (1) hg38: it's the most common
  (2) mm39: similar reasoning
  (3) hs1-t2t: it has different amounts of repeats and Ns
  (4) TAIR10: it's smaller and might behave differently
  (5) tRex1: it's very small

  Reads:
  (1) human: Hodges2011 HSPC / SRR342517
  (2) mouse: SRR5115694
  (4) TAIR10: SRR24436012
  (3) simreads: 10000000 reads, default parameters

  Target regions:
  (1) Promoters
  (2) Randomly shuffled promoters
*/

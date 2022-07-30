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

int
abismalidx(int argc, const char **argv) {

  try {
    bool VERBOSE = false;
    size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "build abismal index",
                           "<genome-fasta> <abismal-index-file>", 2);
    opt_parse.set_show_defaults();
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

    /****************** START BUILDING INDEX *************/
    AbismalIndex abismal_index;
    abismal_index.create_index(genome_file);

    if (VERBOSE)
      cerr << "[writing abismal index to: " << outfile << "]\n";

    abismal_index.write(outfile);
    if (VERBOSE)
      cerr << "[total indexing time: " << omp_get_wtime() - start_time << "]" << endl;
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

#ifndef NO_MAIN
int
main(int argc, const char **argv) {
  return abismalidx(argc, argv);
}
#endif

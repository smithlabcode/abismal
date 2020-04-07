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

#include "AbismalSeed.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

using std::vector;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;
using std::cout;


int main(int argc, const char **argv) {
  try {

    string output_file;

    string pattern_arg;
    uint32_t key_weight;
    uint32_t n_solid_positions;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "");
    opt_parse.add_opt("pattern", 'p', "the pattern", true, pattern_arg);
    opt_parse.add_opt("weight", 'w', "number of hash bits", true, key_weight);
    opt_parse.add_opt("solid", 'S', "number of solid positions",
                      true, n_solid_positions);
    opt_parse.add_opt("output", 'o', "output file", true, output_file);
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
    if (!leftover_args.empty()) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    /****************** END COMMAND LINE OPTIONS *****************/

    const AbismalSeed as(n_solid_positions, key_weight);

    cerr << as << endl;

    as.write(output_file);

    AbismalSeed as2;
    as2.read(output_file);

    if (!(as == as2))
      cerr << as2 << endl;
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  catch (std::bad_alloc &ba) {
    cerr << "ERROR: could not allocate memory" << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

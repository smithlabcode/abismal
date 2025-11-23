/* Copyright (C) 2018-2025 Andrew D. Smith and Guilherme de Sena Brandine
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
#include "smithlab_os.hpp"

#include <config.h>

#include <chrono>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int
abismalidx(int argc, char *argv[]) {

  try {

    const std::string version_str =
      std::string{"(v"} + VERSION + std::string{")"};
    const std::string description = "build abismal index " + version_str;

    std::string target_regions_file;
    bool verbose = false;
    std::size_t n_threads = 1;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), description,
                           "<genome-fasta> <abismal-index-file>", 2);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("targets", 'A', "target regions", false,
                      target_regions_file);
    opt_parse.add_opt("threads", 't', "number of threads", false, n_threads);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << std::endl;
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << std::endl;
      return EXIT_SUCCESS;
    }
    if (leftover_args.size() != 2) {
      std::cerr << opt_parse.help_message() << std::endl;
      return EXIT_SUCCESS;
    }
    const std::string genome_file = leftover_args.front();
    const std::string outfile = leftover_args.back();
    /****************** END COMMAND LINE OPTIONS *****************/

    const auto start_time{std::chrono::high_resolution_clock::now()};
    AbismalIndex::set_n_threads(n_threads);
    AbismalIndex::VERBOSE = verbose;

    if (!std::filesystem::exists(genome_file))
      throw std::runtime_error("file not found: " + genome_file);

    if (!std::filesystem::is_regular_file(genome_file))
      throw std::runtime_error("not regular file: " + genome_file);

    /****************** START BUILDING INDEX *************/
    AbismalIndex abismal_index;
    if (!target_regions_file.empty())
      abismal_index.create_index(target_regions_file, genome_file);
    else
      abismal_index.create_index(genome_file);

    if (verbose)
      std::cerr << "[writing abismal index to: " << outfile << "]\n";

    abismal_index.write(outfile);
    if (verbose) {
      const auto stop_time{std::chrono::high_resolution_clock::now()};
      const auto d = stop_time - start_time;
      std::cerr
        << "[total indexing time: "
        << std::chrono::duration_cast<std::chrono::duration<double>>(d).count()
        << "]\n";
    }
    /****************** END BUILDING INDEX *************/
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

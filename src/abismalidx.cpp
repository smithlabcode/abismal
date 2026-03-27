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
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

static constexpr auto about = "index for another bisulfite mapping algorithm";
static constexpr auto description = R"(
Examples:

abismal idx -t 32 hg38.fa hg38.idx
)";

#include "abismalidx.hpp"

#include "AbismalIndex.hpp"

#include "CLI11/CLI11.hpp"

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
abismalidx(int argc, char *argv[]) {  // NOLINT(*-c-arrays)
  try {
    std::string target_regions_file;
    bool verbose{};
    std::size_t n_threads = 1;
    std::string genome_file;
    std::string outfile;

    CLI::App app{about};
    argv = app.ensure_utf8(argv);
    app.usage("Usage: abismal idx [OPTIONS] input output");
    if (argc >= 2)
      app.footer(description);
    app.get_formatter()->label("REQUIRED", "REQ");
    app.get_formatter()->long_option_alignment_ratio(0.2);
    app.get_formatter()->enable_footer_formatting(false);

    // clang-format off
    app.add_option("input", genome_file,
                   "reference genome in FASTA format (gzip ok)")
      ->option_text("FILE REQUIRED")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("output", outfile,
                   "abismal index file")
      ->option_text("FILE REQUIRED")
      ->required();
    // ->check(CLI::NonexistentPath);
    app.add_option("-A,--targets", target_regions_file,
                   "target regions (BED format)")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-t,--threads", n_threads, "number of threads")
      ->option_text("UINT")
      ->check(CLI::PositiveNumber);
    app.add_flag("-v,--verbose", verbose, "print more run info");
    // clang-format on
    if (argc < 2) {
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const auto start_time{std::chrono::high_resolution_clock::now()};
    AbismalIndex::set_n_threads(n_threads);
    AbismalIndex::VERBOSE = verbose;

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
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

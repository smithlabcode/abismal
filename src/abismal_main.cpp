/* Copyright (C) 2018-2025 Andrew D. Smith and Guilherme Sena
 *
 * Authors: Andrew D. Smith and Guilherme Sena
 *
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 */

#include <config.h>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static const std::string PROGRAM_NAME = "abismal";

struct abismal_command {
  std::string tag;
  std::string description;
  std::function<int(const int, const char **)> fun;

  auto
  operator()(const int argc, const char **argv) const -> int {
    return fun(argc - 1, argv + 1);
  }
};

auto
operator<<(std::ostream &out, const abismal_command &cmd) -> std::ostream & {
  static const std::size_t pad_size = 4;
  static const std::size_t offset = 8;
  static const std::string pad(pad_size, ' ');
  return out << pad << std::left << std::setw(offset) << (cmd.tag + ":")
             << cmd.description;
}

// ADS: not sure of best way to acquire these below beyond simply
// declaring them here
int
abismal(int argc, const char **argv);
int
abismalidx(int argc, const char **argv);
int
simreads(int argc, const char **argv);

void
print_help(const std::vector<abismal_command> &commands) {
  std::cout << "Program: " << PROGRAM_NAME << "\n"
            << "Version: " << VERSION << "\n"
            << "Usage: " << PROGRAM_NAME << " <command> [options]\n"
            << "Commands:" << std::endl;
  for (const auto &c : commands)
    std::cout << c << std::endl;
}

int
main(int argc, const char **argv) {
  try {
    // clang-format off
    std::vector<abismal_command> commands = {
      {"map", "map FASTQ reads to an index or a FASTA reference genome", abismal},
      {"idx", "make an index for a FASTA reference genome", abismalidx},
      {"sim", "simulate WGBS reads for a FASTA reference genome", simreads}
    };
    // clang-format on
    if (argc < 2) {
      print_help(commands);
      return EXIT_SUCCESS;
    }
    const auto has_tag = [&](const abismal_command &a) {
      return a.tag == argv[1];
    };
    const auto the_cmd =
      std::find_if(std::cbegin(commands), std::cend(commands), has_tag);
    if (the_cmd != std::cend(commands))
      return (*the_cmd)(argc, argv);
    std::cerr << "ERROR: invalid command " << argv[1] << std::endl;
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

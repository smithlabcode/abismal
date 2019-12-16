/*
 * Copyright (C) 2018 Andrew D. Smith
 *
 * Author: Andrew D. Smith
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

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <random>
#include <functional>

using std::vector;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;
using std::function;
using std::bind;

static void
sim_frag(const string &genome, const size_t frag_len,
         string &the_frag, size_t &pos) {

  const size_t lim = genome.length() - frag_len + 1;

  do {
    pos = rand() % lim;
    the_frag = string(begin(genome) + pos,
                      begin(genome) + pos + frag_len);
  }
  while (count(begin(the_frag), end(the_frag), 'N') > 0);
}

static void
mutate_frag(const double mutation_rate, const size_t max_mutations,
            function<double()> &unif,
            string &the_frag,
            size_t &n_mutations) {
  for (size_t i = 0; n_mutations < max_mutations && i < the_frag.length(); ++i)
    if (unif() < mutation_rate) {
      char mut = int2base(static_cast<unsigned>(100*unif()) % 4);
      while (mut == the_frag[i])
        mut = int2base(static_cast<unsigned>(100*unif()) % 4);
      the_frag[i] = mut;
      ++n_mutations ;
    }
}


static void
sim_reads(string the_frag, const size_t frag_pos,
          size_t read_len,
          string &read1, string &read2, size_t &read2_pos) {

  read_len = std::min(read_len, the_frag.length());

  read1 = the_frag.substr(0, read_len);

  revcomp_inplace(the_frag);
  read2 = the_frag.substr(0, read_len);
  read2_pos = frag_pos + (the_frag.length() - read_len);
}

static void
sim_conversion(const bool pbat, string &read1, string &read2) {
  if (pbat) {
    replace(begin(read1), end(read1), 'G', 'A');
    replace(begin(read2), end(read2), 'C', 'T');
  }
  else {
    replace(begin(read1), end(read1), 'C', 'T');
    replace(begin(read2), end(read2), 'G', 'A');
  }
}


int main(int argc, const char **argv) {
  srand (time(NULL));
  try {
    string chrom_file;
    string outfile1, outfile2;
    bool VERBOSE = false;

    size_t n_reads = 100;
    size_t read_len = 100;
    size_t frag_len = 250;

    bool pbat = false;
    char strand_arg = 'f';

    size_t rng_seed = std::numeric_limits<size_t>::max();

    double mutation_rate = 0.0;
    size_t max_mutations = std::numeric_limits<size_t>::max();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate reads for "
                           "testing abismal", "<reference-genome-fasta>", 1);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out1", 'o', "output file end1", true, outfile1);
    opt_parse.add_opt("out2", 'p', "output file end2", true, outfile2);
    opt_parse.add_opt("len", 'l', "read length", false, read_len);
    opt_parse.add_opt("frag", 'f', "fragment length", false, frag_len);
    opt_parse.add_opt("n-reads", 'n', "number of reads", false, n_reads);
    opt_parse.add_opt("mut", 'm', "mutation rate", false, mutation_rate);
    opt_parse.add_opt("max-mut", 'M', "max mutations", false, max_mutations);
    opt_parse.add_opt("pbat", 'a', "pbat", false, pbat);
    opt_parse.add_opt("strand", 's', "strand {f, r, b}", false, strand_arg);
    opt_parse.add_opt("seed", '\0', "rng seed (default: from system)", false, rng_seed);
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
    if (leftover_args.size() != 1) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }
    const string genome_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      std::random_device rd;
      rng_seed = rd();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    std::mt19937 gen(rng_seed);
    std::uniform_real_distribution<double> unif(0.0, 1.0);
    function<double()> distr(bind(unif, std::ref(gen)));

    std::ifstream in(genome_file);
    if (!in)
      throw runtime_error("bad genome file: " + genome_file);

    if (VERBOSE)
      cerr << "[loading genome]" << endl;
    string genome;
    ChromLookup cl;
    load_genome(genome_file, genome, cl);

    std::transform(begin(genome), end(genome),
                   begin(genome), toupper);

    std::ofstream out1(outfile1);
    std::ofstream out2(outfile2);

    ProgressBar progress(n_reads);
    if (VERBOSE)
      progress.report(cerr, 0);
    for (size_t i = 0; i < n_reads; ++i) {

      if (VERBOSE && progress.time_to_report(i))
        progress.report(cerr, i);

      // simulate a fragment and a position
      string the_frag;
      size_t frag_pos = 0;
      sim_frag(genome, frag_len, the_frag, frag_pos);

      size_t n_mutations = 0;
      if (mutation_rate > 0.0)
        mutate_frag(mutation_rate, max_mutations, distr,
                    the_frag, n_mutations);

      // flip the strand of the fragment if appropriate
      char strand = '+';
      if (strand_arg == 'r' ||
          (strand_arg == 'b' && distr() > 0.5)) {
        revcomp_inplace(the_frag);
        strand = '-';
      }

      // extract the reads from the ends of the frag
      string read1, read2;
      size_t read2_pos = 0;
      sim_reads(the_frag, frag_pos, read_len,
                read1, read2, read2_pos);

      // simulate the conversion
      sim_conversion(false, read1, read2);

      size_t read1_pos = frag_pos;
      if (strand == '-')
        std::swap(read1_pos, read2_pos);

      uint32_t offset = 0, chrom_idx = 0;
      cl.get_chrom_idx_and_offset(read1_pos, chrom_idx, offset);

      out1 << "@read"
           << i << ":"
           << cl.names[chrom_idx] << ':'
           << offset << ":"
           << strand << ":"
           << n_mutations << endl
           << read1 << endl
           << "+" << endl
           << string(read_len, 'B') << endl;

      const char strand_two = (strand == '+' ? '-' : '+');
      cl.get_chrom_idx_and_offset(read2_pos, chrom_idx, offset);
      out2 << "@read" << i << ":"
           << cl.names[chrom_idx] << ':'
           << offset << ":"
           << strand_two << ":"
           << n_mutations << endl
           << read2 << endl
           << "+" << endl
           << string(read_len, 'B') << endl;
    }
    if (VERBOSE)
      progress.report(cerr, n_reads);
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

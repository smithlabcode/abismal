/* Copyright (C) 2018-2025 Andrew D. Smith
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

#include "simreads.hpp"

#include "AbismalIndex.hpp"
#include "OptionParser.hpp"
#include "cigar_utils.hpp"
#include "sam_record.hpp"
#include "smithlab_utils.hpp"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>  // IWYU pragma: keep
#include <vector>

#include <unistd.h>  // getpid()

template <typename T> using num_lim = std::numeric_limits<T>;

namespace simreads_random {
// ADS: I made this namespace and functions because different
// implementations of rand() on different OS meant that even with
// the same seed, the results could be different. This meant testing
// didn't work.

std::mt19937 e;                                   // NOLINT
bool initialized = false;                         // NOLINT
std::uniform_real_distribution<double> dr;        // NOLINT
std::uniform_int_distribution<std::uint64_t> di;  // NOLINT

void
initialize(const std::size_t the_seed) {
  e = std::mt19937(the_seed);
  initialized = true;
}

auto
rand() -> std::uint64_t {
  assert(initialized);
  // ADS: should have same range as ordinary rand() by properties of
  // std::uniform_int_distribution default constructor.
  return di(e);
}
auto
rand_double() -> double {  // ADS: in the interval [0, 1]
  assert(initialized);
  // ADS: default constructor for std::uniform_real_distribution
  // sets a range of [0,1)
  return dr(e);
}
}  // namespace simreads_random

static inline auto
format_fastq_record(const std::string &name,
                    const std::string &read) -> std::string {
  assert(!name.empty());
  std::string s;
  s += '@';
  s += name;
  s += '\n';
  s += read;
  s += "\n+\n";
  s += std::string(std::size(read), 'B');
  return s;
}

static inline auto
format_fasta_record(const std::string &name,
                    const std::string &read) -> std::string {
  assert(!name.empty());
  std::string s;
  s += '>';
  s += name;
  s += '\n';
  s += read;
  return s;
}

struct FragInfo {
  void
  set_sequential_name() {
    name = "read" + std::to_string(frag_count++);
  }

  [[nodiscard]] auto
  read1() const -> std::string {
    assert(!name.empty());
    std::string read = seq.substr(0, read_length);
    for (std::size_t i = 0; i < read_length - std::size(read); ++i)
      read += random_base();
    return fasta_format ? format_fasta_record(name + ".1", read)
                        : format_fastq_record(name + ".1", read);
  }

  [[nodiscard]] auto
  read2() const -> std::string {
    assert(!name.empty());
    std::string read(seq);
    revcomp_inplace(read);
    read = read.substr(0, read_length);  // cppcheck-suppress uselessCallsSubstr
    for (std::size_t i = 0; i < read_length - std::size(read); ++i)
      read += random_base();
    return fasta_format ? format_fasta_record(name + ".2", read)
                        : format_fastq_record(name + ".2", read);
  }

  void
  erase_info_through_insert() {
    const std::size_t orig_ref_len = end_pos - start_pos;
    if (2 * read_length < std::size(seq)) {
      std::string cigar2(cigar);
      truncate_cigar_q(cigar, read_length);
      reverse_cigar(begin(cigar2), end(cigar2));
      truncate_cigar_q(cigar2, read_length);
      reverse_cigar(begin(cigar2), end(cigar2));
      const std::size_t rseq_ops =
        cigar_rseq_ops(cigar) + cigar_rseq_ops(cigar2);
      cigar = cigar + std::to_string(orig_ref_len - rseq_ops) + "N" + cigar2;
      seq =
        seq.substr(0, read_length) +  // cppcheck-suppress uselessCallsSubstr
        std::string(orig_ref_len - rseq_ops, 'N') +
        seq.substr(std::size(seq) - read_length, read_length);
    }
  }

  void
  remove_cigar_match_symbols() {
    replace(begin(cigar), end(cigar), '=', 'M');
    merge_equal_neighbor_cigar_ops(cigar);
  }

  void
  bisulfite_conversion(const bool random_pbat, const double bs_conv) {
    constexpr auto coin_flip = 0.5;
    if (pbat || (random_pbat && simreads_random::rand_double() < coin_flip)) {
      for (char &it : seq) {
        if (it == 'G' && (simreads_random::rand_double() < bs_conv))
          it = 'A';  // cppcheck-suppress useStlAlgorithm
      }
    }
    else {
      for (char &it : seq) {
        if (it == 'C' && (simreads_random::rand_double() < bs_conv))
          it = 'T';  // cppcheck-suppress useStlAlgorithm
      }
    }
  }

  [[nodiscard]] auto
  rc() const -> bool {
    return strand == '-';
  }

  std::string chrom;
  std::size_t start_pos{};
  std::size_t end_pos{};
  std::string name;
  double score{};
  char strand{};
  std::string seq;
  std::string cigar;

  static bool pbat;
  static bool fasta_format;
  static std::size_t frag_count;
  static std::size_t read_length;
};

static constexpr auto read_length_default = 100;
bool FragInfo::pbat = false;
bool FragInfo::fasta_format = false;
std::size_t FragInfo::frag_count = 0;
std::size_t FragInfo::read_length = read_length_default;

static auto
operator<<(std::ostream &out, const FragInfo &the_info) -> std::ostream & {
  const bool rc = the_info.rc();
  std::uint16_t flags_read = 0;
  std::uint16_t flags_mate = 0;

  samflags::set(flags_read, samflags::read_paired);
  samflags::set(flags_read, samflags::read_pair_mapped);
  samflags::set(flags_read, samflags::template_first);
  samflags::set(flags_read,
                the_info.rc() ? samflags::read_rc : samflags::mate_rc);

  samflags::set(flags_mate, samflags::read_paired);
  samflags::set(flags_mate, samflags::read_pair_mapped);
  samflags::set(flags_mate, samflags::template_last);
  samflags::set(flags_mate,
                the_info.rc() ? samflags::mate_rc : samflags::read_rc);

  const std::size_t read_pos = the_info.start_pos + 1;
  const std::size_t mate_pos = the_info.end_pos - FragInfo::read_length + 1;

  const int tlen = rc ? -static_cast<int>(std::size(the_info.seq))
                      : static_cast<int>(std::size(the_info.seq));

  std::string cigar1 = the_info.cigar;
  std::string cigar2 = the_info.cigar;

  truncate_cigar_q(cigar1, FragInfo::read_length);
  reverse_cigar(cigar2);
  truncate_cigar_q(cigar2, FragInfo::read_length);

  if (rc) {
    reverse_cigar(cigar1);
  }
  else
    reverse_cigar(cigar2);

  const std::string seq1 = the_info.seq.substr(0, FragInfo::read_length);
  const std::string read_rc = revcomp(the_info.seq);
  const std::string seq2 = read_rc.substr(0, FragInfo::read_length);
  const std::size_t pos1 = rc ? mate_pos : read_pos;
  const std::size_t pos2 = rc ? read_pos : mate_pos;

  // clang-format off
  return out << the_info.name << ".1" << '\t'
             << flags_read << '\t'
             << the_info.chrom << '\t'
             << pos1 << '\t'
             << "255" << '\t'
             << cigar1 << '\t'
             << "=" << '\t'
             << pos2 << '\t'
             << tlen << '\t'
             << seq1 << '\t'
             << "*" << '\n'
             << the_info.name << ".2" << '\t'
             << flags_mate << '\t'
             << the_info.chrom << '\t'
             << pos2 << '\t'
             << "255" << '\t'
             << cigar2 << '\t'
             << "=" << '\t'
             << pos1 << '\t'
             << -tlen << '\t'
             << seq2 << '\t'
             << "*";
  // clang-format on
}

// extract the position of the fragment checking all bases are valid
static void
sim_frag_position(const std::string &genome, const std::size_t frag_len,
                  std::string &the_frag, std::size_t &the_posn,
                  const bool require_valid) {
  static auto is_valid = [](const char c) { return valid_base(c); };

  const auto g_beg = std::cbegin(genome);

  const std::size_t lim = std::size(genome) - frag_len + 1;
  the_posn = simreads_random::rand() % lim;
  // NOLINTBEGIN(*-narrowing-conversions)
  the_frag = std::string(g_beg + the_posn, g_beg + the_posn + frag_len);
  while (require_valid &&
         std::all_of(std::cbegin(the_frag), std::cend(the_frag), is_valid)) {
    the_posn = simreads_random::rand() % lim;
    the_frag = std::string(g_beg + the_posn, g_beg + the_posn + frag_len);
  }
  // NOLINTEND(*-narrowing-conversions)
}

// simulate from a uniform distribution in a range
static auto
sim_frag_length(const std::size_t min_length,
                const std::size_t max_length) -> std::size_t {
  assert(max_length >= min_length);
  if (min_length == max_length)
    return min_length;
  const std::size_t diff = max_length - min_length;
  return min_length + (simreads_random::rand() % diff);
}

struct FragSampler {
  FragSampler(const std::string &g, ChromLookup c, const char sc,
              const std::size_t milen, const std::size_t malen,
              const bool require_valid) :
    genome(g), cl(std::move(c)), strand_code(sc), min_length(milen),
    max_length(malen), require_valid(require_valid) {}
  void
  sample_fragment(FragInfo &the_info) const {
    const std::size_t frag_len = sim_frag_length(min_length, max_length);
    sim_frag_position(genome, frag_len, the_info.seq, the_info.start_pos,
                      require_valid);

    uint32_t offset = 0;
    std::int32_t chrom_idx = 0;
    cl.get_chrom_idx_and_offset(the_info.start_pos, chrom_idx, offset);
    the_info.chrom = cl.names[chrom_idx];
    the_info.start_pos = offset;

    the_info.end_pos = the_info.start_pos + frag_len;
    the_info.set_sequential_name();  // default
    the_info.strand = sim_strand();  // based on frag code
    if (the_info.strand == '-')
      revcomp_inplace(the_info.seq);
    the_info.cigar = std::to_string(frag_len) + "M";  // default, no muts
  }
  [[nodiscard]] auto
  sim_strand() const -> char {
    switch (strand_code) {
    case 'f':
      return '+';
    case 'r':
      return '-';
    case 'b':
      return (simreads_random::rand() & 1) ? '+' : '-';
    default:
      std::abort();
    }
  }
  const std::string &genome;  // NOLINT(*-avoid-const-*-data-members)
  ChromLookup cl;
  char strand_code{};
  std::size_t min_length{};
  std::size_t max_length{};
  bool require_valid{};
};

struct FragMutator {
  FragMutator(const double m, const double s, const double i, const double d) :
    mutation_rate(m), substitution_rate(s), insertion_rate(i),
    deletion_rate(d) {
    const double total =
      std::max(substitution_rate + insertion_rate + deletion_rate,
               num_lim<double>::min());
    substitution_rate /= total;
    insertion_rate /= total;
    deletion_rate /= total;
    insertion_rate += substitution_rate;
    deletion_rate += insertion_rate;
  }
  void
  mutate(FragInfo &the_info) const {
    std::string seq, cigar;
    std::size_t i = 0;
    the_info.score = 0;
    while (i < std::size(the_info.seq)) {
      // select a mutation or not
      const char mut = sample_mutation();
      if (mut == 'I') {
        cigar += "I";
        seq += random_base();
        ++the_info.score;
      }
      else if (mut == 'D') {
        cigar += "D";
        ++i;
        ++the_info.score;
      }
      else if (mut == 'M') {
        cigar += "M";
        seq += random_base();
        ++the_info.score;
        ++i;
      }
      else {  // if (mut == '=') {
        cigar += "=";
        seq += the_info.seq[i];
        ++i;
      }
    }
    the_info.cigar.resize(2 * std::size(cigar));
    compress_cigar(std::cbegin(cigar), std::cend(cigar), the_info.cigar);
    std::swap(seq, the_info.seq);
  }
  [[nodiscard]] auto
  sample_mutation() const -> char {
    const double x = simreads_random::rand_double();
    if (x > mutation_rate)
      return '=';
    else {
      const double y = simreads_random::rand_double();
      if (y < substitution_rate)
        return 'M';
      else if (y < insertion_rate)
        return 'I';
      else
        return 'D';
    }
  }
  [[nodiscard]] auto
  tostring() const -> std::string {
    std::ostringstream oss;
    oss << "mutation_rate=" << mutation_rate << '\n'
        << "substitution_rate=" << substitution_rate << '\n'
        << "insertion_rate=" << insertion_rate << '\n'
        << "deletion_rate=" << deletion_rate;
    return oss.str();
  }
  double mutation_rate{};
  double substitution_rate{};
  double insertion_rate{};
  double deletion_rate{};
};

static void
extract_change_type_vals(const std::string &change_type_vals,
                         double &substitution_rate, double &insertion_rate,
                         double &deletion_rate) {
  if (!change_type_vals.empty()) {
    std::istringstream iss(change_type_vals);
    char x{};
    iss >> substitution_rate;
    iss >> x;
    iss >> insertion_rate;
    iss >> x;
    iss >> deletion_rate;
  }
}

int
simreads(int argc, char *argv[]) {  // NOLINT(*-c-arrays)
  static constexpr auto n_reads_default = 100;
  static constexpr auto min_frag_len_default = 100;
  static constexpr auto max_frag_len_defeault = 250;
  try {
    std::string output_prefix;
    std::string locations_file;

    bool VERBOSE = false;
    bool single_end = false;
    bool show_cigar_matches = true;
    bool random_pbat = false;
    bool require_valid = false;

    std::size_t n_reads{n_reads_default};
    std::size_t min_frag_len{min_frag_len_default};
    std::size_t max_frag_len{max_frag_len_defeault};

    char strand_arg = 'b';

    std::size_t rng_seed = num_lim<std::size_t>::max();

    double mutation_rate = 0.0;
    std::string change_type_vals;
    double substitution_rate = 1.0;
    double insertion_rate = 1.0;
    double deletion_rate = 1.0;

    double bs_conv = 1.0;

    std::size_t max_mutations = num_lim<std::size_t>::max();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],  // NOLINT(*-pointer-arithmetic)
                           "simulate reads for "
                           "testing walt2",
                           "<reference-genome-fasta>", 1);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out", 'o', "output file prefix", true, output_prefix);
    opt_parse.add_opt("single", '\0', "output single end", false, single_end);
    opt_parse.add_opt("loc", '\0', "write locations here", false,
                      locations_file);
    opt_parse.add_opt("read-len", 'l', "read length", false,
                      FragInfo::read_length);
    opt_parse.add_opt("min-fraglen", '\0', "min fragment length", false,
                      min_frag_len);
    opt_parse.add_opt("max-fraglen", '\0', "max fragment length", false,
                      max_frag_len);
    opt_parse.add_opt("n-reads", 'n', "number of reads", false, n_reads);
    opt_parse.add_opt("mut", 'm', "mutation rate", false, mutation_rate);
    opt_parse.add_opt("bis", 'b', "bisulfite conversion rate", false, bs_conv);
    opt_parse.add_opt("show-matches", '\0', "show match symbols in cigar",
                      false, show_cigar_matches);
    opt_parse.add_opt("changes", 'c', "change types (comma sep relative vals)",
                      false, change_type_vals);
    opt_parse.add_opt("max-mut", 'M', "max mutations", false, max_mutations);
    opt_parse.add_opt("pbat", 'a', "pbat", false, FragInfo::pbat);
    opt_parse.add_opt("random-pbat", 'R', "random pbat", false, random_pbat);
    opt_parse.add_opt("strand", 's', "strand {f, r, b}", false, strand_arg);
    opt_parse.add_opt("fasta", '\0', "output fasta format (no quality scores)",
                      false, FragInfo::fasta_format);
    opt_parse.add_opt("seed", '\0', "rng seed (default: from system)", false,
                      rng_seed);
    opt_parse.add_opt("require-valid", '\0', "require valid bases in fragments",
                      false, require_valid);
    opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);
    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (std::size(leftover_args) != 1) {
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }
    const std::string genome_file(leftover_args.front());
    /****************** END COMMAND LINE OPTIONS *****************/

    extract_change_type_vals(change_type_vals, substitution_rate,
                             insertion_rate, deletion_rate);

    if (rng_seed == num_lim<std::size_t>::max())
      rng_seed = time(nullptr) + getpid();

    if (VERBOSE)
      std::cerr << "rng seed: " << rng_seed << '\n';
    simreads_random::initialize(rng_seed);

    if (VERBOSE)
      std::cerr << "[loading genome]\n";
    std::ifstream in(genome_file);
    if (!in)
      throw std::runtime_error("bad genome file: " + genome_file);
    std::string genome;
    ChromLookup cl;
    load_genome(genome_file, genome, cl);
    std::transform(std::cbegin(genome), std::cend(genome), std::begin(genome),
                   [](unsigned char c) { return toupper(c); });

    std::ofstream loc_out;
    if (!locations_file.empty()) {
      if (VERBOSE)
        std::cerr << "[opening frag locations file: " << locations_file << "]"
                  << '\n';
      loc_out.open(locations_file);
      if (!loc_out)
        throw std::runtime_error("bad locations output file: " +
                                 locations_file);
    }

    const std::string read1_outfile =
      output_prefix + (FragInfo::fasta_format ? "_1.fa" : "_1.fq");
    if (VERBOSE) {
      if (FragInfo::fasta_format)
        std::cerr << "[opening read1 fastq: " << read1_outfile << "]\n";
      else
        std::cerr << "[opening read1 fasta: " << read1_outfile << "]\n";
    }
    std::ofstream read1_out(read1_outfile);
    if (!read1_out)
      throw std::runtime_error("bad output file: " + read1_outfile);

    std::ofstream read2_out;
    if (!single_end) {
      const std::string read2_outfile =
        output_prefix + (FragInfo::fasta_format ? "_2.fa" : "_2.fq");
      if (VERBOSE) {
        if (FragInfo::fasta_format)
          std::cerr << "[opening read2 fastq: " << read2_outfile << "]\n";
        else
          std::cerr << "[opening read2 fasta: " << read2_outfile << "]\n";
      }
      read2_out.open(read2_outfile);
      if (!read2_out)
        throw std::runtime_error("bad output file: " + read2_outfile);
    }

    FragSampler frag_samp(genome, cl, strand_arg, min_frag_len, max_frag_len,
                          require_valid);
    if (VERBOSE)
      std::cerr << "[constructed fragment sampler]\n";

    FragMutator frag_mut(mutation_rate, substitution_rate, insertion_rate,
                         deletion_rate);
    if (VERBOSE)
      std::cerr << "[constructed mutator]\n"
                << "[simulating frags]\n";

    for (std::size_t i = 0; i < n_reads; ++i) {
      FragInfo info;
      frag_samp.sample_fragment(info);
      frag_mut.mutate(info);
      info.bisulfite_conversion(random_pbat, bs_conv);
      if (!show_cigar_matches)
        info.remove_cigar_match_symbols();
      if (!locations_file.empty())
        loc_out << info << '\n';
      read1_out << info.read1() << '\n';
      if (!single_end)
        read2_out << info.read2() << '\n';
    }
  }
  catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}

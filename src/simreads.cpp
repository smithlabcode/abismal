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
#include "cigar_utils.hpp"
#include "sam_record.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <random>
#include <functional>
#include <algorithm>

#include <unistd.h>

using std::vector;
using std::runtime_error;
using std::string;
using std::cerr;
using std::endl;
using std::function;
using std::istream;
using std::ostream;
using std::ifstream;
using std::ofstream;
using std::ostringstream;
using std::istringstream;
using std::to_string;

double rand_double() { return ((double) rand() / (RAND_MAX)); }

static string
format_fastq_record(const string &name, const string &read) {
  assert(!name.empty());
  ostringstream oss;
  oss << '@' << name << endl << read << endl
      << '+' << endl << string(read.length(), 'B');
  return oss.str();
}


struct FragInfo {
  void set_sequential_name() {
    name = "read" + to_string(frag_count++);
  }
  string
  read1() const {
    assert(!name.empty());
    string read = seq.substr(0, read_length);
    for (size_t i = 0; i < read_length - read.length(); ++i)
      read += random_base();
    return format_fastq_record(name + ".1", read);
  }
  string
  read2() const {
    assert(!name.empty());
    string read(seq);
    revcomp_inplace(read);
    read = read.substr(0, read_length);
    for (size_t i = 0; i < read_length - read.length(); ++i)
      read += random_base();
    return format_fastq_record(name + ".2", read);
  }
  void
  erase_info_through_insert() {
    const size_t orig_ref_len = end_pos - start_pos;
    if (2*read_length < seq.length()) {
      string cigar2(cigar);
      truncate_cigar_q(cigar, read_length);
      reverse_cigar(begin(cigar2), end(cigar2));
      truncate_cigar_q(cigar2, read_length);
      reverse_cigar(begin(cigar2), end(cigar2));
      const size_t rseq_ops = cigar_rseq_ops(cigar) + cigar_rseq_ops(cigar2);
      cigar = cigar + to_string(orig_ref_len - rseq_ops) + "N" + cigar2;
      seq = seq.substr(0, read_length) +
        string(orig_ref_len - rseq_ops, 'N') +
        seq.substr(seq.length() - read_length, read_length);

    }
  }
  void
  remove_cigar_match_symbols() {
    replace(begin(cigar), end(cigar), '=', 'M');
    merge_equal_neighbor_cigar_ops(cigar);
  }
  void
  bisulfite_conversion(const bool random_pbat) {
    if (pbat || (random_pbat && rand_double() < 0.5))
      replace(begin(seq), end(seq), 'G', 'A');
    else replace(begin(seq), end(seq), 'C', 'T');
  }

  bool rc() const { return strand == '-'; }

  string chrom;
  size_t start_pos;
  size_t end_pos;
  string name;
  double score;
  char strand;
  string seq;
  string cigar;

  static bool pbat;
  static size_t frag_count;
  static size_t read_length;
};

bool FragInfo::pbat = false;
size_t FragInfo::frag_count = 0;
size_t FragInfo::read_length = 100;


ostream &
operator<<(ostream &out, FragInfo &the_info) {
  const bool rc = the_info.rc();
  uint16_t flags_read = 0;
  uint16_t flags_mate = 0;

  samflags::set(flags_read, samflags::read_paired);
  samflags::set(flags_read, samflags::read_pair_mapped);
  samflags::set(flags_read, samflags::template_first);
  samflags::set(flags_read, the_info.rc() ? samflags::read_rc : samflags::mate_rc);

  samflags::set(flags_mate, samflags::read_paired);
  samflags::set(flags_mate, samflags::read_pair_mapped);
  samflags::set(flags_mate, samflags::template_last);
  samflags::set(flags_mate, the_info.rc() ? samflags::mate_rc : samflags::read_rc);

  const size_t read_pos = the_info.start_pos + 1;
  const size_t mate_pos = the_info.end_pos - FragInfo::read_length + 1;
  const int tlen = rc ? (-the_info.seq.size()) : (the_info.seq.size());
  string cigar1 = the_info.cigar;
  string cigar2 = the_info.cigar;

  truncate_cigar_q(cigar1, FragInfo::read_length);
  reverse_cigar(cigar2);
  truncate_cigar_q(cigar2, FragInfo::read_length);

  if (rc) {
    reverse_cigar(cigar1);
  }
  else
    reverse_cigar(cigar2);

  const string seq1 = the_info.seq.substr(0, FragInfo::read_length);
  const string read_rc = revcomp(the_info.seq);
  const string seq2 = read_rc.substr(0, FragInfo::read_length);
  const size_t pos1 = rc ? mate_pos : read_pos;
  const size_t pos2 = rc ? read_pos : mate_pos;

  return out << the_info.name << ".1\t"
             << flags_read << '\t'
             << the_info.chrom << '\t'
             << pos1 << '\t'
             << "255\t"
             << cigar1 << '\t'
             << "=\t"
             << pos2 << "\t"
             << tlen << '\t'
             << seq1 << "\t"
             << "*" << endl

             << the_info.name << ".2\t"
             << flags_mate << '\t'
             << the_info.chrom << '\t'
             << pos2 << '\t'
             << "255\t"
             << cigar2 << '\t'
             << "=\t"
             << pos1 << "\t"
             << -tlen << '\t'
             << seq2 << "\t"
             << "*";
}

istream &
operator>>(istream &in, FragInfo &the_info) {
  string line;
  getline(in, line);
  istringstream iss;
  iss.rdbuf()->pubsetbuf(&line[0], line.length());
  return iss >> the_info.chrom >> the_info.start_pos >> the_info.end_pos
             >> the_info.name >> the_info.score >> the_info.strand
             >> the_info.seq >> the_info.cigar;
}


// extract the position of the fragment checking all bases are valid
static void
sim_frag_position(const string &genome, const size_t frag_len,
                  string &the_frag, size_t &the_position) {
  static auto is_invalid = [](const char c) {return !valid_base(c);};

  const size_t lim = genome.length() - frag_len + 1;
  do {
    the_position = rand() % lim;
    the_frag = string(begin(genome) + the_position,
                      begin(genome) + the_position + frag_len);
  }
  while (find_if(begin(the_frag), end(the_frag), is_invalid) != end(the_frag));
}


// simulate from a uniform distribution in a range
static size_t
sim_frag_length(const size_t min_length, const size_t max_length) {
  assert(max_length >= min_length);
  if (min_length == max_length) return min_length;
  const size_t diff = max_length - min_length;
  return min_length + (rand() % diff);
}


struct FragSampler {
  FragSampler(const string &g, const ChromLookup c, const char sc,
              const size_t milen, const size_t malen) :
    genome(g), cl(c), strand_code(sc), min_length(milen), max_length(malen) {}
  void
  sample_fragment(FragInfo &the_info, const bool random_pbat) const {
    const size_t frag_len = sim_frag_length(min_length, max_length);
    sim_frag_position(genome, frag_len, the_info.seq, the_info.start_pos);

    uint32_t offset = 0, chrom_idx = 0;
    cl.get_chrom_idx_and_offset(the_info.start_pos, chrom_idx, offset);
    the_info.chrom = cl.names[chrom_idx];
    the_info.start_pos = offset;

    the_info.end_pos = the_info.start_pos + frag_len;
    the_info.set_sequential_name(); // default
    the_info.strand = sim_strand(); // based on frag code
    if (the_info.strand == '-')
      revcomp_inplace(the_info.seq);
    the_info.cigar = to_string(frag_len) + "M"; // default, no muts

    the_info.bisulfite_conversion(random_pbat);
  }
  char sim_strand() const {
    if (strand_code == 'f') return '+';
    else if (strand_code == 'r') return '-';
    else if (strand_code == 'b') return (rand() & 1) ? '+' : '-';
    else throw runtime_error("bad strand code: " + to_string(strand_code));
    return '\0';
  }
  const string &genome;
  ChromLookup cl;
  char strand_code;
  size_t min_length;
  size_t max_length;
};


ostream &
operator<<(ostream &out, FragSampler &the_samp) {
  return out << "strand_code=" << the_samp.strand_code << '\n'
             << "min_length=" << the_samp.min_length << '\n'
             << "max_length=" << the_samp.max_length << '\n'
             << "cl=" << endl
             << the_samp.cl.tostring();
}


struct FragMutator {
  FragMutator(const double m, const double s, const double i, const double d) :
    mutation_rate(m), substitution_rate(s),
    insertion_rate(i), deletion_rate(d) {
    const double total = substitution_rate + insertion_rate + deletion_rate;
    substitution_rate /= total;
    insertion_rate /= total;
    deletion_rate /= total;
    insertion_rate += substitution_rate;
    deletion_rate += insertion_rate;
  }
  void mutate(FragInfo &the_info) const {
    string seq, cigar;
    size_t i = 0;
    the_info.score = 0;
    while (i < the_info.seq.length()) {
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
      else { //if (mut == '=') {
        cigar += "=";
        seq += the_info.seq[i];
        ++i;
      }
    }
    the_info.cigar.resize(2*cigar.size());
    compress_cigar(begin(cigar), end(cigar), the_info.cigar);
    swap(seq, the_info.seq);
  }
  char sample_mutation() const {
    double x  = (double)rand()/(double)RAND_MAX;
    if (x > mutation_rate) return '=';
    else {
      double y  = (double)rand()/(double)RAND_MAX;
      if (y < substitution_rate) return 'M';
      else if (y < insertion_rate) return 'I';
      else return 'D';
    }
  }
  string tostring() const {
    ostringstream oss;
    oss << "mutation_rate=" << mutation_rate << endl
        << "substitution_rate=" << substitution_rate << endl
        << "insertion_rate=" << insertion_rate << endl
        << "deletion_rate=" << deletion_rate;
    return oss.str();
  }
  double mutation_rate;
  double substitution_rate;
  double insertion_rate;
  double deletion_rate;
};


ostream &
operator<<(ostream &out, FragMutator &the_mut) {
  return out << the_mut.tostring();
}


static void
extract_change_type_vals(const string &change_type_vals,
                         double &substitution_rate,
                         double &insertion_rate,
                         double &deletion_rate) {
  if (!change_type_vals.empty()) {
    istringstream iss(change_type_vals);
    char x;
    iss >> substitution_rate;
    iss >> x;
    iss >> insertion_rate;
    iss >> x;
    iss >> deletion_rate;
  }
}


int main(int argc, const char **argv) {

  try {
    string chrom_file;
    string output_prefix;
    string locations_file;

    bool VERBOSE = false;
    bool write_locations = false;
    bool single_end = false;
    bool show_cigar_matches = false;
    bool random_pbat = false;

    size_t n_reads = 100;
    size_t min_frag_len = 100;
    size_t max_frag_len = 250;

    char strand_arg = 'b';

    size_t rng_seed = std::numeric_limits<size_t>::max();

    double mutation_rate = 0.0;
    string change_type_vals;
    double substitution_rate = 1.0;
    double insertion_rate = 1.0;
    double deletion_rate = 1.0;

    size_t max_mutations = std::numeric_limits<size_t>::max();

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(strip_path(argv[0]), "simulate reads for "
                           "testing walt2", "<reference-genome-fasta>", 1);
    opt_parse.set_show_defaults();
    opt_parse.add_opt("out", 'o', "output file prefix", true, output_prefix);
    opt_parse.add_opt("single", '\0', "output single end", false, single_end);
    opt_parse.add_opt("loc", '\0', "write locations", false, write_locations);
    opt_parse.add_opt("read-len", 'l', "read length", false, FragInfo::read_length);
    opt_parse.add_opt("min-fraglen", '\0', "min fragment length",
                      false, min_frag_len);
    opt_parse.add_opt("max-fraglen", '\0', "max fragment length",
                      false, max_frag_len);
    opt_parse.add_opt("n-reads", 'n', "number of reads", false, n_reads);
    opt_parse.add_opt("mut", 'm', "mutation rate", false, mutation_rate);
    opt_parse.add_opt("show-matches", '\0', "show match symbols in cigar",
                      false, show_cigar_matches);
    opt_parse.add_opt("changes", 'c',
                      "change types (comma sep relative vals)",
                      false, change_type_vals);
    opt_parse.add_opt("max-mut", 'M', "max mutations", false, max_mutations);
    opt_parse.add_opt("pbat", 'a', "pbat", false, FragInfo::pbat);
    opt_parse.add_opt("random-pbat", 'R', "random pbat", false, random_pbat);
    opt_parse.add_opt("strand", 's', "strand {f, r, b}", false, strand_arg);
    opt_parse.add_opt("seed", '\0', "rng seed (default: from system)",
                      false, rng_seed);
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

    extract_change_type_vals(change_type_vals,
                             substitution_rate, insertion_rate, deletion_rate);

    /* standard mersenne_twister_engine seeded with rd()*/
    if (rng_seed == std::numeric_limits<size_t>::max()) {
      rng_seed = time(0) + getpid();
    }
    if (VERBOSE)
      cerr << "rng seed: " << rng_seed << endl;
    srand(rng_seed);


    if (VERBOSE)
      cerr << "[loading genome]" << endl;
    std::ifstream in(genome_file);
    if (!in)
      throw runtime_error("bad genome file: " + genome_file);
    string genome;
    ChromLookup cl;
    load_genome(genome_file, genome, cl);
    transform(begin(genome), end(genome), begin(genome),
              [](unsigned char c){return toupper(c);});


    FragSampler frag_samp(genome, cl, strand_arg, min_frag_len, max_frag_len);
    if (VERBOSE)
      cerr << "[constructed fragment sampler]" << endl
           << frag_samp << endl;
    if (VERBOSE)
      cerr << "[simulating clean frags]" << endl;
    vector<FragInfo> the_info(n_reads);
    for (size_t i = 0; i < n_reads; ++i)
      frag_samp.sample_fragment(the_info[i], random_pbat);


    if (VERBOSE)
      cerr << "[mutating the frags]" << endl;
    FragMutator frag_mut(mutation_rate, substitution_rate,
                         insertion_rate, deletion_rate);
    if (VERBOSE)
      cerr << frag_mut << endl;
    for (size_t i = 0; i < the_info.size(); ++i)
      frag_mut.mutate(the_info[i]);


    if (!show_cigar_matches)
      for (size_t i = 0; i < the_info.size(); ++i)
        the_info[i].remove_cigar_match_symbols();


    if (write_locations) {
      const string locations_file = output_prefix + ".sam";
      if (VERBOSE)
        cerr << "[writing frag locations: " << locations_file << "]" << endl;
      ofstream loc_out(locations_file);
      if (!loc_out)
        throw runtime_error("bad locations file: " + locations_file);
      for (size_t i = 0; i < the_info.size(); ++i) {
        //the_info[i].erase_info_through_insert();
        loc_out << the_info[i] << endl;
      }
    }

    const string read1_outfile = output_prefix + "_1.fq";
    if (VERBOSE)
      cerr << "[writing read1 fastq: " << read1_outfile << "]" << endl;
    ofstream read1_out(read1_outfile);
    if (!read1_out)
      throw runtime_error("bad output file: " + read1_outfile);
    for (size_t i = 0; i < the_info.size(); ++i)
      read1_out << the_info[i].read1() << endl;


    if (!single_end) {
      const string read2_outfile = output_prefix + "_2.fq";
      if (VERBOSE)
        cerr << "[writing read2 fastq: " << read2_outfile << "]" << endl;
      ofstream read2_out(read2_outfile);
      if (!read2_out)
        throw runtime_error("bad output file: " + read2_outfile);
      for (size_t i = 0; i < the_info.size(); ++i)
        read2_out << the_info[i].read2() << endl;
    }
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

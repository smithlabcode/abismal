#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "smithlab_os.cpp"
#include "smithlab_utils.cpp"
#include "OptionParser.cpp"
#include "cigar_utils.hpp"
#include "dna_four_bit.hpp"

using std::cerr;
using std::cout;
using std::ifstream;
using std::ofstream;
using std::runtime_error;
using std::endl;
using std::string;
using std::vector;

static void
apply_cigar(string &seq,
            const string &cigar) {
  string new_seq;
  size_t n;
  char op;
  size_t i = 0;
  std::istringstream iss(cigar);
  while (iss >> n >> op) {
    switch (op)
      {
      case 'M':
        new_seq += seq.substr(i, n);
        i += n;
        break;
      case 'I':
        i += n;
        break;
      case 'D':
        new_seq += string(n, 'N');
        break;
      case 'S':
        i += n;
        break;
      case 'N':
        new_seq += string(n, 'N');
        i += n;
        break;
      }
  }
  // Sum of lengths of the M/I/S/=/X/N operations
  // shall equal the length of seq.
  if (i != seq.length())
    throw runtime_error("inconsistent number of qseq ops! " +
                        seq + " "  + cigar + " " + std::to_string(i) + " " +
                        std::to_string(seq.length()));

  seq = new_seq;
}
int
main(int argc, const char **argv) {
  try {
    bool verbose;
    OptionParser opt_parse(strip_path(argv[0]),
                         "convert abismal mr to regular mr",
                         "<input-mr> <output-mr>");

    opt_parse.set_show_defaults();
    opt_parse.add_opt("verbose", 'v', "print more run info", false, verbose);
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

    const string infile = leftover_args.front(),
                 outfile = leftover_args.back();
    ifstream in(infile);

    if (!in) throw runtime_error("bad in mr file: " + infile);
    ofstream out(outfile);

    string chrom, name, cigar, seq;
    size_t start, end, mapq;
    char strand;
    while(in >> chrom >> start >> end >> name >> mapq
             >> strand >> seq >> cigar) {
      apply_cigar(seq, cigar);
      if (seq.size() != end - start) {
        cerr << "inconsistent number of rseq ops! "  << name << "\n";
        cerr << "got " << std::to_string(seq.size()) << ", expected "
             << std::to_string(end - start) << "\n";
      }
      out << chrom << '\t' << start << '\t' << name << '\t' << mapq << '\t'
          << strand << '\t' << seq << '\n';
    }
    in.close();
    out.close();
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

/* GS: this is an adaptation of the original std::lower_bound
 * implementation that predicts the future cache load of binary
 * search */
template<class ForwardIt, class T, class Compare>
ForwardIt my_lower_bound(ForwardIt first, const ForwardIt last,
                         const T& value, Compare comp) {
  ForwardIt it;
  typename std::iterator_traits<ForwardIt>::difference_type count, step;
  count = std::distance(first, last);

  while (count > 0) {
    it = first;
    step = (count >> 1);

    // GS prefetches both future iterations of the binary search with
    // moderate locality
    __builtin_prefetch(&(*(first + ((count + step) >> 1))), 0, 1);
    __builtin_prefetch(&(*(first + ((step - 1) >> 1))), 0, 1);
    std::advance(it, step);
    if (comp(*it, value)) {
      first = ++it;
      count -= step + 1;
    }
    else
			count = step;
  }
  return first;
}


bool
format_pe(const bool allow_ambig,
          const pe_result &res, const ChromLookup &cl,
          const string &read1, const string &read2,
          const string &name1, const string &name2,
          const string &cig1,  const string &cig2,
          ostream &out) {

  uint32_t r_s1 = 0, r_e1 = 0, chr1 = 0; // positions in chroms (0-based)
  uint32_t r_s2 = 0, r_e2 = 0, chr2 = 0;
  const pe_element p = res.best;

  // does not depend on allow_ambig
  if (!res.should_report(read1, read2))
    return false;

  // PE chromosomes differ or couldn't be found, treat read as unmapped
  if (!chrom_and_posn(cl, cig1, p.r1.pos, r_s1, r_e1, chr1) ||
      !chrom_and_posn(cl, cig2, p.r2.pos, r_s2, r_e2, chr2) || chr1 != chr2)
    return false;

  const bool rc = p.rc();

  // ADS: will this always evaluate correctly with unsigned
  // intermediate vals?
  const int tlen = rc ?
    (static_cast<int>(r_s1) - static_cast<int>(r_e2)) :
    (static_cast<int>(r_e2) - static_cast<int>(r_s1));

  // ADS: +1 to POS & PNEXT; "=" for RNEXT; "*" for QUAL
  sam_rec sr1(name1, 0, cl.names[chr1], r_s1 + 1, 255,
              cig1, "=", r_s2 + 1, tlen, read1, "*");

  sr1.add_tag("NM:i:" + to_string(p.r1.diffs));
  sr1.add_tag(p.r1.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
  set_flag(sr1, samflags::read_paired);
  set_flag(sr1, samflags::read_pair_mapped);
  set_flag(sr1, samflags::template_first);

  sam_rec sr2(name2, 0, cl.names[chr2], r_s2 + 1, 255,
              cig2, "=", r_s1 + 1, -tlen, read2, "*");

  sr2.add_tag("NM:i:" + to_string(p.r2.diffs));
  sr2.add_tag(p.r2.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
  set_flag(sr2, samflags::read_paired);
  set_flag(sr2, samflags::read_pair_mapped);
  set_flag(sr2, samflags::template_last);

  // second mate is reverse strand and richness of 1st mate
  if (rc) {
    set_flag(sr1, samflags::read_rc);
    set_flag(sr2, samflags::mate_rc);
  }
  else {
    set_flag(sr1, samflags::mate_rc);
    set_flag(sr2, samflags::read_rc);
  }

  out << sr1 << '\n'
      << sr2 << '\n';

  return true;
}

static bool
format_se(const bool allow_ambig, se_result res, const ChromLookup &cl,
          string &read, const string &read_name, string &cigar,
          ostream &out) {

  const se_element s = res.best;

  // GS: not the same as "should_report" because we need to account
  // for the allow_ambig flag
  if (!s.valid(read.size()) || (!allow_ambig && res.ambig()))
    return false;

  uint32_t ref_s = 0, ref_e = 0, chrom_idx = 0;
  if (chrom_and_posn(cl, cigar, s.pos, ref_s, ref_e, chrom_idx)) {
    sam_rec sr(read_name, 0, cl.names[chrom_idx], ref_s + 1,
               255, cigar, "*", 0, 0, read, "*");
    if (s.rc())
      set_flag(sr, samflags::read_rc);

    // GS checking allow_ambig to avoid the costly s.ambig()
    // when user has not requested it
    if (allow_ambig && res.ambig())
      set_flag(sr, samflags::secondary_aln);

    sr.add_tag("NM:i:" + to_string(s.diffs));
    sr.add_tag(s.elem_is_a_rich() ? "CV:A:A" : "CV:A:T");
    out << sr << '\n';
    return true;
  }
  return false;
}


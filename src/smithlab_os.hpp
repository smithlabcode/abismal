/* Part of SMITHLAB software
 *
 * Copyright (C) 2008-2018 University of Southern California and
 *                         Andrew D. Smith
 *
 * Authors: Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef SMITHLAB_OS_HPP
#define SMITHLAB_OS_HPP

#include <string>
#include <vector>
#include <unordered_map>

bool
isdir(const char *filename);

bool
is_fastq(const std::string &filename);

bool
is_valid_filename(const std::string &name,
                  const std::string &filename_suffix);

std::string
path_join(const std::string &a, const std::string &b);

void
identify_chromosomes(const std::string &chrom_file,
                     const std::string fasta_suffix, // could be const char*
                     std::unordered_map<std::string,
                     std::string> &chrom_files);

void
identify_and_read_chromosomes(const std::string &chrom_file,
                              const std::string fasta_suffix,
                              std::unordered_map<std::string,
                              std::string> &chrom_files);

void
read_dir(const std::string &dirname,
         std::string filename_suffix, // could be const char*
         std::vector<std::string> &filenames);

void
read_fasta_file(const std::string &filename,
                std::vector<std::string> &names,
                std::vector<std::string> &sequences);
void
read_fasta_file_short_names(const std::string &filename,
                            std::vector<std::string> &names,
                            std::vector<std::string> &sequences);

// This verstion looks for the sequence matching a particular name
void
read_fasta_file(const std::string &filename,
                const std::string &target,
                std::string &sequence);

void
read_fastq_file(const char *filename,
                std::vector<std::string> &names,
                std::vector<std::string> &sequences,
                std::vector<std::vector<double> > &scores);

void
read_fastq_file(const char *filename,
                std::vector<std::string> &names,
                std::vector<std::string> &sequences,
                std::vector<std::string> &scores);

void
read_filename_file(const char *filename,
                   std::vector<std::string> &filenames);

size_t
get_filesize(std::string filename);

std::string
basename(std::string filename);

void
parse_dir_basename_suffix(const std::string &full_path,
                          std::string &dirname,
                          std::string &base_name,
                          std::string &suffix);

void
format_dir_basename_suffix(const std::string &directory_name,
                           const std::string &file_base_name,
                           const std::string &file_name_suffix,
                           std::string &full_path);

std::string
strip_path(std::string full_path);

std::string
strip_path_and_suffix(std::string full_path);

void
read_dir(const std::string& dirname, std::vector<std::string> &filenames);

bool
is_valid_output_file(const std::string &filename);

#endif

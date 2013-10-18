/* Copyright 2012 Tobias Marschall
 *
 * This file is part of HaploClique.
 *
 * HaploClique is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * HaploClique is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with HaploClique.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FASTAREADER_H_
#define FASTAREADER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/unordered_map.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include "NamedDnaSequence.h"

class FastaReader {
private:
	FastaReader() {}
public:
	typedef boost::unordered_map<std::string,NamedDnaSequence*> reference_map_t;
	
	/** Reads input in FASTA format from the given stream and returns a map of sequence names
	 *  to sequences.
	 */
	template<typename Source>
	static std::auto_ptr<boost::unordered_map<std::string,NamedDnaSequence*> > parse(Source& is) {
		std::auto_ptr<boost::unordered_map<std::string,NamedDnaSequence*> > result(new boost::unordered_map<std::string,NamedDnaSequence*>());
		std::string line;
		int n = 1;
		NamedDnaSequence* current = 0;
		while (std::getline(is,line)) {
			if (line.size()==0) continue;
			// strip line from leading or trailing whitespaces
			size_t start = line.find_first_not_of(" \t\r\n");
			size_t end = line.find_last_not_of(" \t\r\n");
			line = line.substr(start, end-start+1);
			if (line.size()==0) continue;
			if (line[0] == '>') {
				// if (current!=0) {
				// 	std::cout << "Finished " << current->getName() << " " << current->size() << std::endl;
				// }
				start = line.find_first_not_of(" \t", 1);
				end = line.find_first_of(" ", start);
				if (end == std::string::npos) end = line.size();
				std::string name = line.substr(start,end-start);
				std::cerr << "Reading reference sequence \"" << name << "\"" << std::endl;
				if (result->find(name)!=result->end()) {
					delete result->at(name);
				}
				current = new NamedDnaSequence(name);
				(*result)[name] = current;
			} else {
				if (current == 0) {
					std::ostringstream oss;
					oss << "Error parsing FASTA input. Offending line: " << n << ": \"" << line << "\"";
					// TODO: free all named sequences
					throw std::runtime_error(oss.str());
				} else {
					current->append(line);
				}
			}
			n += 1;
		}
		// std::cout << "Finished " << current->getName() << " " << current->size() << std::endl;
		return result;
	}

	/** Read FASTA input from file, unzipping it if filename ends on ".gz". */
	inline static std::auto_ptr<boost::unordered_map<std::string,NamedDnaSequence*> > parseFromFile(const std::string& filename) {
		std::ifstream reference_istream(filename.c_str());
		bool input_zipped = filename.substr(filename.size()-3,3).compare(".gz") == 0;
		if (reference_istream.fail()) {
			throw std::runtime_error("Error opening file \"" + filename + "\".");
		}
		boost::iostreams::filtering_istream in;
		if (input_zipped) {
			in.push(boost::iostreams::gzip_decompressor());
		}
		in.push(reference_istream);
		std::auto_ptr<reference_map_t> result = FastaReader::parse(in);
		if (result->size() == 0) {
			if (input_zipped) {
				throw std::runtime_error("Error: references sequences empty or not properly gzipped.");
			} else {
				throw std::runtime_error("Error: references sequences empty.");
			}
		}
		return result;
	}

};

#endif /* FASTAREADER_H_ */

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

#include <fstream>
#include <sstream>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include "PositionSet.h"

using namespace std;

PositionSet::PositionSet() {
}

PositionSet::~PositionSet() {
	position_set_t::const_iterator it = position_set.begin();
	for (; it != position_set.end(); ++it) {
		delete it->second;
	}
}

void PositionSet::loadFromFile(std::string filename, const BamTools::RefVector& ref_data) {
	typedef boost::unordered_map<string,int> map_t;
	map_t map;
	for (size_t i=0; i<ref_data.size(); ++i) {
		map[ref_data[i].RefName] = i;
	}
	ifstream in(filename.c_str());
	if (in.fail()) {
		ostringstream oss;
		oss << "Could not open file \"" << filename << "\".";
		throw std::runtime_error(oss.str());
	}
	typedef boost::tokenizer<boost::char_separator<char> > tokenizer_t;
	boost::char_separator<char> whitespace_separator(" \t");
	string line;
	int linenr = 0;
	while (getline(in,line)) {
		linenr += 1;
		tokenizer_t tokenizer(line, whitespace_separator);
		vector<string> tokens(tokenizer.begin(), tokenizer.end());
		if ((tokens.size() < 2) || (tokens.size() > 3)) {
			ostringstream oss;
			oss << "Error parsing \"" << filename << "\". Offending line: " << linenr;
			throw std::runtime_error(oss.str());
		}
		try {
			map_t::const_iterator it = map.find(tokens[0]);
			if (it == map.end()) {
				ostringstream oss;
				oss << "Unknown chromosome name in \"" << filename << ", line: " << linenr;
				throw std::runtime_error(oss.str());
			}
			int chromosome_id = it->second;
			string chromosome = tokens[0];
			int position = boost::lexical_cast<int>(tokens[1]) - 1;
			set(chromosome_id, position);
		} catch(boost::bad_lexical_cast &){
			ostringstream oss;
			oss << "Error parsing \"" << filename << "\". Offending line: " << linenr;
			throw std::runtime_error(oss.str());
		}
	}
}

bool PositionSet::get(int chromosome_id, size_t position) {
	position_set_t::const_iterator it = position_set.find(chromosome_id);
	if (it == position_set.end()) return false;
	chromosome_bitset_t* c = it->second;
	if (position >= c->size()) return false;
	return (*c)[position];
}

void PositionSet::set(int chromosome_id, size_t position, bool value) {
	position_set_t::const_iterator it = position_set.find(chromosome_id);
	if (it == position_set.end()) {
		chromosome_bitset_t* c = new chromosome_bitset_t(position + 1);
		c->set(position, value);
		position_set[chromosome_id] = c;
	} else {
		chromosome_bitset_t* c = it->second;
		if (position >= c->size()) {
			c->resize(position + 1);
		}
		c->set(position, value);
	}
}

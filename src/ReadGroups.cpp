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

#include <sstream>

#include "ReadGroups.h"

using namespace std;

ReadGroups::ReadGroups(const BamTools::SamReadGroupDictionary& dict, bool merge_samplewise) : merge_samplewise(merge_samplewise) {
	typedef boost::unordered_map<std::string,int> sample_name_map_t;
	sample_name_map_t sample_name_map;
	for (BamTools::SamReadGroupConstIterator it = dict.Begin(); it != dict.End(); ++it) {
		if (!it->HasID()) throw std::runtime_error("Read group without ID");
		if (!it->HasSample()) throw std::runtime_error("Read group without sample field (SM)");
		if (id_map.find(it->ID) != id_map.end()) {
			ostringstream oss;
			oss << "Duplicate read group: " << it->ID;
			throw std::runtime_error(oss.str());
		}
		
		int sample_index = -1;
		sample_name_map_t::const_iterator sample_it = sample_name_map.find(it->Sample);
		if (sample_it == sample_name_map.end()) {
			sample_names.push_back(it->Sample);
			sample_index = sample_names.size() - 1;
			sample_name_map[it->Sample] = sample_index;
		} else {
			sample_index = sample_it->second;
		}
		int index = id_map.size();
		id_map[it->ID] = index;
		sample_map[index] = sample_index;
		names.push_back(it->ID);
	}
}

ReadGroups::~ReadGroups() {
}

int ReadGroups::get_index(const std::string& name) const {
	id_map_t::const_iterator it = id_map.find(name);
	if (it == id_map.end()) {
		return -1;
	} else {
		return it->second;
	}
}

int ReadGroups::getIndex(const std::string& name) const {
	if (merge_samplewise) {
		return getSampleIndex(name);
	} else {
		return get_index(name);
	}
}

int ReadGroups::getIndex(const BamTools::BamAlignment& aln) const {
	string rg = "";
	if (!aln.GetTag("RG", rg)) {
		return -1;
	}
	return getIndex(rg);
}

const string& ReadGroups::getName(int index) const {
	if (merge_samplewise) {
		return getSampleName(index);
	} else {
		assert(index >= 0);
		assert(index < names.size());
		return names[index];
	}
}

size_t ReadGroups::size() const {
	if (merge_samplewise) {
		return sample_names.size();
	} else {
		return names.size();
	}
}

int ReadGroups::getSampleIndex(const std::string& read_group_name) const {
	int index = get_index(read_group_name);
	if (index == -1) return -1;
	sample_map_t::const_iterator it = sample_map.find(index);
	assert(it != sample_map.end());
	return it->second;
}

int ReadGroups::getSampleIndex(const BamTools::BamAlignment& aln) const {
	string rg = "";
	if (!aln.GetTag("RG", rg)) {
		return -1;
	}
	return getSampleIndex(rg);
}

const string& ReadGroups::getSampleName(int index) const {
	assert(index >= 0);
	assert(index < sample_names.size());
	return sample_names[index];
}

size_t ReadGroups::sampleCount() const {
	return sample_names.size();
}

int ReadGroups::readGroupIndexToSampleIndex(int index) const {
	if (merge_samplewise) {
		assert(index >= 0);
		assert(index < sample_names.size());
		return index;
	} else {
		assert(index >= 0);
		assert(index < names.size());
		sample_map_t::const_iterator it = sample_map.find(index);
		assert(it != sample_map.end());
		return it->second;
	}
}

std::ostream& operator<<(std::ostream& os, const ReadGroups& rg) {
	for (size_t i=0; i<rg.size(); ++i) {
		os << "Read group " << i << ": " << rg.getName(i) << endl;
	}
	for (size_t i=0; i<rg.sampleCount(); ++i) {
		os << "Sample " << i << ": " << rg.getSampleName(i) << endl;
	}
	return os;
}

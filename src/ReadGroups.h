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


#ifndef READGROUPS_H
#define READGROUPS_H

#include <ostream>
#include <vector>

#include <boost/unordered_map.hpp>

#include <api/SamReadGroup.h>
#include <api/SamReadGroupDictionary.h>
#include <api/BamAlignment.h>

class ReadGroups {
private:
	// mapping read group names (IDs) to their index
	typedef boost::unordered_map<std::string,int> id_map_t;
	id_map_t id_map;
	// maps read group index to sample index
	typedef boost::unordered_map<int,int> sample_map_t;
	sample_map_t sample_map;
	// names (IDs) of read groups
	std::vector<std::string> names;
	// sample names
	std::vector<std::string> sample_names;
	bool merge_samplewise;
	int get_index(const std::string& name) const;
public:
	/** @param merge_samplewise If true, then this class behaves as if readgroups were samples. */
	ReadGroups(const BamTools::SamReadGroupDictionary& dict, bool merge_samplewise = false);
	virtual ~ReadGroups();
	
	/** Get index of a read group by name, i.e. by what is called ID in BAM headers.
	 *  Returns -1 if name is unknown.
	 */
	int getIndex(const std::string& read_group_name) const;
	
	/** Get index for a given alignment. */
	int getIndex(const BamTools::BamAlignment& aln) const;
	
	/** Get read group by index as returned by getIndex(). */
	const std::string& getName(int index) const;

	/** Returns the number of read groups. */
	size_t size() const;

	/** Returns the sample index for a given read group name.
	 *  Returns -1 if name is unknown.
	 */
	int getSampleIndex(const std::string& read_group_name) const;

	int getSampleIndex(const BamTools::BamAlignment& aln) const;
	
	int readGroupIndexToSampleIndex(int index) const;

	const std::string& getSampleName(int index) const;

	size_t sampleCount() const;

};

std::ostream& operator<<(std::ostream& os, const ReadGroups& rg);

#endif // READGROUPS_H

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

#include <math.h>

#include "CoverageMonitor.h"

using namespace std;

CoverageMonitor::CoverageMonitor(const ReadGroups* read_groups) {
	this->read_groups = read_groups;
	if (read_groups != 0) {
		for (size_t i=0; i<read_groups->size(); ++i) {
			readwise_monitors.push_back(SingleTrackCoverageMonitor());
		}
	}
}

CoverageMonitor::~CoverageMonitor() {
}

void CoverageMonitor::pruneLeftOf(size_t pos) {
	global_monitor.pruneLeftOf(pos);
	if (read_groups != 0) {
		for (size_t i=0; i<read_groups->size(); ++i) {
			readwise_monitors[i].pruneLeftOf(pos);
		}
	}
}


void CoverageMonitor::addAlignment(const AlignmentRecord& ap) {
	global_monitor.addAlignment(ap);
	if (read_groups != 0) {
		assert(ap.getReadGroup() < (int)readwise_monitors.size());
		readwise_monitors[ap.getReadGroup()].addAlignment(ap);
	}
}

size_t CoverageMonitor::probeAlignment(const AlignmentRecord& ap) const {
	return global_monitor.probeAlignment(ap);
}


size_t CoverageMonitor::getCoverage(size_t pos) {
	return global_monitor.getCoverage(pos);
}

bool CoverageMonitor::hasReadGroups() const {
	return read_groups != 0;
}

size_t CoverageMonitor::getReadGroupCoverage(size_t read_group, size_t pos) {
	assert(read_groups != 0);
	assert(read_group < readwise_monitors.size());
	return readwise_monitors[read_group].getCoverage(pos);
}

unique_ptr<vector<size_t> > CoverageMonitor::getReadGroupCoverages(size_t pos) {
	assert(read_groups != 0);
	unique_ptr<vector<size_t> > result(new vector<size_t>());
	for (size_t i=0; i<read_groups->size(); ++i) {
		result->push_back(readwise_monitors[i].getCoverage(pos));
	}
	return result;
}

ostream& operator<<(ostream& os, const CoverageMonitor& cm) {
	os << "Overall coverage: " << cm.global_monitor << endl;
	if (cm.read_groups != 0) {
		for (size_t i=0; i<cm.read_groups->size(); ++i) {
			os << "Read group " << i << " coverage: " << cm.readwise_monitors[i] << endl;
		}
	}
	return os;
}

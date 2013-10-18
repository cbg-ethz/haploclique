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

#include "SingleTrackCoverageMonitor.h"

using namespace std;

const int SingleTrackCoverageMonitor::MAX_LOOKAHEAD = 4096;

SingleTrackCoverageMonitor::SingleTrackCoverageMonitor() {
	offset = 0;
}

SingleTrackCoverageMonitor::~SingleTrackCoverageMonitor() {
}

void SingleTrackCoverageMonitor::extendTo(size_t pos) {
	size_t n = offset+coverage.size();
	while (pos >= n) {
		coverage.push_back(overhanging.size());
		n += 1;
		// discard all overhanging segments whose end is now within [offset..n]
		while (!overhanging.empty()) {
			multiset<size_t>::iterator first = overhanging.begin();
			if (*first < n) {
				overhanging.erase(first);
			} else {
				break;
			}
		}
	}
}

void SingleTrackCoverageMonitor::pruneLeftOf(size_t pos) {
	while (offset < pos) {
		if (coverage.size() > 0) {
			coverage.pop_front();
		}
		offset += 1;
	}
	while (!overhanging.empty()) {
		multiset<size_t>::iterator first = overhanging.begin();
		if (*first < offset) {
			overhanging.erase(first);
		} else {
			break;
		}
	}
}


void SingleTrackCoverageMonitor::addAlignment(const AlignmentRecord& ap) {
	if (ap.getIntervalStart() > ap.getIntervalEnd()) return;
	assert(ap.getIntervalStart() >= offset);
	// pruneLeftOf(ap.getIntervalStart());
	size_t n = min((size_t)ap.getIntervalEnd(), (size_t)(ap.getIntervalStart()+MAX_LOOKAHEAD-1));
	extendTo(n);
	for (size_t i=ap.getIntervalStart(); i<=ap.getIntervalEnd(); ++i) {
		size_t j = i - offset;
		if (j<coverage.size()) {
			coverage[j] += 1;
		} else {
			overhanging.insert(ap.getIntervalEnd());
			break;
		}
	}
}

size_t SingleTrackCoverageMonitor::probeAlignment(const AlignmentRecord& ap) const {
// 	cout << "SingleTrackCoverageMonitor::probeAlignment: "<< ap.getIntervalStart() << ", " << ap.getIntervalEnd() << endl;
	assert(ap.getIntervalStart() >= offset);
	// pruneLeftOf(ap.getIntervalStart());
	size_t result = 0;
	if (ap.getIntervalEnd() < offset+coverage.size()) {
		for (size_t i=ap.getIntervalStart(); i<=ap.getIntervalEnd(); ++i) {
			result = max(result, coverage[i-offset]);
		}
	} else {
		for (size_t i=ap.getIntervalStart()-offset; i<coverage.size(); ++i) {
			result = max(result, coverage[i]);
		}
		size_t n = 0;
		multiset<size_t>::const_iterator it = overhanging.lower_bound(ap.getIntervalStart());
		for (;it!=overhanging.end(); ++it) {
			n += 1;
		}
		result = max(result, n);
	}
	return result + 1;
}


size_t SingleTrackCoverageMonitor::getCoverage(size_t pos) {
	assert(pos >= offset);
	extendTo(pos);
	return coverage[pos - offset];
}

ostream& operator<<(ostream& os, const SingleTrackCoverageMonitor& cm) {
	os << "Region: [" << cm.offset << ".." << cm.offset+cm.coverage.size()-1 << "]: [";
	for (size_t i=0; i<cm.coverage.size(); ++i) {
		if (i>0) os << ",";
		os << cm.coverage[i];
	}
	os << "], overhanging: " << cm.overhanging.size() << ": ";
	multiset<size_t>::const_iterator it = cm.overhanging.begin();
	for (;it!=cm.overhanging.end(); ++it) {
		os << " " << (*it);
	}
	return os;
}


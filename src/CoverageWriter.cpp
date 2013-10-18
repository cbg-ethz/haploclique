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

#include "CoverageWriter.h"

using namespace std;

CoverageWriter::CoverageWriter(const string& filename) : finished(false), ofs(filename.c_str()), rightmost_pos(0), pos(0) {
}

CoverageWriter::~CoverageWriter() {
	finish();
}

void CoverageWriter::addAlignment(const AlignmentRecord& ap) {
	if (ap.getIntervalEnd()>rightmost_pos) {
		rightmost_pos = ap.getIntervalEnd();
	}
	coverage_monitor.addAlignment(ap);
	for (; pos<ap.getIntervalStart(); ++pos) {
		ofs << coverage_monitor.getCoverage(pos) << endl;
	}
	coverage_monitor.pruneLeftOf(pos);
}

void CoverageWriter::finish() {
	if (finished) return;
	for (; pos<=rightmost_pos; ++pos) {
		ofs << coverage_monitor.getCoverage(pos) << endl;
	}
	ofs.close();
}

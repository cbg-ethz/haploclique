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

#ifndef READSETSIGNIFICANCETESTER_H_
#define READSETSIGNIFICANCETESTER_H_

#include <vector>

/** Class to perform a statistical test on a set of read alignments and return a p-value.
 *  It is tested whether a given sum of insert sizes is compatible with the null hypothesis.
 */
class ReadSetSignificanceTester {
public:
	/** Returns the two sided p-value for the event that "count" alignments have the given sum of insert sizes. 
	 *  @param read_group_counts[i] gives the number of reads from read group i.
	 */
	virtual double computeSignificance(int insert_size_sum, int count, const std::vector<int>& read_group_counts) const = 0;
	
	virtual ~ReadSetSignificanceTester() {}
};

#endif /* READSETSIGNIFICANCETESTER_H_ */

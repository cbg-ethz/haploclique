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


#ifndef INDELLENGTHDISTRIBUTION_H
#define INDELLENGTHDISTRIBUTION_H

#include <vector>
#include <iostream>

#include "HistogramBasedDistribution.h"

class IndelLengthDistribution {
private:
	std::vector<int> phred_costs;
public:
	IndelLengthDistribution(const HistogramBasedDistribution& distribution);
	IndelLengthDistribution(const std::vector<int>& phred_costs);

	/** Returns PHRED-like costs of an indel of the given length. */
	int getPhredCost(int length) const;
	
	friend std::ostream& operator<<(std::ostream& os, const IndelLengthDistribution& ild);
};

#endif // INDELLENGTHDISTRIBUTION_H

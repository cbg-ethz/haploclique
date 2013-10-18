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

#ifndef TYPES_H_
#define TYPES_H_

#include <boost/dynamic_bitset.hpp>

typedef unsigned int alignment_id_t;
typedef boost::dynamic_bitset<> alignment_set_t;

typedef struct mean_and_stddev_t {
	double mean;
	double stddev;
	mean_and_stddev_t() : mean(0.0), stddev(1.0) {}
	mean_and_stddev_t(double mean, double stddev) : mean(mean), stddev(stddev) {}
} mean_and_stddev_t;

#endif /* TYPES_H_ */

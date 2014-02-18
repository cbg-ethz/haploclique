/* Copyright 2014 Armin Töpfer
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace boost;

int main() {
    set<string> reads;
    ifstream tsv_stream("data_clique_to_reads.tsv");
    vector<string>::iterator it;
    string line;
    vector<string> split;

    while (getline(tsv_stream, line)) {
        trim_right(line);
        boost::split(split, line, boost::is_any_of("\t,"));
        it=++split.begin();
        std::copy(it, split.end(), std::inserter( reads, reads.begin()));
    }
    tsv_stream.close();

    ifstream ap_stream2("alignment.prior");
    ofstream singles;
    singles.open ("singles.prior", ios::out);
    while (getline(ap_stream2, line)) {
        if (boost::starts_with(line, "Clique")) continue;
        boost::split(split, line, boost::is_any_of(" "));
        if (reads.find(split[0]) == reads.end()) {
            singles << line << endl;
        }
    }
    ap_stream2.close();
    singles.close();

    return 0;
}

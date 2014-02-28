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
#include <map>

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

using namespace std;
using namespace boost;

int main() {
    map<string,int> read_map;
    ifstream tsv_stream("data_clique_to_reads.tsv");
    vector<string>::iterator it;
    string tsv_stream_line;
    vector<string> clique;
    vector<string> words;
    while (getline(tsv_stream, tsv_stream_line)) {
        trim_right(tsv_stream_line);
        boost::split(clique, tsv_stream_line, boost::is_any_of("\t"));
        boost::split(words, clique[1], boost::is_any_of(","));
        for (it=words.begin(); it!=words.end(); ++it) {
            if (boost::starts_with(*it, "HCOUNT|")) continue;
            if (read_map.find(*it) == read_map.end()) {
                read_map[*it] = 1;
            } else {
                read_map[*it]++;
            }
        }
    }
    tsv_stream.close();

    ifstream tsv_stream2("data_clique_to_reads.tsv");
    stringstream ss;
    int count_global = 0;
    while (getline(tsv_stream2, tsv_stream_line)) {
        trim_right(tsv_stream_line);
        boost::split(clique, tsv_stream_line, boost::is_any_of("\t"));
        boost::split(words, clique[1], boost::is_any_of(","));
        double count = 0;
        for (it=words.begin(); it!=words.end(); ++it) {
            if (boost::starts_with(*it, "HCOUNT|")) {
                vector<string> split;
                boost::split(split, *it, boost::is_any_of("|"));
                count += boost::lexical_cast<int>(split[1]);
            } else if (read_map[*it] == 1) {
                count++;
            } else {
                count += 1.0/double(read_map[*it]);
            }
        }
        count_global += count;
        ss << clique[0] << "\t" << count << endl;
    }
    tsv_stream2.close();
    ofstream out;
    out.open ("clique_count.tsv", ios::out);
    out << ss.str();
    out.close();

    // ifstream tsv_stream2("data_clique_to_reads.tsv");
    // stringstream ss;
    // int count_global = 0;
    // while (getline(tsv_stream2, tsv_stream_line)) {
    //     trim_right(tsv_stream_line);
    //     boost::split(clique, tsv_stream_line, boost::is_any_of("\t"));
    //     ss << clique[0] << "\t";
    //     boost::split(words, clique[1], boost::is_any_of(","));
    //     int count = 0;
    //     for (it=words.begin(); it!=words.end(); ++it) {
    //         if (boost::starts_with(*it, "HCOUNT|")) {
    //             vector<string> split;
    //             boost::split(split, *it, boost::is_any_of("|"));
    //             count += boost::lexical_cast<int>(split[1]);
    //             cerr << "XLR: " << *it << "\t" << boost::lexical_cast<int>(split[1]) << endl;
    //         } else if (read_map[*it] == 1) {
    //             read_map.erase(*it);
    //             count++;
    //         } else {
    //             ss << *it << ",";
    //         }
    //     }
    //     count_global += count;
    //     ss << "HCOUNT|" << count << endl;
    // }
    // tsv_stream2.close();
    // ofstream out;
    // out.open ("data_clique_to_reads.tsv", ios::out);
    // out << ss.str();
    // out.close();

    // cerr << read_map.size() << "\t" << count_global << "\t" << read_map.size()+count_global << endl;
    cerr << read_map.size() << endl;

    // map<string,int>::iterator read_map_it;

    // for (read_map_it = read_map.begin(); read_map_it != read_map.end(); ++read_map_it) {
    //     cerr << read_map_it->first << "\t" << read_map_it->second << endl;
    // }

    return 0;
}

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
    int sum = 0;
    map<string,int> read_map;
    ifstream tsv_stream("data_clique_to_reads.tsv");
    vector<string>::iterator it;
    string tmpstring;
    vector<string> clique;
    vector<string> words;
    while (getline(tsv_stream, tmpstring)) {
        trim_right(tmpstring);
        boost::split(clique, tmpstring, boost::is_any_of("\t"));
        boost::split(words, clique[1], boost::is_any_of(","));
        for (it=words.begin(); it!=words.end(); ++it) {
            if (boost::starts_with(*it, "HCOUNT|")) continue;
            if (read_map.find(*it) == read_map.end()) {
                read_map[*it] = 1;
                sum++;
            } else {
                read_map[*it]++;
            }
        }
    }
    tsv_stream.close();

    ifstream tsv_stream2("data_clique_to_reads.tsv");
    stringstream ss;
    map<string,double> clique_count_map;
    while (getline(tsv_stream2, tmpstring)) {
        trim_right(tmpstring);
        boost::split(clique, tmpstring, boost::is_any_of("\t"));
        boost::split(words, clique[1], boost::is_any_of(","));
        double count = 0;
        for (it=words.begin(); it!=words.end(); ++it) {
            if (boost::starts_with(*it, "HCOUNT|")) {
                vector<string> split;
                boost::split(split, *it, boost::is_any_of("|"));
                count += boost::lexical_cast<int>(split[1]);
            } else if (read_map[*it] == 1) {
                count += 1.0/sum;
            } else {
                count += 1.0/double(read_map[*it])/sum;
            }
        }
        clique_count_map[clique[0]] = count;
    }
    tsv_stream2.close();

    ifstream s_stream("data_cliques_single.fastq");
    int line_number = 0;
    string tmp_id;
    map<string, string> clique_id_map_single;
    while (getline(s_stream, tmpstring)) {
        trim_right(tmpstring);
        switch (line_number++ % 4) {
            case 0:
                tmp_id = tmpstring.substr(1);
                //cerr << tmp_id << "\t";
                break;
            case 1:
                clique_id_map_single[tmp_id] = tmpstring;
                //cerr << tmpstring << endl;
                break;
            default:break;
        }
    }

    ifstream p1_stream("data_cliques_paired_R1.fastq");
    line_number = 0;
    map<string, string> clique_id_map_p1;
    while (getline(p1_stream, tmpstring)) {
        trim_right(tmpstring);
        switch (line_number++ % 4) {
            case 0:
                tmp_id = tmpstring.substr(1);
                //cerr << tmp_id << "\t";
                break;
            case 1:
                clique_id_map_p1[tmp_id] = tmpstring;
                //cerr << tmpstring << endl;
                break;
            default:break;
        }
    }

    ifstream p2_stream("data_cliques_paired_R1.fastq");
    line_number = 0;
    map<string, string> clique_id_map_p2;
    while (getline(p2_stream, tmpstring)) {
        trim_right(tmpstring);
        switch (line_number++ % 4) {
            case 0:
                tmp_id = tmpstring.substr(1);
                //cerr << tmp_id << "\t";
                break;
            case 1:
                clique_id_map_p2[tmp_id] = tmpstring;
                //cerr << tmpstring << endl;
                break;
            default:break;
        }
    }

    map<string,double>::iterator key_it = clique_count_map.begin();
    for (;key_it != clique_count_map.end(); ++key_it) {
        if (clique_id_map_single.find(key_it->first) != clique_id_map_single.end()) {
            cout << ">" << key_it->first << "x " << key_it->second << " ";
            cout << clique_id_map_single[key_it->first] << endl;
        } else if (clique_id_map_p1.find(key_it->first) != clique_id_map_p1.end()) {
            cout << ">" << key_it->first << "/1x " << key_it->second << " ";
            cout << clique_id_map_p1[key_it->first] << endl;
            cout << ">" << key_it->first << "/2x " << key_it->second << " ";
            cout << clique_id_map_p2[key_it->first] << endl;
        }
    }

    return 0;
}

#ifndef LOGWRITER_H_
#define LOGWRITER_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
#include <map>

class LogWriter {
private:
    std::map<unsigned int, std::list<unsigned int>> vertices_;
    std::map<unsigned int, std::list<unsigned int>> cliques_;
    //keeps track of how many of original reads contribute to each clique
    std::map<unsigned int, unsigned int> reads_in_cliques_;
    //keep track of how many times a read contributes to a clique
    std::vector<unsigned int> read_has_cliques_;
    std::ofstream logfile_;
    std::string filename_;
    unsigned int iteration_;
public:
    LogWriter(std::string fn, std::vector<unsigned int>& read_clique_counter) : vertices_(), logfile_(fn, std::ofstream::out), filename_(fn), iteration_(0), read_has_cliques_(read_clique_counter) {}
    virtual ~LogWriter() {}

    void reportEdge(unsigned int v1, unsigned int v2) {
        if (vertices_.count(v1) == 0) vertices_[v1] = std::list<unsigned int>();
        if (vertices_.count(v2) == 0) vertices_[v2] = std::list<unsigned int>();

        vertices_[v1].push_back(v2);
        vertices_[v2].push_back(v1);
    };

    void reportClique(unsigned int id, std::list<unsigned int> clique) {cliques_[id] = clique;}
    void reportReadsInCliques(unsigned int id, unsigned int reads) { reads_in_cliques_[id] = reads;}
    void reportReadsHasCliques(const std::set<int>& readNames) {
        for (const auto& i : readNames){
            read_has_cliques_[i]++;
        }

    };
    void initialize() {vertices_.clear(); cliques_.clear(); reads_in_cliques_.clear(); std::fill(read_has_cliques_.begin(), read_has_cliques_.end(), 0);}
    void finish() {
        //logs edges of every node to other nodes
        logfile_ << ">--" << iteration_++ << "--<" << std::endl;

        for (const auto& adj : vertices_) {
            logfile_ << adj.first << " -- {";
            for (const auto& i : adj.second) {
                logfile_ << " " << i;
            }
            logfile_ << "};"<< std::endl;
        }

        logfile_ << "-" << iteration_-1 << "-" << std::endl;

        //logs which clique contains which nodes (note: after merging of Alignment Records, these number does not represent the original alignment but of the merged ones)
        for (const auto& cl : cliques_) {
            logfile_ << cl.first << ":";
            for (const auto& i : cl.second) {
                logfile_ << " " << i;
            }
            logfile_ << " (" << reads_in_cliques_[cl.first] << ")"<< std::endl;
        }
        //logs number of cliques in which each read appears
        logfile_ << "---abundance of reads in cliques---" << std::endl;
        for (int i = 0; i < read_has_cliques_.size(); i++){
            logfile_ << i << ": " << read_has_cliques_[i] << std::endl;
        }
    }
};

#endif //LOGWRITER_H_

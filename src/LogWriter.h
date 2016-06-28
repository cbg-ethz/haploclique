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
    std::ofstream logfile_;
    std::string filename_;
    unsigned int iteration_;
public:
    LogWriter(std::string fn) : vertices_(), logfile_(fn, std::ofstream::out), filename_(fn), iteration_(0) {};
    virtual ~LogWriter() {};

    void reportEdge(unsigned int v1, unsigned int v2) {
        if (vertices_.count(v1) == 0) vertices_[v1] = std::list<unsigned int>();
        if (vertices_.count(v2) == 0) vertices_[v2] = std::list<unsigned int>();

        vertices_[v1].push_back(v2);
        vertices_[v2].push_back(v1);
    };

    void reportClique(unsigned int id, std::list<unsigned int> clique) {cliques_[id] = clique;};
    void initialize() {vertices_.clear(); cliques_.clear();};
    void finish() {
        logfile_ << ">--" << iteration_++ << "--<" << std::endl;

        for (const auto& adj : vertices_) {
            logfile_ << adj.first << " -- {";
            for (const auto& i : adj.second) {
                logfile_ << " " << i;
            }
            logfile_ << "};"<< std::endl;
        }

        logfile_ << "---" << std::endl;

        for (const auto& cl : cliques_) {
            logfile_ << cl.first << ":";
            for (const auto& i : cl.second) {
                logfile_ << " " << i;
            }
            logfile_ << std::endl;
        }

    };
};

#endif //LOGWRITER_H_

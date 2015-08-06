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

#ifndef CLIQUECOLLECTOR_H_
#define CLIQUECOLLECTOR_H_

#include "Clique.h"
#include <deque>

class CliqueCollector {
private:
    std::deque<AlignmentRecord*>* superReads;
    int id;
public:
    CliqueCollector() : id(0) {
        superReads = new std::deque<AlignmentRecord*>;
    };

    virtual ~CliqueCollector() {
        assert(superReads->empty());
        delete superReads;
    };

	void add(std::unique_ptr<Clique> clique) {
        assert(clique.get() != nullptr);

//        std::cerr << "Clique " << id << std::endl;
        std::unique_ptr<std::vector<const AlignmentRecord*>> alignments = clique->getAllAlignments();
        AlignmentRecord* al;

        if (alignments->size() > 1) {
            al = new AlignmentRecord(alignments, this->id++);
        } else {
            al = new AlignmentRecord(*(alignments->front()));
        }

        superReads->push_back(al);
    };

    std::deque<AlignmentRecord*>* finish()
    {
        auto retVal = superReads;
        superReads = new std::deque<AlignmentRecord*>;
        return retVal;
    };
};

#endif /* CLIQUECOLLECTOR_H_ */

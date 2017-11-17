/* Copyright 2012-2014 Tobias Marschall and Armin Töpfer
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

#ifndef NEWEDGECALCULATOR_H_
#define NEWEDGECALCULATOR_H_

#include <set>

#include "AlignmentRecord.h"
#include "EdgeCalculator.h"

using namespace std;

class NewEdgeCalculator : public EdgeCalculator {
private:

    double MIN_OVERLAP_CLIQUES;
    double MIN_OVERLAP_SINGLE;
    static const double FRAME_SHIFT_WEIGHT;
    double Q;
    double EDGE_QUASI_CUTOFF_MIXED;
    double EDGE_QUASI_CUTOFF_SINGLE;
    double EDGE_QUASI_CUTOFF;
    bool FRAMESHIFT_MERGE;
    bool NOPROB0;
    //std::unordered_map<int, double> SIMPSON_MAP;
    std::vector<double> SIMPSON_MAP;
    mutable std::vector<std::pair<int,int>> aub;

    void calculateProbM(const AlignmentRecord::mapValue &val1, const AlignmentRecord::mapValue &val2, double &res) const;
    void calculateProb0(const AlignmentRecord::mapValue &val1, double &res) const;
    bool checkGaps(const std::vector<AlignmentRecord::mapValue> & cov_ap1, const std::vector<AlignmentRecord::mapValue> & cov_ap2, const std::vector<std::pair<int,int>>& aub) const;
    bool similarityCriterion(const AlignmentRecord & a1, const std::vector<AlignmentRecord::mapValue> & cov_ap1, const AlignmentRecord & a2, const std::vector<AlignmentRecord::mapValue> & cov_ap2, double probM, double prob0, int tc, int cc) const;
    void iterateCovAp(bool pe1, unsigned int& pos2, double& probM, unsigned int& pos1, unsigned int& equalBase, bool pe2, const std::vector<AlignmentRecord::mapValue>& cov_ap2, int& tc, double& prob0, const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap1, const AlignmentRecord& ap1, int& cc) const;
    void iterateRemainingCovAp(int& tc, bool pe1, const std::vector<AlignmentRecord::mapValue>& cov_ap2, bool pe2, double& prob0, const AlignmentRecord& ap1, const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap1, unsigned int& pos1) const;
    void iterateRemainingCovAp2(const AlignmentRecord& ap2, const std::vector<AlignmentRecord::mapValue>& cov_ap2, bool pe1, const AlignmentRecord& ap1, int& tc, double& prob0, const std::vector<AlignmentRecord::mapValue>& cov_ap1, bool pe2, unsigned int& pos2) const;

    void computeProbM(const char& base1, const char& qual1, const char& base2, const char& qual2, double &res ) const;
    void computeProb0(int ref_pos, double &res ) const;
    

public:
    NewEdgeCalculator(double Q, double edge_quasi_cutoff, double overlap, bool frameshift_merge, unordered_map<int, double>& simpson_map, double edge_quasi_cutoff_single, double overlap_single, double edge_quasi_cutoff_mixed, unsigned int maxPosition, bool noProb0);
    virtual ~NewEdgeCalculator();

    /** Decides whether an edge is to be drawn between the two given nodes. */
    virtual bool edgeBetween(const AlignmentRecord& ap1, const AlignmentRecord& ap2) const;

    /** Compute a length range. An alignment pair with a length outside this range is
     *  guaranteed not to have an edge to the given pair ap. */
    virtual void getPartnerLengthRange(const AlignmentRecord& ap, unsigned int* min, unsigned int* max) const;

    void setOverlapCliques(double d);
    
    double getOverlapCliques() const;
    
    bool checkGapsCigar(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0 ,int& cc, int& tc) const;
    bool checkGapsSingle(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0, int& cc, int& tc, int i, int j) const;
    bool checkGapsPaired(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0, int& cc, int& tc) const;
    bool checkGapsMixed(const AlignmentRecord& ap1, const AlignmentRecord& ap2, double& probM, double& prob0,int& cc, int& tc) const;
    void noOverlapCheckgap(const AlignmentRecord& ap, int& c_pos, int& q_pos, int& ref_pos, double& prob0, int& tc, int i) const;
    bool overlapCheckgap(const AlignmentRecord& ap1, const AlignmentRecord& ap2, int& c_pos1, int& c_pos2, int& q_pos1, int& q_pos2, int& ref_pos, double& probM,int& cc, int i, int j) const;


};
#endif /* NEWEDGECALCULATOR_H_ */

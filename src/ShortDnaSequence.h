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

#ifndef SHORTDNASEQUENCE_H_
#define SHORTDNASEQUENCE_H_

#include <iostream>
#include <string>
#include <boost/smart_ptr.hpp>

/** A DNA sequence stored in an 8-bit encoding. Can compute and cache its reverse complement. */
class ShortDnaSequence {
private:
    typedef boost::shared_ptr<std::string> string_ptr_t;
    mutable string_ptr_t forward;
    mutable string_ptr_t forward_qualities;
    mutable string_ptr_t backward;
    mutable string_ptr_t backward_qualities;
    ShortDnaSequence(string_ptr_t& forward, string_ptr_t& forward_qualities, string_ptr_t& backward, string_ptr_t& backward_qualities);
public:
    ShortDnaSequence();
    ShortDnaSequence(const std::string& dna, const std::string& qualities);
    virtual ~ShortDnaSequence();

    ShortDnaSequence reverseComplement() const;
    size_t size() const;
    char operator[](size_t pos) const;
    char qualityChar(size_t pos) const;
    double qualityCorrect(size_t pos) const;
    const std::string& toString() const;
    const std::string& qualityString() const;
    double qualityCorrectLog(size_t pos) const;
    friend std::ostream& operator<<(std::ostream& os, const ShortDnaSequence& s);
};

#endif /* SHORTDNASEQUENCE_H_ */

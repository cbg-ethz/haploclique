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

#ifndef NAMEDDNASEQUENCE_H_
#define NAMEDDNASEQUENCE_H_

#include <string>
#include <vector>

/** Stores a DNA sequence (over the alphabet {A,C,G,T,N}) along with its name. Uses four bits per character.
 *  The representation is not case-sensitive unless requested at construction time. */
class NamedDnaSequence {
private:
	bool case_sensitive;
	std::string name;
	size_t length;
	std::vector<unsigned char> sequence;

	/** Convert from 4-bit code to ASCII. */
	char decode(unsigned char c) const;
	/** Convert from ASCII to 4-bit code. */
	unsigned char encode(char c) const;
public:
	explicit NamedDnaSequence(const std::string& name);
	NamedDnaSequence(const std::string& name, const std::string& sequence);
	NamedDnaSequence(const std::string& name, bool case_sensitive);
	NamedDnaSequence(const std::string& name, const std::string& sequence, bool case_sensitive);
	virtual ~NamedDnaSequence() {}

	/** Appends a given sequence. */
	void append(const std::string& s);

	const std::string& getName() const { return name; }
	size_t size() const { return length; }

	char operator[](size_t pos) const;

	std::string substr(size_t pos, size_t length) const;
};

#endif /* NAMEDDNASEQUENCE_H_ */

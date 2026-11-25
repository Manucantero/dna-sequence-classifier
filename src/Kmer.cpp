/*
 * Metodología de la Programación: Kmer0
 * Curso 2023/2024
 */

/* 
 * File:   main.cpp
 * @author Manuel Cantero Ruiz  <manucantero@correo.ugr.es>
 * Created on 26 February 2024
 */

#include <iostream>
#include <cstring>
#include "Kmer.h"

using namespace std;

Kmer::Kmer(int k) {

    if (k <= 0) { // k is less or equal than zero

        throw invalid_argument(// composed string
                string("Kmer (const string& text: ") +
                " k is less or equal than zero ");
    }
    for (int i = 0; i < k; i++) {
        _text += MISSING_NUCLEOTIDE;
    }
}

Kmer::Kmer(const string&text) {

    if (text.length() == 0) { // text is empty
        throw invalid_argument(// composed string
                string("Kmer (const string& text: ") +
                " text is an empty string ");
    }
    _text = text;
}

int Kmer::getK() const {

    return (_text.length());
}

int Kmer::size() const {

    return (_text.length());
}

string Kmer::toString() const {

    return (_text);
}

const char& Kmer::at(int index) const {

    if (index < 0 || index >= (size())) { // out_of_range
        throw out_of_range(// composed string
                string("const char& Kmer::at(int index) const: ") +
                " index is not in the range from 0 to k-1 (both included)");
    }
    return (_text.at(index));
}

char& Kmer::at(int index) {

    if (index < 0 || index >= getK()) { // out_of_range
        throw out_of_range(// composed string
                string("char& Kmer::at(int index): ") +
                " index is not in the range from 0 to k-1");
    }
    return (_text.at(index));
}

void Kmer::normalize(const string& validNucleotides) {

    toUpper();
    for (int i = 0; i < getK(); i++) {
        if (!IsValidNucleotide(_text.at(i), validNucleotides)) {
            _text.at(i) = MISSING_NUCLEOTIDE;
        }
    }
}

Kmer Kmer::complementary(const string& nucleotides,
        const string& complementaryNucleotides) const {

    if ((int) nucleotides.length() != (int) complementaryNucleotides.length()) { // invalid_argument
        throw invalid_argument(// composed string
                string("Kmer Kmer::complementary(const string& nucleotides, const string& complementaryNucleotides) const : ") +
                " the sizes of nucleotides and complementaryNucleotides " +
                "are not the same");
    }
    string complementary;

    for (int i = 0; i < (int) _text.size(); i++) {

        char x = _text.at(i);
        size_t pos = nucleotides.find(x);

        if (pos != string::npos) {
            complementary += complementaryNucleotides.at(pos);
        } else {
            complementary += x;
        }
    }
    return ( Kmer(complementary));
}

bool IsValidNucleotide(char nucleotide, const string& validNucleotides) {

    bool isvalid = false;
    int i = 0;
    int max = validNucleotides.length();

    while (!isvalid && i < max) {

        if (nucleotide == validNucleotides.at(i)) {
            isvalid = true;
        } else {
            i++;
        }
    }
    return isvalid;
}

void Kmer::toLower() {

    for (int i = 0; i < getK(); i++) {
        _text.at(i) = tolower(_text.at(i));
    }
}

void Kmer::toUpper() {

    for (int i = 0; i < getK(); i++) {
        _text.at(i) = toupper(_text.at(i));
    }
}

void ToLower(Kmer& kmer) {

    for (int i = 0; i < kmer.getK(); i++) {
        kmer.at(i) = tolower(kmer.at(i));
    }
}

void ToUpper(Kmer& kmer) {

    for (int i = 0; i < kmer.getK(); i++) {
        kmer.at(i) = toupper(kmer.at(i));
    }
}

/**
 * @brief Writes this object to the given output stream. All the characters in
 * the string of this Kmer (including '\0') are sent to the ouput stream.
 * Query method
 * @param outputStream An output stream where this object will be written
 */
void Kmer::write(ostream& outputStream) const {

     outputStream.write(_text.c_str(), _text.size() + 1);
}

/**
 * @brief Reads this object from the given input stream. It reads characters
 * from the given input stream and put them in the text of this Kmer
 * Modifier method
 * @param inputSstream An input stream from which this object will be read
 */
void Kmer::read(istream& inputStream) {
    
    getline(inputStream, _text, '\0');
}



char& Kmer::operator[](int index) {

    return _text.at(index);
}

const char& Kmer::operator[](int index) const {

    return _text.at(index);
}

ostream &operator<<(std::ostream& os, const Kmer &kmer) {

    os << kmer._text;
    return os;
}

istream& operator>>(std::istream& is, Kmer& kmer) {

    string text;
    is >> text;
    kmer._text = text;
    return is;
}

bool operator>(const Kmer& Kmer1, const Kmer &Kmer2) {

    return (Kmer1._text < Kmer2._text);
}

bool operator==(const Kmer& Kmer1, const Kmer &Kmer2) {

    return (Kmer1._text == Kmer2._text);
}

bool operator!=(const Kmer& Kmer1, const Kmer &Kmer2) {
    return (!(Kmer1 == Kmer2));
}

bool operator<(const Kmer& Kmer1, const Kmer &Kmer2) {

    return (Kmer2 > Kmer1);
}

bool operator<=(const Kmer& Kmer1, const Kmer &Kmer2) {
    return (!(Kmer1 > Kmer2));
}

bool operator>=(const Kmer& Kmer1, const Kmer &Kmer2) {
    return (!(Kmer1 < Kmer2));
}




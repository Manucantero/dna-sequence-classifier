/*
 * Metodología de la Programación: Kmer1
 * Curso 2023/2024
 */

/** 
 * @file KmerFreq.cpp
 * @author Manuel Cantero Ruiz <manucantero@correo.ugr.es>
 * 
 * Created on 11 March 2024
 */

#include "KmerFreq.h"
#include <iostream>

using namespace std;

KmerFreq::KmerFreq() {

    _kmer = Kmer();
    _frequency = 0;
}

const Kmer& KmerFreq::getKmer() const {

    return (_kmer);
}

int KmerFreq::getFrequency() const {

    return (_frequency);
}

void KmerFreq::setKmer(const Kmer& kmer) {

    _kmer = kmer;
}

void KmerFreq::setFrequency(int frequency) {

    if (frequency < 0) { // k is less or equal than zero
        throw out_of_range(// composed string
                string("void KmerFreq::setFrequency(int frequency): ") +
                " frequency is less or equal than zero ");
    }
    _frequency = frequency;
}

string KmerFreq::toString() const {

    return (_kmer.toString() + " " + to_string(_frequency));
}

std::istream& operator>>(std::istream& is, KmerFreq& kmerFreq) {

    Kmer kmer;
    int freq;

    is >> kmer >> freq;

    kmerFreq._kmer = kmer;
    kmerFreq._frequency = freq;

    return is;
}

std::ostream &operator<<(std::ostream& os, const KmerFreq &kmerFreq) {

    os << kmerFreq._kmer << endl << kmerFreq._frequency;
    return os;
}

bool operator>(const KmerFreq& kmerFreq1, const KmerFreq &kmerFreq2) {

    bool mayor = false;

    if (kmerFreq1._frequency > kmerFreq2._frequency) {
        mayor = true;
    } else {

        bool iguales = kmerFreq1._frequency == kmerFreq2._frequency;
        if (iguales) {
            mayor = kmerFreq1._kmer < kmerFreq2._kmer;
        }
    }
    return mayor;
}

bool operator==(const KmerFreq& kmerFreq1, const KmerFreq &kmerFreq2) {

    return ((kmerFreq1._frequency == kmerFreq2._frequency) &&
            (kmerFreq1._kmer == kmerFreq2._kmer));
}

bool operator<(const KmerFreq& kmerFreq1, const KmerFreq &kmerFreq2) {

    return (kmerFreq2 > kmerFreq1);
}

bool operator!=(const KmerFreq& kmerFreq1, const KmerFreq &kmerFreq2) {

    return (!(kmerFreq1 == kmerFreq2));
}

bool operator<=(const KmerFreq& kmerFreq1, const KmerFreq &kmerFreq2) {

    return (!(kmerFreq1 > kmerFreq2));
}

bool operator>=(const KmerFreq& kmerFreq1, const KmerFreq &kmerFreq2) {

    return (!(kmerFreq1 < kmerFreq2));
}

/**
 * @brief Writes this object to the given output stream. It first writes
 * the kmer of this object (using method Kmer::write(ostream&)) 
 * and them the bytes of the frequency (an int value) in binary format 
 * (using method ostream::write(const char* s, streamsize n))
 * Query method
 * @param outputStream An output stream where this object will be written
 */
void KmerFreq::write( std::ostream& outputStream) const {

    this->_kmer.write(outputStream);

    // Escribir la frecuencia en formato binario
    outputStream.write(reinterpret_cast<const char*> (&_frequency), sizeof (int));
}

/**
 * @brief Reads this object from the given input stream. It first reads
 * the Kmer of this object (using method Kmer::read(std::istream&) and 
 * then the bytes of the frequency (an int value) in binary format (using 
 * method istream::read(char* s, streamsize n))
 * Modifier method
 * @param inputStream An input stream from which this object will be read
 */
void KmerFreq::read(std::istream& inputStream) {

    this->_kmer.read(inputStream);

    // Leer la frecuencia en formato binario
    inputStream.read(reinterpret_cast<char*> (&_frequency), sizeof (int));
}

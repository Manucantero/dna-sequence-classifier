/*
 * Metodología de la Programación: Kmer2
 * Curso 2023/2024
 */

/** 
 * @file Profile.cpp
 * @author Manuel Cantero Ruiz <manucantero@correo.ugr.es>
 * 
 * Created on 25 March 2024, 11:00
 */


#include <fstream>         // Para trabajar con ficheros
#include <cstdlib>         // Para exit(), rand(), srand()
#include <ctime>           // Para time()
#include <fstream>

#include "Profile.h"

using namespace std;

const string Profile::MAGIC_STRING_T = "MP-KMER-T-1.0";
const string Profile::MAGIC_STRING_B = "MP-KMER-B-1.0";

Profile::Profile() {

    _profileId = "unknown";
    _size = 0;
    _capacity = INITIAL_CAPACITY;
    reservar(INITIAL_CAPACITY);
}

Profile::Profile(int size) {

    if (size < 0) { // size is not valid
        throw out_of_range(// composed string
                string("Profile::Profile(int size): size is not valid "));
    }
    _profileId = "unknown";
    _capacity = _size = size;
    reservar(size);
}

Profile::Profile(const Profile &orig) {

    liberar();
    reservar(orig._capacity);
    this->_size = orig._size;
    this->_capacity = orig._capacity;
    this->_profileId = orig._profileId;
    copiar(_vectorKmerFreq, orig._vectorKmerFreq, orig._capacity);
}

Profile::~Profile() {
    liberar();
}

Profile& Profile::operator=(Profile orig) {

    if (&orig != this) {
        liberar();
        reservar(orig._capacity);
        this->_size = orig._size;
        this->_capacity = orig._capacity;
        this->_profileId = orig._profileId;
        copiar(_vectorKmerFreq, orig._vectorKmerFreq, orig._size);
    }
    return *this;
}

const std::string& Profile::getProfileId() const {

    return _profileId;
}

void Profile::setProfileId(const std::string& id) {

    _profileId = id;
}

const KmerFreq& Profile::at(int index) const {

    if ((index < 0) || index >= _size) { // k is less or equal than zero
        throw out_of_range(// composed string
                string("const KmerFreq& Profile::at(int index) const: size is not valid) "));
    }
    return _vectorKmerFreq[index];
}

KmerFreq& Profile::at(int index) {

    if ((index < 0) || index >= _size) {
        throw out_of_range(// composed string
                string("KmerFreq& Profile::at(int index):size is not valid) "));
    }

    return _vectorKmerFreq[index];
}

int Profile::getSize()const {

    return _size;
}

int Profile::getCapacity() const {

    return _capacity;
}

double Profile::getDistance(const Profile& otherProfile) const {

    if (_size == 0 || otherProfile.getSize() == 0) {

        throw invalid_argument(// composed string
                string("double Profile::getDistance(const Profile& otherProfile) const: Kmer is empty"));
    }

    double sum = 0.0;
    double denom = _size * otherProfile.getSize();

    for (int i = 0; i < _size; i++) {

        int pos = otherProfile.findKmer(_vectorKmerFreq[i].getKmer());

        if (pos == -1) {
            pos = otherProfile.getSize();
        }

        sum += abs(i - pos);
    }
    return (sum / denom);
}

int Profile::findKmer(const Kmer& kmer, int initialPos, int finalPos) const {

    int pos = -1;
    if ((initialPos >= 0 && initialPos < _size)&& (finalPos >= 0 && finalPos < _size)) {

        if (initialPos > finalPos) { // Check if positions are correct
            int aux = finalPos;
            finalPos = initialPos;
            initialPos = aux;
        }
        bool found = false;
        for (int i = initialPos; i <= finalPos && !found; i++) {
            if (_vectorKmerFreq[i].getKmer() == kmer) {
                pos = i;
                found = true;
            }
        }
    }
    return pos;
}

int Profile::findKmer(const Kmer& kmer)const {

    return (findKmer(kmer, 0, _size - 1));
}

std::string Profile::toString() const {

    string lines;
    if (_size > 0) {

        for (int i = 0; i < _size; i++) {
            lines += _vectorKmerFreq[i].toString() + "\n";
        }
    }

    return (_profileId + "\n" + to_string(_size) + "\n" + lines);
}

void Profile::sort() {

    int max;
    int pos_max = -1;

    for (int left = 0; left < _size - 1; left++) {
        max = _vectorKmerFreq[left].getFrequency();
        pos_max = left;
        for (int i = left + 1; i <= _size - 1; i++) {
            if ((_vectorKmerFreq[i].getFrequency() > max) ||
                    (_vectorKmerFreq[i].getFrequency() == max &&
                    _vectorKmerFreq[i].getKmer() >
                    _vectorKmerFreq[pos_max].getKmer())) {

                max = _vectorKmerFreq[i].getFrequency();
                pos_max = i;
            }
        }
        SwapElementsArrayKmerFreq(left, pos_max);

    }
}

void Profile::save(const char fileName[], char mode) const {

    if (mode != 't' && mode != 'b') {
        throw invalid_argument(// composed string
                string("void Profile::save(const char fileName[], char mode) const: Mode is not valid"));
    }

    if (mode == 'b') {

        ofstream fout("perfil.dat", ios::out | ios::binary); // Apertura del fichero (escritura)
        if (!fout) {
            throw ios_base::failure(// composed string
                    string("void Profile::save(const char fileName[]) const : Given file cannot be ")+
                    ("opened or an error has ocurred while reading from the file"));
        }

        fout << MAGIC_STRING_B << endl;
        fout << _profileId << endl;
        fout << _size << endl;

        for (int i = 0; i < _size; i++) {
            this->_vectorKmerFreq[i].write(fout);
            fout << endl;
        }
        fout.close(); // Cerramos fichero

    } else {

        ofstream fout(fileName); // Apertura del fichero (escritura)
        if (!fout) {
            throw ios_base::failure(// composed string
                    string("void Profile::save(const char fileName[]) const : Given file cannot be ")+
                    ("opened or an error has ocurred while reading from the file"));
        }

        fout << MAGIC_STRING_T << endl << toString(); // Escribimos en el fichero
        fout.close(); // Cerramos fichero
    }
}

void Profile::load(const char fileName[]) { //Cambiar excepcion y poner redimensionar

    _size = 0;
    ifstream finp(fileName);
    liberar();

    if (!finp) { // Apertura con éxito
        throw ios_base::failure(// composed string
                string("void Profile::load(const char fileName[]): ")+
                ("Given file cannot be opened or an error has ocurred ")+
                ("while reading from the file"));
    }
    //Read magic string
    string magic_string;
    finp >> magic_string;

    if (magic_string != MAGIC_STRING_T && magic_string != MAGIC_STRING_B) {
        throw invalid_argument(// composed string
                string("void Profile::load(const char fileName[]): Invalid magic string ")+
                ("is found in the given file)"));
    }
    finp.ignore();

    // Reservamos memoria
    ////

    if (magic_string == MAGIC_STRING_T) {
        finp >> *this;
    } else {
        int n_kmers;
        string id;
        getline(finp, id);
        _profileId = id;
        finp >> n_kmers;

        liberar();
        reservar(n_kmers);

        for (int i = 0; i < n_kmers; i++) {
            this->_vectorKmerFreq[i].read(finp);
        }
    }
    finp.close(); // Cerramos el fichero
}

void Profile::append(const KmerFreq& kmerFreq) {

    int pos = findKmer(kmerFreq.getKmer());

    if (pos == -1) {

        if (_size == _capacity) {
            redimensionar(BLOCK_SIZE);
        }

        _vectorKmerFreq[_size] = kmerFreq;
        _size++;

    } else {
        _vectorKmerFreq[pos].setFrequency(_vectorKmerFreq[pos].getFrequency()
                + kmerFreq.getFrequency());
    }
}

void Profile::normalize(const std::string& validNucleotides) {

    for (int i = 0; i < _size; i++) {

        Kmer aux(_vectorKmerFreq[i].getKmer());
        aux.normalize(validNucleotides);
        _vectorKmerFreq[i].setKmer(aux);
    }
    // Loop to traverse the kmers in array from position 1 to position _size-1
    int j = 1;
    while (j < _size) {

        // index = Position of array[i].getKmer() in the subarray that begins
        //         at position 0 and ends at position i-1

        int index = findKmer(_vectorKmerFreq[j].getKmer(), 0, j - 1);

        // If array[i].getKmer() was found in the the subarray from 0 to i-1 
        // Accumulate the frequencies of the kmers at positions 
        //    index and i in the kmer at position index
        // Delete from the arrhay, the kmer at position i 

        if (index != -1) {
            int freq = _vectorKmerFreq[j].getFrequency() + _vectorKmerFreq[index].getFrequency();
            _vectorKmerFreq[index].setFrequency(freq);
            deletePos(j);
        } else {
            j++;
        }
    }
}

void Profile::deletePos(int pos) {

    if ((pos < 0 || pos >= _size)) { // out_of_range
        throw out_of_range(// composed string
                string("void Profile::deletePos(int pos):Pos is out of range)"));
    }
    for (int i = pos; i < _size - 1; i++) {
        _vectorKmerFreq[i] = _vectorKmerFreq[i + 1];
    }
    _size--;
}

void Profile::zip(bool deleteMissing, int lowerBound) {

    int i = 0;
    while (i < _size) {

        if ((_vectorKmerFreq[i].getFrequency() <= lowerBound) ||
                ((deleteMissing) && ((int) _vectorKmerFreq[i].getKmer().toString().find(Kmer::MISSING_NUCLEOTIDE) != -1))) {

            deletePos(i);
        } else {
            i++;
        }
    }
}

void Profile::join(const Profile& profile) {

    for (int i = 0; i < profile.getSize(); i++) {

        append(profile.at(i));
    }

}

void Profile::SwapElementsArrayKmerFreq(int first, int second) {


    if (((first < 0 || first >= _size) || ((second < 0 || second >= _size)))) { // out_of_range

        throw out_of_range(// composed string
                string("void Profile::SwapElementsArrayKmerFreq(int first, int second): ") +
                " Positions are out of range  ");
    }

    KmerFreq aux;

    aux = _vectorKmerFreq[first];
    _vectorKmerFreq[first] = _vectorKmerFreq[second];
    _vectorKmerFreq[second] = aux;
}

void Profile::reservar(int size) {

    _vectorKmerFreq = nullptr;
    _vectorKmerFreq = new KmerFreq[size];

}

void Profile::redimensionar(int aumento) {

    KmerFreq * aux = new KmerFreq [_capacity + aumento];
    copiar(aux, _vectorKmerFreq, _size);
    delete[] _vectorKmerFreq;
    _vectorKmerFreq = aux;
    _capacity += aumento;
}

void Profile::liberar() {

    delete[] _vectorKmerFreq;
    _vectorKmerFreq = nullptr;
    _size = 0;
    _capacity = 0;
}

void Profile::copiar(KmerFreq * destino, const KmerFreq * v, int util) {

    for (int i = 0; i < util; i++) {
        destino[i] = v[i];
    }
}

std::ostream& operator<<(std::ostream& os, const Profile &profile) {


    os << profile._profileId << endl << profile._size << endl;

    for (int i = 0; i < profile._size; i++) {

        os << profile._vectorKmerFreq[i];
    }
    return os;
}

/**
 * @brief Overloading of the stream extraction operator for Profile class.
 * Note that this operator should remove any Kmer-frequency pairs that 
 * the argument Profile object previously contained.
 * @throw std::out_of_range Throws a std::out_of_range if the number of kmers
 * read from the file is negative.
 * @param is The input stream to be used. Output parameter
 * @param profile the Profile object. Output parameter
 * @return @p is A reference to the input stream
 */
std::istream& operator>>(std::istream& is, Profile& profile) {

    string id;
    int n_kmeros;
    getline(is, id);

    is >> n_kmeros;
    if (n_kmeros < 0) {
        throw out_of_range(// composed string
                string("void Profile::load(const char fileName[]): ") +
                " The number of kmers read from the given file is invalid ");
    }
    profile._profileId = id;
    KmerFreq k;

    for (int i = 0; i < n_kmeros; i++) {
        is >> k;
        profile += k;
    }
    return is;
}

const KmerFreq& Profile::operator[](int index) const {

    return (this->_vectorKmerFreq[index]);
}

KmerFreq& Profile::operator[](int index) {

    return (this->_vectorKmerFreq[index]);
}

Profile& Profile::operator+=(const KmerFreq& kmerFreq) {

    append(kmerFreq);
    return *this;
}

Profile& Profile::operator+=(const Profile& profile) {

    for (int i = 0; i < profile.getSize(); i++) {

        *this += (profile._vectorKmerFreq[i]);
    }
    return *this;
}

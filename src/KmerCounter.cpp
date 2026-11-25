/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/** 
 * @file KmerCounter.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include <fstream>

#include "KmerCounter.h"

using namespace std;

/**
 * DEFAULT_VALID_NUCLEOTIDES is a c-string that contains the set of characters
 * that will be considered as valid nucleotides. 

 * The constructor of the class KmerCounter uses this c-string as a 
 * default parameter. It is possible to use a different c-string if that
 * constructor is used with a different c-string
 */
const char* const KmerCounter::DEFAULT_VALID_NUCLEOTIDES = "ACGT";

/**
 * @brief Constructor of the class. This object uses a matrix of 
 * integers that contains the frequency for each kmer. Each element of the 
 * matrix will be initialized with 0. The constructor also initializes the
 * private data _k, _validNucleotides and _allNucleotides. _allNucleotides
 * is initialized with the characters 
 * Kmer::MISSING_NUCLEOTIDE + validNucleotides
 * @param k The number of nucleotides in each kmer. Input parameter
 * @param validNucleotides The set of nucleotides (characters) that are 
 * considered as part of a kmer. Input parameter
 */
KmerCounter::KmerCounter(int k, const string &validNucleotides) {

    _k = k;
    _validNucleotides = validNucleotides;
    _allNucleotides = Kmer::MISSING_NUCLEOTIDE + validNucleotides;

    allocate(getNumRows(), getNumCols());
    initFrequencies();
}

/**
 * @brief Copy constructor
 * @param orig the KmerCounter object used as source for the copy. 
 * Input parameter
 */
KmerCounter::KmerCounter(const KmerCounter& orig) {

    this->_k = orig._k;
    this-> _validNucleotides = orig._validNucleotides;
    this-> _allNucleotides = orig._allNucleotides;

    allocate(orig.getNumRows(), orig.getNumCols());
    this->copy(orig);
}

/**
 * @brief Destructor
 */
KmerCounter::~KmerCounter() {

    deallocate();
}

/**
 * @brief Returns the number of nucleotides that can be part of a
 * kmer, that is, the number of characters in the private data
 * _allNucleotides. For example, if "_ACGT" are the set of all nucleotides, 
 * then this method will return 5.
 * Query method
 * @return The number of nucleotides that can be part of a kmer
 */
int KmerCounter::getNumNucleotides() const {

    return _allNucleotides.size();
}

/**
 * @brief Returns the number of nucleotides in each kmer
 * Query method
 * @return The number of nucleotides in each kmer
 */
int KmerCounter::getK() const {

    return _k;
}

/**
 * @brief Returns the number of different kmers that can be built using
 * @p _k nucleotides (including the missing nucleotide)
 * Query method
 * @return The number of different kmers that can be built using
 * @p _k nucleotides 
 */
int KmerCounter::getNumKmers() const {

    // n^k
    int n = getNumNucleotides();

    for (int i = 1; i < getK(); i++) {

        n *= n;
    }
    return n;
}

/**
 * @brief Gets the number of kmers with a frequency greater than 0
 * Query method
 * @return the number of kmers with a frequency greater than 0
 */
int KmerCounter::getNumberActiveKmers() const {

    int active = 0;

    for (int f = 0; f < getNumRows(); f++) {
        for (int c = 0; c < getNumCols(); c++) {
            if ((*this)(f, c) > 0) {
                active++;
            }
        }
    }
    return active;
}

/**
 * @brief Obtains a string with the following content:
 * - In the first line, the content of the private data _allNucleotides and
 * the value of _k (number of nucleotides in each kmer) separated by a 
 * whitespace.
 * - In the following lines, each one of the rows in the frequency matrix
 * (frequencies separated by a whitespace).
 * Query method
 * @return A string with the content of this object
 */
std::string KmerCounter::KmerCounter::toString() const {
    string outputString = _allNucleotides + " " + to_string(_k) + "\n";

    for (int row = 0; row<this->getNumRows(); row++) {
        for (int col = 0; col<this->getNumCols(); col++) {
            outputString += to_string((*this)(row, col)) + " ";
        }
        outputString += "\n";
    }

    return outputString;
}

/**
 * @brief Increases the current frequency of the given kmer using the value 
 * provided by @p frequency. If the argument @p frequency is not 
 * provided, then 1 is added to the current frequency of the kmer.
 * Modifier method
 * @throw std::invalid_argument This method throws an 
 * std::invalid_argument exception if the given kmer contains any
 * invalid nucleotide.
 * @param kmer The kmer in which the frequency will be modified. 
 * Input parameter
 * @param frequency The quantity that will be added to the current 
 * frequency. Input parameter
 */
void KmerCounter::increaseFrequency(const Kmer &kmer, int frequency) {

    int row = 0;
    int column = 0;
    getRowColumn(kmer, row, column);


    if (row == -1 || column == -1) { // size is not valid
        throw invalid_argument(// composed string
                string(" void KmerCounter::increaseFrequency(const Kmer &kmer, int frequency): kmer contains any invalid nucleotide."));
    }

    (*this)(row, column) += frequency;

}

/**
 * @brief Overloading of the assignment operator.
 * Modifier method
 * @param orig the KmerCounter object used as source for the assignment.
 * Input parameter
 * @return A reference to this object
 */
KmerCounter& KmerCounter::operator=(const KmerCounter& orig) {

    if (&orig != this) {

        this->deallocate();

        this->_k = orig._k;
        this->_validNucleotides = orig._validNucleotides;
        this->_allNucleotides = orig._allNucleotides;

        // Reservamos memoria
        allocate(orig.getNumRows(), orig.getNumCols());
        this->copy(orig);
    }
    return *this;
}

/**
 * @brief Overloading of the operator +=. It increases the current 
 * frequencies of the kmers of this object with the frequencies of the 
 * kmers of the given object.
 * Modifier method
 * @param kc a KmerCounter object. Input parameter 
 * @throw std::invalid_argument This method throws an 
 * std::invalid_argument exception if the given argument @p kc has a different
 * set of nucleotides or a different K (number of nucleotides in
 * kmers).
 * @return A reference to this object
 */
KmerCounter& KmerCounter::operator+=(const KmerCounter& kc) {

    if ((this->_validNucleotides != kc._validNucleotides) || (this->_k != kc._k)) { // size is not valid
        throw invalid_argument(// composed string
                string(" KmerCounter& KmerCounter::operator+=(const KmerCounter& kc)") +
                "given argumentkc has a different set of nucleotides or a different K ");
    }

    for (int row = 0; row < kc.getNumRows(); row++) {
        for (int col = 0; col < kc.getNumCols(); col++) {
            (*this)(row, col) += kc(row, col);
        }
    }
    return *this;
}

/**
 * @brief Reads the given text file and calculates the frequencies of each 
 * kmer in that file. This method normalizes each found Kmer and then
 * sum 1 at the corresponding element of the frequency matrix.
 * This method sets to zero the frequency of each kmer before starting to 
 * calculate frequencies. In this way, if this method is called twice 
 * consecutively, then this KmerCounter will contain only the frecuencies
 * calculated in the last call.
 * Modifier method
 * @param fileName The name of the file to process. Input parameter
 * @throw std::ios_base::failure Throws a std::ios_base::failure if the
 * given file cannot be opened
 */
void KmerCounter::calculateFrequencies(const char* fileName) {

    initFrequencies();
    ifstream input;
    input.open(fileName);
    if (!input) {
        throw ios_base::failure(// composed string
                string("void KmerCounter::calculateFrequencies(const char* fileName): ")+
                ("Given file cannot be opened or an error has ocurred ")+
                ("while reading from the file"));
    }


    string geneticsequence;
    input >> geneticsequence;
    int num = geneticsequence.size() - _k + 1;

    if (num < 0) {
        num = 0;
    }

    for (int i = 0; i < num; i++) {
        Kmer kmero = Kmer(geneticsequence.substr(i, _k));
        kmero.normalize(_validNucleotides);
        increaseFrequency(kmero);
    }
}

Profile KmerCounter::toProfile() const {

    //Calculamos el numero de frecuencias positivas para saber cuando parar el bucle
    int n_active = getNumberActiveKmers();
    Profile profile;

    // Creamos una variable auxiliar para ayudarnos a inicializar Profile
    KmerFreq tmp;

    bool seguir = true;
   
    for (int r = 0; r < this->getNumRows() && seguir; r++) {
        for (int c = 0; c < this->getNumCols(); c++) {
           
            if ((*this)(r, c) > 0) {
                tmp.setKmer(getKmer(r, c).toString());
                tmp.setFrequency((*this)(r, c));
                profile += tmp;
            }
            
            if(profile.getSize() == n_active) {
                seguir = false;
            }
        }
    }
 

    return profile;
}

/**
 * @brief Returns the numbers of rows of the matrix in this object.
 * Query method
 * @return The numbers of rows of the matrix in this object
 */
int KmerCounter::getNumRows() const {

    // n^((k+1)/2)

    int pow = (_k + 1) / 2;
    int l_f = getNumNucleotides();

    for (int i = 1; i < pow; i++) {
        l_f *= l_f;
    }
    return l_f;

}

/**
 * @brief Returns the numbers of columns of the matrix in this object.
 * Query method
 * @return The numbers of columns of the matrix in this object
 */
int KmerCounter::getNumCols() const {

    //  n^(k-(k+1)/2)

    int pow = _k - ((_k + 1) / 2);
    int l_c = getNumNucleotides();
    for (int i = 1; i < pow; i++) {
        l_c *= l_c;
    }

    return l_c;
}

int KmerCounter::getIndex(const std::string& kmer) const {
    int index = 0;
    int base = 1;

    for (size_t i = 0; i < kmer.size(); i++) {
        size_t pos = _allNucleotides.find(kmer[kmer.size() - i - 1]);
        if (pos == string::npos)
            return -1;
        index += pos * base;
        base *= _allNucleotides.size();
    }
    return index;
}

string KmerCounter::getInvertedIndex(int index, int nCharacters) const {
    string result(nCharacters, Kmer::MISSING_NUCLEOTIDE);

    for (int i = result.size(); i > 0; i--) {
        result[i - 1] = _allNucleotides[index % _allNucleotides.size()];
        index = index / _allNucleotides.size();
    }
    return result;
}

/**
 * @brief Obtains the row and column for the frequency of the given 
 * argument Kmer. -1 will be assigned to row or column if some of the
 * nucleotides in @p kmer does not belong to the set _allNucleotides
 * Query method
 * @param kmer A Kmer. Input parameter
 * @param row The row in the frequency matrix corresponding to Kmer @p kmer.
 * Output parameter
 * @param column The column in the frequency matrix corresponding to 
 * Kmer @p kmer. Output parameter
 */
void KmerCounter::getRowColumn(const Kmer& kmer, int &row, int &column) const {

    /*
    int pos = getIndex(kmer.toString());
   // cout << " LA POSICION ES " << pos << endl << endl;
    int ncol = getNumCols();

    // Dividimos la posicion en la que está Kmer hasta que salga cero.
    // En el momento que sale cero, estaríamos en esa fila ( restamos uno porque
    // empezamos en cero
    int i = 0;
    while ((pos / (ncol * i)) != 0) {
        i++;
    }

    // Restamos para obtener en que posicion de la fila está y restamos uno 
    // porque empieza en cero
    column = pos - (ncol * i) - 1;
    row = i - 1;
     */


    string kmero = kmer.toString();
    string str_row = kmero.substr(0, (_k + 1) / 2);
    string str_col = kmero.substr((_k + 1) / 2, _k - (_k + 1) / 2);
    row = getIndex(str_row);
    column = getIndex(str_col);
}

/**
 * @brief Returns the Kmer that is defined by the provided row and column 
 * Query method
 * @param row A row in the frequency matrix (int value from 0 to nRows-1).
 * Input parameter
 * @param column A column in the frequency matrix (int value from 0 to 
 * nColumns-1). Input parameter
 * @throw std::invalid_argument This method throws an 
 * std::invalid_a rgument exception if the given row or column are out
 * of the correct bounds
 * @return the Kmer that is defined by the provided row and column 
 */
Kmer KmerCounter::getKmer(int row, int column) const {

    if ((row < 0 || row > getNumRows() || column < 0 || column > getNumCols())) { // size is not valid
        throw invalid_argument(// composed string
                string(" void KmerCounter::increaseFrequency(const Kmer &kmer, int frequency): kmer contains any invalid nucleotide."));
    }
    int If = (_k + 1) / 2;
    int Ic = _k - If;

    string kmero = getInvertedIndex(row, If) + getInvertedIndex(column, Ic);
    return ( Kmer(kmero));
}

/**
 * @brief Sets the frequency of each kmer to 0, that is, it fills with 0 the 
 * matrix of frequencies
 * Modifier method
 */
void KmerCounter::initFrequencies() {

    for (int f = 0; f < getNumRows(); f++) {
        for (int c = 0; c < getNumCols(); c++) {
            (*this)(f, c) = 0;
        }
    }
}

/**
 * @brief Overloading of the () operator to access to the element at a 
 * given position.
 * Query method
 * @param row Row of the element. Input parameter
 * @param column Column of the element. Input parameter
 * @return A const reference to the element at the given position
 */
const int &KmerCounter::operator()(int row, int column) const {
    return _frequency[row][column];
}

/**
 * @brief Overloading of the () operator to access to the element at a 
 * given position.
 * Query and modifier method
 * @param row Row of the element. Input parameter
 * @param column Column of the element. Input parameter
 * @return A reference to the element at the given position
 */
int & KmerCounter::operator()(int row, int column) {
    return _frequency[row][column];
}

void KmerCounter::allocate(int nrows, int ncols) {

    _frequency = new int*[nrows];
    for (int i = 0; i < nrows; ++i) {
        _frequency[i] = new int[ncols];
    }
}

void KmerCounter::deallocate() {

    if (_frequency != nullptr) {
        for (int i = 0; i < this->getNumRows(); ++i) {
            delete[] _frequency[i];
        }

        delete[] _frequency;
    }
}

// Pregunatr si está bien

void KmerCounter::copy(const KmerCounter& orig) {

    if (&orig != this) {

        for (int row = 0; row < orig.getNumRows(); row++) {

            for (int col = 0; col < orig.getNumCols(); col++) {

                (*this)(row, col) = orig(row, col);
            }
        }
    }
}

/*
bool KmerCounter::kmer_verifier (const Kmer& kmer) {
    
    bool valid = true;
    int i = 0;
    while ((i < kmer.getK()) && verifica) {

        if (!IsValidNucleotide(kmer[i], _allNucleotides)) {
            valid = false;
        } else {
            i++;
        }
    }
    
    return valid;
}
 * 
 * */
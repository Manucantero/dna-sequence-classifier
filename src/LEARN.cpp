/*
 * Metodología de la Programación: Kmer5
 * Curso 2023/2024
 */

/*
 * @file LEARN.cpp
 * @author Silvia Acid Carrillo <acid@decsai.ugr.es>
 * @author Andrés Cano Utrera <acu@decsai.ugr.es>
 * @author Luis Castillo Vidal <L.Castillo@decsai.ugr.es>
 * @author Javier Martínez Baena <jbaena@ugr.es>
 * 
 * Created on 22 December 2023, 10:00
 */

#include <fstream>
#include <cstring>
#include <iostream>
#include "KmerCounter.h"
using namespace std;

/**
 * Shows help about the use of this program in the given output stream
 * @param outputStream The output stream where the help will be shown (for example,
 * cout, cerr, etc) 
 */
void showEnglishHelp(ostream& outputStream) {
    outputStream << "ERROR in LEARN parameters" << endl;
    outputStream << "Run with the following parameters:" << endl;
    outputStream << "LEARN [-t|-b] [-k kValue] [-n nucleotidesSet] [-p profileId] [-o outputFilename] <file1.dna> [<file2.dna> <file3.dna> .... ]" << endl;
    outputStream << endl;
    outputStream << "Parameters:" << endl;
    outputStream << "-t|-b: text mode or binary mode for the output file (-t by default)" << endl;
    outputStream << "-k kValue: number of nucleotides in a kmer (5 by default)" << endl;
    outputStream << "-n nucleotidesSet: set of possible nucleotides in a kmer (ACGT by default). "
            << "Note that the characters should be provided in uppercase" << endl;
    outputStream << "-p profileId: profile identifier (unknown by default)" << endl;
    outputStream << "-o outputFilename: name of the output file (output.prf by default)" << endl;
    outputStream << "<file1.dna> <file2.dna> <file3.dna> ....: names of the input files (at least one is mandatory)" << endl;
    outputStream << endl;
    outputStream << "This program learns a profile model from a set of " <<
            "input DNA files <file1.dna> <file2.dna> <file3.dna> ...." << endl;
    outputStream << endl;
}

/**
 * This program learns a Profile model from a set of input DNA files (file1.dna,
 * file2.dna, ...). The learned Profile object is then zipped (kmers with any 
 * missing nucleotide or with frequency equals to zero will be removed) 
 * and ordered by frequency and saved in 
 * the file outputFilename (or output.prf if the output file is not provided).
 * 
 * Running sintax:
 * > LEARN [-t|-b] [-k kValue] [-n nucleotidesSet] [-p profileId] [-o outputFilename] <file1.dna> [<file2.dna> <file3.dna> ....]
 * 
 * Running example:
 * > LEARN -k 2 -p bug -o /tmp/unknownACGT.prf ../Genomes/unknownACGT.dna
 * 
 * > cat /tmp/unknownACGT.prf
MP-KMER-T-1.0
bug
7
GG 2
AC 1
AG 1
AT 1
CC 1
GA 1
TA 1
 * 
 * @param argc The number of command line parameters
 * @param argv The vector of command line parameters (cstrings)
 * @return 0 If there is no error; a value > 0 if error
 */
int main(int argc, char *argv[]) {

    // Process the main() arguments
    if (argc < 2) {
        showEnglishHelp(cerr);
        exit(EXIT_FAILURE);
    }

    char mode = 't';
    string id = "unknown";
    int k = 5;
    string nucl = "ACGT";
    string otp = "output.prf";
    int pos = 1;
    bool fin = false;

    while ((pos < argc) && (!fin)) {
        if (argv[pos][0] == '-') {
            if (strlen(argv[pos]) == 2) { // Se necesita una letra
                if (argv[pos][1] == 't') {
                    mode = 't';
                    pos++;
                } else if (argv[pos][1] == 'b') {
                    mode = 'b';
                    pos++;
                } else if (argv[pos][1] == 'p' && (pos + 1 < argc)) {
                    id = argv[pos + 1];
                    pos += 2;
                } else if (argv[pos][1] == 'k' && (pos + 1 < argc)) {
                    k = atoi(argv[pos + 1]);
                    pos += 2;
                } else if (argv[pos][1] == 'n' && (pos + 1 < argc)) {
                    nucl = argv[pos + 1];
                    pos += 2;
                } else if (argv[pos][1] == 'o' && (pos + 1 < argc)) {
                    otp = argv[pos + 1];
                    pos += 2;

                } else { // option after - is not valid 
                    showEnglishHelp(cerr);
                    exit(EXIT_FAILURE);
                    fin = true;
                }

            } else { // - was only introduced

                showEnglishHelp(cerr);
                exit(EXIT_FAILURE);
            }
        } else { // Empezamos a leer ficheros
            fin = true;
        }
    }
    
    // Si son iguales, no habría ficheros para leer
    if (pos == argc) {
        showEnglishHelp(cerr);
        exit(EXIT_FAILURE);
    }
    
    // Loop to calculate the kmer frecuencies of the input genome files using 
    // a KmerCounter object
    KmerCounter base(k, nucl);
    KmerCounter aux(k, nucl);
    for (int i = pos; i < argc; i++) {
        aux.calculateFrequencies(argv[i]);
        base += aux;
    }
    // Obtain a Profile object from the KmerCounter object
    Profile profile = base.toProfile();
    
    profile.setProfileId(id);
    
    // Hasta aqui todo ok (. y me tiro)
    
    // Zip the Profile object
    profile.zip(true);

    // Sort the Profile object
    profile.sort();

   // cout << profile.toString();
    // Save the Profile object in the output file
    profile.save(otp.c_str(), mode);

    
    
    
}

/*
 if (pos != -1) { 

        if (pos < (_k - (_k + 1) / 2)) { // Can't divide by zero
            column = 0;
        } else {
            i = 1;
            while ((pos / (ncol * i)) != 0) {
                i++;
            }

            // Restamos para obtener en que posicion de la fila está y restamos uno 
            // porque empieza en cero
            column = pos - (ncol * i) - 1;
        }

        row = i - 1;
    }
    else { // Not found

        row = column = -1;
    }
 
 
 */
/*
 
 
    int i = 0;
    bool valid = true;

    while ((i < kmer.getK()) && valid) {

        if (!IsValidNucleotide(kmer[i], _allNucleotides)) {
            valid = false;
        } else {
            i++;
        }
    }

    // Si contiene un nucleotido invalido, devuelve -1
    if (valid) {

        int pos = getIndex(kmer.toString());
        int ncol = getNumCols();

        // Dividimos la posicion en la que está Kmer hasta que salga cero.
        // En el momento que sale cero, estaríamos en esa fila ( restamos uno porque
        // empezamos en cero
        i = 0;
        
        int len_c = _k - ((_k + 1) / 2);
        if (pos < len_c) {
            
            row = 0;
            column = pos;
        } else {
            i = 1;
             while ((pos / (ncol * i)) != 0) {
                i++;
            }

            // Restamos para obtener en que posicion de la fila está y restamos uno 
            // porque empieza en cero
            row = pos - (ncol * i) - 1;
             column = i - 1;
        }
        
  
    } else {

        row = column = -1;
    }

 
 */
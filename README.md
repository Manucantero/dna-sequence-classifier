# DNA Sequence Classifier — K-mer Frequency Profiling in C++

Command-line bioinformatics tool developed in C++ for analyzing genetic sequences, building normalized *k*-mer frequency profiles, and classifying unknown DNA samples against reference models using distance-based similarity metrics.

## Key Features

- **K-mer Extraction & Normalization:** Parses raw DNA sequence datasets (`A`, `C`, `G`, `T`), extracts contiguous substrings of length *k*, and handles canonical/invalid nucleotide filtering.
- **Frequency Profile Modeling:** Dynamically manages sorted *k*-mer frequency distributions using custom dynamic memory structures and operator overloading in C++.
- **Profile Training & Classification:** Generates reference statistical profiles (`.dat`) from known genetic datasets and classifies unknown DNA sequences by computing profile distance metrics.
- **Automated Testing & Validation:** Includes dedicated unit/integration tests (`/tests`) and automation scripts (`/scripts`) to verify memory safety and classification accuracy.

## Repository Structure

- **`/src` & `/include`:** Core C++ implementation (`.cpp`) and header files (`.h`) for *k*-mer representation, frequency profiles, file I/O, and classification logic.
- **`/tests`:** Test suite validating data structures, sorting, and classification outputs.
- **`/scripts`:** Automation scripts for batch training and evaluation.
- **`Makefile` & `perfil.dat`:** Build configuration and sample generated profile data.

## Tech Stack

- **Language:** C++ (C++11/14)
- **Build System:** GNU Make (`Makefile`) / `g++`
- **Core Concepts:** Dynamic Memory Management, Pointers & References, File I/O, Object-Oriented Design, Bioinformatics (*k*-mer profiling)

## How to Build & Run

1. Clone the repository:
   ```bash
   git clone [https://github.com/Manucantero/dna-sequence-classifier.git](https://github.com/Manucantero/dna-sequence-classifier.git)
   cd dna-sequence-classifier
   ```

2. Build the project using `make` (or compile directly with `g++`):
   ```bash
   make
   # Alternative direct compilation if Makefile relies on IDE paths:
   # g++ -O2 -std=c++14 -Iinclude src/*.cpp -o dna_classifier
   ```

3. Run the classifier with sample data or execute the test scripts:
   ```bash
   ./dna_classifier
   ```
   *(Note: Check `/scripts` or `Makefile` for the exact binary names if split into separate training and classification executables).*

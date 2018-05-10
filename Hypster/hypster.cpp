#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <vector>

using namespace std;


/*************************************
*  Forward function declarations
**************************************/

string parseFasta(string fastaFile);
string table_exclusion;
string upper;
string lower;
stringstream out;
int mapSize;
string fh;

map<string, long long int> initializeMap(string &genome, int nmerSize, string table_exclusion, int mapSize);
void getNmers(string &genome, int nmerSize, string table_exclusion, int mapSize, string fh);

/***********************************
*                MAIN
***********************************/
int main(int argc, char* argv[]) {


	if (argc < 4) {
		cerr << "\nUSAGE: getNmerFreqs.exe fasta_file Kmer_size <output_file>\n";
		return 0;
	}

	string fastaFile = argv[1];
	int nmerSize = atoi(argv[2]);
	fh = argv[3];


	int mapSize;


	table_exclusion = "+?/()-,[]#*1234567890BJOUXZbjouxz";  // Removable sequence digits, not of AA alphabet
	lower = "acdefghiklmnpqrstvwy";
	upper = "ACDEFGHIKLMNPQRSTVWY";

	mapSize = 20;


	cerr << "\n" << upper;
	cerr << "\nInput: " << fastaFile << endl
		<< "Nmer:  " << nmerSize << endl;

	string genome = parseFasta(fastaFile);

	cerr << "\nGenome length: " << genome.length() << endl << endl;


	getNmers(genome, nmerSize, table_exclusion, mapSize, fh);

	return 0;
}


/*******************************************
* Function reads in a FASTA file and stores
*  the genome in it as one long string
* The string is returned by this function.
********************************************/

string parseFasta(string F) {

	ifstream inputF;
	string line;
	string ret; // what we will return
	string end_char = "e";

	long long int N = 0;

	inputF.open(F.c_str(), ios::in); // open the for reading only, that's what ios::in means

	if (!inputF.is_open()) { // error checking
		cerr << "\nERROR: Unable to open '" << F << "'\nExiting now...\n\n";
		exit(0);
	}

	ret = ""; // prep variable for new sequence data to comes

	int i = 0;

	while (!inputF.eof()) {
		line = "";
		getline(inputF, line); // read in 1 line and assign it to the 'line' variable

		N = (signed)line.length(); // get the length of the current line.

		if (N <= 0) {
			continue;
		} // skip blank lines

		if (line.at(0) == '>') { // check to see if this is sequence or just the header line
								 //ret = ""; // prep variable for new sequence data to come
			ret += "*";
			continue; // skip header line, we won't be doing anything with it.
		}

		else { // should be sequence

			ret += line;
			i++;
		}
	}

	inputF.close();
	return ret;

}


/*=====================================================================
* Function records Kmers found in 'proteome' string and prints them out.
* Proteome is passed by reference to reduce memory usage.
* Kmers are defined using sliding window.
*=====================================================================*/
void getNmers(string &genome, int nmerSize, string table_exclusion, int mapSize, string fh) {

	ofstream file;
	
	map<string, long long int> nmerHash;
	map<string, long long int>::iterator iter;
	string curNmer;
	string findX;
	string find_s = "";

	long long int N = (signed)genome.length();

	clock_t t1;
	t1 = clock();
	nmerHash = initializeMap(genome, nmerSize, table_exclusion, mapSize);
	t1 = clock() - t1;

	long double maxHash = (long double)nmerHash.max_size();
	cerr << "\n";
	long double nHash = (long double)pow(20.0, nmerSize);
	long double lHash = (long double)nmerHash.size();

	cerr << nHash << "\n";
	cerr << maxHash << "\n";
	cerr << lHash << "\n";

	cerr << "\n\nCounting kmers\n";
	clock_t t2;
	t2 = clock();
	for (long long int i = 0; i < (N - nmerSize); i++) {

		curNmer = genome.substr(i, nmerSize);

		if (curNmer.find_first_of(table_exclusion) == string::npos) { //
			nmerHash[curNmer] += 1;
	
		}

	}

	// Determine the data set benchmark.
	t2 = clock() - t2;

	// final output to be printed out from here on
	out << "kmer\tfreq\n"; // header line
	for (iter = nmerHash.begin(); iter != nmerHash.end(); iter++) {
		out << iter->first << "\t" << iter->second << "\n";
	}
	out << endl;



	file.open(fh);
	file << "# Processing Time Step 1- finding words (sec): " << ((float)t1) / CLOCKS_PER_SEC << endl;
	file << "# Processing Time Step 2- counting words(sec): " << ((float)t2) / CLOCKS_PER_SEC << endl;
	file << "# Total AA in Proteome: " << N << endl;
	file << "# K-mer Length: " << nmerSize << endl;
	file << "# Total Word Count: " << lHash << endl;
	file << "# Total Possible Word Count: " << nHash << endl;
	file << "# Max Hash-MAP length: " << maxHash << endl;
	file << out.str();


	file.close();


	cerr << "seconds to process: " << ((float)t1) / CLOCKS_PER_SEC + ((float)t2) / CLOCKS_PER_SEC;
	cerr << "\n\nDone!\n";
}


/*=====================================================================
** Initialize map to Nmers
*======================================================================*/
map<string, long long int> initializeMap(string &genome, int nmer, string table_exclusion, int mapSize) {


	cerr << "\nProcessing kmers to load hash map\n";
	long long int N = (signed)genome.length() - nmer;

	string curStr = "";
	string find_s = "s";


	map<string, long long int> ret; // what we return

	long double LIMIT = pow(mapSize, nmer); // Compute how big the map should be (ie: length of kmers)
											// There are 4 nucleotide characters 4^(nmer) = number
											// of elements in the map
											//   for amino acids:20
	for (long long int i = 0; i <= N; i++) {
		//cerr << "\r" << i;
		curStr = genome.substr(i, nmer);

		if (N <= 0) { //SCOTT   from prior change in code. probabliy not needed anymore.
			continue;
		} // skip blank lines from multi fasta file

		if (curStr.find_first_of(table_exclusion) == string::npos) { 
			ret[curStr] = 0;
		}

		if ((unsigned)ret.size() == LIMIT) break; // map is ready, break out of loop
	}
	return ret;
}


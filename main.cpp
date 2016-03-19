#include "utils.h"
#include <fstream>
#include <chrono>


using namespace std;


int main(int argc, char ** argv){
	srand (time(NULL));
	string refFile("ecoliref.fa"),seq,ref;
	string readFile("colisimul.fa");
	ifstream in(refFile);
	getline(in,ref);
	getline(in,ref);
	uint k(31);
	cout<<"Counting kmers"<<endl;
	unordered_map<kmer,uint64_t> count(kmerCounting(ref, k));
	cout<<"filling"<<endl;
	unordered_map<kmer,vector<position>> kmer2pos;
	fillIndex(refFile, k, kmer2pos);
	cout<<"mapping"<<endl;
	auto startChrono=chrono::system_clock::now();
	mapReadFile(readFile,k,kmer2pos, ref);
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"mapping took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
	return 0;
}

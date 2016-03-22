#include "utils.h"
#include <fstream>
#include <chrono>


using namespace std;


int main(int argc, char ** argv){
	// readCounting("../Unitigs/SRR959239.clean.fa");
	srand(time(NULL));
	if(argc==1){
		cout<<"Usage : "<<endl;
		cout<<"-u : read file (minicoli.fa)"<<endl;
		cout<<"-x : reference file (ecoliref.fa)"<<endl;
		cout<<"-k : size of anchor (31)"<<endl;
		return 0;
	}
	string refFile("ecoliref.fa"),seq,ref,readFile("minicoli.fa");
	uint k(31);
	char c;
	while ((c = getopt (argc, argv, "u:k:x:")) != -1){
		switch(c){
			case 'u':
				readFile=optarg;
				break;
			case 'x':
				refFile=optarg;
				break;
			case 'k':
				k=stoi(optarg);
			break;
		}
	}
	ifstream in(refFile);
	getline(in,ref);
	getline(in,ref);
	cout<<"Counting "<<k<<"mers "<<endl;
	unordered_map<kmer,uint64_t> count(kmerCounting(ref, k));
	cout<<"Filling index of "<<refFile<<endl;
	unordered_map<kmer,vector<position>> kmer2pos;
	fillIndex(refFile, k, kmer2pos);
	cout<<"Mapping "<<readFile<<endl;
	auto startChrono=chrono::system_clock::now();
	mapReadFile(readFile,k,kmer2pos, ref);
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Mapping took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
	return 0;
}

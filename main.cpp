#include "utils.h"
#include "mapping.h"
#include <fstream>
#include <chrono>

#include <getopt.h>

using namespace std;


int main(int argc, char ** argv){
	// readCounting("../Unitigs/SRR959239.clean.fa");
	srand(time(NULL));
	if(argc==1){
		cout<<"Usage : "<<endl;
		cout<<"-u : read file (minicoli.fa)"<<endl;
		cout<<"-x : reference file (ecoliref.fa)"<<endl;
		cout<<"-k : size of anchor (31)"<<endl;
		cout<<"-m : max missmatch (2)"<<endl;
		return 0;
	}
	string refFile("ecoliref.fa"),seq,ref,readFile("minicoli.fa");
	uint k(31);
	uint maxMiss(2);
	char c;
	while ((c = getopt (argc, argv, "u:k:x:m:")) != -1){
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
			case 'm':
				maxMiss=stoi(optarg);
			break;
		}
	}
	ifstream in(refFile);
	getline(in,ref);
	getline(in,ref);
	// cout<<"Counting "<<k<<"mers "<<endl;
	// unordered_map<kmer,uint8_t> count(kmerCounting(ref, k));
	cout<<"Filling index of "<<refFile<<endl;
	unordered_map<kmer,vector<position>> kmer2pos;
	fillIndex(refFile, k, kmer2pos);
	cout<<"Mapping "<<readFile<<endl;
	auto startChrono=chrono::system_clock::now();
	uint nbread(mapReadFile(readFile,k,kmer2pos, ref,maxMiss));
	auto end=chrono::system_clock::now();auto waitedFor=end-startChrono;
	cout<<"Mapping took : "<<(chrono::duration_cast<chrono::seconds>(waitedFor).count())<<" sec"<<endl;
	cout<<"Throughout: "<< nbread/(1000*(chrono::duration_cast<chrono::seconds>(waitedFor).count()))<<"k read by second or "
	<< (nbread*3600/(1000000*(chrono::duration_cast<chrono::seconds>(waitedFor).count())))<<"M by hour"<<endl;
	return 0;
}

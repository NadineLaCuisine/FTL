#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include "utils.h"
#include "distances.h"


using namespace std;


void fillIndex(const string& refFile, const uint64_t k, unordered_map<kmer,vector<position>>& kmer2pos){
	string seq;
	ifstream readS(refFile);
	getline(readS,seq);
	getline(readS,seq);
	uint64_t i(0);
	minimizer kmerS(seq2intStranded((seq.substr(0,k))));
	minimizer kmerRC(rc(kmerS,k));
	minimizer kmer(min(kmerRC,kmerS));
	bool end(false);
	do{
		kmer2pos[kmer].push_back(i);
		if(seq[i+k]==':'){
			i+=k;
			do{++i;}while(seq[i]==':');
			kmerS=(seq2intStranded((seq.substr(i,k))));
			kmerRC=(rc(kmerS,k));
			kmer=(min(kmerRC,kmerS));
		}else if(i+k<seq.size()){
			updateMinimizer(kmerS, seq[i+k], k);
			updateMinimizerRC(kmerRC, seq[i+k], k);
			kmer=min(kmerRC,kmerS);
			++i;
		}else{
			end=true;
		}
	}while(!end);
}


uint mapRead(const  string& read,const uint64_t k, unordered_map<kmer,vector<position>>& kmer2pos, const string& ref,uint maxMiss,string& corrected){
	minimizer kmerS(seq2intStranded((read.substr(0,k))));
	minimizer kmerRC(rc(kmerS,k));
	minimizer kmer(min(kmerRC,kmerS));
	bool end(false),mapped(false);
	uint i(0);
    uint bestScore(6);
	do{
		if(kmer2pos.count(kmer)!=0){
			vector<position> positions(kmer2pos[kmer]);
			for(uint j(0);j<positions.size() and not mapped;++j){
				int possrt(positions[j]-i);
				if(possrt>=0){
                    uint score(distHammingIndel(read,ref.substr(possrt,read.size()+10),maxMiss));  // hamming
                    //~ uint score(nbMismatchesSW(read, ref.substr(possrt,read.size())));  // smith waterman
					if(score<maxMiss){
                        corrected=ref.substr(possrt,read.size());
					    bestScore=score;
					}
				}
			}
		}
		if(i+k<read.size()){
			updateMinimizer(kmerS, read[i+k], k);
			updateMinimizerRC(kmerRC, read[i+k], k);
			kmer=min(kmerRC,kmerS);
			++i;
		}else{
			end=true;
		}
	}while(not end);
	return bestScore;
}


uint mapReadFile(const string& readFile,const uint64_t k, unordered_map<kmer,vector<position>>& kmer2pos, const string& ref,uint maxMiss){
    ofstream mapped("mapped.fa"),notMapped("notMapped.fa");
	ifstream readS(readFile);
	string read,useless,comp,corrected,correctedRC;
	uint mappedRead(0),readNumber(0);
	while(!readS.eof()){
		getline(readS,useless);
		getline(readS,read);
		if(read.empty()){
			break;
		}
		++readNumber;
        uint score(mapRead(read,  k, kmer2pos, ref,maxMiss,corrected));
        uint scorerc(mapRead(reversecomplement(read),  k, kmer2pos, ref,maxMiss,correctedRC));
		if(min(score,scorerc)<maxMiss){
            ++mappedRead;
			if(score<scorerc){
                mapped<<useless<<endl<<corrected<<endl;
            }else{
                mapped<<useless<<endl<<correctedRC<<endl;
            }
		}else{
			notMapped<<useless<<endl<<correctedRC<<endl;
		}

	}
	cout<<"Reads: "<<readNumber<<endl;
	cout<<"Reads mapped: "<<mappedRead<<endl;
	cout<<"Percent Read mapped: "<<((10000*(double)mappedRead)/readNumber)/100<<"%"<<endl;
    return readNumber;
}

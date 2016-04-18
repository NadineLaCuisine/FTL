#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include <thread>

#include "utils.h"
#include "distances.h"
#include "mapping.h"


using namespace std;


ofstream mapped("mapped.fa"),notMapped("notMapped.fa");
ifstream readFile;
atomic<uint> mappedRead(0),readNumber(0),notMappedRead(0);
mutex lockOutFile,lockReadFile;


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


uint mapRead(const  string& read,const uint64_t k, const unordered_map<kmer,vector<position>>& kmer2pos, const string& ref,uint maxMiss,string& corrected){
	minimizer kmerS(seq2intStranded((read.substr(0,k))));
    minimizer kmerRC(rc(kmerS,k));
    minimizer kmer(min(kmerRC,kmerS));
    bool end(false),mapped(false);
    uint i(0);
    uint bestScore(maxMiss);
	vector<position> positions;
    do{
        if(kmer2pos.count(kmer)!=0){
            positions=(kmer2pos.at(kmer));
            for(uint j(0);j<positions.size() and not mapped;++j){
                int possrt(positions[j]-i);
                if(possrt>=0){
                    uint score(distHamming(read,ref.substr(possrt,read.size()+10),maxMiss));//hamming
    				// uint score(distHammingIndel(read,ref.substr(possrt,read.size()+10),maxMiss));  // hammingindel
    				//~ uint score(nbMismatchesSW(read, ref.substr(possrt,read.size())));  // smith waterman
                    if(score<bestScore){
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


void treatRead(const uint64_t k, const unordered_map<kmer,vector<position>>& kmer2pos, const string& ref,uint maxMiss,bool notAlignedSequence){
	string corrected,correctedRC,read,useless,header;
	while(not readFile.eof()){
		lockReadFile.lock();
		getline(readFile,header);
		getline(readFile,read);
		lockReadFile.unlock();
		if(not header.empty()){
			++readNumber;
			uint score(mapRead(read,  k, kmer2pos, ref,maxMiss,corrected));
			uint scorerc(mapRead(reversecomplement(read),  k, kmer2pos, ref,maxMiss,correctedRC));
			if(min(score,scorerc)<maxMiss){
				if(score<scorerc){
					lockOutFile.lock();
					mapped<<useless<<endl<<corrected<<endl;
					lockOutFile.unlock();
				}else{
					lockOutFile.lock();
					mapped<<useless<<endl<<reversecomplement(correctedRC)<<endl;
					lockOutFile.unlock();
				}
				mappedRead++;
			}else{
				if(notAlignedSequence){
					lockOutFile.lock();
					notMapped<<useless<<endl<<read<<endl;
					lockOutFile.unlock();
				} else {
					lockOutFile.lock();
					notMapped<<useless<<endl<<"not_aligned"<<endl;
					lockOutFile.unlock();
				}
				notMappedRead++;
			}
		}else{
			break;
		}
	}
}


uint mapReadFile(const string& readFileName,const uint64_t k, const unordered_map<kmer,vector<position>>& kmer2pos, const string& ref,uint maxMiss, bool notAlignedSequence,uint coreNumber){
    readFile.open(readFileName,ios::in);
    string read,useless,comp,corrected,correctedRC;
	vector<thread> threads;
	// treatRead(k, kmer2pos, ref,  maxMiss,  notAlignedSequence);
	for (size_t i(0); i<coreNumber; ++i){
		threads.push_back(thread(treatRead,k, cref(kmer2pos), cref(ref), maxMiss,notAlignedSequence));
	}
	for(auto &t : threads){t.join();}
    cout<<"Reads: "<<readNumber<<endl;
    cout<<"Reads mapped: "<<mappedRead<<endl;
    cout<<"Percent Read mapped: "<<((10000*(double)mappedRead)/readNumber)/100<<"%"<<endl;
	cout<<"Reads not mapped: "<<notMappedRead<<endl;
	return readNumber;
}

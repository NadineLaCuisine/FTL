#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iterator>
#include <ctime>
#include <unordered_map>
#include <algorithm>
#include "utils.h"


using namespace std;


bool isMappedHamming( string& read, const string& reference, uint maxMissmatch){
    uint error(0);
    for(uint i(0);i<read.size();++i){
        if(read[i]!=reference[i]){
            if(++error>maxMissmatch){
                return false;
            }
        }
    }
    return true;
}

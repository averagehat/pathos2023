/*
PRICE was written by Graham Ruby.
This file is part of PRICE.

PRICE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PRICE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PRICE.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef SEQFILTERLISTENER_CPP
#define SEQFILTERLISTENER_CPP

#include "SeqFilterListener.h"
#include "AssemblyException.h"
#include "SeqFilterListenerNull.h"
#include "SeqFilterListenerStd.h"
				    //#include "SeqFilterListenerVerbose.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <stdio.h>
using namespace::std;


SeqFilterListener* SeqFilterListener::makeNullListener(){ return new SeqFilterListenerNull(); }
SeqFilterListener* SeqFilterListener::makeStdListener(){ return new SeqFilterListenerStd(); }
//SeqFilterListener* SeqFilterListener::makeVerboseListener(){ return new SeqFilterListenerVerbose(); }

SeqFilterListener::~SeqFilterListener() {}

long SeqFilterListener::getTotalSize(set<ScoredSeq*>* contigs){
  long count = 0;
  for (set<ScoredSeq*>::iterator itS = contigs->begin(); itS != contigs->end(); ++itS){
    count += (*itS)->size();
  }
  return count;
}

long SeqFilterListener::getMaxSize(set<ScoredSeq*>* contigs){
  long maxSize = 0;
  for (set<ScoredSeq*>::iterator itS = contigs->begin(); itS != contigs->end(); ++itS){
    if ( (*itS)->size() > maxSize ){ maxSize = (*itS)->size(); }
  }
  return maxSize;
}


long SeqFilterListener::getN50(set<ScoredSeq*>* contigs){
  vector<long> contigLengths;
  for (set<ScoredSeq*>::iterator seqIt = contigs->begin(); seqIt != contigs->end(); ++seqIt){
    contigLengths.push_back( (*seqIt)->size() );
  }
  sort( contigLengths.begin(), contigLengths.end() );
  long halfSize = getTotalSize(contigs) / 2;

  long n50 = 0;
  vector<long>::iterator sizeIt = contigLengths.begin();
  while ( halfSize > 0 and sizeIt != contigLengths.end() ){
    n50 = (*sizeIt); // current n50; will only be updated if half the total size is not reached
    halfSize -= (*sizeIt);
    ++sizeIt;
  }
  return n50;
}



#endif

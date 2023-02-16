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

#ifndef EXTENDJOBMAPPER_CPP
#define EXTENDJOBMAPPER_CPP

#include "ExtendJobMapper.h"
#include <typeinfo>
#include <cstring>
#include <omp.h>
#include <stdio.h>
#include "ScoredSeqSubseq.h"
#include "ScoredSeqMonoScore.h"
using namespace::std;


ExtendJobMapper::~ExtendJobMapper(){}

long ExtendJobMapper::readToLong(char* seq, long length){
  long seqSum = 0;
  for (long n = 0; n < length; ++n){
    switch (seq[n]) {
    case 'N': break; // no info, nothing added
    case 'A': ++seqSum; break;
    case 'C': seqSum += 2; break;
    case 'G': seqSum += 3; break;
    case 'T': seqSum += 4; break;
    default:
      char errorMessage[200];
      char part1[] = "EJM::readToLong, invalid nucleotide:";
      sprintf(errorMessage, "%s %c", part1, seq[n]);
      throw AssemblyException::ArgError(errorMessage);
    }
  }
  return seqSum;
}



//ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(){}

ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId) :
  _fractMapId(fractMapId){
  _edgeSize = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    if ( (*it)->size() > _edgeSize ){ _edgeSize = (*it)->size(); }
  }
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestNull();
}
ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize) :
  _edgeSize(edgeSize),
  _fractMapId(fractMapId){
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestNull();
}
ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(long edgeTest, set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId) :
  _fractMapId(fractMapId){
  _edgeSize = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    if ( (*it)->size() > _edgeSize ){ _edgeSize = (*it)->size(); }
  }
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestEdgeProximity(edgeTest, fractMapId);
}
ExtendJobMapperWithLinks::ExtendJobMapperWithLinks(long edgeTest, set<ScoredSeqWithCarriers*>* contigsWithCarriers,
						   float fractMapId, long edgeSize) :
  _edgeSize(edgeSize),
  _fractMapId(fractMapId){
  constructorHelper(contigsWithCarriers);  
  _alTest = new AlTestEdgeProximity(edgeTest, fractMapId);
}


void ExtendJobMapperWithLinks::constructorHelper(set<ScoredSeqWithCarriers*>* contigsWithCarriers){

  _shortArray = new ScoredSeq*[ contigsWithCarriers->size() * 2 + 1 ];
  _flipArray = new ScoredSeq*[ contigsWithCarriers->size() * 2 + 1 ];

  _plusIndex = ExtendCycle::plusIndex();
  _minusIndex = ExtendCycle::minusIndex();
  _frontIndex = ExtendCycle::frontIndex();
  _backIndex = ExtendCycle::backIndex();

  // use the contigsWithCarriers to construct the front and back BWT subseq collections
  // minContigSplitLen -> edgeSize

  // first the front edges
  set<ScoredSeq*> inputSet;
  _disposableIndex = 0;

  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    ScoredSeqWithCarriers* carrier = (*it);
    ScoredSeq* keySeqFront;
    ScoredSeq* keySeqBack;
    if ( _edgeSize >= 0 and carrier->size() >= _edgeSize ){
      keySeqFront = new ScoredSeqSubseq( carrier, 0, _edgeSize);
      keySeqBack = new ScoredSeqSubseq( carrier, carrier->size() - _edgeSize, _edgeSize);
    } else {
      keySeqFront = new ScoredSeqSubseq( carrier, 0, carrier->size() );
      keySeqBack = new ScoredSeqSubseq( carrier, 0, carrier->size());
    }
    ScoredSeq* inSeqFront = ScoredSeqFlip::getFlip(keySeqFront,'-');
    ScoredSeq* inSeqBack = ScoredSeqFlip::getFlip(keySeqBack,'+');
    inputSet.insert(inSeqFront);
    inputSet.insert(inSeqBack);

    _shortArray[_disposableIndex] = keySeqFront;
    _flipArray[_disposableIndex] = inSeqFront;
    ++_disposableIndex;
    _shortArray[_disposableIndex] = keySeqBack;
    _flipArray[_disposableIndex] = inSeqBack;
    ++_disposableIndex;

    keySeqFront->buffer();
    keySeqBack->buffer();
  }
  _mappingCollection = new ScoredSeqCollectionBwt(&inputSet, _fractMapId, false);
  _mappingCollection->disableEdgeScaling();
}



ExtendJobMapperWithLinks::~ExtendJobMapperWithLinks(){
  // get rid of the subseq keys, not the contig values
  for (long n = 0; n < _disposableIndex; ++n){
    delete _shortArray[n];
    delete _flipArray[n];
  }
  delete [] _shortArray;
  delete [] _flipArray;

  // and delete the collections themselves
  delete _mappingCollection;
  delete _alTest;
}


void ExtendJobMapperWithLinks::mapQueriesHelper(vector<Alignment*>* matchesArray, ScoredSeqWithCarriers** queries,
						bool* seekMatch, long numQueries){

  long* queryLengths = new long[ numQueries+1 ];
  ScoredSeq** typeQueries = new ScoredSeq*[ numQueries+1 ];
  vector<Alignment*>** indexedMatches = new vector<Alignment*>*[ numQueries+1 ];
  long* originalIndexes = new long[ numQueries+1 ];

  long numFiltered = 0;
  for (long n = 0; n < numQueries; ++n){
    if (seekMatch[n]){
      typeQueries[numFiltered] = queries[n];
      queryLengths[numFiltered] = queries[n]->size();
      // only create these for the elements that will be searched
      indexedMatches[numFiltered] = &matchesArray[n];
      originalIndexes[numFiltered] = n;
      ++numFiltered;
    }
  }

  // this is just to the '+' strand because the collection was made up to consider what strands things should map
  // to; in the second mapping step, it is all of both strands, but the data structure is set up in the same manner
  // as for the first mapping step, where only distinct portions of each strand are used.
  _mappingCollection->getBestMatches(numFiltered, indexedMatches, typeQueries, '+', queryLengths,
				     ScoredSeqCollectionBwt::notThreaded, ScoredSeqCollectionBwt::softMinOvl);

  delete [] typeQueries;
  delete [] queryLengths;
  delete [] indexedMatches;
  delete [] originalIndexes;
}


long ExtendJobMapperWithLinks::mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
						bool* seekReadMatch, bool* seekPairMatch,
						bool* readMatchFound, bool* pairMatchFound, long numQueries){

  // this should have every value filled in; readPass and pairPass will be sparse
  bool* eitherPass = new bool[ numQueries+1 ];
  vector<Alignment*>* tempMatches = new vector<Alignment*>[ numQueries+1 ];

  vector<ScoredSeqWithCarriers*>* readMatchContigsPlus = new vector<ScoredSeqWithCarriers*>[ numQueries+1 ];
  vector<ScoredSeqWithCarriers*>* readMatchContigsMinus = new vector<ScoredSeqWithCarriers*>[ numQueries+1 ];
  vector<ScoredSeqWithCarriers*>* pairMatchContigsPlus = new vector<ScoredSeqWithCarriers*>[ numQueries+1 ];
  vector<ScoredSeqWithCarriers*>* pairMatchContigsMinus = new vector<ScoredSeqWithCarriers*>[ numQueries+1 ];

  // first, map the reads
  mapQueriesHelper(tempMatches, readQueries, seekReadMatch, numQueries);

  bool* readPass = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    if (seekReadMatch[n]){
      // i need to index all of the hits because i don't know about the PEs yet
      readPass[n] = _alTest->passes( &tempMatches[n], readQueries[n] );
      sortMatchesHelper(&tempMatches[n], &readMatchContigsPlus[n], &readMatchContigsMinus[n]);
      for (vector<Alignment*>::iterator it = tempMatches[n].begin(); it != tempMatches[n].end(); ++it){ delete *it; }
      eitherPass[n] = readPass[n];
    } else { eitherPass[n] = false; }
  }
  delete [] tempMatches;
  tempMatches = new vector<Alignment*>[ numQueries+1 ];

  // second, map the pairs
  mapQueriesHelper(tempMatches, pairQueries, seekPairMatch, numQueries);

  bool* pairPass = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    if (seekPairMatch[n]){
      if (! eitherPass[n]){
	pairPass[n] = _alTest->passes( &tempMatches[n], pairQueries[n] );
	eitherPass[n] = pairPass[n] or eitherPass[n];
      }
      if (eitherPass[n]){
	sortMatchesHelper(&tempMatches[n], &pairMatchContigsPlus[n], &pairMatchContigsMinus[n]);
	// go ahead and shrink the vectors now
	removeReciprocalMatches(&pairMatchContigsPlus[n], &pairMatchContigsMinus[n]);
	if (seekReadMatch[n]){ removeReciprocalMatches(&readMatchContigsPlus[n], &readMatchContigsMinus[n]); }
      }
      for (vector<Alignment*>::iterator it = tempMatches[n].begin(); it != tempMatches[n].end(); ++it){ delete *it; }
    }
  }
  delete [] tempMatches;


  // third, combine the results and create the linkages
  long count = 0;
  for (long n = 0; n < numQueries; ++n){
    if (eitherPass[n]){
      if (seekReadMatch[n] and readMatchContigsPlus[n].size() + readMatchContigsMinus[n].size() > 0){
	addToContigMatchesHelper(readQueries[n], _plusIndex, _backIndex, &readMatchContigsPlus[n]);
	addToContigMatchesHelper(readQueries[n], _minusIndex, _frontIndex, &readMatchContigsMinus[n]);
	readMatchFound[n] = true;
	++count;
      } else {
	readMatchFound[n] = false;
      }
      if (seekPairMatch[n] and pairMatchContigsPlus[n].size() + pairMatchContigsMinus[n].size() > 0){
	addToContigMatchesHelper(pairQueries[n], _plusIndex, _backIndex, &pairMatchContigsPlus[n]);
	addToContigMatchesHelper(pairQueries[n], _minusIndex, _frontIndex, &pairMatchContigsMinus[n]);
	pairMatchFound[n] = true;
	++count;
      } else {
	pairMatchFound[n] = false;
      }
    } else {
      readMatchFound[n] = false;
      pairMatchFound[n] = false;
    }
  }

  // fourth, clean up any leftover data structures
  delete [] readPass;
  delete [] pairPass;
  delete [] eitherPass;
  delete [] readMatchContigsPlus;
  delete [] readMatchContigsMinus;
  delete [] pairMatchContigsPlus;
  delete [] pairMatchContigsMinus;

  return count;
}


/*
long ExtendJobMapperWithLinks::mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
						bool* seekReadMatch, bool* seekPairMatch,
						bool* readMatchFound, bool* pairMatchFound, long numQueries){

  long numWithMatch = 0;
  // these will be filled in by the setup function
  vector<ScoredSeqWithCarriers*>** readMatchCarrier = new vector<ScoredSeqWithCarriers*>*[2];
  vector<ScoredSeqWithCarriers*>** pairMatchCarrier = new vector<ScoredSeqWithCarriers*>*[2];

  // do everything twice: once for reads, once for pairs
  // READS:
  long* qToUniqueReads = new long[numQueries+1];
  bool* uReadPasses = mapQueriesSetup(readQueries, seekReadMatch, numQueries, qToUniqueReads, readMatchCarrier);
  vector<ScoredSeqWithCarriers*>* readMatchContigsPlus = readMatchCarrier[0];
  vector<ScoredSeqWithCarriers*>* readMatchContigsMinus = readMatchCarrier[1];
  // PAIRS:
  long* qToUniquePairs = new long[numQueries+1];
  bool* uPairPasses = mapQueriesSetup(pairQueries, seekPairMatch, numQueries, qToUniquePairs, pairMatchCarrier);
  vector<ScoredSeqWithCarriers*>* pairMatchContigsPlus = pairMatchCarrier[0];
  vector<ScoredSeqWithCarriers*>* pairMatchContigsMinus = pairMatchCarrier[1];

  // now create links where appropriate
  for (long n = 0; n < numQueries; ++n){
    long urN = qToUniqueReads[n];
    long upN = qToUniquePairs[n];
    bool readPasses = seekReadMatch[n] and uReadPasses[urN];
    bool pairPasses = seekPairMatch[n] and uPairPasses[upN];
    if (readPasses or pairPasses){
      // this pair is ok, so fill in the values of reads, then pairs, as appropriate
      // READS:
      if ( seekReadMatch[n] ){
	removeReciprocalMatches(&readMatchContigsPlus[urN], &readMatchContigsMinus[urN]);
	if (readMatchContigsPlus[urN].size() + readMatchContigsMinus[urN].size() > 0){
	  addToContigMatchesHelper( readQueries[n], _plusIndex, _backIndex, &readMatchContigsPlus[urN]);
	  addToContigMatchesHelper( readQueries[n], _minusIndex, _frontIndex, &readMatchContigsMinus[urN]);
	  readMatchFound[n] = true;
	  ++numWithMatch;
	} else { readMatchFound[n] = false; }
      } else { readMatchFound[n] = false; }
      // PAIRS:
      if ( seekPairMatch[n] ){
	removeReciprocalMatches(&pairMatchContigsPlus[upN], &pairMatchContigsMinus[upN]);
	if (pairMatchContigsPlus[upN].size() + pairMatchContigsMinus[upN].size() > 0){
	  addToContigMatchesHelper( pairQueries[n], _plusIndex, _backIndex, &pairMatchContigsPlus[upN]);
	  addToContigMatchesHelper( pairQueries[n], _minusIndex, _frontIndex, &pairMatchContigsMinus[upN]);
	  pairMatchFound[n] = true;
	  ++numWithMatch;
	} else { pairMatchFound[n] = false; }
      } else { pairMatchFound[n] = false; }
    } else {
      readMatchFound[n] = false;
      pairMatchFound[n] = false;
    }
    readQueries[n]->deepUnbuffer();
    pairQueries[n]->deepUnbuffer();
  }

  delete [] readMatchCarrier;
  delete [] pairMatchCarrier;
  delete [] qToUniqueReads;
  delete [] qToUniquePairs;
  delete [] uReadPasses;
  delete [] uPairPasses;
  delete [] readMatchContigsPlus;
  delete [] readMatchContigsMinus;
  delete [] pairMatchContigsPlus;
  delete [] pairMatchContigsMinus;

  return numWithMatch;
}
*/


bool* ExtendJobMapperWithLinks::mapQueriesSetup(ScoredSeqWithCarriers** queries, bool* seekMatch,
						long numQueries, long* qToUnique,
						vector<ScoredSeqWithCarriers*>** matchCarrier){

  // applies the test and sorts out the results but doesn't create actual linkages
  ScoredSeq** monoQueries = new ScoredSeq*[numQueries+1];
  vector<Alignment*>** matches = new vector<Alignment*>*[numQueries+1];
  long numUnique = mapQueriesHelper(queries, qToUnique, seekMatch, matches, monoQueries, numQueries);
  vector<ScoredSeqWithCarriers*>* matchContigsPlus = new vector<ScoredSeqWithCarriers*>[numUnique+1];
  vector<ScoredSeqWithCarriers*>* matchContigsMinus = new vector<ScoredSeqWithCarriers*>[numUnique+1];
  bool* uPass = new bool[numUnique+1];

  // i NEED to sort the matches even if the sequence doesn't pass the test because any
  // one of the reads with that sequences' paired ends could pass
  for (long n = 0; n < numUnique; ++n){
    uPass[n] = _alTest->passes( matches[n], monoQueries[n] );
    sortMatchesHelper(matches[n], &matchContigsPlus[n], &matchContigsMinus[n]);
    for (vector<Alignment*>::iterator it = matches[n]->begin(); it != matches[n]->end(); ++it){ delete *it; }
    delete matches[n];
    delete monoQueries[n];
  }
  delete [] matches;
  delete [] monoQueries;
  matchCarrier[0] = matchContigsPlus;
  matchCarrier[1] = matchContigsMinus;
  return uPass;
}


long ExtendJobMapperWithLinks::mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchesFound, long numQueries){
  vector<Alignment*>** matchesArray = new vector<Alignment*>*[numQueries+1];
  ScoredSeq** monoQueries = new ScoredSeq*[numQueries+1];
  long* qToUnique = new long[numQueries+1];
  long numUnique = mapQueriesHelper(queries, qToUnique, seekMatch, matchesArray, monoQueries, numQueries);
  long numThatPass = 0;

  vector<ScoredSeqWithCarriers*>* matchingContigsPlus = new vector<ScoredSeqWithCarriers*>[numUnique+1];
  vector<ScoredSeqWithCarriers*>* matchingContigsMinus = new vector<ScoredSeqWithCarriers*>[numUnique+1];

  bool* uniquePasses = new bool[numUnique+1];
  for (long n = 0; n < numUnique; ++n){
    uniquePasses[n] = _alTest->passes(matchesArray[n], monoQueries[n]);
    if (uniquePasses[n]){
      sortMatchesHelper(matchesArray[n], &matchingContigsPlus[n], &matchingContigsMinus[n]);
    }
    delete monoQueries[n];
  }
  delete [] monoQueries;

  for (long n = 0; n < numQueries; ++n){
    long uN = qToUnique[n];
    if ( seekMatch[n] and uniquePasses[qToUnique[n]] ){
      removeReciprocalMatches(&matchingContigsPlus[uN], &matchingContigsMinus[uN]);
      if (matchingContigsPlus[uN].size() + matchingContigsMinus[uN].size() > 0){
	addToContigMatchesHelper( queries[n], _plusIndex, _backIndex, &matchingContigsPlus[uN]);
	addToContigMatchesHelper( queries[n], _minusIndex, _frontIndex, &matchingContigsMinus[uN]);
	++numThatPass;
	matchesFound[n] = true;
      } else { matchesFound[n] = false; }
    } else { matchesFound[n] = false; }
    queries[n]->deepUnbuffer();
  }

  for (long n = 0; n < numUnique; ++n){
    for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
    delete matchesArray[n];
  }
  delete [] matchesArray;
  delete [] matchingContigsPlus;
  delete [] matchingContigsMinus;
  delete [] uniquePasses;
  delete [] qToUnique;

  return numThatPass;
}

/*
bool ExtendJobMapperWithLinks::passesTest(vector<Alignment*>* matches){
  return matches->size() > 0;
  //return true;
}
*/

ExtendJobMapperWithLinks::AlignmentTest::~AlignmentTest(){}

// the null test always returns true
ExtendJobMapperWithLinks::AlTestNull::AlTestNull(){}
ExtendJobMapperWithLinks::AlTestNull::~AlTestNull(){}
bool ExtendJobMapperWithLinks::AlTestNull::passes(vector<Alignment*>* matches, ScoredSeq* query){
  return matches->size() > 0;
}
ExtendJobMapperWithLinks::AlTestEdgeProximity::AlTestEdgeProximity(long edgeProximity, float fractId) :
  _edgeProximity(edgeProximity),
  _fractId(fractId){}
ExtendJobMapperWithLinks::AlTestEdgeProximity::~AlTestEdgeProximity(){}

bool ExtendJobMapperWithLinks::AlTestEdgeProximity::passes(vector<Alignment*>* matches, ScoredSeq* query){
  bool foundPass = false;
  // this allows the query to overhang the window to a degree acceptable by the percent ID,
  // while also including a -1 that compensates for the zero-indexing of the final aligned pos
  long effQsizeM1 = long(float(query->size()) * _fractId) - 1;
  long aLastPos = query->size() - 1;

  vector<Alignment*>::iterator alIt = matches->begin();
  vector<Alignment*>::iterator alEnd = matches->end();
  while (alIt != alEnd){
    Alignment* al = *alIt;
    long bLastPos;
    if ( al->isLinked(aLastPos, query) ){ bLastPos = al->getLinkage(aLastPos, query); }
    else { bLastPos = al->gapPairedAfter(aLastPos, query); }

    // test if the full sequence of the query would fit in the window given where the
    // alignment ends (it is possible that gaps push the 5p end of the query beyond the
    // window, but i don't want that to penalize the match).
    foundPass = _edgeProximity >= al->seqB()->size() - bLastPos + effQsizeM1;

    // if the test passed, skip to the end
    if (foundPass){ alIt = alEnd; }
    else { ++alIt; }
  }
  return foundPass;
}


// MODIFIES both input vectors
void ExtendJobMapperWithLinks::removeReciprocalMatches(vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
						       vector<ScoredSeqWithCarriers*>* matchingContigsMinus){
  // the smaller of the two vectors will be made into a set (vector A); the
  // order-of-growth is NlogN for size(contigsA) and N for size(contigsB)
  vector<ScoredSeqWithCarriers*>* contigsA;
  vector<ScoredSeqWithCarriers*>* contigsB;
  if (matchingContigsPlus->size() < matchingContigsMinus->size()){
    contigsA = matchingContigsPlus;
    contigsB = matchingContigsMinus;
  } else {
    contigsB = matchingContigsPlus;
    contigsA = matchingContigsMinus;
  }

  // i don't need to do anything if the smaller vector is size zero
  if (contigsA->size() > 0){

    // this keeps the contigs that are ok from B; if not all of the contigs are being kept,
    // then the contents of contigsB can be replaced by the contents of this vector
    vector<ScoredSeqWithCarriers*> keepB;

    // this is the painful part - both of these operations require log entry/search times
    set<ScoredSeqWithCarriers*> nrSetA;
    for (vector<ScoredSeqWithCarriers*>::iterator it = contigsA->begin(); it != contigsA->end(); ++it){ nrSetA.insert(*it); }
    for (vector<ScoredSeqWithCarriers*>::iterator it = contigsB->begin(); it != contigsB->end(); ++it){
      set<ScoredSeqWithCarriers*>::iterator findIt = nrSetA.find(*it);
      if (findIt == nrSetA.end()){ keepB.push_back(*it); }
      else { nrSetA.erase( findIt ); }
    }

    // i don't have to do anything if contigsB and keepB are the same length
    if (keepB.size() < contigsB->size()){
      // replace the contents of B
      contigsB->clear();
      contigsB->insert(contigsB->begin(), keepB.begin(), keepB.end());
      // replace the contents of A (the non-redundant set can be used)
      contigsA->clear();
      for (set<ScoredSeqWithCarriers*>::iterator it = nrSetA.begin(); it != nrSetA.end(); ++it){ contigsA->push_back(*it); }
    }
  }
}


long ExtendJobMapperWithLinks::mapQueriesHelper(ScoredSeqWithCarriers** queries, long* qToUnique, bool* seekMatch,
						vector<Alignment*>** matchesArray, ScoredSeq** monoQueries, long numQueries){

  // make a non-redundant sequence set
  map<string, vector<long>*> seqToQuerySet;
  long numNr = 0;
  for (long n = 0; n < numQueries; ++n){
    if (seekMatch[n]){
      char* querySeq = queries[n]->getSeq('+');
      queries[n]->gatherSeq(querySeq,'+');
      map<string, vector<long>*>::iterator foundSeqIt = seqToQuerySet.find( querySeq );
      if (foundSeqIt == seqToQuerySet.end()){
	vector<long>* newCarrier = new vector<long>;
	newCarrier->push_back(n);
	seqToQuerySet.insert( pair<string, vector<long>*> (querySeq, newCarrier) );
      } else {
	foundSeqIt->second->push_back(n);
      }
      delete [] querySeq;
      queries[n]->deepUnbuffer();
    }
  }

  // set up the arrays
  long numMonoP1 = seqToQuerySet.size() + 1;
  vector<long>** likeSequences = new vector<long>*[ numMonoP1 ];
  //ScoredSeq** monoQueries = new ScoredSeq*[ numMonoP1 ];
  long* queryLengths = new long[ numMonoP1 ];
  // now use this as an index iterator - it will return to the full value
  long numMono = 0;
  for (map<string, vector<long>*>::iterator seqIt = seqToQuerySet.begin(); seqIt != seqToQuerySet.end(); ++seqIt){
    likeSequences[numMono] = seqIt->second;
    queryLengths[numMono] = seqIt->first.size();
    monoQueries[numMono] = new ScoredSeqMonoScore(seqIt->first, 1);
    matchesArray[numMono] = new vector<Alignment*>;
    ++numMono;
  }

  // this is just to the '+' strand because the collection was made up to consider what strands things should map
  // to; in the second mapping step, it is all of both strands, but the data structure is set up in the same manner
  // as for the first mapping step, where only distinct portions of each strand are used.
  _mappingCollection->getBestMatches(numMono, matchesArray, monoQueries, '+', queryLengths,
				     ScoredSeqCollectionBwt::notThreaded, ScoredSeqCollectionBwt::softMinOvl);

  for (long mqN = 0; mqN < numMono; ++mqN){
    for (vector<long>::iterator it = likeSequences[mqN]->begin(); it != likeSequences[mqN]->end(); ++it){
      queries[*it]->deepUnbuffer();
      qToUnique[*it] = mqN;
    }
    //delete monoQueries[mqN];
    delete likeSequences[mqN];
  }
  //delete [] monoQueries;
  delete [] likeSequences;
  delete [] queryLengths;

  return numMono;
}


// MODIFIES: matchingContigs AND alreadyQueriedToHits
void ExtendJobMapperWithLinks::findMatchesHelper(ScoredSeq* currentQuery,
						 vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
						 vector<ScoredSeqWithCarriers*>* matchingContigsMinus){
  vector<Alignment*> matchAlignments;
  _mappingCollection->getBestMatches( &matchAlignments, currentQuery, '+', currentQuery->size(), ScoredSeqCollectionBwt::softMinOvl);

  for (vector<Alignment*>::iterator alIt = matchAlignments.begin(); alIt != matchAlignments.end(); ++alIt){
    // find the full contig
    ScoredSeqFlip* flipContig = dynamic_cast<ScoredSeqFlip*>( (*alIt)->seqB() );
    ScoredSeqSubseq* partialContig = dynamic_cast<ScoredSeqSubseq*>( flipContig->getNested() );
    ScoredSeqWithCarriers* fullContig = dynamic_cast<ScoredSeqWithCarriers*>( partialContig->getNested() );
    switch ( flipContig->getSense() ){
    case '+': matchingContigsPlus->push_back( fullContig ); break;
    case '-': matchingContigsMinus->push_back( fullContig ); break;
    default: throw AssemblyException::LogicError("EJMwLinks::findMatchesHelper got a bad sense char");
    }
    delete *alIt;
  }
}

// MODIFIES: matchingContigs AND alreadyQueriedToHits
void ExtendJobMapperWithLinks::sortMatchesHelper(vector<Alignment*>* matches,
						 vector<ScoredSeqWithCarriers*>* matchingContigsPlus,
						 vector<ScoredSeqWithCarriers*>* matchingContigsMinus){

  for (vector<Alignment*>::iterator alIt = matches->begin(); alIt != matches->end(); ++alIt){
    // find the full contig
    ScoredSeqFlip* flipContig = dynamic_cast<ScoredSeqFlip*>( (*alIt)->seqB() );
    ScoredSeqSubseq* partialContig = dynamic_cast<ScoredSeqSubseq*>( flipContig->getNested() );
    ScoredSeqWithCarriers* fullContig = dynamic_cast<ScoredSeqWithCarriers*>( partialContig->getNested() );
    switch ( flipContig->getSense() ){
    case '+': matchingContigsPlus->push_back( fullContig ); break;
    case '-': matchingContigsMinus->push_back( fullContig ); break;
    default: throw AssemblyException::LogicError("EJMwLinks::findMatchesHelper got a bad sense char");
    }
  }
}


void ExtendJobMapperWithLinks::addToContigMatchesHelper(ScoredSeqWithCarriers* seqToAdd, int side, int frontOrBack,
							vector<ScoredSeqWithCarriers*>* matchingContigs){
  #pragma omp critical (EJMlink)
  {
    for (vector<ScoredSeqWithCarriers*>::iterator mcIt = matchingContigs->begin(); mcIt != matchingContigs->end(); ++mcIt){
      (*mcIt)->addCarriedSeq( seqToAdd, frontOrBack );
      seqToAdd->addCarriedSeq( (*mcIt), side );
    }
  }
}




ExtendJobMapperNoLinks::ExtendJobMapperNoLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId) :
  _fractMapId(fractMapId){
  _edgeSize = 0;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    if ( (*it)->size() > _edgeSize ){ _edgeSize = (*it)->size(); }
  }
  constructorHelper(contigsWithCarriers);  
}
ExtendJobMapperNoLinks::ExtendJobMapperNoLinks(set<ScoredSeqWithCarriers*>* contigsWithCarriers, float fractMapId, long edgeSize) :
  _edgeSize(edgeSize),
  _fractMapId(fractMapId){
  constructorHelper(contigsWithCarriers);  
}
void ExtendJobMapperNoLinks::constructorHelper(set<ScoredSeqWithCarriers*>* contigsWithCarriers){

  // use the contigsWithCarriers to construct the front and back BWT subseq collection (one collection)
  set<ScoredSeq*> shortSet;

  // first the front edges
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    ScoredSeqWithCarriers* carrier = (*it);
    ScoredSeq* keySeq;
    if ( _edgeSize >= 0 and carrier->size() >= _edgeSize ){
      keySeq = new ScoredSeqSubseq( carrier, 0, _edgeSize);
    } else {
      keySeq = new ScoredSeqSubseq( carrier, 0, carrier->size() );
    }
    shortSet.insert(keySeq);
    keySeq->buffer();
  }

  // now the back edges
  set<ScoredSeq*> shortBackSet;
  for (set<ScoredSeqWithCarriers*>::iterator it = contigsWithCarriers->begin(); it != contigsWithCarriers->end(); ++it){
    ScoredSeqWithCarriers* carrier = (*it);
    ScoredSeq* keySeq;
    if ( _edgeSize >= 0 and carrier->size() >= _edgeSize ){
      keySeq = new ScoredSeqSubseq( carrier, carrier->size() - _edgeSize, _edgeSize);
    } else {
      keySeq = new ScoredSeqSubseq( carrier, 0, carrier->size());
    }
    shortSet.insert(keySeq);
    keySeq->buffer();
  }

  _mappingCollection = new ScoredSeqCollectionBwt(&shortSet, _fractMapId, false);
  _mappingCollection->disableEdgeScaling();

  // MORE HERE???
}

ExtendJobMapperNoLinks::~ExtendJobMapperNoLinks(){
  // get rid of the subseq keys, not the contig values
  set<ScoredSeq*> outerToDelete;
  _mappingCollection->getSeqs( &outerToDelete );
  for (set<ScoredSeq*>::iterator it = outerToDelete.begin(); it != outerToDelete.end(); ++it){ delete (*it); }
  delete _mappingCollection;
}


ExtendJobMapper::QueryCarrier::QueryCarrier(char* query, long size) : _query(query), _size(size) {}
ExtendJobMapper::QueryCarrier::~QueryCarrier(){}


long ExtendJobMapperNoLinks::mapPairedQueries(ScoredSeqWithCarriers** readQueries, ScoredSeqWithCarriers** pairQueries,
					      bool* seekReadMatch, bool* seekPairMatch,
					      bool* readMatchFound, bool* pairMatchFound, long numQueries){
  // simply forward the method - for now
  long hitCount = mapQueries(readQueries, seekReadMatch, readMatchFound, numQueries);
  hitCount += mapQueries(pairQueries, seekPairMatch, pairMatchFound, numQueries);
  return hitCount;
}



// this version has no redundancy collapse
long ExtendJobMapperNoLinks::mapQueries(ScoredSeqWithCarriers** queries, bool* seekMatch, bool* matchFound, long numQueries){

  long* queryLengths = new long[ numQueries+1 ];
  ScoredSeq** typeQueries = new ScoredSeq*[ numQueries+1 ];
  long* originalIndexes = new long[ numQueries+1 ];
  long numFiltered = 0;
  for (long n = 0; n < numQueries; ++n){
    if (seekMatch[n]){
      typeQueries[numFiltered] = queries[n];
      queryLengths[numFiltered] = queries[n]->size();
      originalIndexes[numFiltered] = n;
      ++numFiltered;
    } else {
      matchFound[n] = false;
    }
  }

  bool* qHasMatch = _mappingCollection->hasMatch(numFiltered, typeQueries, queryLengths,
						 ScoredSeqCollectionBwt::notThreaded,
						 ScoredSeqCollectionBwt::softMinOvl );

  long countWithHits = 0;
  for (long fn = 0; fn < numFiltered; ++fn){
    matchFound[originalIndexes[fn]] = qHasMatch[fn];
    if (qHasMatch[fn]){ ++countWithHits; }
  }
  delete [] typeQueries;
  delete [] queryLengths;
  delete [] qHasMatch;
  delete [] originalIndexes;
  return countWithHits;
}



#endif

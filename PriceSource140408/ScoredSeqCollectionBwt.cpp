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

#ifndef SCOREDSEQCOLLECTIONBWT_CPP
#define SCOREDSEQCOLLECTIONBWT_CPP

#include "ScoredSeqCollectionBwt.h"

#include "AssemblyException.h"
#include "AlignmentNull.h"
#include "AlignmentUngapped.h"
#include "ScoredSeqFlip.h"
#include <iostream>
#include <omp.h>
#include <algorithm>
using namespace::std;


long ScoredSeqCollectionBwt::_defaultMinScore = 0;

ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, bool penalizeEdgeGaps) : 
  _asif(0),
  _gapMode(ungapped),
  _fractId(fractId),
  _minOverlap(0),
  _penalizeEdgeGaps(penalizeEdgeGaps),
  _usingMaxOverlap(false),
  _maxOverlap(-1) {
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap) : 
  _asif(0),
  _gapMode(ungapped),
  _fractId(fractId),
  _minOverlap(minOverlap),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1) {
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap, long maxOverlap) : 
  _asif(0),
  _gapMode(ungapped),
  _fractId(fractId),
  _minOverlap(minOverlap),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap) {
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DyProAlignerFactory * asif) :
  _fractId(asif->getFractId()),
  _minOverlap(asif->getMinOverlap()),
  _gapMode(gapped),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(false),
  _maxOverlap(-1),
  _asif(asif->copy())
{
  constructorHelper(inputSeqs);
}
ScoredSeqCollectionBwt::ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DyProAlignerFactory * asif, long maxOverlap) :
  _fractId(asif->getFractId()),
  _minOverlap(asif->getMinOverlap()),
  _gapMode(gapped),
  _penalizeEdgeGaps(false),
  _usingMaxOverlap(true),
  _maxOverlap(maxOverlap),
  _asif(asif->copy())
 {
  constructorHelper(inputSeqs);
}


void ScoredSeqCollectionBwt::constructorHelper(set<ScoredSeq*>* inputSeqs){
  // fill _allSeqs now so that the log of the collection doesn't need to be faced again
  _allSeqs.insert( inputSeqs->begin(), inputSeqs->end() );

  // by default, edge scaling will happen.  it can be disabled to save
  // time during read mapping.
  _edgeScalingEnabled = true;

  // minOverlap cannot be zero for this class; sub-string words must be defined
  // and statistically scaled for this overlap minimum, so a hard bottom line
  // needs to be defined here
  long minOverlapLimit = 4;
  if (_minOverlap < minOverlapLimit){ _minOverlap = minOverlapLimit; }

  // fractId must be at least that expected for random sequence matches
  if (_fractId < 0.25){ _fractId = 0.25; }

  // figure out the windows and make the bins
  //long seqCount = 0;
  long seqLenSum = 0;
  // determine the length of the shortest sequence and whether or not it is
  // smaller than _minOverlap
  long minSeqLength = 0;
  long seqCount = inputSeqs->size();

  set<ScoredSeq*>::iterator seqIt = inputSeqs->begin();
  if (seqIt != inputSeqs->end()){
    minSeqLength = (*seqIt)->size();
    seqLenSum = minSeqLength + 1;
    ++seqIt;
    while (seqIt != inputSeqs->end()){
      long seqSize = (*seqIt)->size();
      if (seqSize < minSeqLength){ minSeqLength = seqSize; }
      seqLenSum += seqSize + 1;
      ++seqIt;
    }
  }
  if (minSeqLength < _minOverlap){ _minFragmentSize = _minOverlap; }
  else { _minFragmentSize = long(float(minSeqLength) * _fractId); }

  // some data that will be useful to have pre-computed
  _maxFractMis = 1.0 - _fractId;
  if (_asif == NULL){
    long matchScore = 1;
    long misScore = -1;
    _scoreAsm = new AlignmentScoreMatrix(matchScore,misScore,misScore,misScore);
  } else {
    _scoreAsm = _asif->getScoreMatrix();
  }
  _matchCountAsm = new AlignmentScoreMatrix(1,0,0,0);
  _misCountAsm = new AlignmentScoreMatrix(0,1,1,1); // gaps count as mismatches (terminal gaps will exist)

  // figure out the windows and make the bins
  if (seqCount==0){
    _binSize = 0; // this value is only possible in this case
    _maxBin = 0;
  } else {
    _binSize = (seqLenSum / seqCount) + 1; // +1 to round up
    _maxBin = seqLenSum / _binSize;
  }
  _binToSeqs = new vector<ScoredSeqLocus*>[_maxBin+1];

  // convert the full text to a numeric array, filling in the bins
  if (seqLenSum==0){ seqLenSum = 1; } // at least an end char will be added
  _textSize = seqLenSum;

  // a local array of text length
  long* textAsNums = new long[seqLenSum];

  long currentIndex = 0;
  _numContigLoci = seqCount;
  if (_numContigLoci > 0){
    _orderedContigLoci = new ScoredSeqLocus*[ _numContigLoci ];
  }
  int numThreads = omp_get_max_threads();
  _contigIndexWasHit = new bool*[ numThreads ];
  _localContigIndex = new long*[ numThreads ];
  // to minimize the number of "new" calls
  _contigIndexWasHitMemory = new bool[ _numContigLoci * numThreads + 1 ];
  _localContigIndexMemory = new long[ _numContigLoci * numThreads + 1 ];

  for (int n = 0; n < numThreads; ++n){
    _contigIndexWasHit[n] = &_contigIndexWasHitMemory[ n * _numContigLoci ];
    for (long n2 = 0; n2 < _numContigLoci; ++n2){ _contigIndexWasHit[n][n2] = false; }
    _localContigIndex[n] = &_localContigIndexMemory[ n * _numContigLoci ];
  }

  // the currentOffset number is unique to each locus, and is therefore also
  // unique to each contig, and is continuous and zero-indexed
  long currentOffset = 0;
  // this will be used to separate the sequences
  char separator[2];
  separator[0] = '&';
  separator[1] = '\0';

  long seqStringLen = 0;
  char* seqString = new char[1];

  for (set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it){
    // add the sequence, preceeded by a separator if necessary
    if ( it != _allSeqs.begin() ){
      burrowsWheelerUtilities::makeNumString( &separator[0], long(1), currentOffset, textAsNums );
      currentOffset++;
    }

    // possibly replace seqString to make one that is long enough
    ScoredSeq* seq = (*it);
    if (seq->size() > seqStringLen){
      delete [] seqString;
      seqString = new char[ seq->size() + 1 ];
    }
    seq->gatherSeq(seqString,'+');
    burrowsWheelerUtilities::makeNumString( seqString, seq->size(), currentOffset, textAsNums );

    // fill the overlapping bins with this ScoredSeq in Locus form
    long firstBin = currentOffset / _binSize;
    long lastBinP1 = ( currentOffset + seq->size() ) / _binSize + 1;
    ScoredSeqLocus* seqAsLocus = new ScoredSeqLocus(seq, currentOffset, currentOffset + seq->size() - 1, currentIndex);
    for (long bin = firstBin; bin < lastBinP1; ++bin){ _binToSeqs[bin].push_back( seqAsLocus ); }

    _orderedContigLoci[currentIndex] = seqAsLocus;
    currentIndex++;
    currentOffset += seq->size();
  }
  delete [] seqString;

  // end the sequence with an end char (even if there were no sequence entries)
  char endChar[2];
  endChar[0] = '#';
  endChar[1] = '\0';
  burrowsWheelerUtilities::makeNumString( &endChar[0], long(1), currentOffset, textAsNums );

  // make the burrows-wheeler transform
  _sortedSuffixes = new long[seqLenSum];
  long alphaSize = burrowsWheelerUtilities::getAlphabetSize();
  burrowsWheelerUtilities::sortSuffixes( textAsNums, _sortedSuffixes, seqLenSum, alphaSize );
  _bwTransform = new long[seqLenSum];
  burrowsWheelerUtilities::bwTransform(textAsNums, _sortedSuffixes, seqLenSum, _bwTransform);

  delete [] textAsNums;

  // make the accessory data structures for alignment
  _occCounts = new long[seqLenSum];
  _tableC = burrowsWheelerUtilities::fillBwtHelperArrays(_bwTransform, seqLenSum, alphaSize, _occCounts);
  //OK();
}



void ScoredSeqCollectionBwt::OK(){
  if ( _fractId > 1 ){
    throw AssemblyException::ArgError("fraction ID is a fraction; cannot be greater than 1");
  }
  if ( _allSeqs.empty() ){
    if (_textSize != 1){ throw AssemblyException::LogicError("BWT text len is 1 if there are no seqs"); }
  } else {
    long seqLenSum = 0;
    for (set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it){
      seqLenSum += (*it)->size() + 1; // +1 because of the separator/end char
    }
    if (_textSize != seqLenSum){ throw AssemblyException::LogicError("BWT text len is wrong"); }
  }
}



ScoredSeqCollectionBwt::~ScoredSeqCollectionBwt(){
  delete _asif;
  delete _matchCountAsm;
  delete _misCountAsm;
  delete _scoreAsm;
  delete [] _sortedSuffixes;
  delete [] _bwTransform;
  delete [] _occCounts;
  delete [] _tableC;
  for (long n = 0; n < _numContigLoci; ++n){ delete _orderedContigLoci[n]; }
  if (_numContigLoci > 0){ delete [] _orderedContigLoci; }

  delete [] _binToSeqs;

  delete [] _contigIndexWasHit;
  delete [] _localContigIndex;
  delete [] _contigIndexWasHitMemory;
  delete [] _localContigIndexMemory;
}




ScoredSeqCollectionBwt * ScoredSeqCollectionBwt::copy(){
  // uses a private constructor to efficiently copy stuff
  throw AssemblyException::ImplementationError("not implemented");
}

float ScoredSeqCollectionBwt::getFractId(){ return _fractId; }
long ScoredSeqCollectionBwt::getMinOverlap(){ return _minOverlap; }



// DEFAULT: edge scaling is enabled
void ScoredSeqCollectionBwt::enableEdgeScaling(){ _edgeScalingEnabled = true; }
void ScoredSeqCollectionBwt::disableEdgeScaling(){ _edgeScalingEnabled = false; }



bool ScoredSeqCollectionBwt::contains(ScoredSeq* s){
  //OK();
  return ( _allSeqs.find(s) != _allSeqs.end() );
}



void ScoredSeqCollectionBwt::getSeqs( set<ScoredSeq*>* allSeqs ){
  //OK();
  allSeqs->insert( _allSeqs.begin(), _allSeqs.end() );
}



// UNTHREADED VERSIONS
// min overlap is specified - these four methods cascade down to the bottom one
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq, long minOverlap,
					MinOvlStringency ovlStringency){
  getMatches(matches, seq, '.', minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, long minOverlap,
					MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  getMatches(matches, seq, sense, seqTest, minOverlap, ovlStringency);
  delete seqTest;
}
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq,
					MatchSeqTest* matchTest, long minOverlap,
					MinOvlStringency ovlStringency){
  getMatches(matches, seq, '.', matchTest, minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
					MatchSeqTest* matchTest, long minOverlap,
					MinOvlStringency ovlStringency){
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense,
		       matchTest, minOverlap,
		       semiGlobal, allMatches, ovlStringency, notThreaded);
  matches->insert( matches->end(), tempMatches.begin(), tempMatches.end() );
}



// THREADED VERSIONS
// these four methods cascade down to the one on the bottom
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getMatches(numQueries, matchesArray, seqArray, '.', minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  MatchSeqTest** matchTestArray = new MatchSeqTest*[  numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ matchTestArray[n] = seqTest; }
  getMatches(numQueries, matchesArray, seqArray, sense, matchTestArray, minOvlArray, threadedness, ovlStringency);
  delete seqTest;
  delete [] matchTestArray;
}
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
					MatchSeqTest** matchTestArray, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getMatches(numQueries, matchesArray, seqArray, '.', matchTestArray, minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
					MatchSeqTest** matchTestArray, long* minOvlArray,
					AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){


  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, matchTestArray, minOvlArray,
		     semiGlobal, allMatches, ovlStringency, threadedness);

  for (long n = 0; n < numQueries; ++n){
    matchesArray[n]->insert( matchesArray[n]->end(), tempMatches[n]->begin(), tempMatches[n]->end() );
    delete tempMatches[n];
  }
  delete [] tempMatches;
}





long ScoredSeqCollectionBwt::makeAlignmentHelper(OffsetTracker* tracker, ScoredSeq* seq, vector<Alignment*>* matches,
						 long minOverlap, AlignmentMaker* alMaker, ScoredSeq* seqFlip, char sense,
						 long bestScoreSoFar, bool onlyBestAlignment, GapMode gapMode,
						 MatchesSought matchesSought){

  bool justFirstMatch = (matchesSought == firstMatch);
  ScoredSeq* contig = tracker->targetContig();
  long numOffsets = tracker->numOffsets();
  OffsetAndQueryCarrier** allOffsets = tracker->getOffsets(gapMode == ungapped);

  OffsetAndQueryCarrier** adequateOffsets = new OffsetAndQueryCarrier*[ numOffsets + 1 ];
  long numAdequateOffsets = 0;
  for (long allIndex = 0; allIndex < numOffsets; ++allIndex){
    // determine how long the alignment will be??  i will need to do that to make sure
    // that only legit alignments are kept.
    long alLength = getOverlapFromOffset(allOffsets[allIndex]->_offset, seq->size(), contig->size());
    if (alLength < 0){ throw AssemblyException::LogicError("SSCBwt::getMatches, alLength should be > 0"); }
    if ( (_penalizeEdgeGaps or alLength >= minOverlap) and ( (! _usingMaxOverlap) or alLength <= _maxOverlap) ){
      adequateOffsets[numAdequateOffsets] = allOffsets[allIndex];
      ++numAdequateOffsets;
    }
  }

  // this will be different for gapped versus ungapped alignments
  if (gapMode == ungapped){ // UNGAPPED ALIGNMENT
    vector<Alignment*> tempMatches;
    long inputBestScore = bestScoreSoFar;

    // this will reduce the number of seq replacements and (more time-consuming) RcFlips
    // that need to be performed on alignments (only in the ungapped case; deal with gapped differently?)
    for (long offN = 0; offN < numAdequateOffsets; ++offN){
      Alignment* alignCand = alMaker->makeAlignment(seqFlip, contig, '+', adequateOffsets[offN]);
      if ( alignCand->isNull() ){ delete alignCand; }
      else {
	if (justFirstMatch){ offN = numAdequateOffsets; }
	bool returnMatch = true;
	if (onlyBestAlignment){
	  long newScore = alMaker->getScore(alignCand);
	  if (newScore < bestScoreSoFar){ returnMatch = false; }
	  // if they are equal, then the match is just added
	  else if (newScore > bestScoreSoFar){
	    bestScoreSoFar = newScore;
	    for (vector<Alignment*>::iterator it = tempMatches.begin(); it != tempMatches.end(); ++it){ delete *it; }
	    tempMatches.clear();
	    alMaker->setMinScore(newScore);
	  }
	}
	if (returnMatch){ tempMatches.push_back(alignCand); }
	else { delete alignCand; }
      }
    }
    // now compare to and deal with the input set's contents
    if (onlyBestAlignment and bestScoreSoFar > inputBestScore){
      for (vector<Alignment*>::iterator it = matches->begin(); it != matches->end(); ++it){ delete *it; }
      matches->clear();
    }
    // and modify the output alignments appropriately according to the seq's orientation
    matches->insert(matches->end(), tempMatches.begin(), tempMatches.end());

  } else { // GAPPED ALIGNMENT
    if (numAdequateOffsets > 0){
      DyProAlignerFactory * localAsif = _asif->minScoreCopy(bestScoreSoFar);
      DynamicProgrammingAligner* localAsi = localAsif->makeDyProAligner(seqFlip,contig,'+');
      long* justOffsetNums = new long[ numAdequateOffsets ];
      for (long offN = 0; offN < numAdequateOffsets; ++offN){ justOffsetNums[offN] = adequateOffsets[offN]->_offset; }

      Alignment* alignCand = localAsi->align(justOffsetNums,numAdequateOffsets);

      delete [] justOffsetNums;
      if ( alignCand->isNull() ){ delete alignCand; }
      else {
	bool returnMatch = true;
	if (onlyBestAlignment){
	  long newScore = alignCand->score(_scoreAsm,_penalizeEdgeGaps);
	  if (newScore < bestScoreSoFar){ returnMatch = false; }
	  // if they are equal, then the match is just added
	  else if (newScore > bestScoreSoFar){
	    bestScoreSoFar = newScore;
	    for (vector<Alignment*>::iterator it = matches->begin(); it != matches->end(); ++it){ delete *it; }
	    matches->clear();
	  }
	}
	if (returnMatch){ matches->push_back(alignCand); }
	else { delete alignCand; }
      }
      delete localAsi;
      delete localAsif;
    }
  }

  for (long allIndex = 0; allIndex < numOffsets; ++allIndex){
    for (long wqN = 0; wqN < allOffsets[allIndex]->_numQueries; ++wqN){ delete allOffsets[allIndex]->_queries[wqN]; }
    delete allOffsets[allIndex];
  }
  delete [] allOffsets;
  delete [] adequateOffsets;

  return bestScoreSoFar;
}




// min overlap is specified - these four methods cascade down
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq,
					    long minOverlap, MinOvlStringency ovlStringency){
  getBestMatches(matches, seq, '.', minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
					    long minOverlap, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  getBestMatches(matches, seq, sense, seqTest, minOverlap, ovlStringency);
  delete seqTest;
}
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, MatchSeqTest* seqTest,
					    long minOverlap, MinOvlStringency ovlStringency){
  getBestMatches(matches, seq, '.', seqTest, minOverlap, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* seqTest,
					    long minOverlap, MinOvlStringency ovlStringency){
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, seqTest, minOverlap, semiGlobal, bestMatches, ovlStringency, notThreaded);
  matches->insert( matches->end(), tempMatches.begin(), tempMatches.end() );
}

// THREADED VERSIONS
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getBestMatches(numQueries, matchesArray, seqArray, '.', minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  MatchSeqTest** seqTestArray = new MatchSeqTest*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ seqTestArray[n] = seqTest; }
  getBestMatches(numQueries, matchesArray, seqArray, sense, seqTestArray, minOvlArray, threadedness, ovlStringency);
  delete seqTest;
  delete [] seqTestArray;
}
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, MatchSeqTest** seqTestArray,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  getBestMatches(numQueries, matchesArray, seqArray, '.', seqTestArray, minOvlArray, threadedness, ovlStringency);
}
void ScoredSeqCollectionBwt::getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray,
					    long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, seqTestArray, minOvlArray, semiGlobal, bestMatches, ovlStringency, threadedness);

  for (long n = 0; n < numQueries; ++n){
    matchesArray[n]->insert( matchesArray[n]->end(), tempMatches[n]->begin(), tempMatches[n]->end() );
    delete tempMatches[n];
  }
  delete [] tempMatches;
}



// THE BEGINNING OF HASMATCH
bool ScoredSeqCollectionBwt::hasFullMatch(ScoredSeq* seq, char sense, MatchSeqTest* seqTest){
  MatchSeqTest* useSeqTest;
  if (seqTest == NULL){ useSeqTest = MatchSeqTest::getNullTest(); }
  else { useSeqTest = seqTest; }
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, useSeqTest, _minFragmentSize, fullAlignment, firstMatch, hardMinOvl, notThreaded);
  bool hasMatch = (tempMatches.size() > 0);
  for (vector<Alignment*>::iterator it = tempMatches.begin(); it != tempMatches.end(); ++it){ delete *it; }
  if (seqTest == NULL){ delete useSeqTest; }
  return hasMatch;
}


bool* ScoredSeqCollectionBwt::hasFullMatch(long numQueries, ScoredSeq** seqArray,
					  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray){
  return hasFullMatch(numQueries, seqArray, '.', threadedness, seqTestArray);
}

bool* ScoredSeqCollectionBwt::hasFullMatch(long numQueries, ScoredSeq** seqArray, char sense,
					  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray){
  MatchSeqTest** useSeqTestArray;
  MatchSeqTest* dummyTest = NULL;
  if (seqTestArray == NULL){
    dummyTest = MatchSeqTest::getNullTest();
    useSeqTestArray = new MatchSeqTest*[ numQueries+1 ];
    for (long n = 0; n < numQueries; ++n){ useSeqTestArray[n] = dummyTest; }
  } else { useSeqTestArray = seqTestArray; }
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  long* minOvlArray = new long[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ minOvlArray[n] = _minFragmentSize; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, useSeqTestArray, minOvlArray, fullAlignment, firstMatch, hardMinOvl, threadedness);

  delete [] minOvlArray;

  bool* answers = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    answers[n] = (tempMatches[n]->size() > 0);
    for (vector<Alignment*>::iterator it = tempMatches[n]->begin(); it != tempMatches[n]->end(); ++it){ delete *it; }
  }
  for (long n = 0; n < numQueries; ++n){ delete tempMatches[n]; }
  delete [] tempMatches;
  if (seqTestArray == NULL){
    delete dummyTest;
    delete [] useSeqTestArray;
  }
  return answers;
}



// min overlap is specified - these four methods cascade down
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, long minOverlap, MinOvlStringency ovlStringency){
  return hasMatch(seq, '.', minOverlap, ovlStringency);
}
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, char sense, long minOverlap, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  bool returnVal = hasMatch(seq, sense, seqTest, minOverlap, ovlStringency);
  delete seqTest;
  return returnVal;
}
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency){
  return hasMatch(seq, '.', seqTest, minOverlap, ovlStringency);
}
bool ScoredSeqCollectionBwt::hasMatch(ScoredSeq* seq, char sense, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency){
  vector<Alignment*> tempMatches;
  runSearchesConverter(&tempMatches, seq, sense, seqTest, minOverlap, semiGlobal, firstMatch, ovlStringency, notThreaded);
  bool hasMatch = (tempMatches.size() > 0);
  for (vector<Alignment*>::iterator it = tempMatches.begin(); it != tempMatches.end(); ++it){ delete *it; }
  return hasMatch;
}


// THREADED
// min overlap is specified - these four methods cascade down
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  return hasMatch(numQueries, seqArray, '.', minOvlArray, threadedness, ovlStringency);
}
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, char sense, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  MatchSeqTest* seqTest = MatchSeqTest::getNullTest();
  MatchSeqTest** seqTestArray = new MatchSeqTest*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ seqTestArray[n] = seqTest; }
  bool* returnVal = hasMatch(numQueries, seqArray, sense, seqTestArray, minOvlArray, threadedness, ovlStringency);
  delete [] seqTestArray;
  delete seqTest;
  return returnVal;
}
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, MatchSeqTest** seqTestArray, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  return hasMatch(numQueries, seqArray, '.', seqTestArray, minOvlArray, threadedness, ovlStringency);
}
bool* ScoredSeqCollectionBwt::hasMatch(long numQueries, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray, long* minOvlArray,
				      AlignmentThreadedness threadedness, MinOvlStringency ovlStringency){
  vector<Alignment*>** tempMatches = new vector<Alignment*>*[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){ tempMatches[n] = new vector<Alignment*>; }

  runSearchesPrivate(numQueries, tempMatches, seqArray, sense, seqTestArray, minOvlArray, semiGlobal, firstMatch, ovlStringency, threadedness);

  bool* answers = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    answers[n] = (tempMatches[n]->size() > 0);
    for (vector<Alignment*>::iterator it = tempMatches[n]->begin(); it != tempMatches[n]->end(); ++it){ delete *it; }
  }
  for (long n = 0; n < numQueries; ++n){ delete tempMatches[n]; }
  delete [] tempMatches;
  return answers;
}




// THE END OF HASMATCH


// ASSUMES: matches is empty; otherwise, it will be modified
void ScoredSeqCollectionBwt::runSearchesConverter(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
						  MatchSeqTest* seqTest, long minOverlap,
						  AlignmentType alType, MatchesSought matchesSought,
						  MinOvlStringency ovlStringency, AlignmentThreadedness threadedness){

  // since these are all short arrays, I can just declare them here
  vector<Alignment*>* matchesArray[ 1 ];
  matchesArray[ 0 ] = matches;
  ScoredSeq* seqArray[ 1 ];
  seqArray[0] = seq;
  MatchSeqTest* seqTestArray[ 1 ];
  seqTestArray[0] = seqTest;

  long minOvlArray[ 1 ];
  minOvlArray[0] = minOverlap;

  runSearchesPrivate(1, &matchesArray[0], &seqArray[0], sense,
		     &seqTestArray[0], &minOvlArray[0],
		     alType, matchesSought, ovlStringency, threadedness);
}



// ASSUMES: matches is empty; otherwise, it will be modified
void ScoredSeqCollectionBwt::runSearchesPrivate(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
						MatchSeqTest** seqTestArray, long* minOvlArray,
						AlignmentType alType, MatchesSought matchesSought,
						MinOvlStringency ovlStringency, AlignmentThreadedness threadedness){
  int numSenses;
  char senses[2];
  if (sense == '.'){
    numSenses = 2;
    senses[0] = '+';
    senses[1] = '-';
  } else if (sense == '+' or sense == '-'){
    numSenses = 1;
    senses[0] = sense;
  } else {
    cerr << '"' << sense << '"' << endl;
    throw AssemblyException::ArgError("SSCBwt:getBestMatchesPrivate bad sense char when '.' is allowed.");
  }

  long* bestScoreArray = new long[ numQueries+1 ];
  bool* lookForMatches = new bool[ numQueries+1 ];
  for (long n = 0; n < numQueries; ++n){
    // note: I tried optimizing this score once, it didn't help
    bestScoreArray[n] = _defaultMinScore;
    lookForMatches[n] = true;
  }

  long* correctedMinOvls = new long[ numQueries+1 ];
  if (ovlStringency == softMinOvl){
    for (long n = 0; n < numQueries; ++n){ correctedMinOvls[n] = long(float(minOvlArray[n]) * _fractId); }
  } else {
    for (long n = 0; n < numQueries; ++n){ correctedMinOvls[n] = minOvlArray[n]; }
  }

  // i am going to try and avoid seeking offsets in both senses for sequences that
  // already find matches in one sense if it is a "hasMatch" call
  OffsetTracker*** otArraysBySense[2];
  for (int sn = 0; sn < numSenses; ++sn){
    otArraysBySense[sn] = new OffsetTracker**[ numQueries+1 ];
    if (threadedness == threaded){
      long step = (numQueries / omp_get_max_threads()) + 2;
      #pragma omp parallel for schedule(dynamic)
      for (long start = 0; start < numQueries; start += step){
	long localNQ;
	if (start+step > numQueries){ localNQ = numQueries - start; }
	else { localNQ = step; }
	getOffsets(localNQ, start, seqArray, senses[sn], minOvlArray, seqTestArray, lookForMatches, alType, otArraysBySense[sn]);
      }
    } else {
      getOffsets(numQueries, 0, seqArray, senses[sn], minOvlArray, seqTestArray, lookForMatches, alType, otArraysBySense[sn]);
    }

    if ((matchesSought != allMatches) or (_gapMode == ungapped)){
      runAlignmentsHelper(numQueries, otArraysBySense[sn], matchesArray, seqArray, senses[sn], ungapped, seqTestArray,
			  correctedMinOvls, bestScoreArray, matchesSought, lookForMatches, threadedness);
    }

    // block new offsets from being generated for "hasMatch" if a match was found
    if (matchesSought == firstMatch){
      for (long n2 = 0; n2 < numQueries; ++n2){ lookForMatches[n2] = matchesArray[n2]->empty(); }
    }

    // delete the offsets early if they need not be used again
    if (_gapMode != gapped){
      for (int qn = 0; qn < numQueries; ++qn){
	long hitIndex = 0;
	while (otArraysBySense[sn][qn][hitIndex] != NULL){
	  delete otArraysBySense[sn][qn][hitIndex];
	  hitIndex++;
	}
	delete [] otArraysBySense[sn][qn];
      }
      delete [] otArraysBySense[sn];
    }
  }

  delete [] lookForMatches;

  if (_gapMode == gapped){
    // only get rid of the ungapped alignments if optimal alignments need be found
    if (matchesSought != firstMatch){
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
	matchesArray[n]->clear();
      }
    }
    for (int n = 0; n < numSenses; ++n){
      runAlignmentsHelper(numQueries, otArraysBySense[n], matchesArray, seqArray, senses[n], gapped, seqTestArray,
			  correctedMinOvls, bestScoreArray, matchesSought, lookForMatches, threadedness);
    }

    // delete the offset tracker arrays and their contents (gapped => not deleted before)
    for (int sn = 0; sn < numSenses; ++sn){
      for (int qn = 0; qn < numQueries; ++qn){
	long hitIndex = 0;
	while (otArraysBySense[sn][qn][hitIndex] != NULL){
	  delete otArraysBySense[sn][qn][hitIndex];
	  hitIndex++;
	}
	delete [] otArraysBySense[sn][qn];
      }
      delete [] otArraysBySense[sn];
    }
  }
  //OK();

  delete [] correctedMinOvls;
  delete [] bestScoreArray;
}



// this code looks redundant because it is optimized for speed
inline long ScoredSeqCollectionBwt::getOverlapFromOffset(long offset, long seqSizeA, long seqSizeB){
  if (offset < 0){
    if (offset + seqSizeB >= seqSizeA){
      return seqSizeA;
    } else {
      return offset + seqSizeB;
    }
  } else {
    if (offset + seqSizeB >= seqSizeA){
      return seqSizeA - offset;
    } else {
      return seqSizeB;
    }
  }
}


// REQUIRES: all members of "matches" (if any) have the same score
void ScoredSeqCollectionBwt::runAlignmentsHelper(long numQueries, OffsetTracker*** contigToOffsetsArray,
						 vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, GapMode gapMode,
						 MatchSeqTest** seqTestArray, long* minOvlArray, long* bestScoreArray,
						 MatchesSought matchesSought, bool* seekMatch, AlignmentThreadedness threadedness){

  bool justFirstMatch = (matchesSought == firstMatch);
  bool getAllMatches = (matchesSought == allMatches);

  long* initBestScoreArray = new long[ numQueries+1 ];
  long* otCountByQuery = new long[ numQueries+1 ];
  long totalOtCount = 0;

  bool* matchesStillSought = new bool[ numQueries+1 ];
  ScoredSeq** seqFlipArray = new ScoredSeq*[ numQueries+1 ];

  // this is so that all of the alignments produced here can be modified in the end
  vector<Alignment*>* tempMatchesArray = new vector<Alignment*>[ numQueries+1 ];
  // this is threaded by query, the best I can do at this point, but sets up for threading by hit
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numQueries; ++n){
      initBestScoreArray[n] = bestScoreArray[n];
      matchesStillSought[n] = ( (! justFirstMatch) or matchesArray[n]->empty() );
      otCountByQuery[n] = 0;
      if (matchesStillSought[n]){
	while (contigToOffsetsArray[n][otCountByQuery[n]] != NULL){ ++otCountByQuery[n]; }
      }
      seqFlipArray[n] = ScoredSeqFlip::getFlip(seqArray[n], sense);

      #pragma omp atomic
      totalOtCount += otCountByQuery[n];
    }
  } else {
    for (long n = 0; n < numQueries; ++n){
      initBestScoreArray[n] = bestScoreArray[n];
      matchesStillSought[n] = ( (! justFirstMatch) or matchesArray[n]->empty() );
      otCountByQuery[n] = 0;
      if (matchesStillSought[n]){
	while (contigToOffsetsArray[n][otCountByQuery[n]] != NULL){ ++otCountByQuery[n]; }
      }
      seqFlipArray[n] = ScoredSeqFlip::getFlip(seqArray[n], sense);
      totalOtCount += otCountByQuery[n];
    }
  }

  // this sets things up for threading across multiple hits per query by defining
  // the query and the target for each array element across an array that encompasses
  // all of the hits for all queries
  long totalOtCountP1 = totalOtCount + 1;
  long* useQueryIndex = new long[ totalOtCountP1 ];
  OffsetTracker** contigToOffsets = new OffsetTracker*[ totalOtCountP1 ];
  ScoredSeq** useSeqArray = new ScoredSeq*[ totalOtCountP1 ];
  vector<Alignment*>** useMatchesArray = new vector<Alignment*>*[ totalOtCountP1 ];
  long* useMinOvlArray = new long[ totalOtCountP1 ];
  ScoredSeq** useSeqFlipArray = new ScoredSeq*[ totalOtCountP1 ];
  long totalIndex = 0;

  // i am unsure about which is the better way to order these jobs (n outer or n2 outer) for threading
  for (long n = 0; n < numQueries; ++n){
    for (long n2 = 0; n2 < otCountByQuery[n]; ++n2){
      useQueryIndex[totalIndex] = n;
      contigToOffsets[totalIndex] = contigToOffsetsArray[n][n2];
      useSeqArray[totalIndex] = seqArray[n];
      useMatchesArray[totalIndex] = &tempMatchesArray[n];
      useMinOvlArray[totalIndex] = minOvlArray[n];
      useSeqFlipArray[totalIndex] = seqFlipArray[n];
      ++totalIndex;
    }
  }
  delete [] otCountByQuery;

  // make alignments
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long hitIndex = 0; hitIndex < totalIndex; ++hitIndex){
      // set up the alMaker to reflect the most up-to-date min score
      long oldBestScore;
      bool matchSought;
      #pragma omp critical (SSCBwtScoreStatus)
      {
	oldBestScore = bestScoreArray[useQueryIndex[hitIndex]];
	matchSought = matchesStillSought[useQueryIndex[hitIndex]];
      }
      if (matchSought){
	AlignmentMaker * alMaker = makeAlMaker(gapMode, oldBestScore);
	vector<Alignment*> tempMatches;
	long newBestScore = makeAlignmentHelper(contigToOffsets[hitIndex], useSeqArray[hitIndex], &tempMatches,
						useMinOvlArray[hitIndex], alMaker, useSeqFlipArray[hitIndex],
						sense, oldBestScore, true, gapMode, matchesSought);
        #pragma omp critical (SSCBwtScoreStatus)
	{
	dealWithMatches(&tempMatches, useMatchesArray[hitIndex],
			getAllMatches, justFirstMatch,
			useQueryIndex[hitIndex], bestScoreArray, newBestScore, matchesStillSought);
	}
	delete alMaker;
      }
    }
  } else {
    for (long hitIndex = 0; hitIndex < totalIndex; ++hitIndex){
      // set up the alMaker to reflect the most up-to-date min score
      long oldBestScore = bestScoreArray[useQueryIndex[hitIndex]];
      bool matchSought = matchesStillSought[useQueryIndex[hitIndex]];
      if (matchSought){
	AlignmentMaker * alMaker = makeAlMaker(gapMode, oldBestScore);
	vector<Alignment*> tempMatches;
	long newBestScore = makeAlignmentHelper(contigToOffsets[hitIndex], useSeqArray[hitIndex], &tempMatches,
						useMinOvlArray[hitIndex], alMaker, useSeqFlipArray[hitIndex],
						sense, oldBestScore, true, gapMode, matchesSought);
	dealWithMatches(&tempMatches, useMatchesArray[hitIndex],
			getAllMatches, justFirstMatch,
			useQueryIndex[hitIndex], bestScoreArray, newBestScore, matchesStillSought);
	delete alMaker;
      }
    }
  }

  // replace the sequences in the alignment
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numQueries; ++n){
      if (bestScoreArray[n] > initBestScoreArray[n]){
	for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
	matchesArray[n]->clear();
      }
    }
  } else {
    for (long n = 0; n < numQueries; ++n){
      if (bestScoreArray[n] > initBestScoreArray[n]){
	for (vector<Alignment*>::iterator it = matchesArray[n]->begin(); it != matchesArray[n]->end(); ++it){ delete *it; }
	matchesArray[n]->clear();
      }
    }
  }
  delete [] initBestScoreArray;

  if (sense == '+'){
    if (threadedness == threaded){
      #pragma omp parallel for schedule(dynamic)
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  (*it)->seqReplace(seqArray[n], (*it)->seqB());
	}
	matchesArray[n]->insert(matchesArray[n]->end(), tempMatchesArray[n].begin(), tempMatchesArray[n].end());
      }
    } else {
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  (*it)->seqReplace(seqArray[n], (*it)->seqB());
	}
	matchesArray[n]->insert(matchesArray[n]->end(), tempMatchesArray[n].begin(), tempMatchesArray[n].end());
      }
    }
  } else if (sense == '-'){
    if (threadedness == threaded){
      #pragma omp parallel for schedule(dynamic)
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  matchesArray[n]->push_back( (*it)->copyRcSeqA(seqArray[n]) );
	  delete *it;
	}
      }
    } else {
      for (long n = 0; n < numQueries; ++n){
	for (vector<Alignment*>::iterator it = tempMatchesArray[n].begin(); it != tempMatchesArray[n].end(); ++it){
	  matchesArray[n]->push_back( (*it)->copyRcSeqA(seqArray[n]) );
	  delete *it;
	}
      }
    }
  } else { throw AssemblyException::ArgError("SSCBWT:alignment method under construction, bad sense char"); }

  delete [] tempMatchesArray;
  delete [] matchesStillSought;

  delete [] useQueryIndex;
  delete [] contigToOffsets;
  delete [] useSeqArray;
  delete [] useMatchesArray;
  delete [] useMinOvlArray;
  delete [] useSeqFlipArray;

  for (long n = 0; n < numQueries; ++n){ delete seqFlipArray[n]; }
  delete [] seqFlipArray;
}



ScoredSeqCollectionBwt::AlignmentMaker* ScoredSeqCollectionBwt::makeAlMaker(GapMode gapMode, long bestScoreSoFar){
  if (gapMode == gapped){ return NULL; }
  else if (_penalizeEdgeGaps){ return new AlMakerPenalizeEdgeA(_scoreAsm, bestScoreSoFar, _maxFractMis); }
  else { return new AlMakerNoEdge(_scoreAsm, bestScoreSoFar, _maxFractMis); }
}

void ScoredSeqCollectionBwt::dealWithMatches(vector<Alignment*>* tempMatches, vector<Alignment*>* realMatches,
					     bool getAllMatches, bool justFirstMatch,
					     long index, long* bestScoreArray, long newBestScore, bool* matchSoughtArray){
  if (getAllMatches){
    realMatches->insert(realMatches->end(), tempMatches->begin(), tempMatches->end());
  } else {
    if (newBestScore < bestScoreArray[index]){
      for (vector<Alignment*>::iterator it = tempMatches->begin(); it != tempMatches->end(); ++it){ delete *it; }
    } else {
      if (newBestScore > bestScoreArray[index]){
	for (vector<Alignment*>::iterator it = realMatches->begin(); it != realMatches->end(); ++it){ delete *it; }
	realMatches->clear();
	bestScoreArray[index] = newBestScore;
      }
      realMatches->insert(realMatches->end(), tempMatches->begin(), tempMatches->end());
      if (justFirstMatch and tempMatches->size() > 0){ matchSoughtArray[index] = false; }
    }
  }
}

bool ScoredSeqCollectionBwt::SortWqByLength::operator() (WordQuery* wqA, WordQuery* wqB){
  return wqA->size() < wqB->size();
}




void ScoredSeqCollectionBwt::getOffsets(long numQ, long startQn, ScoredSeq** seqs, char sense, long* minOverlaps,
					MatchSeqTest** matchTests, bool* lookForMatches, AlignmentType alType, OffsetTracker*** resultArray){

  // i need to use this to get the correct modifiable array
  int threadNum = omp_get_thread_num();
  bool* localIndexWasHit = _contigIndexWasHit[threadNum];
  long* localIndex = _localContigIndex[threadNum];
  long* hitIndexes = new long[ _numContigLoci + 1 ];

  // if this is false, then the queries don't need to be saved and the memory load can be reduced more quickly
  // but even if it is true, the immutable form of the OffsetTrackers will have condensed, NR query sets
  bool queriesNeeded = true; //(_gapMode == gapped);

  // these will be used to check whether or not the WQBO needs to be re-generated
  WordQueryBinOrganizer* wqbo = NULL;
  long priorQuerySize;
  long priorMinOverlap;

  // this can also be re-used through the different queries
  OffsetTrackerMutable** hitArrayMutable = new OffsetTrackerMutable*[ _allSeqs.size() + 1 ];
  char* hitArrayMutableMemory = new char[ sizeof(OffsetTrackerMutable) * ( _allSeqs.size() + 1 ) ];

  // this creates a re-usable data structure
  long maxSeqSize = 0;
  for (long qN = startQn; qN < startQn + numQ; ++qN){
    if (seqs[qN]->size() > maxSeqSize){ maxSeqSize = seqs[qN]->size(); }
  }
  char* seqString = new char[ maxSeqSize + 1 ];

  //1
  for (long qN = startQn; qN < startQn + numQ; ++qN){

    // if I am not going to look for matches, just make a null array
    if (! lookForMatches[qN]){
      resultArray[qN] = new OffsetTracker*[1];
      resultArray[qN][0] = NULL;
    } 
    //2
    else {

      ScoredSeq* seq = seqs[qN];
      long minOverlap = minOverlaps[qN];
      MatchSeqTest* matchTest = matchTests[qN];

      // this is what I am collecting now so that I can sort the return values
      long currentIndex = 0;

      long querySize = seq->size();
      seq->gatherSeq(seqString,sense);

      // GET THE WORD QUERIES

      // if these conditions are true, the WQBO can just be reset
      if (wqbo!=NULL and querySize==priorQuerySize and minOverlap==priorMinOverlap){
	wqbo->reset();
      }
      //3
      else {
	if (wqbo != NULL){ wqbo->deleteWithWordQueries(); }
	// get the size and position of the sub-strings for which matches will be sought
	// here is the block data structure that will be used to keep track of redundant queries
	if (alType == semiGlobal){
	  wqbo = getSubtextWordQueries(seq, seqString, sense, minOverlap);
	} else if (alType == fullAlignment){
	  wqbo = getSubtextWordQueriesFullAlignment(seq, seqString, sense, minOverlap);
	} else {
	  throw AssemblyException::ArgError("SSCBwt:getOffsets got a bad alType");
	}
	priorQuerySize = querySize;
	priorMinOverlap = minOverlap;
      }
      //3

      // here is where i get rid of the seeds with Ns
      subtextWordQueriesFilterNs(seq, seqString, wqbo);

      // this is blocked separately so that i can delete the seqString afterwards
      long* queryNumArray;
      bool qnaCreated;

      if (wqbo->numValidQueries() == 0){
	qnaCreated = false;
      }
      //3
      else {
	queryNumArray = new long[querySize];
	burrowsWheelerUtilities::makeNumString( seqString, querySize, 0, queryNumArray );
	qnaCreated = true;


	// FIND THE HITS/OFFSETS
	WordQuery** originalWq = wqbo->getQueries();

	// sort word queries short to long - doing it by sorting the lengths, since those will
	// be highly redundant
	long maxWqLen = 0;
	long numWq = 0;
	while (originalWq[numWq] != NULL){
	  if (originalWq[numWq]->size() > maxWqLen){ maxWqLen = originalWq[numWq]->size(); }
	  ++numWq;
	}

	long* qlenArrayMemory = new long[maxWqLen * 2 + 2];
	long* qlenToStart = &qlenArrayMemory[0];
	long* qlenToCtSparse = &qlenArrayMemory[maxWqLen + 1];

	for (long n = 0; n < numWq; ++n){ qlenToCtSparse[originalWq[n]->size()] = 0; }
	vector<long> sortedNrQlens;
	for (long n = 0; n < numWq; ++n){
	  long qlen = originalWq[n]->size();
	  if (qlenToCtSparse[qlen] == 0){ sortedNrQlens.push_back(qlen); }
	  qlenToCtSparse[qlen]++;
	}
	sort(sortedNrQlens.begin(), sortedNrQlens.end());

	long curStart = 0;
	for (vector<long>::iterator qlIt = sortedNrQlens.begin(); qlIt != sortedNrQlens.end(); ++qlIt){
	  long ql = *qlIt;
	  qlenToStart[ql] = curStart;
	  curStart += qlenToCtSparse[ql];
	}

	WordQuery** sortedWqArray = new WordQuery*[ numWq + 1 ];
	for (long n = 0; n < numWq; ++n){
	  long qlen = originalWq[n]->size();
	  sortedWqArray[qlenToStart[qlen]] = originalWq[n];
	  qlenToStart[qlen]++;
        }
	delete [] qlenArrayMemory;
	delete [] originalWq;
	// now i have the word queries sorted shortest to longest

	// DO NOT optimize by checking that there are still queries - I need to delete these anyway
	//for (vector<WordQuery*>::iterator wqIt = sortedWordQueries.begin(); wqIt != sortedWordQueries.end(); ++wqIt){

	long matchArraySize = 1000;
	long* matchArray = new long[matchArraySize+1];
	matchArray[0] = 0;

	//4
        for (long wqN = 0; wqN < numWq; ++wqN){
	  WordQuery* wordQuery = sortedWqArray[wqN];
	  long queryIndex = wordQuery->getIndex();

	  // make sure that the query is still valid before doing the search
	  if (wqbo->queryStillValid( queryIndex )){

	    // these will get used a lot
	    long stSize = wordQuery->size();
	    long startNum = wordQuery->start();

	    // do the search with the sub-query word
	    // deal with possible redundant matchLoci if multiple words match
	    matchArray = burrowsWheelerUtilities::bwtCountAndFind(queryNumArray, startNum, stSize,
								  _bwTransform, _sortedSuffixes, _textSize,
								  _tableC, _occCounts, matchArray, matchArraySize);

	    // keep track of the starts AND how many non-overlapping windows have matched
	    long numMatches = matchArray[0];
	    //cout << numMatches << " " << matchArraySize << endl;
	    if (numMatches == 0){ wqbo->removeSupersetQueries(queryIndex); }
	    else {
	      if (numMatches > matchArraySize){ matchArraySize = numMatches; }
	      currentIndex = getOffsetsSub2(currentIndex, wordQuery, seq, matchArray, hitArrayMutable, hitArrayMutableMemory, hitIndexes, localIndex, localIndexWasHit, matchTest);
	    }
	  }
	}
	delete [] matchArray;
	delete [] sortedWqArray;
      }
      //3

      OffsetTracker** sortedHitWithOffsets = new OffsetTracker*[ currentIndex + 1 ];
      // the last index is a pair of nulls so that I don't have to return a length
      sortedHitWithOffsets[currentIndex] = NULL;

      // OffsetTrackers will be sorted, but that only needs to be gone through if there ARE OffsetTrackers
      if (currentIndex > 0){
	// BIG operation abstracted out into sub-function
	getOffsetsSub1(currentIndex, hitArrayMutable, sortedHitWithOffsets, queriesNeeded, hitIndexes, localIndexWasHit);
      }
      if (qnaCreated){ delete [] queryNumArray; }
      resultArray[qN] = sortedHitWithOffsets;
    }
    //2
  }
  //1
  delete [] seqString;
  if (wqbo != NULL){ wqbo->deleteWithWordQueries(); }
  delete [] hitIndexes;
  delete [] hitArrayMutable;
  delete [] hitArrayMutableMemory;
}



void ScoredSeqCollectionBwt::getOffsetsSub1(long currentIndex, OffsetTrackerMutable** hitArrayMutable, OffsetTracker** sortedHitWithOffsets,
					    bool queriesNeeded, long* hitIndexes, bool* localIndexWasHit){
  // ORGANIZE THE HITS/OFFSETS

  // this value (hitArrayMutable[n]->maxOffsetCount()) is obtained many times below, so
  // this array allows it to be obtained only once
  long* hamMoc = new long[ currentIndex ];

  // find the highest hit count so i can make an appropriately-sized array
  long maxHitCount = 0;
  for (long n = 0; n < currentIndex; ++n){
    hamMoc[n] = hitArrayMutable[n]->maxOffsetCount();
    if (hamMoc[n] > maxHitCount){ maxHitCount = hamMoc[n]; }
  }
  long* hitCtArrayMemory = new long[ maxHitCount * 2 + 2 ];
  long* hitCtToCtSparse = &hitCtArrayMemory[0];
  long* hitCountsToStart = &hitCtArrayMemory[maxHitCount+1];


  // initialize an array of the number of OffsetTrackers with each hit count, then add up the values
  for (long n = 0; n < currentIndex; ++n){ hitCtToCtSparse[hamMoc[n]] = 0; }

  // get the number of OffsetTrackers with each hit count, and create a non-redundant sorted vector of hit count values
  vector<long> sortedNrHitCts;
  for (long n = 0; n < currentIndex; ++n){
    if (hitCtToCtSparse[hamMoc[n]]==0){ sortedNrHitCts.push_back( hamMoc[n] ); }
    hitCtToCtSparse[hamMoc[n]] += 1;
    // while I'm at it, I need to reset the local array of is-it-found values
    localIndexWasHit[ hitIndexes[n] ] = false;
  }
  sort(sortedNrHitCts.begin(), sortedNrHitCts.end());

  // now get the final-array starting positions for OffsetTrackers with each given # hits
  long currentStart = currentIndex;
  for (vector<long>::iterator hcIt = sortedNrHitCts.begin(); hcIt != sortedNrHitCts.end(); ++hcIt){
    long hcN = *hcIt;
    currentStart -= hitCtToCtSparse[hcN];
    hitCountsToStart[hcN] = currentStart;
  }

  // now, create the new sorted list and fill it
  for (long n = 0; n < currentIndex; ++n){
    // COPY HERE AND DELETE OLD OT
    sortedHitWithOffsets[hitCountsToStart[hamMoc[n]]] = hitArrayMutable[n]->immutableCopy(queriesNeeded);
    hitCountsToStart[hamMoc[n]]++;
    hitArrayMutable[n]->~OffsetTrackerMutable();
  }

  delete [] hitCtArrayMemory;
  delete [] hamMoc;
}


long ScoredSeqCollectionBwt::getOffsetsSub2(long currentIndex, WordQuery* wordQuery, ScoredSeq* seq, long* matchArray, OffsetTrackerMutable** hitArrayMutable, 
					    char* hitArrayMutableMemory, long* hitIndexes, long* localIndex, bool* localIndexWasHit, MatchSeqTest* matchTest){

  long stSizeM1 = wordQuery->size() - 1;
  long startNum = wordQuery->start();

  long numMatches = matchArray[0];

  for (long n = 1; n <= numMatches; ++n){
    long hitStart = matchArray[n];
    ScoredSeqLocus* contigLocus = _orderedContigLoci[ getOverlapContigIndex(hitStart, hitStart + stSizeM1) ];
    ScoredSeq* contig = contigLocus->_seq;

    if (seq != contig and matchTest->seqIsOk(contig) ){ // don't align to self!!!
      // determine the offset for the alignment of the query to this seq;
      // number is the start position of target in seq (i.e. a seq coord)
      long contigStart = hitStart - contigLocus->_start;

      // make sure that the hit is on an appropriate region of the target contig
      if ( wordQuery->acceptableOverlap(contigStart, contig) ){

	long offset = startNum - contigStart;
	long contigIndex;

	long contigN = contigLocus->_index;
	if (localIndexWasHit[contigN]){ contigIndex = localIndex[contigN]; }
	else {
	  localIndexWasHit[contigN] = true;
	  contigIndex = currentIndex;
	  char* hamMem = &hitArrayMutableMemory[ sizeof(OffsetTrackerMutable) * contigIndex ];
	  hitArrayMutable[contigIndex] = new (hamMem) OffsetTrackerMutable(contig);
	  localIndex[contigN] = contigIndex;
	  hitIndexes[contigIndex] = contigN;
	  currentIndex++;
	}
	hitArrayMutable[contigIndex]->addOffset( offset, wordQuery, numMatches );
      }
    }
  }
  return currentIndex;
}


ScoredSeqCollectionBwt::WordQueryBinOrganizer* ScoredSeqCollectionBwt::getSubtextWordQueries(ScoredSeq* querySeq, char* seqSeq, char sense, long minOverlap){
  // local so that filterN's can be done in here but not repeated
  vector<WordQuery*> localQueries;

  // this is the whole-sequence case
  long blockStartPos = 0;
  long blockEndPos = querySeq->size();
  getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, false);

  // these are the edge cases
  if (_edgeScalingEnabled){
    long fragmentLen = minOverlap;
    long qSize = querySeq->size();
    while (fragmentLen < qSize){
      // the front end
      blockStartPos = 0;
      blockEndPos = fragmentLen;
      getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, true);
      // the back end
      blockStartPos = qSize - fragmentLen;
      blockEndPos = qSize;
      getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, true);
      // scale up the window
      fragmentLen *= 2;
    }
  }

  return new WordQueryBinOrganizer( &localQueries );
}


ScoredSeqCollectionBwt::WordQueryBinOrganizer* ScoredSeqCollectionBwt::getSubtextWordQueriesFullAlignment(ScoredSeq* querySeq, char* seqSeq, char sense, long minTargetLength){
  // local so that filterN's can be done in here but not repeated
  vector<WordQuery*> localQueries;

  // this is the whole-sequence case; it will allow hits to the entire query to be aligned
  long blockStartPos = 0;
  long blockEndPos = querySeq->size();
  getSubtextWordQueriesHelper(querySeq, &localQueries, sense, blockStartPos, blockEndPos, true);

  long fragmentLen = minTargetLength / 2;
  while (fragmentLen < querySeq->size()){
    getSubtextWordQueriesHelper(querySeq, &localQueries, sense, fragmentLen);
    // scale up the window
    fragmentLen *= 2;
  }
  return new WordQueryBinOrganizer( &localQueries );
}


void ScoredSeqCollectionBwt::subtextWordQueriesFilterNs(ScoredSeq* querySeq, char* seqSeq, WordQueryBinOrganizer* wqbo){
  if (wqbo->numValidQueries() > 0){
    // now do the filtering
    long numBlocks = wqbo->numBlocks();
    for (long block = 0; block < numBlocks; ++block){
      if (wqbo->blockHasQueries(block)){
	// scan for N
	long pos = wqbo->blockStart(block);
	long blockEnd = wqbo->blockEnd(block);
	while (pos < blockEnd and seqSeq[pos] != 'N'){ ++pos; }
	if (pos != blockEnd){
	  // an N was found, so all overlapping windows can be gotten rid of!
	  long* queriesInBin = wqbo->getBlockQueryIndexes(block);
	  long qibCount = wqbo->getBlockQueryIndexCount(block);
	  for (long qN = 0; qN < qibCount; ++qN){
	    wqbo->removeQuery( queriesInBin[qN], block );
	    //delete wqbo->getQuery( queriesInBin[qN] );
	  }
	  delete [] queriesInBin;
	}
      }
    }
  }
}



ScoredSeqCollectionBwt::WordQueryBinOrganizer::WordQueryBinOrganizer(vector<WordQuery*>* wordQueries){
  _numIndexes = wordQueries->size();
  _numValidQueries = _numIndexes;
  long numQp1 = _numIndexes + 1; // just an optimization; it is used more later in this constructor
  _queries = new WordQuery*[ numQp1 ];
  _queryIsValid = new bool[ numQp1 ];
  _invalidQueries = new long[ numQp1 ];
  vector<WordQuery*>::iterator inputIt = wordQueries->begin();
  for (long n = 0; n < _numIndexes; ++n){
    _queryIsValid[n] = true;
    _queries[n] = *inputIt;
    _queries[n]->setIndex(n);
    ++inputIt;
  }

  // get all the non-redundant edges sorted
  set<long> edgeSetStarts;
  set<long> edgeSetEnds;
  long maxEdge = 0;
  for (long n = 0; n < _numIndexes; ++n){
    edgeSetStarts.insert( _queries[n]->start() );
    long end = _queries[n]->end();
    if (end > maxEdge){ maxEdge = end; }
    edgeSetEnds.insert( end );
  }

  // only the elements that are actually starts will be occupied
  _startToBlock = new long[ maxEdge+1 ];

  // this can save time because there should be a lot of redundancy between
  // starts, and also between ends, but not combining starts and ends
  set<long> edgeSet;
  edgeSet.insert(edgeSetStarts.begin(), edgeSetStarts.end());
  edgeSet.insert(edgeSetEnds.begin(), edgeSetEnds.end());

  // use the edges to construct blocks; they will be ordered in the set
  long edgeSetSizeP1 = edgeSet.size() + 1;
  _blockToStart = new long[ edgeSetSizeP1 ];
  _blockToEnd = new long[ edgeSetSizeP1 ];

  _blockToBinStart = new long[ edgeSetSizeP1 ];
  _blockToBinEnd = new long[ edgeSetSizeP1 ];
  _blockToAllQisInBin = new long[ (edgeSetSizeP1 - 1) * (numQp1 - 1) + 1 ];
  _blockToAllQiBinIndex = new long[ (edgeSetSizeP1 - 1) * (numQp1 - 1) + 1 ];

  _numBlocks = edgeSet.size() - 1;
  if (_numBlocks < 0){ _numBlocks = 0; }

  long block = 0;
  long blockStart = 0;
  for (set<long>::iterator it = edgeSet.begin(); it != edgeSet.end(); ++it){
    _blockToStart[block] = *it;
    if (block > 0){ _blockToEnd[block-1] = *it; }
    _startToBlock[*it] = block;
    _blockToBinStart[block] = blockStart;
    _blockToBinEnd[block] = blockStart;
    ++block;
    blockStart += _numIndexes;
  }

  // insert each word query into all of its blocks
  for (long n = 0; n < _numIndexes; ++n){
    long nb = _startToBlock[ _queries[n]->start() ];
    long endEdge = _startToBlock[ _queries[n]->end() ];
    while (nb < endEdge){
      _blockToAllQisInBin[ _blockToBinEnd[nb] ] = n;
      _blockToAllQiBinIndex[_blockToBinStart[nb] + n] = _blockToBinEnd[nb];
      _blockToBinEnd[nb]++;
      ++nb;
    }
  }
}
ScoredSeqCollectionBwt::WordQueryBinOrganizer::~WordQueryBinOrganizer(){
  delete [] _blockToStart;
  delete [] _blockToEnd;
  delete [] _startToBlock;

  delete [] _blockToBinStart;
  delete [] _blockToBinEnd;
  delete [] _blockToAllQisInBin;
  delete [] _blockToAllQiBinIndex;

  delete [] _queries;
  delete [] _queryIsValid;
  delete [] _invalidQueries;
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::deleteWithWordQueries(){
  for (long n = 0; n < _numIndexes; ++n){ delete _queries[n]; }
  delete this;
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::reset(){
  // make all of the queries valid again
  for (long n = _numIndexes; n > _numValidQueries; --n){ _queryIsValid[_invalidQueries[n]] = true; }
  //for (long n = 0; n < _numIndexes; ++n){ _queryIsValid[n] = true; }
  _numValidQueries = _numIndexes;


  // insert each word query into all of its blocks and reset the counts
  for (long block = 0; block < _numBlocks; ++block){ _blockToBinEnd[block] = _blockToBinStart[block]; }
  for (long n = 0; n < _numIndexes; ++n){
    long nb = _startToBlock[ _queries[n]->start() ];
    long endEdge = _startToBlock[ _queries[n]->end() ];
    while (nb < endEdge){
      _blockToAllQisInBin[ _blockToBinEnd[nb] ] = n;
      _blockToAllQiBinIndex[_blockToBinStart[nb] + n] = _blockToBinEnd[nb];
      _blockToBinEnd[nb]++;
      ++nb;
    }
  }
}




long ScoredSeqCollectionBwt::WordQueryBinOrganizer::numIndexes(){ return _numIndexes; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::numValidQueries(){ return _numValidQueries; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::numBlocks(){ return _numBlocks; }
bool ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockHasQueries(long block){
  return _blockToBinEnd[block] > _blockToBinStart[block];
}
bool ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockHasQuery(long block, long qIndex){
  return _queryIsValid[qIndex] and _blockToStart[block] >= _queries[qIndex]->start() and _blockToEnd[block] <= _queries[qIndex]->end();
}
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockStart(long block){ return _blockToStart[block]; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::blockEnd(long block){ return _blockToEnd[block]; }
long ScoredSeqCollectionBwt::WordQueryBinOrganizer::getBlockQueryIndexCount(long block){
  return _blockToBinEnd[block] - _blockToBinStart[block];
}
long* ScoredSeqCollectionBwt::WordQueryBinOrganizer::getBlockQueryIndexes(long block){
  long count = _blockToBinEnd[block] - _blockToBinStart[block];
  long* indexes = new long[ count + 1 ];
  long* shiftedArray = &_blockToAllQisInBin[ _blockToBinStart[block] ];
  for (long n = 0; n < count; ++n){ indexes[n] = shiftedArray[n]; }
  return indexes;
}
// PRIVATE USE ONLY!!!! EXPOSES REP!!!
long* ScoredSeqCollectionBwt::WordQueryBinOrganizer::viewBlockQueryIndexes(long block){
  return &_blockToAllQisInBin[ _blockToBinStart[block] ];
}


bool ScoredSeqCollectionBwt::WordQueryBinOrganizer::queryStillValid(long qIndex){
  return _queryIsValid[qIndex];
}
ScoredSeqCollectionBwt::WordQuery* ScoredSeqCollectionBwt::WordQueryBinOrganizer::getQuery(long qIndex){
  return _queries[qIndex];
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::removeQuery(long qIndex){
  if (_queryIsValid[qIndex]){
    long nb = _startToBlock[ _queries[qIndex]->start() ];
    long endEdge = _startToBlock[ _queries[qIndex]->end() ];
    long* shiftedArray = &_blockToAllQiBinIndex[qIndex];
    while (nb < endEdge){
      _blockToBinEnd[nb]--;
      if (_blockToBinEnd[nb] > _blockToBinStart[nb]){
	// this line moves the last element in the array of qIndexes to the position of the one
	// being removed, so that all active qIndexes are now found at the front of the array
	long movedQi = _blockToAllQisInBin[ _blockToBinEnd[nb] ];
	//long newIndex = _blockToAllQiBinIndex[_blockToBinStart[nb] + qIndex];
	long newIndex = shiftedArray[_blockToBinStart[nb]];
	_blockToAllQisInBin[newIndex] = movedQi;
	// since the qIndex above was moved, the array that helps find it must also be updated
	_blockToAllQiBinIndex[movedQi] = newIndex;
      }
      ++nb;
    }
    _queryIsValid[ qIndex ] = false;
    _invalidQueries[ _numValidQueries ] = qIndex;
    --_numValidQueries;
  }
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::removeQuery(long qIndex, long currentBlock){
  if (_queryIsValid[ qIndex ]){
    _queryIsValid[ qIndex ] = false;
    _invalidQueries[ _numValidQueries ] = qIndex;
    --_numValidQueries;
    long nb = _startToBlock[ _queries[qIndex]->start() ];
    if (nb <= currentBlock){ nb = currentBlock + 1; }
    long endEdge = _startToBlock[ _queries[qIndex]->end() ];
    long* shiftedArray = &_blockToAllQiBinIndex[qIndex];
    while (nb < endEdge){
      _blockToBinEnd[nb]--;
      if (_blockToBinEnd[nb] > _blockToBinStart[nb]){
	long movedQi = _blockToAllQisInBin[ _blockToBinEnd[nb] ];
	//long newIndex = _blockToAllQiBinIndex[_blockToBinStart[nb] + qIndex];
	long newIndex = shiftedArray[_blockToBinStart[nb]];
	_blockToAllQisInBin[newIndex] = movedQi;
	// since the qIndex above was moved, the array that helps find it must also be updated
	_blockToAllQiBinIndex[movedQi] = newIndex;
      }
      ++nb;
    }
  }
}
void ScoredSeqCollectionBwt::WordQueryBinOrganizer::removeSupersetQueries(long qIndex){
  long firstBin = _startToBlock[ _queries[qIndex]->start() ];
  long lastBin = _startToBlock[ _queries[qIndex]->end() ] - 1;
  //long* inFirstBin = getBlockQueryIndexes(firstBin);
  long* inFirstBin = viewBlockQueryIndexes(firstBin);
  long qiCount = getBlockQueryIndexCount(firstBin);
  if (firstBin==lastBin){
    for (long n = 0; n < qiCount; ++n){ removeQuery( inFirstBin[n] ); }
  } else {
    for (long n = 0; n < qiCount; ++n){
      if ( blockHasQuery(lastBin, inFirstBin[n]) ){ removeQuery( inFirstBin[n] ); }
    }
  }
  //delete [] inFirstBin;
}

// element after last valid entry of returned array is NULL
ScoredSeqCollectionBwt::WordQuery** ScoredSeqCollectionBwt::WordQueryBinOrganizer::getQueries(){
  WordQuery** wqa = new WordQuery*[ _numIndexes + 1 ];
  long aN = 0;
  for (long n = 0; n < _numIndexes; ++n){
    if (_queryIsValid[n]){ wqa[aN] = _queries[n]; ++aN; }
  }
  wqa[aN] = NULL;
  return wqa;
}



int ScoredSeqCollectionBwt::_specialCaseLenDenom[6] =     { 1,  2,  3,  3,  4,  5};//, 6, 7, 8, 9, 10};
int ScoredSeqCollectionBwt::_specialCaseStepDenom[6] =    { 1,  1,  1,  2,  2,  1};//, 1, 1, 1, 1, 1};
// product of the two above
int ScoredSeqCollectionBwt::_specialCaseStartDenom[6] =   { 1,  2,  3,  6,  8,  5};//, 6, 7, 8, 9, 10};
// num allowed errors:                                     0   1   2   3   4   5   6   7   8   9
int ScoredSeqCollectionBwt::_specialCaseWindowNum[170] = {1,   2,  3,  5,  7,  5,  6,  6,  7,  7, // 0s
                                                          8,   8,  8,  9,  9, 10, 10, 11, 11, 11, // 10s
							  12, 12, 13, 13, 13, 14, 14, 14, 15, 15, // 20s
							  15, 16, 16, 16, 17, 17, 17, 18, 18, 18, // 30s
                                                          19, 19, 19, 20, 20, 20, 21, 21, 21, 21, // 40s
                                                          22, 22, 22, 23, 23, 23, 24, 24, 24, 24, // 50s
                                                          25, 25, 25, 26, 26, 26, 27, 27, 27, 27, // 60s
                                                          28, 28, 28, 29, 29, 29, 29, 30, 30, 30, // 70s
                                                          31, 31, 31, 31, 32, 32, 32, 32, 33, 33, // 80s
                                                          33, 34, 34, 34, 34, 35, 35, 35, 35, 36, // 90s
                                                          36, 36, 36, 37, 37, 37, 38, 38, 38, 38, // 100s
                                                          39, 39, 39, 39, 40, 40, 40, 40, 41, 41, // 110s
                                                          41, 41, 42, 42, 42, 42, 43, 43, 43, 43, // 120s
							  44, 44, 44, 45, 45, 45, 45, 46, 46, 46, // 130s
                                                          46, 47, 47, 47, 47, 48, 48, 48, 48, 49, // 140s
                                                          49, 49, 49, 49, 50, 50, 50, 50, 51, 51, // 150s
                                                          51, 51, 52, 52, 52, 52, 53, 53, 53, 53, // 160s
};

inline long ScoredSeqCollectionBwt::getNumWindowsHelper(long numAllowedMis){
  if (numAllowedMis < 170){ return _specialCaseWindowNum[numAllowedMis]; }
  else { return numAllowedMis / 3; }
}


void ScoredSeqCollectionBwt::getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
							 long blockStartPos, long blockEndPos, bool useLimits){
  long seqLen = blockEndPos - blockStartPos;
  bool fiveEdge;
  if (blockStartPos == 0){ fiveEdge = true; }
  else if (blockEndPos == querySeq->size()){ fiveEdge = false; }
  else { throw AssemblyException::ArgError("SSCBwt::getSWQHelper block is at neither the front nor back of the query."); }

  long numAllowedMis = long(float(seqLen) * _maxFractMis);
  // a safety precaution for a rediculous situation
  if (numAllowedMis==seqLen){ numAllowedMis -= 1; }

  // figure these out for the special or general cases
  long startUpperLimit;
  long startFactor;
  long subLen;
  long fragLenT2;

  // special cases
  if (numAllowedMis < 6){
    // OLD WAY - new way is the same, but value is pre-computed by hand
    //long startDenom = _specialCaseLenDenom[numAllowedMis] * _specialCaseStepDenom[numAllowedMis];
    //long startDenom = _specialCaseStartDenom[numAllowedMis];
    subLen = seqLen / _specialCaseLenDenom[numAllowedMis];
    long numWindows = getNumWindowsHelper(numAllowedMis);
    //startFactor = seqLen / startDenom;
    startFactor = seqLen / _specialCaseStartDenom[numAllowedMis];
    startUpperLimit = blockStartPos + (numWindows * startFactor);
  } else { // general cases
    long numWindows = getNumWindowsHelper(numAllowedMis);
    subLen = seqLen / numWindows;
    if (subLen < 2){ subLen = 2; }
    startFactor =  seqLen / numWindows;
    // i am putting the booleans outside the loop to save time
    // the input args below with the equation calculates the start position
    startUpperLimit = blockStartPos + (numWindows * startFactor);
  }

  // now create the word queries using the parameters calculated above
  if (! useLimits){
    for (long start = blockStartPos; start < startUpperLimit; start += startFactor){
      wordQuerySet->push_back(new WordQueryNoEdgeLimit(start, subLen, sense));
    }
  } else if (fiveEdge){
    for (long start = blockStartPos; start < startUpperLimit; start += startFactor){
      wordQuerySet->push_back(new WordQueryFiveEdgeLimit(start, subLen, sense, seqLen));
    }
  } else {
    for (long start = blockStartPos; start < startUpperLimit; start += startFactor){
      wordQuerySet->push_back(new WordQueryThreeEdgeLimit(start, subLen, sense, seqLen));
    }
  }
}



void ScoredSeqCollectionBwt::getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
							 long fragmentLength){

  long numAllowedMis = long(float(fragmentLength) * _maxFractMis);
  // a safety precaution for a rediculous situation
  if (numAllowedMis==fragmentLength){ numAllowedMis -= 1; }
  long querySize = querySeq->size();

  // figure these out for the special or general cases
  long startUpperLimit;
  long startFactor;
  long subLen;
  long fragLenT2;

  // special cases
  if (numAllowedMis < 6){
    //long startDenom = _specialCaseLenDenom[numAllowedMis] * _specialCaseStepDenom[numAllowedMis];
    subLen = fragmentLength / _specialCaseLenDenom[numAllowedMis];
    long numWindows = getNumWindowsHelper(numAllowedMis) * querySize / fragmentLength;
    //startFactor = fragmentLength / startDenom;
    startFactor = fragmentLength / _specialCaseStartDenom[numAllowedMis];
    fragLenT2 = fragmentLength * 2;
    // increment up using the startFactor
    startUpperLimit = numWindows * startFactor;
    // the highest word cannot read off the edge of the query sequence
    if (startUpperLimit + subLen > querySize + 1){ startUpperLimit = querySize + 1 - subLen; }
  } else { // general cases
    long numWindows = getNumWindowsHelper(numAllowedMis);
    // adjust for the query length
    subLen = fragmentLength / numWindows;
    if (subLen < 2){ subLen = 2; }
    numWindows = numWindows * querySize / fragmentLength;
    // other constants that only need be calculated once
    fragLenT2 = fragmentLength * 2;
    startFactor = fragmentLength / numWindows;
    // increment up using the startFactor
    startUpperLimit = numWindows * startFactor;
    // the highest word cannot read off the edge of the query sequence
    if (startUpperLimit + subLen > querySize + 1){ startUpperLimit = querySize + 1 - subLen; }
  }

  // now create the windows using the variables calculated above
  for (long start = 0; start < startUpperLimit; start += startFactor){
    wordQuerySet->push_back(new WordQueryFullTarget(querySize, start, subLen, sense, _fractId, fragLenT2));
  }

}





void ScoredSeqCollectionBwt::bufferSeqs(){
  //OK();
  for ( set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it ){ (*it)->buffer(); }
  //OK();
}

void ScoredSeqCollectionBwt::unbufferSeqs(){
  //OK();
  for ( set<ScoredSeq*>::iterator it = _allSeqs.begin(); it != _allSeqs.end(); ++it ){ (*it)->unbuffer(); }
  //OK();
}



long ScoredSeqCollectionBwt::size(){
  //OK();
  return _allSeqs.size();
}


// overlaps is the return value, locus is the query
// ASSUMES that the locus only overlaps one contig
long ScoredSeqCollectionBwt::getOverlapContigIndex(long start, long end){

  // determine the bin to look in
  long startBin = start / _binSize;
  if (startBin < 0){ startBin = 0; }

  // return the legit match's index
  vector<ScoredSeqLocus*>::iterator firstLocusIt = _binToSeqs[startBin].begin();
  while (firstLocusIt != _binToSeqs[startBin].end() and (! (*firstLocusIt)->overlaps(start,end))){ ++firstLocusIt; }
  if ( firstLocusIt == _binToSeqs[startBin].end() ){ throw AssemblyException::LogicError("SSCBwt::gOCI can't find match."); }
  return (*firstLocusIt)->_index;
}


ScoredSeqCollectionBwt::ScoredSeqLocus::ScoredSeqLocus(){}
ScoredSeqCollectionBwt::ScoredSeqLocus::ScoredSeqLocus(ScoredSeq* seq, long start, long end, long index) :
  _start(start),
  _end(end),
  _seq(seq),
  _index(index) {
  }
bool ScoredSeqCollectionBwt::ScoredSeqLocus::overlaps(long otherStart, long otherEnd){
  //return (! (_start > otherEnd or otherStart > _end) );
  return (_start <= otherEnd and otherStart <= _end);
}
bool ScoredSeqCollectionBwt::ScoredSeqLocus::contains(long otherStart, long otherEnd){
  //return (! (_start > otherStart or otherEnd > _end) );
  return (_start <= otherStart and otherEnd <= _end);
}



ScoredSeqCollectionBwt::AlignmentMaker::~AlignmentMaker(){}

inline int ScoredSeqCollectionBwt::AlignmentMaker::scoreNucMisMatch(char nA, char nB){
  switch(nA) {
    // allowed nucleotide characters
  case 'N': return 0;
  case 'A':
    switch(nB) {
    case 'N': return 0;
    case 'A': return 0;
    default: return 1;
    }
  case 'C':
    switch(nB) {
    case 'N': return 0;
    case 'C': return 0;
    default: return 1;
    }
  case 'G':
    switch(nB) {
    case 'N': return 0;
    case 'G': return 0;
    default: return 1;
    }
  case 'T':
    switch(nB) {
    case 'N': return 0;
    case 'T': return 0;
    default: return 1;
    }
  default:
    throw AssemblyException::ArgError("nucleotide seq contains invalid char");
  }
}


bool ScoredSeqCollectionBwt::AlignmentMaker::makeAlignmentHelperPlus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB,
								     long bStart, long maxNumMis, long exBlockCount,
								     long* examineStarts, long* examineEnds){

  long subBinSize = 25;
  char fragA[subBinSize+1]; 
  char fragB[subBinSize+1]; 

  // gaps are irrelevant
  long misCount = 0;

  // initial block is defined by the last block element
  --exBlockCount;

  long startA = examineStarts[exBlockCount];
  long startB = bStart - aStart + startA;
  long endA = examineEnds[exBlockCount];
  long blockSize = endA - startA;
  long blockSizeMsub = blockSize - subBinSize;
  long startN = 0;
  while (startN < blockSizeMsub and maxNumMis >= misCount){
     seqA->gatherSubseq(&fragA[0], startA+startN, subBinSize, '+');
     seqB->gatherSubseq(&fragB[0], startB+startN, subBinSize, '+');
     for (long n2 = 0; n2 < subBinSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
     startN += subBinSize;
  }
  if (maxNumMis >= misCount){
    long endSize = blockSize % subBinSize;
    seqA->gatherSubseq(&fragA[0], startA+startN, endSize, '+');
    seqB->gatherSubseq(&fragB[0], startB+startN, endSize, '+');
    for (long n2 = 0; n2 < endSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
  }

  while ( (maxNumMis >= misCount) and (exBlockCount > 0) ){
    --exBlockCount;

    if (startA < examineEnds[exBlockCount]){ endA = startA; }
    else { endA = examineEnds[exBlockCount]; }
    startA = examineStarts[exBlockCount];
    blockSize = endA - startA;
    if (blockSize > 0){
      startB = bStart - aStart + startA;
      blockSizeMsub = blockSize - subBinSize;
      startN = 0;
      while (startN < blockSizeMsub and maxNumMis >= misCount){
        seqA->gatherSubseq(&fragA[0], startA+startN, subBinSize, '+');
        seqB->gatherSubseq(&fragB[0], startB+startN, subBinSize, '+');
        for (long n2 = 0; n2 < subBinSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
        startN += subBinSize;
      }
      if (maxNumMis >= misCount){
	long endSize = blockSize % subBinSize;
	seqA->gatherSubseq(&fragA[0], startA+startN, endSize, '+');
	seqB->gatherSubseq(&fragB[0], startB+startN, endSize, '+');
	for (long n2 = 0; n2 < endSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
      }
    }
  }

  return maxNumMis >= misCount;
}
bool ScoredSeqCollectionBwt::AlignmentMaker::makeAlignmentHelperMinus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB,
								      long bStart, long maxNumMis, long exBlockCount,
								      long* examineStarts, long* examineEnds){


  long subBinSize = 25;
  char fragA[subBinSize+1]; 
  char fragB[subBinSize+1]; 

  // gaps are irrelevant
  long misCount = 0;

  // initial block is defined by the last block element
  --exBlockCount;

  long startA = examineStarts[exBlockCount];
  long startB = bStart - aStart + startA;
  long endA = examineEnds[exBlockCount];
  long blockSize = endA - startA;
  long blockSizeMsub = blockSize - subBinSize;
  long startN = 0;
  while (startN < blockSizeMsub and maxNumMis >= misCount){
     seqA->gatherSubseq(&fragA[0], startA+startN, subBinSize, '+');
     seqB->gatherSubseq(&fragB[0], startB+startN, subBinSize, '-');
     for (long n2 = 0; n2 < subBinSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
     startN += subBinSize;
  }
  if (maxNumMis >= misCount){
    long endSize = blockSize % subBinSize;
    seqA->gatherSubseq(&fragA[0], startA+startN, endSize, '+');
    seqB->gatherSubseq(&fragB[0], startB+startN, endSize, '-');
    for (long n2 = 0; n2 < endSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
  }

  while ( (maxNumMis >= misCount) and (exBlockCount > 0) ){
    --exBlockCount;

    if (startA < examineEnds[exBlockCount]){ endA = startA; }
    else { endA = examineEnds[exBlockCount]; }
    startA = examineStarts[exBlockCount];
    blockSize = endA - startA;
    if (blockSize > 0){
      startB = bStart - aStart + startA;
      blockSizeMsub = blockSize - subBinSize;
      startN = 0;
      while (startN < blockSizeMsub and maxNumMis >= misCount){
        seqA->gatherSubseq(&fragA[0], startA+startN, subBinSize, '+');
        seqB->gatherSubseq(&fragB[0], startB+startN, subBinSize, '-');
        for (long n2 = 0; n2 < subBinSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
        startN += subBinSize;
      }
      if (maxNumMis >= misCount){
	long endSize = blockSize % subBinSize;
	seqA->gatherSubseq(&fragA[0], startA+startN, endSize, '+');
	seqB->gatherSubseq(&fragB[0], startB+startN, endSize, '-');
	for (long n2 = 0; n2 < endSize; ++n2){ misCount += scoreNucMisMatch(fragA[n2], fragB[n2]); }
      }
    }
  }

  return maxNumMis >= misCount;
}




//ScoredSeqCollectionBwt::AlMakerNoEdge::AlMakerNoEdge(){}
ScoredSeqCollectionBwt::AlMakerNoEdge::AlMakerNoEdge(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis) :
  _scoreMatrix(scoreMatrix),
  _minScore(minScore),
  _maxFractMis(maxFractMis)
  //_misCountAsm(new AlignmentScoreMatrix(0,1,1,1))
{
  if (_scoreMatrix->_match <= _scoreMatrix->_mismatch){
    throw AssemblyException::ArgError("SSCBwt::AlMakerNoEdge needs a lower mismatch score than match score (i.e. mismatch should be negative)");
  }
}
ScoredSeqCollectionBwt::AlMakerNoEdge::~AlMakerNoEdge(){
  //delete _misCountAsm;
}
//AlignmentScoreMatrix* ScoredSeqCollectionBwt::AlMakerNoEdge::accessMisCountAsm(){ return _misCountAsm; }

Alignment* ScoredSeqCollectionBwt::AlMakerNoEdge::makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense,
								OffsetAndQueryCarrier* offsetCarrier){
  long offset = offsetCarrier->_offset;

  // get the starts and alignment length
  long aStart;
  long bStart;
  if (offset < 0){
    aStart = 0;
    bStart = 0 - offset;
  } else {
    aStart = offset;
    bStart = 0;
  }
  long aEnd;
  if (offset + seqB->size() >= seqA->size()){ aEnd = seqA->size(); }
  else { aEnd = offset + seqB->size(); }
  long alLength = aEnd - aStart;

  // i can maybe take out these tests
  if (aStart != 0 and bStart != 0){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither start is zero");
  }
  if ( aStart + alLength != seqA->size() and bStart + alLength != seqB->size() ){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither end is the end of the sequence");
  }

  // the min score must be beaten, and the fractId must also be satisfied
  long maxNumMisByScore = (_minScore - _scoreMatrix->_match * alLength) / (_scoreMatrix->_mismatch - _scoreMatrix->_match );
  long maxNumMisByFid = long(float(alLength) * _maxFractMis);
  long maxNumMis;
  if (maxNumMisByScore < maxNumMisByFid){ maxNumMis = maxNumMisByScore; }
  else { maxNumMis = maxNumMisByFid; }

  // lowest to highest
  long exBlockCount = 1 + offsetCarrier->_numQueries;
  bool alIsGood;

  // optimization for using stack memory if it won't take up too much space
  if (exBlockCount > 50){
    long* examineStarts = new long[exBlockCount];
    long* examineEnds = new long[exBlockCount];
    examineStarts[0] = aStart;
    long* examineStartsP1 = &examineStarts[1];
    examineEnds[exBlockCount-1] = aStart + alLength;
    for (long n = 0; n < offsetCarrier->_numQueries; ++n){
      examineEnds[n] = offsetCarrier->_queries[n]->start();
      examineStartsP1[n] = offsetCarrier->_queries[n]->end();
    }
    switch (sense){
    case '+': alIsGood = makeAlignmentHelperPlus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
    case '-': alIsGood = makeAlignmentHelperMinus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
    default: throw AssemblyException::ArgError("SSCBwt:AlMakerNoEdge:makeAl sense char is bad");
    }
    delete [] examineStarts;
    delete [] examineEnds;
  } else {
    long examineStarts[exBlockCount];
    long examineEnds[exBlockCount];
    examineStarts[0] = aStart;
    long* examineStartsP1 = &examineStarts[1];
    examineEnds[exBlockCount-1] = aStart + alLength;
    for (long n = 0; n < offsetCarrier->_numQueries; ++n){
      examineEnds[n] = offsetCarrier->_queries[n]->start();
      examineStartsP1[n] = offsetCarrier->_queries[n]->end();
    }
    switch (sense){
    case '+': alIsGood = makeAlignmentHelperPlus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, &examineStarts[0], &examineEnds[0]); break;
    case '-': alIsGood = makeAlignmentHelperMinus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, &examineStarts[0], &examineEnds[0]); break;
    default: throw AssemblyException::ArgError("SSCBwt:AlMakerNoEdge:makeAl sense char is bad");
    }
  }

  if (alIsGood){ return new AlignmentUngapped(seqA, seqB, sense, offset); }
  else { return new AlignmentNull(seqA, seqB, sense); }
}

long ScoredSeqCollectionBwt::AlMakerNoEdge::getScore(Alignment* al){
  return al->score(_scoreMatrix, false);
}
long ScoredSeqCollectionBwt::AlMakerNoEdge::getMinScore(){ return _minScore; }
void ScoredSeqCollectionBwt::AlMakerNoEdge::setMinScore(long minScore){ _minScore = minScore; }
AlignmentScoreMatrix* ScoredSeqCollectionBwt::AlMakerNoEdge::getScoreMatrix(){ return _scoreMatrix; }



//ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::AlMakerPenalizeEdgeA(){}
ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::AlMakerPenalizeEdgeA(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis) :
  _scoreMatrix(scoreMatrix),
  _minScore(minScore),
  _maxFractMis(maxFractMis)
  //_misCountAsm(new AlignmentScoreMatrix(0,1,1,1))
{
  if (_scoreMatrix->_match <= _scoreMatrix->_mismatch){
    throw AssemblyException::ArgError("SSCBwt::AlMakerNoEdge needs a lower mismatch score than match score (i.e. mismatch should be negative)");
  }
}
ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::~AlMakerPenalizeEdgeA(){
  //delete _misCountAsm;
}
//AlignmentScoreMatrix* ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::accessMisCountAsm(){ return _misCountAsm; }


Alignment* ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense,
								       OffsetAndQueryCarrier* offsetCarrier){
  long offset = offsetCarrier->_offset;

  // get the starts and alignment length
  long aStart;
  long bStart;
  if (offset < 0){
    aStart = 0;
    bStart = 0 - offset;
  } else {
    aStart = offset;
    bStart = 0;
  }
  long aEnd;
  if (offset + seqB->size() >= seqA->size()){ aEnd = seqA->size(); }
  else { aEnd = offset + seqB->size(); }
  long alLength = aEnd - aStart;

  // i can maybe take out these tests
  if (aStart != 0 and bStart != 0){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither start is zero");
  }
  if ( aStart + alLength != seqA->size() and bStart + alLength != seqB->size() ){
    throw AssemblyException::ArgError("in SSCBwt::AlMaker, neither end is the end of the sequence");
  }

  // the min score must be beaten, and the fractId must also be satisfied
  long maxNumMisByScore = (_minScore - _scoreMatrix->_match * seqA->size()) / (_scoreMatrix->_mismatch - _scoreMatrix->_match );
  long maxNumMisByFid = long(float(seqA->size()) * _maxFractMis);
  long maxNumMis;
  if (maxNumMisByScore < maxNumMisByFid){ maxNumMis = maxNumMisByScore; }
  else { maxNumMis = maxNumMisByFid; }
  // adjust for the edges
  maxNumMis -= seqA->size() - alLength;

  // lowest to highest
  long exBlockCount = 1 + offsetCarrier->_numQueries;
  bool alIsGood;

  // optimization for using stack memory if it won't take up too much space
  if (exBlockCount > 50){
    long* examineStarts = new long[exBlockCount];
    long* examineEnds = new long[exBlockCount];
    examineStarts[0] = aStart;
    long* examineStartsP1 = &examineStarts[1];
    examineEnds[exBlockCount-1] = aStart + alLength;
    for (long n = 0; n < offsetCarrier->_numQueries; ++n){
      examineEnds[n] = offsetCarrier->_queries[n]->start();
      examineStartsP1[n] = offsetCarrier->_queries[n]->end();
    }
    switch (sense){
    case '+': alIsGood = makeAlignmentHelperPlus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
    case '-': alIsGood = makeAlignmentHelperMinus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, examineStarts, examineEnds); break;
    default: throw AssemblyException::ArgError("SSCBwt:AlMakerNoEdge:makeAl sense char is bad");
    }
    delete [] examineStarts;
    delete [] examineEnds;
  } else {
    long examineStarts[exBlockCount];
    long examineEnds[exBlockCount];
    examineStarts[0] = aStart;
    long* examineStartsP1 = &examineStarts[1];
    examineEnds[exBlockCount-1] = aStart + alLength;
    for (long n = 0; n < offsetCarrier->_numQueries; ++n){
      examineEnds[n] = offsetCarrier->_queries[n]->start();
      examineStartsP1[n] = offsetCarrier->_queries[n]->end();
    }
    switch (sense){
    case '+': alIsGood = makeAlignmentHelperPlus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, &examineStarts[0], &examineEnds[0]); break;
    case '-': alIsGood = makeAlignmentHelperMinus(seqA, aStart, seqB, bStart, maxNumMis, exBlockCount, &examineStarts[0], &examineEnds[0]); break;
    default: throw AssemblyException::ArgError("SSCBwt:AlMakerNoEdge:makeAl sense char is bad");
    }
  }

  if (alIsGood){ return new AlignmentUngapped(seqA, seqB, sense, offset); }
  else { return new AlignmentNull(seqA, seqB, sense); }
}
long ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::getScore(Alignment* al){
  return al->scoreOverhangA(_scoreMatrix);
}
long ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::getMinScore(){ return _minScore; }
void ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::setMinScore(long minScore){ _minScore = minScore; }
AlignmentScoreMatrix* ScoredSeqCollectionBwt::AlMakerPenalizeEdgeA::getScoreMatrix(){ return _scoreMatrix; }



ScoredSeqCollectionBwt::OffsetTracker::~OffsetTracker(){}



ScoredSeqCollectionBwt::OffsetTrackerMutable::OffsetTrackerMutable(ScoredSeq* target) : _targetContig(target){
  _maxOffsetCount = 0;
  // this value is meaningless unless the above values are >0
  _maxCountOffset = 0;
}
ScoredSeqCollectionBwt::OffsetTrackerMutable::~OffsetTrackerMutable(){
  for (map<long,OffsetFields*>::iterator it = _offsetToFields.begin(); it != _offsetToFields.end(); ++it){
    delete it->second;
  }
}

ScoredSeqCollectionBwt::OffsetTrackerMutable::OffsetFields::OffsetFields(float count, WordQuery* maxWord) : 
  _count(count),
  _maxWord(maxWord->copy()){}
ScoredSeqCollectionBwt::OffsetTrackerMutable::OffsetFields::~OffsetFields(){
  delete _maxWord;
  for (vector<WordQuery*>::iterator wqIt = _otherWords.begin(); wqIt != _otherWords.end(); ++wqIt){ delete *wqIt; }
}
void ScoredSeqCollectionBwt::OffsetTrackerMutable::OffsetFields::addOffset(float localAddition, WordQuery* seedingQuery){
  _count += localAddition;
  if ( seedingQuery->notContinuous(_maxWord) ){
    // figure out which word is more important (longer)
    if (seedingQuery->size() > _maxWord->size()){
      _otherWords.push_back(_maxWord);
      _maxWord = seedingQuery->copy();
    } else { _otherWords.push_back(seedingQuery->copy()); }
  } else {
    _maxWord->merge(seedingQuery->start(), seedingQuery->end());
  }
}

void ScoredSeqCollectionBwt::OffsetTrackerMutable::addOffset(long offset, WordQuery* seedingQuery, long normFactor){
  float localCount;
  if (normFactor < 1){ throw AssemblyException::ArgError("SSCBwt::OffsetTrackerMutable::addOffset, normFactor cannot be less than 1"); }
  float localAddition = float(1) / float(normFactor);
  // now find or create the offset entry
  map<long,OffsetFields*>::iterator offsetEntry = _offsetToFields.find(offset);
  if (offsetEntry == _offsetToFields.end()){
    localCount = localAddition;
    _offsetToFields.insert( pair<long,OffsetFields*>(offset, new OffsetFields(localAddition, seedingQuery)) );
  } else {
    offsetEntry->second->addOffset(localAddition, seedingQuery);
    localCount = offsetEntry->second->_count;
  }
  if (localCount > _maxOffsetCount){
    _maxOffsetCount = localCount;
    _maxCountOffset = offset;
  }
}
long ScoredSeqCollectionBwt::OffsetTrackerMutable::maxOffsetCount(){ return _maxOffsetCount; }
long ScoredSeqCollectionBwt::OffsetTrackerMutable::maxCountOffset(){ return _maxCountOffset; }
ScoredSeqCollectionBwt::OffsetAndQueryCarrier** ScoredSeqCollectionBwt::OffsetTrackerMutable::getOffsets(bool includeQueries){

  // figure out how many offsets there are with each count
  map<float,long> countToNum;
  for (map<long,OffsetFields*>::iterator it = _offsetToFields.begin(); it != _offsetToFields.end(); ++it){
    map<float,long>::iterator foundCtm = countToNum.find( it->second->_count );
    if (foundCtm == countToNum.end()){ countToNum.insert(pair<float,long>(it->second->_count,1)); }
    else { foundCtm->second++; }
  }
  // figure out the index in the sorted array at which each offset count will begin
  map<float,long> countToIndex;
  map<float,long>::iterator ctiIt = countToIndex.begin();
  long index = _offsetToFields.size();
  for (map<float,long>::iterator it = countToNum.begin(); it != countToNum.end(); ++it){
    index -= it->second;
    // i can skip to the end because the keys were already sorted in the other map
    ctiIt = countToIndex.insert(ctiIt, pair<float,long>(it->first,index));
  }

  long numOutput = _offsetToFields.size();
  OffsetAndQueryCarrier** output = new OffsetAndQueryCarrier*[ numOutput + 1];
  long outIndex = 0;
  map<long,OffsetFields*>::iterator fieldIt = _offsetToFields.begin();
  while (outIndex < numOutput){
    // this is a test - i can comment it out later

    OffsetAndQueryCarrier* newOut;
    if (includeQueries){
      // determine if there are any other queries to deal with
      if (fieldIt->second->_otherWords.empty()){
	newOut = new OffsetAndQueryCarrier(fieldIt->first, 1);
	newOut->_queries[0] = fieldIt->second->_maxWord->copy();
      } else {
	// sort the queries and make them non-redundant while creating the output
	vector<WordQuery*> sortQueries;
	sortQueries.insert(sortQueries.end(), fieldIt->second->_otherWords.begin(), fieldIt->second->_otherWords.end());
	sortQueries.push_back( fieldIt->second->_maxWord );
	SortWqByStart wqSorter;
	sort( sortQueries.begin(), sortQueries.end(), wqSorter );

	WordQuery** wqArray = new WordQuery*[ sortQueries.size() ];
	vector<WordQuery*>::iterator sortIt = sortQueries.begin();
	long wqIndex = 0;
	wqArray[wqIndex] = (*sortIt)->copy();
	++sortIt;
	while (sortIt != sortQueries.end()){
	  if (wqArray[wqIndex]->notContinuous( *sortIt )){
	    ++wqIndex;
	    wqArray[wqIndex] = (*sortIt)->copy();
	  } else {
	    wqArray[wqIndex]->merge( (*sortIt)->start(), (*sortIt)->end()  );
	  }
	  ++sortIt;
	}
	// advance wqIndex one beyond the end of the array's last active element
	++wqIndex;
	newOut = new OffsetAndQueryCarrier(fieldIt->first, wqIndex);
	for (long n = 0; n < wqIndex; ++n){ newOut->_queries[n] = wqArray[n]; }
	delete [] wqArray;
      }
    } else {
      newOut = new OffsetAndQueryCarrier(fieldIt->first, 0);
    }
    map<float,long>::iterator foundCtm = countToIndex.find( fieldIt->second->_count );
    output[ foundCtm->second ] = newOut;
    foundCtm->second += 1;
    ++outIndex;
    ++fieldIt;
  }
  return output;
}
bool ScoredSeqCollectionBwt::OffsetTrackerMutable::SortWqByStart::operator() (WordQuery* wqA, WordQuery* wqB){
  return wqA->start() < wqB->start();
}
long ScoredSeqCollectionBwt::OffsetTrackerMutable::numOffsets(){ return long(_offsetToFields.size()); }
ScoredSeq* ScoredSeqCollectionBwt::OffsetTrackerMutable::targetContig(){ return _targetContig; }
bool ScoredSeqCollectionBwt::OffsetTrackerMutable::queriesIncluded(){ return true; }

ScoredSeqCollectionBwt::OffsetTrackerImmutable* ScoredSeqCollectionBwt::OffsetTrackerMutable::immutableCopy(bool queryIncluded){
  return new OffsetTrackerImmutable(_targetContig, long(_offsetToFields.size()), getOffsets(queryIncluded),
				    _maxCountOffset, _maxOffsetCount, queryIncluded);
}





ScoredSeqCollectionBwt::OffsetTrackerImmutable::OffsetTrackerImmutable(ScoredSeq* targetContig,
								       long numOffsets, OffsetAndQueryCarrier** offsets,
								       long maxCountOffset, long maxOffsetCount, bool qIncluded) :
  _targetContig(targetContig),
  _numOffsets(numOffsets),
  _offsets(offsets),
  _maxCountOffset(maxCountOffset),
  _maxOffsetCount(maxOffsetCount),
  _queriesIncluded(qIncluded)
{}
ScoredSeqCollectionBwt::OffsetTrackerImmutable::~OffsetTrackerImmutable(){
  for (long n = 0; n < _numOffsets; ++n){ _offsets[n]->deleteWithQueries(); }
  delete [] _offsets;
}
long ScoredSeqCollectionBwt::OffsetTrackerImmutable::numOffsets(){ return _numOffsets; }
ScoredSeqCollectionBwt::OffsetAndQueryCarrier** ScoredSeqCollectionBwt::OffsetTrackerImmutable::getOffsets(bool includeQueries){
  if (includeQueries and (! _queriesIncluded)){
    throw AssemblyException::ArgError("SSCBwt::OTImmutable::getOffsets cannot include queries if they weren't included at construction");
  }
  OffsetAndQueryCarrier** newArray = new OffsetAndQueryCarrier*[_numOffsets+1];
  for (long n = 0; n < _numOffsets; ++n){ newArray[n] = _offsets[n]->copy(); }
  return newArray;
}
ScoredSeq* ScoredSeqCollectionBwt::OffsetTrackerImmutable::targetContig(){ return _targetContig; }
bool ScoredSeqCollectionBwt::OffsetTrackerImmutable::queriesIncluded(){ return _queriesIncluded; }
long ScoredSeqCollectionBwt::OffsetTrackerImmutable::maxOffsetCount(){ return _maxOffsetCount; }
long ScoredSeqCollectionBwt::OffsetTrackerImmutable::maxCountOffset(){ return _maxCountOffset; }
ScoredSeqCollectionBwt::OffsetTrackerImmutable* ScoredSeqCollectionBwt::OffsetTrackerImmutable::immutableCopy(bool queryIncluded){
  // the getOffsets method that copies the offsets will also verify the legit-ness of queryIncluded
  return new OffsetTrackerImmutable(_targetContig, _numOffsets, getOffsets(queryIncluded),
				    _maxCountOffset, _maxOffsetCount, queryIncluded);
}





ScoredSeqCollectionBwt::OffsetAndQueryCarrier::OffsetAndQueryCarrier(long offset, long numQueries) :
  _offset(offset), _numQueries(numQueries)
{
  _queries = new WordQuery*[ numQueries+1 ];
}
ScoredSeqCollectionBwt::OffsetAndQueryCarrier::~OffsetAndQueryCarrier(){
  delete [] _queries;
}
void ScoredSeqCollectionBwt::OffsetAndQueryCarrier::deleteWithQueries(){
  for (long n = 0; n < _numQueries; ++n){ delete _queries[n]; }
  delete this;
}
ScoredSeqCollectionBwt::OffsetAndQueryCarrier* ScoredSeqCollectionBwt::OffsetAndQueryCarrier::copy(){
  OffsetAndQueryCarrier* carrierCopy = new OffsetAndQueryCarrier(_offset, _numQueries);
  for (long n = 0; n < _numQueries; ++n){ carrierCopy->_queries[n] = _queries[n]->copy(); }
  //for (long n = 0; n < _numQueries; ++n){ carrierCopy->_queries[n] = _queries[n]; }
  return carrierCopy;
}





// VERSION 1: BASIC


ScoredSeqCollectionBwt::WordQuery::~WordQuery(){}


ScoredSeqCollectionBwt::WordQuerySimple::WordQuerySimple(long start, long size) :
  _start(start), _size(size), _end(start+size), _index(0)
{}
ScoredSeqCollectionBwt::WordQuerySimple::WordQuerySimple(long start, long end, long size, long index) :
  _start(start), _size(size), _end(end), _index(index)
{}
ScoredSeqCollectionBwt::WordQuerySimple::~WordQuerySimple(){}
void ScoredSeqCollectionBwt::WordQuerySimple::merge(long start, long end){
  if (start < _start){ _start = start; }
  if (end > _end){ _end = end; }
  _size = _end - _start;
}
ScoredSeqCollectionBwt::WordQuerySimple* ScoredSeqCollectionBwt::WordQuerySimple::copy(){
  return new WordQuerySimple(_start, _end, _size, _index);
}

long ScoredSeqCollectionBwt::WordQuerySimple::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQuerySimple::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQuerySimple::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQuerySimple::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQuerySimple::acceptableOverlap(long targetStart, ScoredSeq* target){ return true; }
long ScoredSeqCollectionBwt::WordQuerySimple::getIndex(){ return _index; }
void ScoredSeqCollectionBwt::WordQuerySimple::setIndex(long index){ _index = index; }


ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::WordQueryNoEdgeLimit(long start, long size, char sense) :
  _start(start), _size(size), _sense(sense), _end(start+size), _index(0)
{}
ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::WordQueryNoEdgeLimit(long start, long end, long size, char sense, long index) :
  _start(start), _size(size), _sense(sense), _end(end), _index(index)
{}
ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::~WordQueryNoEdgeLimit(){}
ScoredSeqCollectionBwt::WordQueryNoEdgeLimit* ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::copy(){
  return new WordQueryNoEdgeLimit(_start, _end, _size, _sense, _index);
}
void ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::merge(long start, long end){
  if (start < _start){ _start = start; }
  if (end > _end){ _end = end; }
  _size = _end - _start;
}
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::acceptableOverlap(long targetStart, ScoredSeq* target){ return true; }
long ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::getIndex(){ return _index; }
void ScoredSeqCollectionBwt::WordQueryNoEdgeLimit::setIndex(long index){ _index = index; }


ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::WordQueryFiveEdgeLimit(long start, long size, char sense, long edgeLimit) :
  _start(start), _size(size), _sense(sense), _edgeLimitX2(edgeLimit * 2), _end(start+size), _index(0) {}
ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::WordQueryFiveEdgeLimit(long start, long end, long size, char sense, long edgeLimitX2, long index) :
  _start(start), _size(size), _sense(sense), _edgeLimitX2(edgeLimitX2), _end(end), _index(index) {}
ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::~WordQueryFiveEdgeLimit(){}
ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit* ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::copy(){
  return new WordQueryFiveEdgeLimit(_start, _end, _size, _sense, _edgeLimitX2, _index);
}
void ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::merge(long start, long end){
  if (start < _start){ _start = start; }
  if (end > _end){ _end = end; }
  _size = _end - _start;
}
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::acceptableOverlap(long targetStart, ScoredSeq* target){
  // distance to edge
  return target->size() - targetStart <= _edgeLimitX2;
}
long ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::getIndex(){ return _index; }
void ScoredSeqCollectionBwt::WordQueryFiveEdgeLimit::setIndex(long index){ _index = index; }

ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::WordQueryThreeEdgeLimit(long start, long size, char sense, long edgeLimit) :
  _start(start), _size(size), _sense(sense), _edgeLimitX2(edgeLimit * 2), _end(start+size), _index(0) {}
ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::WordQueryThreeEdgeLimit(long start, long end, long size, char sense, long edgeLimitX2, long index) :
  _start(start), _size(size), _sense(sense), _edgeLimitX2(edgeLimitX2), _end(end), _index(index) {}
ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::~WordQueryThreeEdgeLimit(){}
ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit* ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::copy(){
  return new WordQueryThreeEdgeLimit(_start, _end, _size, _sense, _edgeLimitX2, _index);
}
void ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::merge(long start, long end){
  if (start < _start){ _start = start; }
  if (end > _end){ _end = end; }
  _size = _end - _start;
}
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::acceptableOverlap(long targetStart, ScoredSeq* target){
  return targetStart + _size <= _edgeLimitX2;
}
long ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::getIndex(){ return _index; }
void ScoredSeqCollectionBwt::WordQueryThreeEdgeLimit::setIndex(long index){ _index = index; }


ScoredSeqCollectionBwt::WordQueryFullTarget::WordQueryFullTarget(long sourceSize, long start, long size, char sense, float fractId, long maxTargetSize) :
  _sourceSize(sourceSize), _start(start), _end(start+size), _size(size), _sense(sense),
  _fractId(fractId), _maxTargetSize(maxTargetSize), _index(0)
{}
ScoredSeqCollectionBwt::WordQueryFullTarget::WordQueryFullTarget(long sourceSize, long start, long end, long size, char sense, float fractId, long maxTargetSize, long index) :
  _sourceSize(sourceSize), _start(start), _end(end), _size(size), _sense(sense),
  _fractId(fractId), _maxTargetSize(maxTargetSize), _index(index)
{}
ScoredSeqCollectionBwt::WordQueryFullTarget::~WordQueryFullTarget(){}
ScoredSeqCollectionBwt::WordQueryFullTarget* ScoredSeqCollectionBwt::WordQueryFullTarget::copy(){
  return new WordQueryFullTarget(_sourceSize, _start, _end, _size, _sense, _fractId, _maxTargetSize, _index);
}
void ScoredSeqCollectionBwt::WordQueryFullTarget::merge(long start, long end){
  if (start < _start){ _start = start; }
  if (end > _end){ _end = end; }
  _size = _end - _start;
}
long ScoredSeqCollectionBwt::WordQueryFullTarget::start(){ return _start; }
long ScoredSeqCollectionBwt::WordQueryFullTarget::end(){ return _end; }
long ScoredSeqCollectionBwt::WordQueryFullTarget::size(){ return _size; }
bool ScoredSeqCollectionBwt::WordQueryFullTarget::notContinuous(WordQuery* otherWq){
  return ( _start > otherWq->end() or _end < otherWq->start() );
}
bool ScoredSeqCollectionBwt::WordQueryFullTarget::acceptableOverlap(long targetStart, ScoredSeq* target){
  //the target must be smaller than the source sequence, excpting for the fractId compensation
  if (float(_sourceSize) < float(target->size()) * _fractId){ return false; }
  else { 
    long targetTrim = long((float(1) - _fractId) * float(target->size()) );
    return (targetStart - targetTrim <= _start) and 
      (_start + target->size() - targetStart <= _sourceSize + targetTrim);
  }
}
long ScoredSeqCollectionBwt::WordQueryFullTarget::getIndex(){ return _index; }
void ScoredSeqCollectionBwt::WordQueryFullTarget::setIndex(long index){ _index = index; }





#endif



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

#ifndef ECOFILTER_CPP
#define ECOFILTER_CPP

#include "EcoFilter.h"
#include "AssemblyException.h"
#include "Alignment.h"
#include "ScoredSeqCollectionBwt.h"
#include "ScoredSeqNormalized.h"

EcoFilter::~EcoFilter(){}


EcoFilterNull::EcoFilterNull(){}
EcoFilterNull::~EcoFilterNull(){}
EcoFilterNull* EcoFilterNull::copy(){ return new EcoFilterNull(); }
void EcoFilterNull::filterSeqs(set<ScoredSeq*>* fullSet,
			       set<ScoredSeq*>* retainedSet,
			       set<ScoredSeq*>* removedSet,
			       int cycleNum, FilterThreadedness threadedness){
  retainedSet->insert(fullSet->begin(), fullSet->end());
}


EcoFilterMinCount::EcoFilterMinCount(){}
EcoFilterMinCount::EcoFilterMinCount(float minCount) : _minCount(minCount){}
EcoFilterMinCount::~EcoFilterMinCount(){}
EcoFilterMinCount* EcoFilterMinCount::copy(){ return new EcoFilterMinCount(_minCount); }
void EcoFilterMinCount::filterSeqs(set<ScoredSeq*>* fullSet,
				   set<ScoredSeq*>* retainedSet,
				   set<ScoredSeq*>* removedSet,
				   int cycleNum, FilterThreadedness threadedness){
  long numSeqs = fullSet->size();
  ScoredSeq** seqArray = new ScoredSeq*[ numSeqs + 1 ];
  long arrayIndex = 0;
  for (set<ScoredSeq*>::iterator seqIt = fullSet->begin(); seqIt != fullSet->end(); ++seqIt){
    seqArray[arrayIndex] = *seqIt;
    ++arrayIndex;
  }
  bool* okArray = new bool[ numSeqs + 1 ];
  if (threadedness == threaded){
    #pragma omp parallel for schedule(dynamic)
    for (long n = 0; n < numSeqs; ++n){ okArray[n] = isOk( seqArray[n] ); }
  } else {
    for (long n = 0; n < numSeqs; ++n){ okArray[n] = isOk( seqArray[n] ); }
  }

  for (long n = 0; n < numSeqs; ++n){
    if (okArray[n]){ retainedSet->insert(seqArray[n]); }
    else { removedSet->insert(seqArray[n]); }
  }
  delete [] seqArray;
  delete [] okArray;
}
bool EcoFilterMinCount::isOk(ScoredSeq* seq){
  long seqSizeM1 = seq->size() - 1;
  long n = 0;
  bool isNotOk = true;
  while (n < seqSizeM1 and isNotOk){
    if (seq->scoreAtPosPlus(n) >= _minCount or seq->linkAfterPosPlus(n) >= _minCount){ isNotOk = false; }
    ++n;
  }
  if (isNotOk and seqSizeM1 >= 0){ isNotOk = seq->scoreAtPosPlus(seqSizeM1) < _minCount; }
  return (! isNotOk);
}
EcoFilterMinCount* EcoFilterMinCount::defaultFilter(){
  return new EcoFilterMinCount(1.6);
}




EcoFilterMinLength::EcoFilterMinLength(){
  _maxCycleNum = 0;
  _minLenByCycle = new long[1];
  _minLenByCycle[0] = 0;
  _isTransition =  new bool[1];
  _isTransition[0] = false;
}
EcoFilterMinLength::EcoFilterMinLength(long minLength, int skipCycles){
  _maxCycleNum = skipCycles;
  _minLenByCycle = new long[_maxCycleNum+1];
  _isTransition =  new bool[_maxCycleNum+1];
  for (long n = 0; n < _maxCycleNum; ++n){
    _minLenByCycle[n] = 0;
    _isTransition[n] = false;
  }
  _minLenByCycle[_maxCycleNum] = minLength;
  _isTransition[_maxCycleNum] = true;
}
EcoFilterMinLength::EcoFilterMinLength(EcoFilterMinLength* efml){
  _maxCycleNum = efml->_maxCycleNum;
  _minLenByCycle = new long[_maxCycleNum+1];
  _isTransition =  new bool[_maxCycleNum+1];
  for (long n = 0; n <= _maxCycleNum; ++n){
    _minLenByCycle[n] = efml->_minLenByCycle[n];
    _isTransition[n] = efml->_isTransition[n];
  }
}
EcoFilterMinLength::~EcoFilterMinLength(){
  delete [] _minLenByCycle;
  delete [] _isTransition;
}
void EcoFilterMinLength::addMinLength(long minLength, int skipCycles){
  if (skipCycles > _maxCycleNum){
    long* newMinLenByCycle = new long[skipCycles+1];
    bool* newIsTransition = new bool[skipCycles+1];
    // copy the old values
    for (long n = 0; n <= _maxCycleNum; ++n){
      newMinLenByCycle[n] = _minLenByCycle[n];
      newIsTransition[n] = _isTransition[n];
    }
    // fill in the new values
    for (long n = _maxCycleNum + 1; n < skipCycles; ++n){
      newMinLenByCycle[n] = _minLenByCycle[_maxCycleNum];
      newIsTransition[n] = false;
    }
    newMinLenByCycle[skipCycles] = minLength;
    newIsTransition[skipCycles] = true;
    // replace the old values with the new ones
    _maxCycleNum = skipCycles;
    delete [] _minLenByCycle;
    _minLenByCycle = newMinLenByCycle;
    delete [] _isTransition;
    _isTransition = newIsTransition;
  } else {
    _minLenByCycle[skipCycles] = minLength;
    _isTransition[skipCycles] = true;
    // i don't need to do this if the values are equal
    if (skipCycles < _maxCycleNum){
      long n = skipCycles+1;
      while (! _isTransition[n]){
	_minLenByCycle[n] = minLength;
	++n;
      }
    }
  }
}
EcoFilterMinLength* EcoFilterMinLength::copy(){
  return new EcoFilterMinLength(this);
}
void EcoFilterMinLength::filterSeqs(set<ScoredSeq*>* fullSet,
				    set<ScoredSeq*>* retainedSet,
				    set<ScoredSeq*>* removedSet,
				    int cycleNum, FilterThreadedness threadedness){
  // I am not going to bother threading this
  long minLength;
  if (cycleNum > _maxCycleNum){ minLength = _minLenByCycle[_maxCycleNum]; }
  else { minLength = _minLenByCycle[cycleNum]; }
  // there cannot be zero-length contigs, so anything less than 2 (i.e. 1 or 0) is a non-filter
  if (minLength < 2){
    retainedSet->insert(fullSet->begin(), fullSet->end());
  } else {
    for (set<ScoredSeq*>::iterator seqIt = fullSet->begin(); seqIt != fullSet->end(); ++seqIt){
      if ( (*seqIt)->size() < minLength ){ removedSet->insert( *seqIt ); }
      else { retainedSet->insert( *seqIt ); }
    }
  }
}





EcoFilterInitialContigMatch::EcoFilterInitialContigMatch(){}
EcoFilterInitialContigMatch::EcoFilterInitialContigMatch(float minFractId, int cyclesToSkip, bool fullFile) : 
  _minFractId(minFractId),
  _cyclesToSkip(cyclesToSkip),
  _numFilterCycles(1),
  _numSkipCycles(0),
  _fullFile(fullFile)
{
  _asMatrix = AlignmentScoreMatrix::getDefault();
}
EcoFilterInitialContigMatch::EcoFilterInitialContigMatch(float minFractId, int cyclesToSkip,
							 int numFilterCycles, int numSkipCycles, bool fullFile) :
  _minFractId(minFractId),
  _cyclesToSkip(cyclesToSkip),
  _numFilterCycles(numFilterCycles),
  _numSkipCycles(numSkipCycles),
  _fullFile(fullFile)
{
  _asMatrix = AlignmentScoreMatrix::getDefault();
}
EcoFilterInitialContigMatch::~EcoFilterInitialContigMatch(){
  delete _asMatrix;
}
EcoFilterInitialContigMatch* EcoFilterInitialContigMatch::copy(){ 
  EcoFilterInitialContigMatch* copy = new EcoFilterInitialContigMatch(_minFractId,_cyclesToSkip,
								      _numFilterCycles,_numSkipCycles,_fullFile);
  copy->addFiles(&_initialContigFiles);
  copy->setAlignmentScoreMatrix(_asMatrix);
  return copy;
}

void EcoFilterInitialContigMatch::addFiles(set<ParameterizedInitialFile*>* initialContigFiles){
  _initialContigFiles.insert(initialContigFiles->begin(), initialContigFiles->end());
}
void EcoFilterInitialContigMatch::setAlignmentScoreMatrix(AlignmentScoreMatrix* alSM){
  delete _asMatrix;
  _asMatrix = new AlignmentScoreMatrix(alSM);
}

void EcoFilterInitialContigMatch::filterSeqs(set<ScoredSeq*>* fullSet,
					     set<ScoredSeq*>* retainedSet,
					     set<ScoredSeq*>* removedSet,
					     int cycleNum, FilterThreadedness threadedness){
  // too early in the assembly, so nothing happens
  if (cycleNum < _cyclesToSkip or
      (cycleNum - _cyclesToSkip) % (_numFilterCycles + _numSkipCycles) >= _numFilterCycles ){
    retainedSet->insert(fullSet->begin(), fullSet->end());

  } else {
    // get the input contigs from all of the prior cycles AND the current cycle
    set<ScoredSeq*> anchorContigs;
    long numFiles = _initialContigFiles.size();
    int* numCycles = new int[ numFiles+1 ];
    ParameterizedInitialFile** initFiles = new ParameterizedInitialFile*[ numFiles+1 ];
    long fN = 0;
    for (set<ParameterizedInitialFile*>::iterator it = _initialContigFiles.begin(); it != _initialContigFiles.end(); ++it){
      initFiles[fN] = *it;
      if (_fullFile){ numCycles[fN] = (*it)->totalCycles(); }
      else { numCycles[fN] = cycleNum; }
      ++fN;
    }

    if (threadedness == threaded){
      #pragma omp parallel for schedule(dynamic)
      for (long fN = 0; fN < numFiles; ++fN){
	set<ScoredSeq*> tempInput;
	for (int n = 0; n <= numCycles[fN]; ++n){ initFiles[fN]->getContigs(&tempInput, n); }
        #pragma omp critical (EFICM_fS_getSeqs)
	{ anchorContigs.insert(tempInput.begin(), tempInput.end()); }
      }
    } else {
      for (long fN = 0; fN < numFiles; ++fN){
	for (int n = 0; n <= numCycles[fN]; ++n){ initFiles[fN]->getContigs(&anchorContigs, n); }
      }
    }
    delete [] numCycles;
    delete [] initFiles;

    // this is the same strategy that I use for AJGraph targeting - it is fastest to 
    // generate two collections and use each set of sequences as queries, looking only for
    // matches to bigger sequences
    // THAT STILL NEEDS TO BE DONE!!!!!!!!!!

    DyProAlignerFactory* dpaf = new DyProAlignerFactory(_minFractId, 0, _asMatrix);
    ScoredSeqCollectionBwt* anchorSeqCollection = new ScoredSeqCollectionBwt(&anchorContigs, dpaf);

    // set up for threaededness (make an array)
    ScoredSeq** fullSetArray = new ScoredSeq*[ fullSet->size() + 1 ];
    long* fullSetLengths = new long[ fullSet->size() + 1 ];
    long fullIndex = 0;
    for (set<ScoredSeq*>::iterator it = fullSet->begin(); it != fullSet->end(); ++it){
      fullSetArray[ fullIndex ] = *it;
      fullSetLengths[fullIndex] = (*it)->size();
      ++fullIndex;
    }

    ScoredSeqCollectionBwt::AlignmentThreadedness mapThreads;
    if (threadedness == threaded){ mapThreads = ScoredSeqCollectionBwt::threaded; }
    else { mapThreads = ScoredSeqCollectionBwt::notThreaded; }

    // PART 1: CONTIGS ARE QUERIES

    // search using the contigs as queries - the array tracks the results
    // this time around, the contigs have to be smaller/equal to their matches
    bool* boolArray = anchorSeqCollection->hasMatch(fullIndex, fullSetArray, fullSetLengths, mapThreads,
						    ScoredSeqCollectionBwt::softMinOvl);

    // PART 2: longer searches with queries that didn't already find matches


    delete anchorSeqCollection;

    // sort found vs unfound
    set<ScoredSeq*> foundContigs1;
    set<ScoredSeq*> unfoundContigs1;
    for (long n = 0; n < fullIndex; ++n){
      if ( boolArray[n] ){ foundContigs1.insert(foundContigs1.end(), fullSetArray[n]); }
      else{ unfoundContigs1.insert(unfoundContigs1.end(), fullSetArray[n]); }
    }
    delete [] fullSetArray;
    delete [] fullSetLengths;
    delete [] boolArray;


    // PART 2: CONTIGS ARE DATABASE

    // I only need to proceed if there are contigs without matches
    ///*
    if (unfoundContigs1.size() > 0){

      ScoredSeqCollectionBwt* unfoundCollect1 = new ScoredSeqCollectionBwt(&unfoundContigs1, dpaf);
      set<ScoredSeq*> foundContigs2;
      set<ScoredSeq*> unfoundContigs2;


      set<ScoredSeq*>::iterator anchIt = anchorContigs.begin();

      long anchN = 0;

      vector<Alignment*> matchCollector;
      MatchSeqTestNotInSet* contigFoundTest = new MatchSeqTestNotInSet();

      while(foundContigs1.size() + foundContigs2.size() < fullSet->size() and anchIt != anchorContigs.end()){
	// find the matches
	ScoredSeq* anchSeq = *anchIt;
	++anchIt;

	unfoundCollect1->getMatches(&matchCollector, anchSeq, contigFoundTest, anchSeq->size(),
				    ScoredSeqCollectionBwt::softMinOvl);

	// record the found matches and block them from being re-examined
	for (vector<Alignment*>::iterator alIt = matchCollector.begin(); alIt != matchCollector.end(); ++alIt){
	  // the found contig is seqB - it need not be found again
	  ScoredSeq* foundSeq = (*alIt)->seqB();
	  foundContigs2.insert(foundSeq);
	  contigFoundTest->addToSet(foundSeq);
	  delete *alIt;
	}

	matchCollector.clear();
      }

      delete unfoundCollect1;
      delete contigFoundTest;

      foundContigs1.insert(foundContigs2.begin(), foundContigs2.end());

    }
    //*/
    // now sort out the unfound contigs by running the full & found iterators in parallel
    set<ScoredSeq*>::iterator foundIt = foundContigs1.begin();
    for (set<ScoredSeq*>::iterator it = fullSet->begin(); it != fullSet->end(); ++it){
      if (foundIt==foundContigs1.end() or (*foundIt != *it)){ removedSet->insert(removedSet->end(), *it); }
      else { ++foundIt; }
    }
    retainedSet->insert(foundContigs1.begin(), foundContigs1.end());

    delete dpaf;
    for (set<ScoredSeq*>::iterator it = anchorContigs.begin(); it != anchorContigs.end(); ++it){ (*it)->deepDelete(); }
  }
}




EcoFilterMulti::EcoFilterMulti(){}
EcoFilterMulti::EcoFilterMulti(EcoFilter** testArray, int numTests) : _numTests(numTests){
  _testArray = new EcoFilter*[ _numTests ];
  for (int n = 0; n < _numTests; ++n){ _testArray[n] = testArray[n]->copy(); }
}
EcoFilterMulti::~EcoFilterMulti(){
  for (int testNum = 0; testNum < _numTests; ++testNum){ delete _testArray[testNum]; }
  delete [] _testArray;
}
EcoFilterMulti* EcoFilterMulti::copy(){ return new EcoFilterMulti(_testArray,_numTests); }
void EcoFilterMulti::filterSeqs(set<ScoredSeq*>* fullSet,
				set<ScoredSeq*>* retainedSet,
				set<ScoredSeq*>* removedSet,
				int cycleNum, FilterThreadedness threadedness){
  set<ScoredSeq*> runningSet;
  runningSet.insert(fullSet->begin(), fullSet->end());
  for (int testNum = 0; testNum < _numTests; ++testNum){
    set<ScoredSeq*> tempSuccess;
    _testArray[testNum]->filterSeqs(&runningSet, &tempSuccess, removedSet, cycleNum, threadedness);
    runningSet.clear();
    runningSet.insert(tempSuccess.begin(), tempSuccess.end());
  }
  retainedSet->insert(runningSet.begin(), runningSet.end());
}


#endif


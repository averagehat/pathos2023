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

#ifndef SEQFILTER_CPP
#define SEQFILTER_CPP

#include "SeqFilter.h"


				    //#include "OutputFileMulti.h"
#include "RepeatDetector.h"
#include "ReadFileForWriting.h"
#include "WritableSeq.h"
#include <fstream>
#include <sstream>
#include <string.h>
#include <omp.h>
#include <limits.h>



SeqFilter::SeqFilter(){
  _numPairsPerRound = 1000;

  _filenameA = NULL;
  _filenameB = NULL;
  _areSequencesPaired = false;
  _outfileNameA = NULL;
  _outfileNameB = NULL;
  _eitherVsBoth = true;

  _alScoreMatrix = AlignmentScoreMatrix::getDefault();
  _listener = SeqFilterListener::makeStdListener();

  _numThreads = omp_get_max_threads();
}

void SeqFilter::setEitherFailMode(){ _eitherVsBoth = true; }
void SeqFilter::setBothFailMode(){ _eitherVsBoth = false; }


void SeqFilter::setNucMatchScore(long score){ _alScoreMatrix->_match = score; }
void SeqFilter::setNucMismatchPenalty(long penalty){ _alScoreMatrix->_mismatch = penalty; }
void SeqFilter::setOpenGapPenalty(long penalty){ _alScoreMatrix->_newGap = penalty; }
void SeqFilter::setExtendGapPenalty(long penalty){ _alScoreMatrix->_extendGap = penalty; }


SeqFilter::~SeqFilter(){
  delete _alScoreMatrix;

  if(_filenameA != NULL){ delete [] _filenameA; }
  if(_filenameB != NULL){ delete [] _filenameB; }
  if(_outfileNameA != NULL){ delete [] _outfileNameA; }
  if(_outfileNameB != NULL){ delete [] _outfileNameB; }

  delete _listener;

  for (vector<ReadPairFilter*>::iterator it = _waitingRpf.begin(); it != _waitingRpf.end(); ++it){ delete *it; }
}






void SeqFilter::addInputFile(string filename){
  if (_filenameA != NULL){ throw AssemblyException::ArgError("SeqFilter input files were already specified."); }
  _filenameA = new char[ filename.size()+1 ];
  strcpy (_filenameA, filename.c_str());
  _areSequencesPaired = false;
}
void SeqFilter::addInputFiles(string filenameA, string filenameB){
  if (_filenameA != NULL){ throw AssemblyException::ArgError("SeqFilter input files were already specified."); }
  if (filenameA==filenameB){
    throw AssemblyException::ArgError("SeqFilter does not allow paired-end read files to be the same file.");
  }
  _filenameA = new char[ filenameA.size()+1 ];
  strcpy (_filenameA, filenameA.c_str());
  _filenameB = new char[ filenameB.size()+1 ];
  strcpy (_filenameB, filenameB.c_str());
  _areSequencesPaired = true;
}



void SeqFilter::addOutputFile(string filename){
  if (_outfileNameA != NULL){ throw AssemblyException::ArgError("SeqFilter output files were already specified."); }
  _outfileNameA = new char[ filename.size()+1 ];
  strcpy (_outfileNameA, filename.c_str());
}
void SeqFilter::addOutputFiles(string filenameA, string filenameB){
  if (_outfileNameA != NULL){ throw AssemblyException::ArgError("SeqFilter output files were already specified."); }
  _outfileNameA = new char[ filenameA.size()+1 ];
  strcpy (_outfileNameA, filenameA.c_str());
  _outfileNameB = new char[ filenameB.size()+1 ];
  strcpy (_outfileNameB, filenameB.c_str());
}


// FILTER READS

void SeqFilter::addBadSequenceFilter(string badSeqFile, float fractId){
  ReadFileForWriting* rf = ReadFileForWriting::makeBasicReadFile(badSeqFile);
  set<ScoredSeq*> badSeqs;
  rf->open();
  while ( rf->hasRead() ){ badSeqs.insert( rf->getWritableRead() ); }
  rf->close();
  _waitingRpf.push_back( new ReadPairFilterAvoidSeqs(&badSeqs, fractId) );
  for (set<ScoredSeq*>::iterator badIt = badSeqs.begin(); badIt != badSeqs.end(); ++badIt){ delete *badIt; }
  delete rf;
}
void SeqFilter::addGoodSequenceFilter(string goodSeqFile, float fractId){
  ReadFileForWriting* rf = ReadFileForWriting::makeBasicReadFile(goodSeqFile);
  set<ScoredSeq*> goodSeqs;
  rf->open();
  while ( rf->hasRead() ){ goodSeqs.insert( rf->getWritableRead() ); }
  rf->close();
  _waitingRpf.push_back( new ReadPairFilterRetainSeqs(&goodSeqs, fractId) );
  for (set<ScoredSeq*>::iterator goodIt = goodSeqs.begin(); goodIt != goodSeqs.end(); ++goodIt){ delete *goodIt; }
  delete rf;
}
void SeqFilter::addHomopolymerFilter(long maxLength){
  _waitingRpf.push_back( new ReadPairFilterHomoPolymer(maxLength) );
}
void SeqFilter::addDinucRepeatFilter(long maxLength){
  _waitingRpf.push_back( new ReadPairFilterDinucleotide(maxLength) );
}
void SeqFilter::addReadQualityFilter(float minFractGood, float minProbCorrect){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterQuality(maxFractBad, minProbCorrect) );
}
void SeqFilter::addReadQualityFilter(float minFractGood, float minProbCorrect, int numSkipCycles, int numRunCycles){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterQuality(maxFractBad, minProbCorrect, numSkipCycles, numRunCycles) );
}
void SeqFilter::addReadCalledBasesFilter(float minFractGood){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterUncalledBases(maxFractBad) );
}
void SeqFilter::addReadCalledBasesFilter(float minFractGood, int numSkipCycles, int numRunCycles){
  float maxFractBad = 1.0 - minFractGood;
  _waitingRpf.push_back( new ReadPairFilterUncalledBases(maxFractBad, numSkipCycles, numRunCycles) );
}



// FILTER NEW CONTIGS

void SeqFilter::addLengthFilter(long minLength){
  _waitingRpf.push_back( new ReadPairFilterMinLength(minLength) );
}



void SeqFilter::makeLogNull(){
  delete _listener;
  _listener = SeqFilterListener::makeNullListener();
}





void SeqFilter::runFilter(){

  checkFileSpecs();

  bool usePairFile = _filenameB != NULL;
  ReadFileForWriting* readfileA = ReadFileForWriting::makeBasicReadFile(_filenameA);
  ReadFileForWriting* readfileB;
  if (usePairFile){ readfileB = ReadFileForWriting::makeBasicReadFile(_filenameB); }

  ofstream fileStreamA;
  ofstream fileStreamB;
  fileStreamA.open( _outfileNameA );
  if (usePairFile){ fileStreamB.open( _outfileNameB ); }


  // for building arrays
  int numThreadsP1 = _numThreads + 1;
  long numPairsPerRoundP1 = _numPairsPerRound + 1;

  // keeps track of the number of reads that are filtered out by RPFs
  long readsFiltered = 0;


  // these can be set up ahead of time and re-used
  WritableSeq*** seqsToFilterA = new WritableSeq**[ numThreadsP1 ];
  WritableSeq*** seqsToFilterB = new WritableSeq**[ numThreadsP1 ];
  bool** seqPassedA = new bool*[ numThreadsP1 ];
  bool** seqPassedB = new bool*[ numThreadsP1 ];
  bool** comboPassed = new bool*[ numThreadsP1 ];
  for (int tN = 0; tN < _numThreads; ++tN){
    seqsToFilterA[tN] = new WritableSeq*[ numPairsPerRoundP1 ];
    seqsToFilterB[tN] = new WritableSeq*[ numPairsPerRoundP1 ];
    seqPassedA[tN] = new bool[ numPairsPerRoundP1 ];
    seqPassedB[tN] = new bool[ numPairsPerRoundP1 ];
    comboPassed[tN] = new bool[ numPairsPerRoundP1 ];
  }


  // set up an array of filters
  //throw AssemblyException::ImplementationError("filter array needs to be implemented.");


  ReadPairFilter* rpf;
  if (_waitingRpf.empty() ){ rpf = new ReadPairFilterNull(); }
  else if (_waitingRpf.size() == 1){ rpf = *(_waitingRpf.begin()); }
  else {
    ReadPairFilter** rpfArray = new ReadPairFilter*[ _waitingRpf.size() ];
    vector<ReadPairFilter*>::iterator rpfIt = _waitingRpf.begin();
    for (int n = 0; n < _waitingRpf.size(); ++n){
      rpfArray[n] = *rpfIt;
      ++rpfIt;
    }
    rpf = new ReadPairFilterMulti(rpfArray, _waitingRpf.size());
  }


  //open the files
  readfileA->open();
  if (usePairFile){ readfileB->open(); }

  _listener->startingFilter(usePairFile);


  while (readfileA->hasRead()){


    // prepare a subset of reads for mapping - fill these in below
    long numPairs[ numThreadsP1 ];
    for (int tN = 0; tN < _numThreads; ++tN){
      long sN = 0;
      while (readfileA->hasRead() and sN < _numPairsPerRound){
	seqsToFilterA[tN][sN] = readfileA->getWritableRead();
	if (usePairFile){ seqsToFilterB[tN][sN] = readfileB->getWritableRead(); }
	++sN;
      }
      numPairs[tN] = sN;
    }
      

    #pragma omp parallel for schedule(dynamic)
    for (int tN = 0; tN < _numThreads; ++tN){
      if (_areSequencesPaired){
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  seqPassedA[tN][sN] = rpf->isReadOk(seqsToFilterA[tN][sN]);
	  seqPassedB[tN][sN] = rpf->isReadOk(seqsToFilterB[tN][sN]);
	}
      } else {
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  seqPassedA[tN][sN] = rpf->isReadOk(seqsToFilterA[tN][sN]);
	}
      }
    }


    // figure out the right combo
    if (_eitherVsBoth and _areSequencesPaired){
      for (int tN = 0; tN < _numThreads; ++tN){
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  comboPassed[tN][sN] = seqPassedA[tN][sN] && seqPassedB[tN][sN];
	}
      }
    } else if (! _areSequencesPaired){
      for (int tN = 0; tN < _numThreads; ++tN){
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  comboPassed[tN][sN] = seqPassedA[tN][sN];
	}
      }
    } else {
      for (int tN = 0; tN < _numThreads; ++tN){
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  comboPassed[tN][sN] = seqPassedA[tN][sN] || seqPassedB[tN][sN];
	}
      }
    }



    long numProcessed = 0;
    long numRemoved = 0;
    char* tempSeqString;
    if (_areSequencesPaired){
      for (int tN = 0; tN < _numThreads; ++tN){
	numProcessed += numPairs[tN];
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  if (comboPassed[tN][sN]){
	    tempSeqString = seqsToFilterA[tN][sN]->getFileString();
	    fileStreamA << tempSeqString;
	    delete [] tempSeqString;
	    tempSeqString = seqsToFilterB[tN][sN]->getFileString();
	    fileStreamB << tempSeqString;
	    delete [] tempSeqString;
	  } else { ++numRemoved; }
	  delete seqsToFilterA[tN][sN];
	  delete seqsToFilterB[tN][sN];
	}
      }
    } else {
      for (int tN = 0; tN < _numThreads; ++tN){
	numProcessed += numPairs[tN];
	for (long sN = 0; sN < numPairs[tN]; ++sN){
	  if (comboPassed[tN][sN]){
	    tempSeqString = seqsToFilterA[tN][sN]->getFileString();
	    fileStreamA << tempSeqString;
	    delete [] tempSeqString;
	  } else { ++numRemoved; }
	  delete seqsToFilterA[tN][sN];
	}
      }
    }

    _listener->updateProgress(numProcessed,numRemoved);
  }

  _listener->finishingFilter();

  readfileA->close();
  fileStreamA.close();
  delete readfileA;

  if (usePairFile){
    readfileB->close();
    fileStreamB.close();
    delete readfileB;
  }


  if (_waitingRpf.size() != 1){ delete rpf; }


  for (int tN = 0; tN < _numThreads; ++tN){
    delete [] seqsToFilterA[tN];
    delete [] seqsToFilterB[tN];
    delete [] seqPassedA[tN];
    delete [] seqPassedB[tN];
    delete [] comboPassed[tN];
  }
  delete [] seqsToFilterA;
  delete [] seqsToFilterB;
  delete [] seqPassedA;
  delete [] seqPassedB;
  delete [] comboPassed;


}






void SeqFilter::checkFileSpecs(){
  // check specifications
  if (_filenameA==NULL){ throw AssemblyException::ArgError("SeqFilter input files were not specified."); }
  if (_outfileNameA==NULL){ throw AssemblyException::ArgError("SeqFilter output files were not specified."); }
  if (_filenameB==NULL){
    if (_outfileNameB!=NULL){ throw AssemblyException::ArgError("SeqFilter same num of input and output files required."); }
    if (_filenameA==_outfileNameA){ throw AssemblyException::ArgError("Input and output files must have different names."); }
  } else {
    if (_outfileNameB==NULL){ throw AssemblyException::ArgError("SeqFilter same num of input and output files required."); }
    if (_filenameA==_outfileNameA or _filenameB==_outfileNameB or
	_filenameA==_outfileNameB or _filenameB==_outfileNameA){ 
      throw AssemblyException::ArgError("Input and output files must have different names.");
    }
  }
  // check file types
  if (fileUtilities::getFileType(string(_filenameA)) != fileUtilities::getFileType(string(_outfileNameA))){
    cerr << _filenameA <<endl;
    cerr << _outfileNameA <<endl;
    throw AssemblyException::ArgError("input and output file format tags did not match");
  }
  if (_filenameB!=NULL){
    if (fileUtilities::getFileType(string(_filenameA)) != fileUtilities::getFileType(string(_filenameB))){
      cerr << _filenameA <<endl;
      cerr << _filenameB <<endl;
      throw AssemblyException::ArgError("input file format tags did not match");
    } else if (fileUtilities::getFileType(string(_outfileNameA)) != fileUtilities::getFileType(string(_outfileNameB))){
      cerr << _outfileNameA <<endl;
      cerr << _outfileNameB <<endl;
      throw AssemblyException::ArgError("output file format tags did not match");
    }
  }
}








#endif



/*
THE AUTHOR: generated by www.degraeve.com/img2txt.php
::::;::::;::::i:,,,tiiiitt;;,;,,,;,,,,,,
::::;::::;:,i:i;i;;jfGLfji:,,;,,,;,,,,,,
::::;::::,:;,jGGfLGLLGLffLtit;,,,;,,,,,,
::::,::::;;jLGGGEEEDGGGLjfLLt;::,;,,,,,,
::::,::::ijLDGDGLKWKEDGGLffLt;,::;:,,,,,
::::,:::ifDGDDGGDEKWKEDLLLGGft:::;::,,,:
::::,::;tiLDDEDDDDDEKKKKEDLLLLi,:;,::,,,
:..:,::jitEDKKDGLLLLGGDEEDDfGLf,:;::,,,,
:::::::L,LDEWDLLffffLLGDDDDfLGLi:,:::,,:
..::::;jiDDKKLLLfffffffLGDEGLLLf:,:.::::
...:,.;:tfEWDGLLffffjjjjjLEDLDLf::...:.:
:...:..;iEEKGGLLfffjfjtttjLDDDLf::...::.
....:t:.LEEKLGLLfffjjjttiitDEEGLL:....:.
:...:LGLEDDDGLLffLLffjttiiijEKGj;:......
:...:jKKKKWELDDDKEDLfjjjjjiiLEEj::......
..:.tiEEKKWDWWKKWKKDffDKDGj;DEDG;:....:.
..:;iGKEKKKLWKEKKWKELEWEEGG#EEDEKtj;,;;,
ii;iiLEEDKDGKDEEEKEDjEEEGfGt#EDDEDi,,,,,
i;;iiiKKKKDLGGGDEEELtWLEEGjtGGGDDL,,,,,,
iiiiiiGfGEDLLLLGDDELiGfLLftfLGEDLfi;;iii
iiiiiiitKDDLLGEEEEGfiiGLLji;LGEGijiiiiii
iiiitiitLfDLGGDDDLGLiitELfi;jDGiijiiiiii
iiiiiiiittGLGDEEEEEDfLfDDLtiGLiiijiiiitt
tttttiiiitEGGGEDDEDDGffLGDtiGLi;itiiiiii
fffffffLGfEDGGDDWEEGGLLKGLjLjfjiitiiiiii
ttttttttitLEGGLGDEGLjDftftiDjiiiitiiiiii
#########WL#GGGGGDDGGjitftDKjiiiitiiiitt
WWWWWWWWKKfWLGGGGGGGLjtttftiiiiiijtttttt
WWWWWDWWWfEKLGGDLGLLfjtittiiitiiijtttttt
WWWWWKKKfGEDGGGDDGLLfttjitiiitiiitiiiitt
WWWWWKLfLfLDGGGDEDGGLjfjjiittttttftttttj
WWWLLGGLLLGEGGGGDDDDDLjtKG;;tttttftttjjj
jjfjjLLfLDKKGGGGGDGGjtttWD;jt;tttftjjjjj
ttjLLLffLLDKDGDGGGGfjtttKf;Lfi,;tjjjjjjj
tjDLDfGLtjGEDGDGDDDfjjjGWL;;tL;,;:tjjjjf
fGLDLLtjjtfDDDDDDDGLjfLEED,iGfLt,,;::;tj
jfLDfjGGGLtLEDGDDGGLLGELfGt;ijGD;ti;;,,j
LjfDfLGKWWWjfDGGDDDDGDjGjft;;fDLij,ti,:D
fDLLLDEKKWKGLjGLGGDGGtLGjfD:itfjfjt,j;tD
LLDGGEDWDEDGDjjtijLitfDjfjKKitDDjff,fiiL
DjLDfDGGGGEGGGjjttjtGfffGDKWKfGDiGL;tiji
jDLjDfGLLLGDGEftiiijtfLLfGEEWttLfLLf;jit
fLEGLjGjffGGGGEftitjfjGjtjjLKiLGttGftjit
GjfDjffjjfLLLLjDfjfffGLti;jjtiGL;fELfijj
LfjGfffjjfLLfffjDLLfGGjji;ttfLLLtEDLLjGj
fLGGffjjjjfGjjjjfGjLGLttiijtttjGjjfDiLLi
LffLjfjjjffGjjjfLGfGDGjtttttjtjjtfLLfff;
LfffjjjjfffGjjffLtfEDDjtfjjtjjffjffjtfjt
GLffjjfjfffGjjfftiGDDGLfffjtjjfLjfGGjjLj
*/

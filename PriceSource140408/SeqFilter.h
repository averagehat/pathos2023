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

/* This is the programmatic interface.
 */


#ifndef SEQFILTER_H
#define SEQFILTER_H

#include <set>
#include <vector>
#include "ScoredSeq.h"
#include "fileUtilities.h"
#include "ParamsMapping.h"
#include "ParamsCycleManagement.h"
#include "ParamsAlignment.h"
#include "OutputFile.h"
#include "EcoFilter.h"
#include "SeqFilterListener.h"
#include "ReadPairFilter.h"
using namespace::std;

class SeqFilter {

 public:
  SeqFilter();
  ~SeqFilter();

  void runFilter();

  // PAIRED ENDS
  // single files, alternating pairs
  void addInputFile(string filename);
  // two files
  void addInputFiles(string filenameA, string filenameB);


  void addOutputFile(string filename);
  void addOutputFiles(string filenameA, string filenameB);

  // FILTERS
  void addLengthFilter(long minLength);
  void addReadQualityFilter(float minFractGood, float minProbCorrect);
  void addReadQualityFilter(float minFractGood, float minProbCorrect, int numSkipCycles, int numRunCycles);
  void addReadCalledBasesFilter(float minFractGood);
  void addReadCalledBasesFilter(float minFractGood, int numSkipCycles, int numRunCycles);
  void addBadSequenceFilter(string badSeqFile, float fractId);
  void addGoodSequenceFilter(string goodSeqFile, float fractId);
  void addHomopolymerFilter(long maxLength);
  void addDinucRepeatFilter(long maxLength);

  void findRepsBeforeCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
			   long minRepeatSize, float matchFractId);
  void findRepsBeforeCycle(int cycleNum, float numStdDevs, float minFoldIncrease,
			   long minRepeatSize, float matchFractId, string outfileName);


  void makeLogNull();




  // sequences are eliminated if either fail or if both fail
  void setEitherFailMode();
  void setBothFailMode();


  // ...and alignment score matrix values
  void setNucMatchScore(long score);
  void setNucMismatchPenalty(long penalty);
  void setOpenGapPenalty(long penalty);
  void setExtendGapPenalty(long penalty);

 private:

  long _numPairsPerRound;

  // only A will be used if file is only single-direction
  char* _filenameA;
  char* _filenameB;
  bool _areSequencesPaired;
  char* _outfileNameA;
  char* _outfileNameB;

  // true = eliminate if either fails, false = eliminate if both fail
  bool _eitherVsBoth;

  void checkFileSpecs();

  int _numThreads;

  // output files as they are being collected
  set<OutputFile*> _outfileSet;
  // the one that will be used; should be generated when assembly is run and deleted at the end of the run
  OutputFile* _outfile;

  // these two are used for set-up
  SeqFilterListener* _listener;


  // FILTERING OF READS
  // collects the input
  vector<ReadPairFilter*> _waitingRpf;

  // the alignment score parameters are the same for both assembly steps
  AlignmentScoreMatrix* _alScoreMatrix;

};

#endif

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

/*
 */

#ifndef WRITABLESEQFASTQ_H
#define WRITABLESEQFASTQ_H

#include "WritableSeq.h"
#include "AssemblyException.h"
#include "ReadFileFastqSingle.h"

class WritableSeqFastq : public WritableSeq {
 public:
  // the third line of fastq entries is not required
  WritableSeqFastq(string line1, string line2, string line4, fileUtilities::FileType fileType, string errString);
  void OK();
  ~WritableSeqFastq();


  /* see parent abstract class WritableSeq */
  char* getFileString();
  fileUtilities::FileType fileType();

  /* see grandparent abstract class ScoredSeq */
  float scoreAtPosition(long position, char sense);
  float scoreAtPosPlus(long position);
  float scoreAtPosMinus(long position);
  float linkAfterPosition(long position, char sense);
  float linkAfterPosPlus(long position);
  float linkAfterPosMinus(long position);
  char nucAtPosition(long position, char sense);
  char nucAtPosPlus(long position);
  char nucAtPosMinus(long position);
  char* getSeq(char sense);
  float* getScores(char sense);
  float* getLinks(char sense);
  char* getSubseq(long position, long length, char sense);
  float* getSubScores(long position, long length, char sense);
  float* getSubLinks(long position, long length, char sense);
  void gatherSeq(char* collector, char sense);
  void gatherScores(float* collector, char sense);
  void gatherLinks(float* collector, char sense);
  void gatherSubseq(char* collector, long position, long length, char sense);
  void gatherSubScores(float* collector, long position, long length, char sense);
  void gatherSubLinks(float* collector, long position, long length, char sense);
  void buffer();
  void unbuffer();
  void bottomBuffer();
  long size();
  bool hasPairedEnd();
  ScoredSeq * getPairedEnd();
  ScoredSeq * getTempPairedEnd();
  ScoredSeq * shallowCopy();
  ScoredSeq * flipCopy();
  bool isNested();
  void deepDelete();
  void deepUnbuffer();

 private:
  char* _seq;
  long _size;
  char* _rcSeq;
  float* _scores;
  char* _nameLine;
  long _nameSize;
  char* _scoreChars;
  fileUtilities::FileType _fileType;

  // note: this should always be 1!!!!
  float _countFactor;
};

#endif

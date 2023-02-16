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


#ifndef WRITABLESEQFASTQ_CPP
#define WRITABLESEQFASTQ_CPP

#include "WritableSeqFastq.h"
#include <cstring>
using namespace::std;
using namespace::fileUtilities;


void WritableSeqFastq::OK() {
}

// OVERLOADED CONSTRUCTORS (this goes on for a while)

WritableSeqFastq::~WritableSeqFastq(){
  delete [] _seq;
  delete [] _rcSeq;
  delete [] _scores;
  delete [] _scoreChars;
  delete [] _nameLine;
}



WritableSeqFastq::WritableSeqFastq(string line1, string line2, string line4, fileUtilities::FileType fileType, string errString) :
  _fileType(fileType),
  _countFactor(1.0){
  bool invert = false;

  // make sure the file type and encoding agree
  if (! hasAsciiScoreEncoding(fileType)){
    throw AssemblyException::FileError("RFFQ encoding requested with inappropriate file type, in WSFastq.");
  }

  long line2Size = line2.size();
  char* line2temp = new char[line2Size+1];
  strcpy(line2temp, line2.c_str());

  // interpret the sequence line
  set<char> okToSkip;
  getOkToSkipChars(&okToSkip, FASTQFILE);
  ReadFileStringInterpreter* si = new ReadFileStringInterpreter(line2temp, line2Size, line2Size, &okToSkip, invert);
  _seq = si->_seqString;
  _size = si->_seqLen;
  delete [] line2temp;
  delete si;

  // calculate the scores
  ScoreConverter* scoreConvo = getScoreConverter(fileType, _countFactor, errString);

  _scoreChars = new char[line4.size()+1];
  strcpy(_scoreChars, line4.c_str());
  _scoreChars[_size] = '\0';
  _scores = scoreConvo->cstringToScores(_scoreChars, _size, invert, _seq);
  delete scoreConvo;

  // determine the RC, now that the seq has been modified base on scores
  _rcSeq = reverseComplement(_seq,_size);

  // save the name line and its size with end-of-line whitespace removed
  _nameLine = new char[line1.size()+1];
  strcpy(_nameLine, line1.c_str());
  _nameSize = line1.size() - 1;
  while (_nameSize >=0 and okToSkip.count(_nameLine[_nameSize]) > 0){ --_nameSize; }
  _nameSize += 1;

  // the other name line will just get the '+' (meets specs)
}




char* WritableSeqFastq::getFileString(){
  // seq/score, '+' line, title, 4 newlines
  long entrySize = _size * 2 + 1 + _nameSize + 4;

  char* fileString = new char[entrySize+1];
  for (long n = 0; n < _nameSize; ++n){ fileString[n] = _nameLine[n]; }
  fileString[_nameSize] = '\n';
  long sectionStart = _nameSize + 1;
  for (long n = 0; n < _size; ++n){ fileString[n+sectionStart] = _seq[n]; }
  fileString[sectionStart + _size] = '\n';
  fileString[sectionStart + _size + 1] = '+';
  fileString[sectionStart + _size + 2] = '\n';
  sectionStart += _size + 3;
  for (long n = 0; n < _size; ++n){ fileString[n+sectionStart] = _scoreChars[n]; }
  fileString[sectionStart + _size] = '\n';
  sectionStart += _size + 1;
  fileString[sectionStart] = '\0';

  return fileString;
}
fileUtilities::FileType WritableSeqFastq::fileType(){ return _fileType; }



float WritableSeqFastq::scoreAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _scores[position];
  case '-': return _scores[_size - position - 1];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float WritableSeqFastq::scoreAtPosPlus(long position) { return _scores[position]; }
float WritableSeqFastq::scoreAtPosMinus(long position) { return _scores[_size - position - 1]; }

float WritableSeqFastq::linkAfterPosition(long position, char sense) { return _countFactor;}
float WritableSeqFastq::linkAfterPosPlus(long position){ return _countFactor; }
float WritableSeqFastq::linkAfterPosMinus(long position){ return _countFactor; }


char WritableSeqFastq::nucAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _seq[ position ];
  case '-': return _rcSeq[ position ];
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
char WritableSeqFastq::nucAtPosPlus(long position) { return _seq[ position ]; }
char WritableSeqFastq::nucAtPosMinus(long position) { return _rcSeq[ position ]; }




char* WritableSeqFastq::getSeq(char sense) {
  char* seq = new char[_size + 1];
  gatherSeq(seq,sense);
  return seq;
}
float* WritableSeqFastq::getScores(char sense) {
  float* scores = new float[_size + 1];
  gatherScores(scores,sense);
  return scores;
}
float* WritableSeqFastq::getLinks(char sense) {
  float* links = new float[_size + 1];
  gatherLinks(links,sense);
  return links;
}

char* WritableSeqFastq::getSubseq(long position, long length, char sense){
  char* subseq = new char[ length + 1 ];
  gatherSubseq(subseq,position,length,sense);
  return subseq;
}
float* WritableSeqFastq::getSubScores(long position, long length, char sense){
  float* scores = new float[length + 1];
  gatherSubScores(scores,position,length,sense);
  return scores;
}
float* WritableSeqFastq::getSubLinks(long position, long length, char sense){
  float* links = new float[length + 1];
  gatherSubLinks(links,position,length,sense);
  return links;
}


void WritableSeqFastq::gatherSeq(char* collector, char sense) {
  switch (sense){
  case '+': strcpy(collector, _seq); break;
  case '-': strcpy(collector, _rcSeq); break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFastq::gatherScores(float* collector, char sense) {
  long sizeM1 = _size - 1;
  switch (sense){
  case '+': for (long n = 0; n < _size; ++n){ collector[n] = _scores[n]; } break;
  case '-': for (long n = 0; n < _size; ++n){ collector[n] = _scores[sizeM1 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFastq::gatherLinks(float* collector, char sense) {
  long sizeM1 = _size - 1;
  for (long n = 0; n < sizeM1; ++n){ collector[n] = _countFactor; }
}



void WritableSeqFastq::gatherSubseq(char* collector, long position, long length, char sense){
  collector[length] = '\0';
  char* seqOs;
  switch (sense){
  case '+': seqOs = &_seq[position]; for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; } break;
  case '-': seqOs = &_rcSeq[position]; for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFastq::gatherSubScores(float* collector, long position, long length, char sense){
  float* plusOs;
  long sizeM1mPos = length - 1 - position;
  switch (sense){
  case '+': plusOs = &_scores[position]; for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _scores[sizeM1mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFastq::gatherSubLinks(float* collector, long position, long length, char sense){
  for (long n = 0; n < length; ++n){ collector[n] = _countFactor; }
}



long WritableSeqFastq::size() {
  //OK();
  return _size;
}

void WritableSeqFastq::buffer(){}
void WritableSeqFastq::unbuffer(){}
void WritableSeqFastq::bottomBuffer(){}
void WritableSeqFastq::deepDelete(){ delete this; }
void WritableSeqFastq::deepUnbuffer(){} // unbuffer doesn't do anything, so no need to call it
bool WritableSeqFastq::hasPairedEnd() { return false; }
ScoredSeq * WritableSeqFastq::getPairedEnd() {
  throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
}
ScoredSeq * WritableSeqFastq::getTempPairedEnd() {
  throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
}
ScoredSeq * WritableSeqFastq::shallowCopy(){ return new ScoredSeqShallow(this,'+'); }
ScoredSeq * WritableSeqFastq::flipCopy(){ return new ScoredSeqShallow(this,'-'); }
bool WritableSeqFastq::isNested(){ return false; }


#endif

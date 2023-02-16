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


#ifndef WRITABLESEQPRICEQ_CPP
#define WRITABLESEQPRICEQ_CPP

#include "WritableSeqPriceq.h"
#include <cstring>
using namespace::std;
using namespace::fileUtilities;

void WritableSeqPriceq::OK() {
}

// OVERLOADED CONSTRUCTORS (this goes on for a while)

WritableSeqPriceq::~WritableSeqPriceq(){
  delete [] _seq;
  delete [] _rcSeq;
  delete [] _scores;
  delete [] _scoreChars;
  delete [] _links;
  delete [] _linkChars;
  delete [] _nameLine;
}



WritableSeqPriceq::WritableSeqPriceq(string line1, string line2, string line4, string line6, string errString){

  bool invert = false;
  set<char> okToSkip;
  getOkToSkipChars(&okToSkip, PRICEQFILE);

  // save the name line and its size with end-of-line whitespace removed
  _nameLine = new char[line1.size()+1];
  strcpy(_nameLine, line1.c_str());
  _nameSize = line1.size() - 1;
  while (_nameSize >=0 and okToSkip.count(_nameLine[_nameSize]) > 0){ --_nameSize; }
  _nameSize += 1;
  // the other name lines will just get the '+' and '~' (meets specs)

  long line2Size = line2.size();
  char* line2temp = new char[line2Size+1];
  strcpy(line2temp, line2.c_str());

  // interpret the sequence line
  ReadFileStringInterpreter* si = new ReadFileStringInterpreter(line2temp, line2Size, line2Size, &okToSkip, invert);
  _seq = si->_seqString;
  _size = si->_seqLen;
  delete [] line2temp;
  delete si;

  // the scores
  _scoreChars = new char[line4.size()+1];
  strcpy(_scoreChars, line4.c_str());
  _scoreChars[_size] = '\0';
  _scores = cstringToScoresPriceq(_scoreChars, _size, invert);

  // determine the RC, now that the seq has been modified base on scores
  _rcSeq = reverseComplement(_seq,_size);

  // the links
  _linkChars = new char[line6.size()+1];
  strcpy(_linkChars, line6.c_str());
  _linkChars[_size-1] = '\0';
  _links = cstringToScoresPriceq(_linkChars, _size-1, invert);
}




char* WritableSeqPriceq::getFileString(){
  // seq/score/links (last is 1 shorter), '+' and '~' lines, title, 6 newlines
  long entrySize = _size * 3 - 1 + 2 + _nameSize + 6;

  // name and seq lines
  char* fileString = new char[entrySize+1];
  for (long n = 0; n < _nameSize; ++n){ fileString[n] = _nameLine[n]; }
  fileString[_nameSize] = '\n';
  long sectionStart = _nameSize + 1;
  for (long n = 0; n < _size; ++n){ fileString[n+sectionStart] = _seq[n]; }
  fileString[sectionStart + _size] = '\n';
  sectionStart += _size + 1;

  // nuc score lines
  fileString[sectionStart] = '+';
  fileString[sectionStart + 1] = '\n';
  sectionStart += 2;

  for (long n = 0; n < _size; ++n){ fileString[n+sectionStart] = _scoreChars[n]; }
  fileString[sectionStart + _size] = '\n';
  sectionStart += _size + 1;

  // link score lines
  fileString[sectionStart] = '~';
  fileString[sectionStart + 1] = '\n';
  sectionStart += 2;

  long sizeM1 = _size - 1;
  for (long n = 0; n < sizeM1; ++n){ fileString[n+sectionStart] = _linkChars[n]; }
  fileString[sectionStart + sizeM1] = '\n';
  sectionStart += sizeM1 + 1;

  fileString[sectionStart] = '\0';
  return fileString;
}
fileUtilities::FileType WritableSeqPriceq::fileType(){ return fileUtilities::PRICEQFILE; }



float WritableSeqPriceq::scoreAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _scores[position];
  case '-': return _scores[_size - position - 1];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float WritableSeqPriceq::scoreAtPosPlus(long position) { return _scores[position]; }
float WritableSeqPriceq::scoreAtPosMinus(long position) { return _scores[_size - position - 1]; }

float WritableSeqPriceq::linkAfterPosition(long position, char sense) { 
  switch (sense){
  case '+': return _links[position];
  case '-': return _links[_size - position - 2];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float WritableSeqPriceq::linkAfterPosPlus(long position){ return _links[position]; }
float WritableSeqPriceq::linkAfterPosMinus(long position){ return _links[_size - position - 2]; }


char WritableSeqPriceq::nucAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _seq[ position ];
  case '-': return _rcSeq[ position ];
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
char WritableSeqPriceq::nucAtPosPlus(long position) { return _seq[ position ]; }
char WritableSeqPriceq::nucAtPosMinus(long position) { return _rcSeq[ position ]; }



char* WritableSeqPriceq::getSeq(char sense) {
  char* seq = new char[_size + 1];
  gatherSeq(seq,sense);
  return seq;
}
float* WritableSeqPriceq::getScores(char sense) {
  float* scores = new float[_size + 1];
  gatherScores(scores,sense);
  return scores;
}
float* WritableSeqPriceq::getLinks(char sense) {
  float* links = new float[_size + 1];
  gatherLinks(links,sense);
  return links;
}

char* WritableSeqPriceq::getSubseq(long position, long length, char sense){
  char* subseq = new char[ length + 1 ];
  gatherSubseq(subseq,position,length,sense);
  return subseq;
}
float* WritableSeqPriceq::getSubScores(long position, long length, char sense){
  float* scores = new float[length + 1];
  gatherSubScores(scores,position,length,sense);
  return scores;
}
float* WritableSeqPriceq::getSubLinks(long position, long length, char sense){
  float* links = new float[length + 1];
  gatherSubLinks(links,position,length,sense);
  return links;
}


void WritableSeqPriceq::gatherSeq(char* collector, char sense) {
  switch (sense){
  case '+': strcpy(collector, _seq); break;
  case '-': strcpy(collector, _rcSeq); break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqPriceq::gatherScores(float* collector, char sense) {
  long sizeM1 = _size - 1;
  switch (sense){
  case '+': for (long n = 0; n < _size; ++n){ collector[n] = _scores[n]; } break;
  case '-': for (long n = 0; n < _size; ++n){ collector[n] = _scores[sizeM1 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqPriceq::gatherLinks(float* collector, char sense) {
  long sizeM1 = _size - 1;
  long sizeM2 = _size - 2;
  switch (sense){
  case '+': for (long n = 0; n < sizeM1; ++n){ collector[n] = _links[n]; } break;
  case '-': for (long n = 0; n < sizeM1; ++n){ collector[n] = _links[sizeM2 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}

void WritableSeqPriceq::gatherSubseq(char* collector, long position, long length, char sense){
  collector[length] = '\0';
  char* seqOs;
  switch (sense){
  case '+': seqOs = &_seq[position]; for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; } break;
  case '-': seqOs = &_rcSeq[position]; for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqPriceq::gatherSubScores(float* collector, long position, long length, char sense){
  float* plusOs;
  long sizeM1mPos = length - 1 - position;
  switch (sense){
  case '+': plusOs = &_scores[position]; for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _scores[sizeM1mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqPriceq::gatherSubLinks(float* collector, long position, long length, char sense){
  float* plusOs;
  long sizeM2mPos = length - 2 - position;
  switch (sense){
  case '+': plusOs = &_links[position]; for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _links[sizeM2mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}



long WritableSeqPriceq::size() {
  //OK();
  return _size;
}

void WritableSeqPriceq::buffer(){}
void WritableSeqPriceq::unbuffer(){}
void WritableSeqPriceq::bottomBuffer(){}
void WritableSeqPriceq::deepDelete(){ delete this; }
void WritableSeqPriceq::deepUnbuffer(){} // unbuffer doesn't do anything, so no need to call it
bool WritableSeqPriceq::hasPairedEnd() { return false; }
ScoredSeq * WritableSeqPriceq::getPairedEnd() {
  throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
}
ScoredSeq * WritableSeqPriceq::getTempPairedEnd() {
  throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
}
ScoredSeq * WritableSeqPriceq::shallowCopy(){ return new ScoredSeqShallow(this,'+'); }
ScoredSeq * WritableSeqPriceq::flipCopy(){ return new ScoredSeqShallow(this,'-'); }
bool WritableSeqPriceq::isNested(){ return false; }


#endif

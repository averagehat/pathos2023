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


#ifndef WRITABLESEQFASTA_CPP
#define WRITABLESEQFASTA_CPP

#include "WritableSeqFasta.h"
#include <cstring>
using namespace::std;
using namespace::fileUtilities;


void WritableSeqFasta::OK() {
}

// OVERLOADED CONSTRUCTORS (this goes on for a while)

WritableSeqFasta::~WritableSeqFasta(){
  delete [] _seq;
  delete [] _rcSeq;
  delete [] _scores;
  delete [] _nameLine;
}



WritableSeqFasta::WritableSeqFasta(string nameLine, vector<string>* seqLines) :
  _countFactor(1.0),
  _lineSize(60){
  bool invert = false;

  set<char> whitespace;
  set<char> okToSkip;
  getWhitespaceChars(&whitespace, FASTAFILE);
  getOkToSkipChars(&okToSkip, FASTAFILE);

  // save the name line and its size with end-of-line whitespace removed
  _nameLine = new char[nameLine.size()+1];
  strcpy(_nameLine, nameLine.c_str());
  _nameSize = nameLine.size() - 1;
  while (_nameSize >=0 and whitespace.count(_nameLine[_nameSize]) > 0){ --_nameSize; }
  _nameSize += 1;

  // convert the strings to char arrays and get the length stats
  long maxLen = 0;
  char** allLines = new char*[ seqLines->size() + 1 ];
  long* lineSizes = new long[ seqLines->size() + 1 ];
  long lineNum = 0;
  for (vector<string>::iterator it =  seqLines->begin(); it != seqLines->end(); ++it){ 
    string oneLine = *it;
    long lineSize = oneLine.size();
    maxLen += lineSize;
    lineSizes[lineNum] = lineSize;
    allLines[lineNum] = new char[ lineSize + 1 ];
    strcpy(allLines[lineNum], oneLine.c_str());
    lineNum++;
  }

  // actual conversion into a sequence
  if (lineNum == 0){
    _seq = new char[1];
    _seq[0] = '\0';
    _rcSeq = new char[1];
    _rcSeq[0] = '\0';
    _size = 0;
    _scores = new float[1];
  } else {
    ReadFileStringInterpreter* si = new ReadFileStringInterpreter(allLines[0], lineSizes[0], maxLen, &okToSkip, false);
    for (long n = 1; n < lineNum; ++n){ si->addSeq(allLines[n], lineSizes[n], &okToSkip); }
    _seq = si->_seqString;
    _size = si->_seqLen;
    _rcSeq = reverseComplement(_seq,_size);
    _scores = new float[_size + 1];
    for (long n = 0; n < _size; ++n){
      if (_seq[n]=='N'){ _scores[n] = 0.0; }
      else { _scores[n] = _countFactor; }
    }
    delete si;
  }

  // clean up
  for (long n = 0; n < lineNum; ++n){ delete [] allLines[n]; }
  delete [] allLines;
  delete [] lineSizes;
}




char* WritableSeqFasta::getFileString(){
  long numLines = _size / _lineSize;
  if (_size % _lineSize != 0){ numLines += 1; }

  // seq length, seq newlines, title, t newline
  long entrySize = _size + numLines + _nameSize + 1;
  char* fileString = new char[entrySize+1];

  for (long n = 0; n < _nameSize; ++n){ fileString[n] = _nameLine[n]; }
  fileString[_nameSize] = '\n';
  long sectionStart = _nameSize + 1;

  // add only the specified number of characters at a time MAX
  long startPos = 0;
  while (startPos < _size){
    int lineLen;
    if (startPos + _lineSize > _size){ lineLen = _size - startPos; }
    else { lineLen = _lineSize; }
    for (long n = 0; n < lineLen; ++n){ fileString[n+sectionStart] = _seq[n+startPos]; }
    fileString[lineLen+sectionStart] = '\n';
    sectionStart += lineLen + 1;
    startPos += lineLen;
  }
  fileString[sectionStart] = '\0';

  return fileString;
}
fileUtilities::FileType WritableSeqFasta::fileType(){ return fileUtilities::FASTAFILE; }



float WritableSeqFasta::scoreAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _scores[position];
  case '-': return _scores[_size - position - 1];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float WritableSeqFasta::scoreAtPosPlus(long position) { return _scores[position]; }
float WritableSeqFasta::scoreAtPosMinus(long position) { return _scores[_size - position - 1]; }

float WritableSeqFasta::linkAfterPosition(long position, char sense) { return _countFactor;}
float WritableSeqFasta::linkAfterPosPlus(long position){ return _countFactor; }
float WritableSeqFasta::linkAfterPosMinus(long position){ return _countFactor; }


char WritableSeqFasta::nucAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _seq[ position ];
  case '-': return _rcSeq[ position ];
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
char WritableSeqFasta::nucAtPosPlus(long position) { return _seq[ position ]; }
char WritableSeqFasta::nucAtPosMinus(long position) { return _rcSeq[ position ]; }

char* WritableSeqFasta::getSeq(char sense) {
  char* seq = new char[_size + 1];
  gatherSeq(seq,sense);
  return seq;
}
float* WritableSeqFasta::getScores(char sense) {
  float* scores = new float[_size + 1];
  gatherScores(scores,sense);
  return scores;
}
float* WritableSeqFasta::getLinks(char sense) {
  float* links = new float[_size + 1];
  gatherLinks(links,sense);
  return links;
}

char* WritableSeqFasta::getSubseq(long position, long length, char sense){
  char* subseq = new char[ length + 1 ];
  gatherSubseq(subseq,position,length,sense);
  return subseq;
}
float* WritableSeqFasta::getSubScores(long position, long length, char sense){
  float* scores = new float[length + 1];
  gatherSubScores(scores,position,length,sense);
  return scores;
}
float* WritableSeqFasta::getSubLinks(long position, long length, char sense){
  float* links = new float[length + 1];
  gatherSubLinks(links,position,length,sense);
  return links;
}



void WritableSeqFasta::gatherSeq(char* collector, char sense) {
  //OK();
  switch (sense){
  case '+': strcpy(collector, _seq); break;
  case '-': strcpy(collector, _rcSeq); break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFasta::gatherScores(float* collector, char sense) {
  long sizeM1 = _size - 1;
  switch (sense){
  case '+': for (long n = 0; n < _size; ++n){ collector[n] = _scores[n]; } break;
  case '-': for (long n = 0; n < _size; ++n){ collector[n] = _scores[sizeM1 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFasta::gatherLinks(float* collector, char sense) {
  long sizeM1 = _size - 1;
  for (long n = 0; n < sizeM1; ++n){ collector[n] = _countFactor; }
}

void WritableSeqFasta::gatherSubseq(char* collector, long position, long length, char sense){
  collector[length] = '\0';
  char* seqOs;
  switch (sense){
  case '+': seqOs = &_seq[position]; for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; } break;
  case '-': seqOs = &_rcSeq[position]; for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFasta::gatherSubScores(float* collector, long position, long length, char sense){
  float* plusOs;
  long sizeM1mPos = length - 1 - position;
  switch (sense){
  case '+': plusOs = &_scores[position]; for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _scores[sizeM1mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void WritableSeqFasta::gatherSubLinks(float* collector, long position, long length, char sense){
  for (long n = 0; n < length; ++n){ collector[n] = _countFactor; }
}



long WritableSeqFasta::size() {
  //OK();
  return _size;
}

void WritableSeqFasta::buffer(){}
void WritableSeqFasta::unbuffer(){}
void WritableSeqFasta::bottomBuffer(){}
void WritableSeqFasta::deepDelete(){ delete this; }
void WritableSeqFasta::deepUnbuffer(){} // unbuffer doesn't do anything, so no need to call it
bool WritableSeqFasta::hasPairedEnd() { return false; }
ScoredSeq * WritableSeqFasta::getPairedEnd() {
  throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
}
ScoredSeq * WritableSeqFasta::getTempPairedEnd() {
  throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
}
ScoredSeq * WritableSeqFasta::shallowCopy(){ return new ScoredSeqShallow(this,'+'); }
ScoredSeq * WritableSeqFasta::flipCopy(){ return new ScoredSeqShallow(this,'-'); }
bool WritableSeqFasta::isNested(){ return false; }


#endif

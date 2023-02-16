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


#ifndef SCOREDSEQSHALLOW_CPP
#define SCOREDSEQSHALLOW_CPP

#include "ScoredSeqShallow.h"
#include <cstring>
#include <iostream>
using namespace::std;

long ScoredSeqShallow::_readCount = 0;

ScoredSeqShallow::ScoredSeqShallow(){}; //default constructor

void ScoredSeqShallow::OK() {
}

// OVERLOADED CONSTRUCTORS (this goes on for a while)

ScoredSeqShallow::~ScoredSeqShallow(){
  delete [] _seq;
  delete [] _rcSeq;
  delete [] _scores;
  delete [] _links;
  //_readCount--;
}


void ScoredSeqShallow::constructorSeqHelper(string seq){
  _size = seq.size();
  _seq = new char [_size+1];
  strcpy (_seq, seq.c_str());
  _rcSeq = reverseComplement(_seq,_size);
}


ScoredSeqShallow::ScoredSeqShallow(string seq, vector<float>* scores) :
  _pe(0)
{
  constructorSeqHelper(seq);
  long seqSizeMin1 = _size - 1;
  _links = new float[ seqSizeMin1 ];
  _scores = new float[ seqSizeMin1 + 1 ];
  for (long n = 0; n < seqSizeMin1; ++n){
    _links[n] = 1.0;
    _scores[n] = (*scores)[n];
  }
  _scores[ seqSizeMin1 ] = (*scores)[ seqSizeMin1 ];
}


ScoredSeqShallow::ScoredSeqShallow(string seq, vector<float>* scores, vector<float>* links) :
  _pe(0)
{
  constructorSeqHelper(seq);
  long seqSizeMin1 = _size - 1;
  _links = new float[ seqSizeMin1 ];
  _scores = new float[ seqSizeMin1 + 1 ];
  for (long n = 0; n < seqSizeMin1; ++n){
    _links[n] = (*links)[n];
    _scores[n] = (*scores)[n];
  }
  _scores[ seqSizeMin1 ] = (*scores)[ seqSizeMin1 ];
}

// STREAMLINED
ScoredSeqShallow::ScoredSeqShallow(string seq, float* scores) :
  _pe(0)
{
  constructorSeqHelper(seq);
  long seqSizeMin1 = _size - 1;
  _links = new float[ seqSizeMin1 ];
  _scores = new float[ seqSizeMin1 + 1 ];
  for (long n = 0; n < seqSizeMin1; ++n){
    _links[n] = 1.0;
    _scores[n] = scores[n];
  }
  _scores[ seqSizeMin1 ] = scores[ seqSizeMin1 ];
}
ScoredSeqShallow::ScoredSeqShallow(string seq, float* scores, float* links) :
  _pe(0)
{
  constructorSeqHelper(seq);
  long seqSizeMin1 = _size - 1;
  _links = new float[ seqSizeMin1 ];
  _scores = new float[ seqSizeMin1 + 1 ];
  for (long n = 0; n < seqSizeMin1; ++n){
    _links[n] = links[n];
    _scores[n] = scores[n];
  }
  _scores[ seqSizeMin1 ] = scores[ seqSizeMin1 ];
}



// CHAR ARRAYS
ScoredSeqShallow::ScoredSeqShallow(char* seq, float* scores, long size) :
  _size(size),
  _pe(0)
{
  long seqSizeMin1 = _size - 1;
  _seq = new char[_size + 1];
  _seq[_size] = '\0';
  _links = new float[seqSizeMin1];
  _scores = new float[_size];
  for (long n = 0; n < seqSizeMin1; ++n){
    _links[n] = 1.0;
    _scores[n] = scores[n];
    _seq[n] = seq[n];
  }
  _scores[ seqSizeMin1 ] = scores[ seqSizeMin1 ];
  _seq[ seqSizeMin1 ] = seq[ seqSizeMin1 ];
  _rcSeq = reverseComplement(_seq,_size);
}
ScoredSeqShallow::ScoredSeqShallow(char* seq, float* scores, float* links, long size) :
  _size(size),
  _pe(0)
{
  long seqSizeMin1 = _size - 1;
  _seq = new char[_size + 1];
  _seq[_size] = '\0';
  _links = new float[seqSizeMin1];
  _scores = new float[_size];
  for (long n = 0; n < seqSizeMin1; ++n){
    _links[n] = links[n];
    _scores[n] = scores[n];
    _seq[n] = seq[n];
  }
  _scores[ seqSizeMin1 ] = scores[ seqSizeMin1 ];
  _seq[ seqSizeMin1 ] = seq[ seqSizeMin1 ];
  _rcSeq = reverseComplement(_seq,_size);
}


// REP EXPOSURE
ScoredSeqShallow::ScoredSeqShallow(bool expRep, char* seq, float* scores, float* links, long size) :
  _size(size),
  _pe(0)
{
  if (expRep){
    _seq = seq;
    _scores = scores;
    _links = links;
  } else {
    long seqSizeMin1 = _size - 1;
    _seq = new char[_size + 1];
    _seq[_size] = '\0';
    _links = new float[seqSizeMin1];
    _scores = new float[_size];
    for (long n = 0; n < seqSizeMin1; ++n){
      _links[n] = links[n];
      _scores[n] = scores[n];
      _seq[n] = seq[n];
    }
    _scores[ seqSizeMin1 ] = scores[ seqSizeMin1 ];
    _seq[ seqSizeMin1 ] = seq[ seqSizeMin1 ];
  }
  _rcSeq = reverseComplement(_seq,_size);
}




ScoredSeqShallow::ScoredSeqShallow(ScoredSeq *seqToCopy, char sense) : _size(seqToCopy->size()) {
  _seq = seqToCopy->getSeq(sense);
  _scores = seqToCopy->getScores(sense); //new float[ seqSizeMin1 + 1 ];
  _links = seqToCopy->getLinks(sense); //new float[ seqSizeMin1 ];
  _rcSeq = reverseComplement(_seq,_size);
  _pe = 0;
  //_readCount++;
}


float ScoredSeqShallow::scoreAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _scores[position];
  case '-': return _scores[_size - position - 1];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float ScoredSeqShallow::scoreAtPosPlus(long position) { return _scores[position]; }
float ScoredSeqShallow::scoreAtPosMinus(long position) { return _scores[_size - position - 1]; }

float ScoredSeqShallow::linkAfterPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _links[position];
  case '-': return _links[_size - position - 2];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float ScoredSeqShallow::linkAfterPosPlus(long position){ return _links[position]; }
float ScoredSeqShallow::linkAfterPosMinus(long position){ return _links[_size - position - 2]; }


char ScoredSeqShallow::nucAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _seq[ position ];
  case '-': return _rcSeq[ position ];
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
char ScoredSeqShallow::nucAtPosPlus(long position) { return _seq[ position ]; }
char ScoredSeqShallow::nucAtPosMinus(long position) { return _rcSeq[ position ]; }

char* ScoredSeqShallow::getSeq(char sense) {
  //OK();
  char* seq = new char[_size + 1];
  gatherSeq(seq,sense);
  return seq;
}
float* ScoredSeqShallow::getScores(char sense) {
  float* scores = new float[_size + 1];
  gatherScores(scores,sense);
  return scores;
}
float* ScoredSeqShallow::getLinks(char sense) {
  float* links = new float[_size + 1];
  gatherLinks(links,sense);
  return links;
}


char* ScoredSeqShallow::getSubseq(long position, long length, char sense){
  char* subseq = new char[ length + 1 ];
  gatherSubseq(subseq,position,length,sense);
  return subseq;
}
float* ScoredSeqShallow::getSubScores(long position, long length, char sense){
  float* scores = new float[length + 1];
  gatherSubScores(scores,position,length,sense);
  return scores;
}
float* ScoredSeqShallow::getSubLinks(long position, long length, char sense){
  float* links = new float[length + 1];
  gatherSubLinks(links,position,length,sense);
  return links;
}



void ScoredSeqShallow::gatherSeq(char* collector, char sense) {
  //OK();
  switch (sense){
  case '+': strcpy(collector, _seq); break;
  case '-': strcpy(collector, _rcSeq); break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqShallow::gatherScores(float* collector, char sense) {
  long sizeM1 = _size - 1;
  switch (sense){
  case '+': for (long n = 0; n < _size; ++n){ collector[n] = _scores[n]; } break;
  case '-': for (long n = 0; n < _size; ++n){ collector[n] = _scores[sizeM1 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqShallow::gatherLinks(float* collector, char sense) {
  long sizeM1 = _size - 1;
  long sizeM2 = _size - 2;
  switch (sense){
  case '+': for (long n = 0; n < sizeM1; ++n){ collector[n] = _links[n]; } break;
  case '-': for (long n = 0; n < sizeM1; ++n){ collector[n] = _links[sizeM2 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}


void ScoredSeqShallow::gatherSubseq(char* collector, long position, long length, char sense){
  //OK();
  collector[length] = '\0';
  char* seqOs;
  switch (sense){
  case '+':
    seqOs = &_seq[position];
    for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; }
    break;
  case '-':
    seqOs = &_rcSeq[position];
    for (long n = 0; n < length; ++n){ collector[n] = seqOs[n]; }
    break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqShallow::gatherSubScores(float* collector, long position, long length, char sense){
  long sizeM1mPos = length - 1 - position;
  float* plusOs = &_scores[position];
  switch (sense){
  case '+': for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _scores[sizeM1mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqShallow::gatherSubLinks(float* collector, long position, long length, char sense){
  long sizeM2mPos = length - 2 - position;
  float* plusOs = &_links[position];
  switch (sense){
  case '+': for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _links[sizeM2mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}



long ScoredSeqShallow::size() {
  //OK();
  return _size;
}

void ScoredSeqShallow::buffer(){}
void ScoredSeqShallow::unbuffer(){}
void ScoredSeqShallow::bottomBuffer(){}
void ScoredSeqShallow::deepDelete(){ delete this; }
void ScoredSeqShallow::deepUnbuffer(){} // unbuffer doesn't do anything, so no need to call it


bool ScoredSeqShallow::hasPairedEnd() {
  //OK();
  return _pe != 0;
}

ScoredSeq * ScoredSeqShallow::getPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}
ScoredSeq * ScoredSeqShallow::getTempPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}

ScoredSeq * ScoredSeqShallow::shallowCopy(){
  //OK();
  return new ScoredSeqShallow(this,'+');
}
ScoredSeq * ScoredSeqShallow::flipCopy(){
  //OK();
  return new ScoredSeqShallow(this,'-');
}


bool ScoredSeqShallow::isNested(){ return false; }


#endif

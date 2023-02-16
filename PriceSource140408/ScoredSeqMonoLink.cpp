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


#ifndef SCOREDSEQMONOLINK_CPP
#define SCOREDSEQMONOLINK_CPP

#include "ScoredSeqMonoLink.h"
#include <cstring>
#include <iostream>
using namespace::std;

//long ScoredSeqMonoLink::_readCount = 0;

void ScoredSeqMonoLink::OK() {
}

// OVERLOADED CONSTRUCTORS (this goes on for a while)

ScoredSeqMonoLink::~ScoredSeqMonoLink(){
  delete [] _seq;
  delete [] _rcSeq;
  delete [] _scores;
  //_readCount--;
}


void ScoredSeqMonoLink::constructorSeqHelper(string seq){
  _size = seq.size();
  _seq = new char [_size+1];
  strcpy (_seq, seq.c_str());
  _rcSeq = reverseComplement(_seq,_size);
}


ScoredSeqMonoLink::ScoredSeqMonoLink(string seq, vector<float>* scores) :
  _pe(0)
{
  constructorSeqHelper(seq);
  long seqSizeMin1 = _size - 1;
  _link = 1.0;
  _scores = new float[ seqSizeMin1 + 1 ];
  for (long n = 0; n < seqSizeMin1; ++n){
    _scores[n] = (*scores)[n];
  }
  _scores[ seqSizeMin1 ] = (*scores)[ seqSizeMin1 ];
}


ScoredSeqMonoLink::ScoredSeqMonoLink(string seq, vector<float>* scores, float link) :
  _pe(0),
  _link(link)
{
  constructorSeqHelper(seq);
  _scores = new float[ _size + 1 ];
  for (long n = 0; n < _size; ++n){ _scores[n] = (*scores)[n]; }
}

// STREAMLINED
ScoredSeqMonoLink::ScoredSeqMonoLink(string seq, float* scores) :
  _pe(0)
{
  constructorSeqHelper(seq);
  _link = 1.0;
  _scores = new float[ _size + 1 ];
  for (long n = 0; n < _size; ++n){ _scores[n] = scores[n]; }
}
ScoredSeqMonoLink::ScoredSeqMonoLink(string seq, float* scores, float link) :
  _pe(0),
  _link(link)
{
  constructorSeqHelper(seq);
  _scores = new float[ _size + 1 ];
  for (long n = 0; n < _size; ++n){ _scores[n] = scores[n]; }
}



// CHAR ARRAYS
ScoredSeqMonoLink::ScoredSeqMonoLink(char* seq, float* scores, long size) :
  _size(size),
  _pe(0)
{
  _seq = new char[_size + 1];
  _seq[_size] = '\0';
  _link = 1.0;
  _scores = new float[_size+1];
  for (long n = 0; n < _size; ++n){
    _scores[n] = scores[n];
    _seq[n] = seq[n];
  }
  _rcSeq = reverseComplement(_seq,_size);
}
ScoredSeqMonoLink::ScoredSeqMonoLink(char* seq, float* scores, float link, long size) :
  _size(size),
  _link(link),
  _pe(0)
{
  _seq = new char[_size + 1];
  _seq[_size] = '\0';
  _scores = new float[_size+1];
  for (long n = 0; n < _size; ++n){
    _scores[n] = scores[n];
    _seq[n] = seq[n];
  }
  _rcSeq = reverseComplement(_seq,_size);
}


// REP EXPOSURE
ScoredSeqMonoLink::ScoredSeqMonoLink(bool expRep, char* seq, float* scores, float link, long size) :
  _size(size),
  _link(link),
  _pe(0)
{
  if (expRep){
    _seq = seq;
    _scores = scores;
  } else {
    _seq = new char[_size + 1];
    _seq[_size] = '\0';
    _scores = new float[_size+1];
    for (long n = 0; n < _size; ++n){
      _scores[n] = scores[n];
      _seq[n] = seq[n];
    }
  }
  _rcSeq = reverseComplement(_seq,_size);
}





float ScoredSeqMonoLink::scoreAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _scores[position];
  case '-': return _scores[_size - position - 1];
  default: throw AssemblyException::ArgError("sense must be (+) or (-)");
  }
}
float ScoredSeqMonoLink::scoreAtPosPlus(long position) { return _scores[position]; }
float ScoredSeqMonoLink::scoreAtPosMinus(long position) { return _scores[_size - position - 1]; }

float ScoredSeqMonoLink::linkAfterPosition(long position, char sense) { return _link; }
float ScoredSeqMonoLink::linkAfterPosPlus(long position){ return _link; }
float ScoredSeqMonoLink::linkAfterPosMinus(long position){ return _link; }


char ScoredSeqMonoLink::nucAtPosition(long position, char sense) {
  //OK();
  switch (sense){
  case '+': return _seq[ position ];
  case '-': return _rcSeq[ position ];
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
char ScoredSeqMonoLink::nucAtPosPlus(long position) { return _seq[ position ]; }
char ScoredSeqMonoLink::nucAtPosMinus(long position) { return _rcSeq[ position ]; }

char* ScoredSeqMonoLink::getSeq(char sense) {
  //OK();
  char* seq = new char[_size + 1];
  gatherSeq(seq,sense);
  return seq;
}
float* ScoredSeqMonoLink::getScores(char sense) {
  float* scores = new float[_size + 1];
  gatherScores(scores,sense);
  return scores;
}
float* ScoredSeqMonoLink::getLinks(char sense) {
  float* links = new float[_size + 1];
  gatherLinks(links,sense);
  return links;
}


char* ScoredSeqMonoLink::getSubseq(long position, long length, char sense){
  char* subseq = new char[ length + 1 ];
  gatherSubseq(subseq,position,length,sense);
  return subseq;
}
float* ScoredSeqMonoLink::getSubScores(long position, long length, char sense){
  float* scores = new float[length + 1];
  gatherSubScores(scores,position,length,sense);
  return scores;
}
float* ScoredSeqMonoLink::getSubLinks(long position, long length, char sense){
  float* links = new float[length + 1];
  gatherSubLinks(links,position,length,sense);
  return links;
}



void ScoredSeqMonoLink::gatherSeq(char* collector, char sense) {
  //OK();
  switch (sense){
  case '+': strcpy(collector, _seq); break;
  case '-': strcpy(collector, _rcSeq); break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqMonoLink::gatherScores(float* collector, char sense) {
  long sizeM1 = _size - 1;
  switch (sense){
  case '+': for (long n = 0; n < _size; ++n){ collector[n] = _scores[n]; } break;
  case '-': for (long n = 0; n < _size; ++n){ collector[n] = _scores[sizeM1 - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqMonoLink::gatherLinks(float* collector, char sense) {
  long sizeM1 = _size - 1;
  for (long n = 0; n < sizeM1; ++n){ collector[n] = _link; }
}


void ScoredSeqMonoLink::gatherSubseq(char* collector, long position, long length, char sense){
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
void ScoredSeqMonoLink::gatherSubScores(float* collector, long position, long length, char sense){
  long sizeM1mPos = length - 1 - position;
  float* plusOs = &_scores[position];
  switch (sense){
  case '+': for (long n = 0; n < length; ++n){ collector[n] = plusOs[n]; } break;
  case '-': for (long n = 0; n < length; ++n){ collector[n] = _scores[sizeM1mPos - n]; } break;
  default: throw AssemblyException::ArgError("sense must be '+' or '-'");
  }
}
void ScoredSeqMonoLink::gatherSubLinks(float* collector, long position, long length, char sense){
  for (long n = 0; n < length; ++n){ collector[n] = _link; }
}



long ScoredSeqMonoLink::size() {
  //OK();
  return _size;
}

void ScoredSeqMonoLink::buffer(){}
void ScoredSeqMonoLink::unbuffer(){}
void ScoredSeqMonoLink::bottomBuffer(){}
void ScoredSeqMonoLink::deepDelete(){ delete this; }
void ScoredSeqMonoLink::deepUnbuffer(){} // unbuffer doesn't do anything, so no need to call it


bool ScoredSeqMonoLink::hasPairedEnd() {
  //OK();
  return _pe != 0;
}

ScoredSeq * ScoredSeqMonoLink::getPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}
ScoredSeq * ScoredSeqMonoLink::getTempPairedEnd() {
  if (! hasPairedEnd()) {
    throw AssemblyException::CallingError("This ScoredSeq doesn't have a paired end.");
  }
  return _pe;
}

ScoredSeq * ScoredSeqMonoLink::shallowCopy(){
  //OK();
  return new ScoredSeqMonoLink(_seq, _scores, _link, _size);
}
ScoredSeq * ScoredSeqMonoLink::flipCopy(){
  char* seq = getSeq('-');
  float* scores = getScores('-');
  return new ScoredSeqMonoLink(true, seq, scores, _link, _size);
}

bool ScoredSeqMonoLink::isNested(){ return false; }


#endif

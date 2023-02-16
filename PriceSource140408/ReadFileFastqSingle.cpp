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


#ifndef READFILEFASTQSINGLE_CPP
#define READFILEFASTQSINGLE_CPP


#include "ReadFileFastqSingle.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <algorithm>
#include "ScoredSeqWithCarriers.h"
#include <vector>
#include <limits.h>
#include <cstring>
#include <cmath>
#include "fileUtilities.h"

using namespace::std;
using namespace::fileUtilities;

long ReadFileFastqSingle::_maxBlockSize = LONG_MAX;
long ReadFileFastqSingle::_posBlockSizeInc = 500;



ReadFileFastqSingle::ReadFileFastqSingle(){}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, FileType encoding) :
  _splitMode(false),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( 1.0 ),
  _invert(false) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, FileType encoding, float countFactor) :
  _splitMode(false),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( countFactor ),
  _invert(false) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}

ReadFileFastqSingle::ReadFileFastqSingle(bool invert, string filename, FileType encoding, float countFactor) :
  _splitMode(false),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( countFactor ),
  _invert(invert) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, FileType encoding, long readStart, long readLength) :
  _splitMode(false),
  _readStart( readStart ),
  _readLength( readLength ),
  _countFactor( 1.0 ),
  _invert(false) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}

ReadFileFastqSingle::ReadFileFastqSingle(string filename, FileType encoding, long readStart, long readLength, float countFactor) :
  _splitMode(false),
  _readStart( readStart ),
  _readLength( readLength ),
  _countFactor( countFactor),
  _invert(false) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}


ReadFileFastqSingle::ReadFileFastqSingle(bool invert, string filename, FileType encoding, long readStart, long readLength, float countFactor) :
  _splitMode(false),
  _readStart( readStart ),
  _readLength( readLength ),
  _countFactor( countFactor),
  _invert(invert) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}


ReadFileFastqSingle::ReadFileFastqSingle(bool invert, long readLength, string filename, FileType encoding, float countFactor) :
  _splitMode(true),
  _maxReadSize( readLength ),
  _readStart( 0 ),
  _readLength( 0 ),
  _countFactor( countFactor ),
  _invert(invert) {
  _encodingConverter = getScoreConverter(encoding, _countFactor, filename);
   constructorHelper(filename,encoding);
}




void ReadFileFastqSingle::constructorHelper(string filename, FileType fileType){
  if (fileType != FASTQFILE and fileType != ILLUMINAFILE){
    throw AssemblyException::ArgError("illegal filetype used for RFFqS constructor");
  }
  getOkToSkipChars(&_okToSkip, fileType);

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::FileError("RFFQ fastq file does not exist."); }

  if ( _readLength < 0 ){
    throw AssemblyException::ArgError("read length cannot be less than zero.");
  }

  _numReads = 0; // this initial value will be replaced if the true num immediately
  _numReadsDetermined = false;

  _initialBlock = 0;
  _initialPos = 0;

  _posBlockSize = _posBlockSizeInc;
  _posBlocks = new long[ _posBlockSize ];
  _posBlocks[0] = 0; // the first block starts with the zero coord
  _numPosBlocks = 1;

  _isOpen = false;
}





// PRIVATE CONTSTRUCTOR
ReadFileFastqSingle::ReadFileFastqSingle(ReadFileFastqSingle* precursor, long initialBlock, long initialPos, long numReads, bool numReadsDetermined) :
  _initialBlock( initialBlock ),
  _initialPos( initialPos ),
  _numReads( numReads ),
  _numReadsDetermined( numReadsDetermined )
 {
   getOkToSkipChars(&_okToSkip, FASTQFILE);
   _readStart = precursor->_readStart;
   _readLength = precursor->_readLength;
   _countFactor = precursor->_countFactor;
   _invert = precursor->_invert;
   _splitMode = precursor->_splitMode;
   _encodingConverter = getScoreConverter(precursor->_encodingConverter->getEncoding(), _countFactor, precursor->getName());
   if (_splitMode){ _maxReadSize = precursor->_maxReadSize; }

  _hasReadWaiting = false;

  // run checks on filename
  string filename = precursor->getName();
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  _posBlockSize = precursor->_posBlockSize;
  _posBlocks = new long[ _posBlockSize ];
  _numPosBlocks = precursor->_numPosBlocks;
  for (long n = 0; n < _numPosBlocks; ++n){ _posBlocks[n] = precursor->_posBlocks[n]; }
  _isOpen = false;
}




ReadFileFastqSingle::~ReadFileFastqSingle(){
  if (_isOpen){ close(); }
  delete [] _posBlocks;
  delete [] _filename;
  delete _encodingConverter;
}





void ReadFileFastqSingle::splitFile(set<ReadFile*>* putFilesHere, long numFiles){
  // determine the number of reads to include in each file
  long readsPerFile[numFiles];
  long currentTally = 0;
  long totalReads = numReads();
  for (long fileNum = 0; fileNum < numFiles; ++fileNum){
    long priorTally = currentTally;
    currentTally = (fileNum + 1) * totalReads / numFiles;
    readsPerFile[fileNum] = currentTally - priorTally;
  }
  // now create the files
  open();
  for (long fileNum = 0; fileNum < numFiles; ++fileNum){
    if (readsPerFile[fileNum] > 0){
      ReadFile* subFile = new ReadFileFastqSingle(this, _currentBlock, _currentPos, readsPerFile[fileNum], true);
      putFilesHere->insert( subFile );
      skipReads( readsPerFile[fileNum] );
    }
  }
  close();
}


ReadFile* ReadFileFastqSingle::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  open();
  skipReads( numReadsSkip );
  ReadFile* subFile = new ReadFileFastqSingle(this, _currentBlock, _currentPos, numReadsKeep, true);
  close();
  return subFile;
}

ReadFileFastqSingle* ReadFileFastqSingle::copy(){
  ReadFileFastqSingle* subFile = new ReadFileFastqSingle(this, _initialBlock, _initialPos, _numReads, _numReadsDetermined);
  return subFile;
}



string ReadFileFastqSingle::getName(){ return string(_filename); }


long ReadFileFastqSingle::numReads(){
  // if the count is zero, then there is no guarantee that a count has taken place;
  // if zero is the true count, then re-counting will take little time to finish.
  #pragma omp critical (RFFQS)
  {
    if (! _numReadsDetermined){
      open();
      while ( hasRead() ){
	skipRead();
	_numReads++;
      }
      close();
      _numReadsDetermined = true;
    }
  }
  return _numReads;
}


long ReadFileFastqSingle::ampliconSize(){
  throw AssemblyException::CallingError("A single-read fastq file does not have an amplicon size by definition.");
}

// open the file
void ReadFileFastqSingle::open(){
  if (_isOpen){
    throw AssemblyException::CallingError("don't open RFFSingle, it is already open.");
  }
  if (_getReadsFile.is_open()){
    throw AssemblyException::CallingError("don't open RFFSingle, it is already open.");
  }

  _getReadsFile.open(_filename,ifstream::in);
  _openWithCarriers = false;
  _openNormalized = false;
  _hasReadWaiting = false;
  _currentBlock = 0;
  _currentPos = 0;
  _numReadsRead = 0;
  // move to the specified initial position
  while (_currentBlock < _initialBlock){
    _currentBlock++;
    _getReadsFile.seekg( _posBlocks[_currentBlock], ios::cur );
  }
  _getReadsFile.seekg( _initialPos, ios::cur );
  _currentPos = _initialPos;
  _isOpen = true;
}
// PE reads are not returned by this file type, so the bool is meaningless
// and the method is just forwarded.
void ReadFileFastqSingle::open(bool getBothEnds){ open(); }

void ReadFileFastqSingle::openWithCarriers(long numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFileFastqSingle::openWithCarriers(long numSets, bool getBothEnds){ openWithCarriers(numSets); }

void ReadFileFastqSingle::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFileFastqSingle::openNormalized(bool getBothEnds){ openNormalized(); }

void ReadFileFastqSingle::openNormWithCarriers(long numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFileFastqSingle::openNormWithCarriers(long numSets, bool getBothEnds){ openWithCarriers(numSets); }


bool ReadFileFastqSingle::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFSingle while it is closed.");
  }
  if (! _hasReadWaiting){ readFromFileHelper(); }
  return _hasReadWaiting;
}

ScoredSeq* ReadFileFastqSingle::getRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't try to get reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){
    ScoredSeq* newRead = getReadHelper();
    if (newRead == NULL){
      throw AssemblyException::LogicError("cannot get read if there is no read to get.");
    }
    return newRead;
  } else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFileFastqSingle::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){ _hasReadWaiting = false; }
  else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFileFastqSingle::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  for (long n = 0; n < numToSkip; ++n){
    skipRead();
  }
}



void ReadFileFastqSingle::close(){
  // these haven't been accessed, so this is the only way to delete them
  if (_getReadsFile.is_open() ){
    _getReadsFile.close();
    _getReadsFile.clear();
  }
  while ( hasRead() ){ delete getRead(); }
  _hasReadWaiting = false;
  _isOpen = false;
}



void ReadFileFastqSingle::readFromFileHelper(){
  if ( _getReadsFile.is_open() ){
    long localPos = 0; // all will be defined by this and then adjusted for the block pos
    while ( (! _getReadsFile.eof() ) and _getReadsFile.peek() != '@' ){
      string uselessLine;
      getline(_getReadsFile,uselessLine);
      localPos += uselessLine.size() + 1;
    }


    if ( (! _getReadsFile.eof() ) and _getReadsFile.peek() == '@' and ( (! _numReadsDetermined) or _numReadsRead < _numReads)  ){
      _numReadsRead++;

      long stepSize;

      string uselessLine; // just something that getline needs to write to; gets discarded.
      getline(_getReadsFile,uselessLine); // the first name line
      stepSize = uselessLine.size() + 1;
      _seqStart = localPos + stepSize;
      localPos = _seqStart;

      getline(_getReadsFile,uselessLine); // the sequence line
      _seqLen = uselessLine.size();
      localPos += _seqLen + 1;

      getline(_getReadsFile,uselessLine); // the second name line
      _scoreStart = localPos + uselessLine.size() + 1;
      localPos = _scoreStart;

      getline(_getReadsFile,uselessLine); // the quality score line
      _scoreLen = uselessLine.size();
      localPos += _scoreLen + 1;

      // check if the block has run out of space; if so
      // create or move on to a new block
      long blockSpaceLeft = _maxBlockSize - _currentPos;
      if (localPos >= blockSpaceLeft){
	_currentBlock++;
	if (_currentBlock == _posBlockSize){
	  _posBlockSize += _posBlockSizeInc;
	  long* newPB = new long[ _posBlockSize ];
	  for (long n = 0; n < _numPosBlocks; ++n){ newPB[n] = _posBlocks[n]; }
	  delete [] _posBlocks;
	  _posBlocks = newPB;
	}
	if (_currentBlock == _numPosBlocks){
	  _posBlocks[_currentBlock] = _currentPos;
	  _numPosBlocks++;
	}
	_currentPos = 0;
      }
      // adjust the positions to global block positions
      _seqStart += _currentPos;
      _scoreStart += _currentPos;
      _currentPos += localPos;
      _hasReadWaiting = true;

    } else {
      _getReadsFile.close();
      _hasReadWaiting = false;
    }
  }
}



ScoredSeq* ReadFileFastqSingle::getReadHelper(){
  if (! _hasReadWaiting){ readFromFileHelper(); }
  if ( _hasReadWaiting ){

    // create the reads and add them to the carrier
    if ( _seqLen != _scoreLen ){
      throw AssemblyException::LogicError("The sequence and score list are different lengths.");
    } else if ( _seqLen < _readStart + _readLength ){
      throw AssemblyException::LogicError("The read's specified dimensions extend off the edge of the available sequence.");
    }

    // create the read itself
    ReadFileIndex * newIndex;
    if (_readLength == 0){
      newIndex = new ReadFileIndex2start1len(_seqStart + _readStart, _scoreStart + _readStart, _seqLen - _readStart, _currentBlock);
    } else {
      newIndex = new ReadFileIndex2start1len(_seqStart + _readStart, _scoreStart + _readStart, _readLength, _currentBlock);
    }
    ScoredSeq* seq = ReadFileCommunicator::makeRead(this, newIndex);
    //ScoredSeq* seq = new Read(this, newIndex);

    // make the type adjustment if appropriate
    if (_openNormalized){ seq = new ScoredSeqNormalized( seq ); }
    if (_openWithCarriers){ seq = new ScoredSeqWithCarriers( seq, _numCarrierSets ); }
    _hasReadWaiting = false;

    return seq;

  } else { return 0; }
}


void ReadFileFastqSingle::bufferQueue(Read * r){
  #pragma omp critical (RFFQSbuffer)
  { _bufferSet.insert(r); }
}
void ReadFileFastqSingle::bufferUnqueue(Read * r){
  #pragma omp critical (RFFQSbuffer)
  { _bufferSet.erase(r); }
}


void ReadFileFastqSingle::bufferFill(){
  bufferFill(bufferNotThreaded);
}


void ReadFileFastqSingle::bufferFill(BufferThreadedness threadedness){

  if (threadedness == bufferThreaded){

      // get reads sorted by file position
      long numReads = _bufferSet.size();
      Read** readArray = ReadFileFastqSingle::bufferedReadsSortedByFilePosition();
      long bufSizeP1 = numReads + 1;

      // create structures for gathering data
      ReadFileIndex** rfiRefs = new ReadFileIndex*[bufSizeP1];
      long* bufferSizes = new long[bufSizeP1];
      long* bufferStarts = new long[bufSizeP1];
      long totalBuffer = 0;
      for (long n = 0; n < numReads; ++n){
	rfiRefs[n] = ReadFileCommunicator::getRfiRef(readArray[n]);
	bufferSizes[n] = rfiRefs[n]->maxBufferSize() + 1;
	bufferStarts[n] = totalBuffer;
	totalBuffer += bufferSizes[n];
      }
      char** rawSeqArray = new char*[ bufSizeP1 ];
      char** rawScoresArray = new char*[ bufSizeP1 ];
      char* rawSeqArrayMemory = new char[ totalBuffer + 1 ];
      char* rawScoresArrayMemory = new char[ totalBuffer + 1 ];
      for (long n = 0; n < numReads; ++n){
	rawSeqArray[n] = &rawSeqArrayMemory[bufferStarts[n]];
	rawScoresArray[n] = &rawScoresArrayMemory[bufferStarts[n]];
      }

      // gather the data
      gatherRawDataForBufferReads(numReads, rfiRefs, rawSeqArray, rawScoresArray);
      _bufferSet.clear();

      // interpret the data
      ScoredSeq** bufferArray = new ScoredSeq*[bufSizeP1];
      if ( _splitMode ){
      #pragma omp parallel for schedule(static)
	for (long seqN = 0; seqN < numReads; ++seqN){
	  bufferArray[seqN] = bufferFillInterpretSplitHelper(rawSeqArray[seqN], rawScoresArray[seqN], rfiRefs[seqN]);
	}
      } else {
      #pragma omp parallel for schedule(static)
	for (long seqN = 0; seqN < numReads; ++seqN){
	  bufferArray[seqN] = bufferFillInterpretHelper(rawSeqArray[seqN], rawScoresArray[seqN], rfiRefs[seqN]);
	}
      }
      for (long seqN = 0; seqN < numReads; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
      delete [] bufferArray;

      delete [] rawSeqArray;
      delete [] rawScoresArray;
      delete [] rawSeqArrayMemory;
      delete [] rawScoresArrayMemory;

      delete [] readArray;
      delete [] rfiRefs;
      delete [] bufferSizes;
      delete [] bufferStarts;

  } else {

    #pragma omp critical (RFFQSbuffer)
    {
      // get reads sorted by file position
      long numReads = _bufferSet.size();
      Read** readArray = ReadFileFastqSingle::bufferedReadsSortedByFilePosition();
      long bufSizeP1 = numReads + 1;

      // create structures for gathering data
      ReadFileIndex** rfiRefs = new ReadFileIndex*[bufSizeP1];
      long* bufferSizes = new long[bufSizeP1];
      long* bufferStarts = new long[bufSizeP1];
      long totalBuffer = 0;
      for (long n = 0; n < numReads; ++n){
	rfiRefs[n] = ReadFileCommunicator::getRfiRef(readArray[n]);
	bufferSizes[n] = rfiRefs[n]->maxBufferSize() + 1;
	bufferStarts[n] = totalBuffer;
	totalBuffer += bufferSizes[n];
      }
      char** rawSeqArray = new char*[ bufSizeP1 ];
      char** rawScoresArray = new char*[ bufSizeP1 ];
      char* rawSeqArrayMemory = new char[ totalBuffer + 1 ];
      char* rawScoresArrayMemory = new char[ totalBuffer + 1 ];
      for (long n = 0; n < numReads; ++n){
	rawSeqArray[n] = &rawSeqArrayMemory[bufferStarts[n]];
	rawScoresArray[n] = &rawScoresArrayMemory[bufferStarts[n]];
      }

      // gather the data
      gatherRawDataForBufferReads(numReads, rfiRefs, rawSeqArray, rawScoresArray);
      _bufferSet.clear();

      // interpret the data
      ScoredSeq** bufferArray = new ScoredSeq*[bufSizeP1];
      if ( _splitMode ){
	for (long seqN = 0; seqN < numReads; ++seqN){
	  bufferArray[seqN] = bufferFillInterpretSplitHelper(rawSeqArray[seqN], rawScoresArray[seqN], rfiRefs[seqN]);
	}
      } else {
	for (long seqN = 0; seqN < numReads; ++seqN){
	  bufferArray[seqN] = bufferFillInterpretHelper(rawSeqArray[seqN], rawScoresArray[seqN], rfiRefs[seqN]);
	}
      }
      for (long seqN = 0; seqN < numReads; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
      delete [] bufferArray;

      delete [] rawSeqArray;
      delete [] rawScoresArray;
      delete [] rawSeqArrayMemory;
      delete [] rawScoresArrayMemory;

      delete [] readArray;
      delete [] rfiRefs;
      delete [] bufferSizes;
      delete [] bufferStarts;
    }
  }
}


Read** ReadFileFastqSingle::bufferedReadsSortedByFilePosition(){
  long numReads = _bufferSet.size();
  long currentBlockP1 = _currentBlock + 1;
  vector<Read*> vectorByBlock[ currentBlockP1 ];
  for (set<Read*>::iterator buffIt = _bufferSet.begin(); buffIt != _bufferSet.end(); ++buffIt){
    Read* read = (*buffIt);
    long block = ReadFileCommunicator::getRfiRef(read)->blockNum();
    vectorByBlock[block].push_back(read);
  }

  Read** sortedReads = new Read*[numReads+1];
  long rN = 0;
  SortReadByRfiStart sorter;
  for (long block = 0; block <= _currentBlock; ++block){
    sort( vectorByBlock[block].begin(), vectorByBlock[block].end(), sorter );
    for (vector<Read*>::iterator readIt = vectorByBlock[block].begin(); readIt != vectorByBlock[block].end(); ++readIt){
      sortedReads[rN] = *readIt;
      ++rN;
    }
  }
  return sortedReads;
}


void ReadFileFastqSingle::gatherRawDataForBufferReads(long numReads, ReadFileIndex** rfiArray, char** readData, char** scoreData){
  BlockPosFileCarrier* bpc = new BlockPosFileCarrier(this);
  for (long n = 0; n < numReads; ++n){
    // get the raw sequence
    bpc->seekBlockPos(rfiArray[n]->blockNum(), rfiArray[n]->readStart());
    bpc->getString(readData[n], streamsize( rfiArray[n]->readSize() + 1 ) );
    // get the raw scores
    bpc->seekBlockPos(rfiArray[n]->blockNum(), rfiArray[n]->scoreStart());
    bpc->getString(scoreData[n], streamsize( rfiArray[n]->scoreSize() + 1 ) );
  }
  delete bpc;
}


ScoredSeq* ReadFileFastqSingle::bufferFillInterpretHelper(char* rawSeqArray, char* rawScoresArray, ReadFileIndex* rfi){

  ReadFileStringInterpreter* si = new ReadFileStringInterpreter(rawSeqArray, rfi->readSize(), rfi->readSize(), &_okToSkip, _invert);
  // use this count for read size; the other strings may include bogus chars;
  // this is especially true if the file has carriage returns.
  float* scores = _encodingConverter->cstringToScores(rawScoresArray, si->_seqLen, _invert, si->_seqString);
  // make the ScoredSeq
  ScoredSeq* bufSeq = ScoredSeq::repExposedSeq( si->_seqString, scores, _countFactor, si->_seqLen );
  delete si;
  return bufSeq;
}

ScoredSeq* ReadFileFastqSingle::bufferFillInterpretSplitHelper(char* rawSeqArray, char* rawScoresArray, ReadFileIndex* rfi){

  ReadFileStringInterpreter* si = new ReadFileStringInterpreter(rawSeqArray, rfi->readSize(), rfi->readSize(), &_okToSkip, _invert);

  // now adjsut the parameters to get only part of the sequence and adjusted scores
  long seqLen = si->_seqLen;
  if (seqLen > _maxReadSize){
    // null-char-terminate at the end of the legit seq
    si->_seqString[ _maxReadSize ] = '\0';
    seqLen = _maxReadSize;
  }

  // figure out where the scores would change to half their value and adjust as necessary
  long scoreChange = si->_seqLen - seqLen;
  if (scoreChange < 0){ scoreChange = 0; }
  // i will iterate to scoreChange below, so make sure it isn't bigger than seqLen
  if (scoreChange > seqLen){ scoreChange = seqLen; }

  // get the scores; the ones after scoreChange will still need to be halved
  float* scores = _encodingConverter->cstringToScores(rawScoresArray, seqLen, _invert, si->_seqString);

  ScoredSeq* bufSeq;
  // no scores need be cut in half, so the link can be provided as a mono score
  if (scoreChange >= seqLen){
    bufSeq = ScoredSeq::repExposedSeq( si->_seqString, scores, _countFactor, seqLen );
  } else {
    // links has an extra element just because that is faster below than excluding
    // the last element and treating it differently
    float* links = new float[ seqLen ];
    // first the full values...
    float halfCount = _countFactor / 2;
    for (long n = 0; n < scoreChange; ++n){ links[n] = _countFactor; }
    // ...now the half values
    for (long n = scoreChange; n < seqLen; ++n){
      scores[n] /= 2;
      links[n] = halfCount;
    }
    // make the ScoredSeqShallow
    bufSeq = ScoredSeq::repExposedSeq( si->_seqString, scores, links, seqLen );
  }

  delete si;
  return bufSeq;
}




bool ReadFileFastqSingle::SortReadByRfiStart::operator() (Read* readA, Read* readB){
  return ReadFileCommunicator::getRfiRef(readA)->readStart() > ReadFileCommunicator::getRfiRef(readB)->readStart();
}




void ReadFileFastqSingle::updateBuffer(){}
void ReadFileFastqSingle::emptyBuffer(){}

// set up some initial conditions
ReadFileFastqSingle::BlockPosFileCarrier::BlockPosFileCarrier(){}
ReadFileFastqSingle::BlockPosFileCarrier::BlockPosFileCarrier(ReadFileFastqSingle* host) :
  _host(host){
  if (_file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is already open."); }
  _file.open(host->_filename,ifstream::in);
  if (! _file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is still closed."); }
  _block = 0;
  _pos = 0;
}
ReadFileFastqSingle::BlockPosFileCarrier::~BlockPosFileCarrier(){
  _file.close();
  _file.clear();
}


void ReadFileFastqSingle::BlockPosFileCarrier::getString(char* cstringToFill, streamsize length){
  _file.get(cstringToFill, length );
  _pos += length - 1;
}

void ReadFileFastqSingle::BlockPosFileCarrier::seekBlockPos(long newBlock, long newPos){
  // increment up or down until the correct block is reached
  _file.clear();
  while (_block > newBlock){
    _file.seekg( 0 - _pos, ios::cur );
    _pos = _host->_posBlocks[_block];
    _block--;
  }
  while (_block < newBlock){
    _block++;
    _file.seekg( _host->_posBlocks[_block] - _pos, ios::cur );
    _pos = 0;
  }
  // go to the correct spot in the correct block
  _file.seekg( newPos - _pos, ios::cur );
  _pos = newPos;
}




// THE RFI CLASS
ReadFileFastqSingle::ReadFileIndex2start1len::ReadFileIndex2start1len(){}
ReadFileFastqSingle::ReadFileIndex2start1len::ReadFileIndex2start1len(long readStart, long scoreStart, long readSize, long blockNum) :
  _readStart(readStart),
  _readSize(readSize),
  _scoreStart(scoreStart),
  _blockNum(blockNum) {
}
ReadFileFastqSingle::ReadFileIndex2start1len::~ReadFileIndex2start1len(){}
long ReadFileFastqSingle::ReadFileIndex2start1len::readStart(){ return _readStart; }
long ReadFileFastqSingle::ReadFileIndex2start1len::readSize(){ return _readSize; }
long ReadFileFastqSingle::ReadFileIndex2start1len::scoreStart(){ return _scoreStart; }
long ReadFileFastqSingle::ReadFileIndex2start1len::scoreSize(){ return _readSize; }
long ReadFileFastqSingle::ReadFileIndex2start1len::linkStart(){
  throw AssemblyException::CallingError("linkStart is irrelevant to the RFI imp for RFFastqSingle.");
}
long ReadFileFastqSingle::ReadFileIndex2start1len::linkSize(){
  throw AssemblyException::CallingError("linkSize is irrelevant to the RFI imp for RFFastqSingle.");
}
long ReadFileFastqSingle::ReadFileIndex2start1len::blockNum(){ return _blockNum; }
long ReadFileFastqSingle::ReadFileIndex2start1len::maxBufferSize(){ return _readSize; }





#endif

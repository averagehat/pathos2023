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


#ifndef READFILEPRICEQSINGLE_CPP
#define READFILEPRICEQSINGLE_CPP

#include "ReadFilePriceqSingle.h"
#include "fileUtilities.h"

#include <iostream>
#include <fstream>
#include <cstddef>
#include <string>
#include <algorithm>
#include "ScoredSeqWithCarriers.h"
#include "OutputFilePriceq.h"
#include <vector>
#include <limits.h>
#include <cstring>
using namespace::std;
using namespace::fileUtilities;

long ReadFilePriceqSingle::_maxBlockSize = LONG_MAX;
long ReadFilePriceqSingle::_posBlockSizeInc = 100;


ReadFilePriceqSingle::ReadFilePriceqSingle(string filename) :
  _readStart( 0 ),
  _readLength( 0 ),
  _invert(false) {
  constructorHelper(filename);
}
ReadFilePriceqSingle::ReadFilePriceqSingle(bool invert, string filename) :
  _readStart( 0 ),
  _readLength( 0 ),
  _invert(invert) {
  constructorHelper(filename);
}

ReadFilePriceqSingle::ReadFilePriceqSingle(string filename, long readStart, long readLength) :
  _readStart( readStart ),
  _readLength( readLength ),
  _invert(false) {
  constructorHelper(filename);
}
ReadFilePriceqSingle::ReadFilePriceqSingle(bool invert, string filename, long readStart, long readLength) :
  _readStart( readStart ),
  _readLength( readLength ),
  _invert(invert) {
  constructorHelper(filename);
}


void ReadFilePriceqSingle::constructorHelper(string filename){
  getOkToSkipChars(&_okToSkip, PRICEQFILE);

  _hasReadWaiting = false;

  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());

  // test that the file exists
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  if(inp.fail()){ cout << _filename << endl; throw AssemblyException::ArgError("RFPQ priceq file does not exist."); }

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
ReadFilePriceqSingle::ReadFilePriceqSingle(ReadFilePriceqSingle* precursor, long initialBlock, long initialPos, long numReads, bool numReadsDetermined) :
  _initialBlock( initialBlock ),
  _initialPos( initialPos ),
  _numReads( numReads ),
  _numReadsDetermined( numReadsDetermined )
 {
   _readStart = precursor->_readStart;
   _readLength = precursor->_readLength;
   _invert = precursor->_invert;

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




ReadFilePriceqSingle::~ReadFilePriceqSingle(){
  if (_isOpen){ close(); }
  delete[] _posBlocks;
  delete[] _filename;
}





void ReadFilePriceqSingle::splitFile(set<ReadFile*>* putFilesHere, long numFiles){
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
      ReadFile* subFile = new ReadFilePriceqSingle(this, _currentBlock, _currentPos, readsPerFile[fileNum], true);
      putFilesHere->insert( subFile );
      skipReads( readsPerFile[fileNum] );
    }
  }
  close();
}


ReadFile* ReadFilePriceqSingle::subFile(long numReadsSkip, long numReadsKeep){
  if (numReadsSkip + numReadsKeep > numReads()){
    throw AssemblyException::LogicError("RFFSingle was asked for a subset file with more reads than it has.");
  }
  open();
  skipReads( numReadsSkip );
  ReadFile* subFile = new ReadFilePriceqSingle(this, _currentBlock, _currentPos, numReadsKeep, true);
  close();
  return subFile;
}

ReadFilePriceqSingle* ReadFilePriceqSingle::copy(){
  ReadFilePriceqSingle* subFile = new ReadFilePriceqSingle(this, _initialBlock, _initialPos, _numReads, _numReadsDetermined);
  return subFile;
}



string ReadFilePriceqSingle::getName(){ return string(_filename); }


long ReadFilePriceqSingle::numReads(){
  // if the count is zero, then there is no guarantee that a count has taken place;
  // if zero is the true count, then re-counting will take little time to finish.
  #pragma omp critical (RFPQS)
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


long ReadFilePriceqSingle::ampliconSize(){
  throw AssemblyException::CallingError("A single-read priceq file does not have an amplicon size by definition.");
}

// open the file
void ReadFilePriceqSingle::open(){
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
void ReadFilePriceqSingle::open(bool getBothEnds){ open(); }

void ReadFilePriceqSingle::openWithCarriers(long numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
}
void ReadFilePriceqSingle::openWithCarriers(long numSets, bool getBothEnds){ openWithCarriers(numSets); }

void ReadFilePriceqSingle::openNormalized(){
  open();
  _openNormalized = true;
}
void ReadFilePriceqSingle::openNormalized(bool getBothEnds){ openNormalized(); }

void ReadFilePriceqSingle::openNormWithCarriers(long numSets){
  open();
  _openWithCarriers = true;
  _numCarrierSets = numSets;
  _openNormalized = true;
}
void ReadFilePriceqSingle::openNormWithCarriers(long numSets, bool getBothEnds){ openWithCarriers(numSets); }


bool ReadFilePriceqSingle::hasRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't check for reads from RFFSingle while it is closed.");
  }
  if (! _hasReadWaiting){ readFromFileHelper(); }
  return _hasReadWaiting;
}

ScoredSeq* ReadFilePriceqSingle::getRead(){
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

void ReadFilePriceqSingle::skipRead(){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  if ( hasRead() ){ _hasReadWaiting = false; }
  else {
    throw AssemblyException::LogicError("cannot get read if there is no read to get.");
  }
}

void ReadFilePriceqSingle::skipReads(long numToSkip){
  if (! _isOpen){
    throw AssemblyException::CallingError("don't skip reads from RFFSingle while it is closed.");
  }
  for (long n = 0; n < numToSkip; ++n){
    skipRead();
  }
}



void ReadFilePriceqSingle::close(){
  // these haven't been accessed, so this is the only way to delete them
  if (_getReadsFile.is_open() ){
    _getReadsFile.close();
    _getReadsFile.clear();
  }
  while ( hasRead() ){ delete getRead(); }
  _hasReadWaiting = false;
  _isOpen = false;
}



void ReadFilePriceqSingle::readFromFileHelper(){
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

      getline(_getReadsFile,uselessLine); // the third name line
      _linkStart = localPos + uselessLine.size() + 1;
      localPos = _linkStart;

      getline(_getReadsFile,uselessLine); // the link score line
      _linkLen = uselessLine.size();
      localPos += _linkLen + 1;

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
      _linkStart += _currentPos;
      _currentPos += localPos;
      _hasReadWaiting = true;

    } else {
      _getReadsFile.close();
      _hasReadWaiting = false;
    }
  }
}



ScoredSeq* ReadFilePriceqSingle::getReadHelper(){
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
      newIndex = new ReadFileIndex3start1len(_seqStart + _readStart,
					     _scoreStart + _readStart,
					     _linkStart + _readStart,
					     _seqLen - _readStart,
					     _currentBlock);
    } else {
      newIndex = new ReadFileIndex3start1len(_seqStart + _readStart,
					     _scoreStart + _readStart,
					     _linkStart + _readStart,
					     _readLength,
					     _currentBlock);
    }
    //ScoredSeq* seq = new Read(this, newIndex);
    ScoredSeq* seq = ReadFileCommunicator::makeRead(this, newIndex);

    // make the type adjustment if appropriate
    if (_openNormalized){ seq = new ScoredSeqNormalized( seq ); }
    if (_openWithCarriers){ seq = new ScoredSeqWithCarriers( seq, _numCarrierSets ); }
    _hasReadWaiting = false;

    return seq;

  } else { return 0; }
}


void ReadFilePriceqSingle::bufferQueue(Read * r){
  #pragma omp critical (RFPQSbufferFill)
  { _bufferSet.insert(r); }
}
void ReadFilePriceqSingle::bufferUnqueue(Read * r){
  #pragma omp critical (RFPQSbufferFill)
  { _bufferSet.erase(r); }
}


void ReadFilePriceqSingle::bufferFill(){
  bufferFill(bufferNotThreaded);
}



void ReadFilePriceqSingle::bufferFill(BufferThreadedness threadedness){

  if (threadedness == bufferThreaded){

      // get reads sorted by file position
      long numReads = _bufferSet.size();
      Read** readArray = ReadFilePriceqSingle::bufferedReadsSortedByFilePosition();
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
      char** rawLinksArray = new char*[ bufSizeP1 ];
      char* rawSeqArrayMemory = new char[ totalBuffer + 1 ];
      char* rawScoresArrayMemory = new char[ totalBuffer + 1 ];
      char* rawLinksArrayMemory = new char[ totalBuffer + 1 ];
      for (long n = 0; n < numReads; ++n){
	rawSeqArray[n] = &rawSeqArrayMemory[bufferStarts[n]];
	rawScoresArray[n] = &rawScoresArrayMemory[bufferStarts[n]];
	rawLinksArray[n] = &rawLinksArrayMemory[bufferStarts[n]];
      }

      // gather the data
      gatherRawDataForBufferReads(numReads, rfiRefs, rawSeqArray, rawScoresArray, rawLinksArray);
      _bufferSet.clear();

      // interpret the data
      ScoredSeq** bufferArray = new ScoredSeq*[bufSizeP1];
      #pragma omp parallel for schedule(static)
      for (long seqN = 0; seqN < numReads; ++seqN){
	bufferArray[seqN] = bufferFillInterpretHelper(rawSeqArray[seqN], rawScoresArray[seqN], rawLinksArray[seqN], rfiRefs[seqN]);
      }
      for (long seqN = 0; seqN < numReads; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
      delete [] bufferArray;

      delete [] rawSeqArray;
      delete [] rawScoresArray;
      delete [] rawLinksArray;
      delete [] rawSeqArrayMemory;
      delete [] rawScoresArrayMemory;
      delete [] rawLinksArrayMemory;

      delete [] readArray;
      delete [] rfiRefs;
      delete [] bufferSizes;
      delete [] bufferStarts;

  } else {

    #pragma omp critical (RFFQSbuffer)
    {
      // get reads sorted by file position
      long numReads = _bufferSet.size();
      Read** readArray = ReadFilePriceqSingle::bufferedReadsSortedByFilePosition();
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
      char** rawLinksArray = new char*[ bufSizeP1 ];
      char* rawSeqArrayMemory = new char[ totalBuffer + 1 ];
      char* rawScoresArrayMemory = new char[ totalBuffer + 1 ];
      char* rawLinksArrayMemory = new char[ totalBuffer + 1 ];
      for (long n = 0; n < numReads; ++n){
	rawSeqArray[n] = &rawSeqArrayMemory[bufferStarts[n]];
	rawScoresArray[n] = &rawScoresArrayMemory[bufferStarts[n]];
	rawLinksArray[n] = &rawLinksArrayMemory[bufferStarts[n]];
      }

      // gather the data
      gatherRawDataForBufferReads(numReads, rfiRefs, rawSeqArray, rawScoresArray, rawLinksArray);
      _bufferSet.clear();

      // interpret the data
      ScoredSeq** bufferArray = new ScoredSeq*[bufSizeP1];
      for (long seqN = 0; seqN < numReads; ++seqN){
	bufferArray[seqN] = bufferFillInterpretHelper(rawSeqArray[seqN], rawScoresArray[seqN], rawLinksArray[seqN], rfiRefs[seqN]);
      }
      for (long seqN = 0; seqN < numReads; ++seqN){
	ReadFileCommunicator::addScoredSeqImp(readArray[seqN], bufferArray[seqN]);
      }
      delete [] bufferArray;

      delete [] rawSeqArray;
      delete [] rawScoresArray;
      delete [] rawLinksArray;
      delete [] rawSeqArrayMemory;
      delete [] rawScoresArrayMemory;
      delete [] rawLinksArrayMemory;

      delete [] readArray;
      delete [] rfiRefs;
      delete [] bufferSizes;
      delete [] bufferStarts;
    }
  }
}


bool ReadFilePriceqSingle::SortReadByRfiStart::operator() (Read* readA, Read* readB){
  return ReadFileCommunicator::getRfiRef(readA)->readStart() > ReadFileCommunicator::getRfiRef(readB)->readStart();
}






Read** ReadFilePriceqSingle::bufferedReadsSortedByFilePosition(){
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


void ReadFilePriceqSingle::gatherRawDataForBufferReads(long numReads, ReadFileIndex** rfiArray, char** readData, char** scoreData, char** linkData){
  BlockPosFileCarrier* bpc = new BlockPosFileCarrier(this);
  for (long n = 0; n < numReads; ++n){
    // get the raw sequence
    bpc->seekBlockPos(rfiArray[n]->blockNum(), rfiArray[n]->readStart());
    bpc->getString(readData[n], streamsize( rfiArray[n]->readSize() + 1 ) );
    // get the raw scores
    bpc->seekBlockPos(rfiArray[n]->blockNum(), rfiArray[n]->scoreStart());
    bpc->getString(scoreData[n], streamsize( rfiArray[n]->scoreSize() + 1 ) );
    // get the raw links
    bpc->seekBlockPos(rfiArray[n]->blockNum(), rfiArray[n]->linkStart());
    bpc->getString(linkData[n], streamsize( rfiArray[n]->linkSize() + 1 ) );
  }
  delete bpc;
}


ScoredSeq* ReadFilePriceqSingle::bufferFillInterpretHelper(char* rawSeqArray, char* rawScoresArray, char* rawLinksArray, ReadFileIndex* rfi){

  ReadFileStringInterpreter* si = new ReadFileStringInterpreter(rawSeqArray, rfi->readSize(), rfi->readSize(), &_okToSkip, _invert);
  // use this count for read size; the other strings may include bogus chars;
  // this is especially true if the file has carriage returns.
  float* scores = cstringToScoresPriceq(rawScoresArray, si->_seqLen, _invert);
  float* links = cstringToScoresPriceq(rawLinksArray, si->_seqLen - 1, _invert);
  // make the ScoredSeqShallow
  ScoredSeq* bufSeq = ScoredSeq::repExposedSeq( si->_seqString, scores, links, si->_seqLen );
  delete si;
  return bufSeq;
}



void ReadFilePriceqSingle::updateBuffer(){}
void ReadFilePriceqSingle::emptyBuffer(){}

// set up some initial conditions
ReadFilePriceqSingle::BlockPosFileCarrier::BlockPosFileCarrier(){}
ReadFilePriceqSingle::BlockPosFileCarrier::BlockPosFileCarrier(ReadFilePriceqSingle* host) :
  _host(host){
  if (_file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is already open."); }
  _file.open(host->_filename,ifstream::in);
  if (! _file.is_open()) { throw AssemblyException::ArgError("RFFS::BPFC file is still closed."); }
  _block = 0;
  _pos = 0;
}
ReadFilePriceqSingle::BlockPosFileCarrier::~BlockPosFileCarrier(){
  _file.close();
  _file.clear();
}


void ReadFilePriceqSingle::BlockPosFileCarrier::getString(char* cstringToFill, streamsize length){
  _file.get(cstringToFill, length );
  _pos += length - 1;
}

void ReadFilePriceqSingle::BlockPosFileCarrier::seekBlockPos(long newBlock, long newPos){
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
ReadFilePriceqSingle::ReadFileIndex3start1len::ReadFileIndex3start1len(){}
ReadFilePriceqSingle::ReadFileIndex3start1len::ReadFileIndex3start1len(long readStart, long scoreStart, long linkStart, long readSize, long blockNum) :
  _readStart(readStart),
  _scoreStart(scoreStart),
  _linkStart(linkStart),
  _readSize(readSize),
  _blockNum(blockNum) {}
ReadFilePriceqSingle::ReadFileIndex3start1len::~ReadFileIndex3start1len(){}
long ReadFilePriceqSingle::ReadFileIndex3start1len::readStart(){ return _readStart; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::readSize(){ return _readSize; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::scoreStart(){ return _scoreStart; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::scoreSize(){ return _readSize; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::linkStart(){ return _linkStart; }
// there are one fewer links than nucleotides
long ReadFilePriceqSingle::ReadFileIndex3start1len::linkSize(){ return _readSize - 1; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::blockNum(){ return _blockNum; }
long ReadFilePriceqSingle::ReadFileIndex3start1len::maxBufferSize(){ return _readSize; }



#endif

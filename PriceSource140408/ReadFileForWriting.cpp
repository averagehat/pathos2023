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


#ifndef READFILEFORWRITING_CPP
#define READFILEFORWRITING_CPP

#include <iostream>
#include <string.h>
#include "ReadFileForWriting.h"
#include "AssemblyException.h"
#include "WritableSeqFasta.h"
#include "WritableSeqFastq.h"
#include "WritableSeqPriceq.h"
using namespace::std;


ReadFileForWriting* ReadFileForWriting::makeBasicReadFile(string filename){
  FileType fileType = fileUtilities::getFileType(filename);
  ReadFileForWriting* rf;
  if (fileType == FASTQFILE){ rf = new ReadFileForWritingFastq(filename,fileType); }
  else if (fileType == PRICEQFILE){ rf = new ReadFileForWritingPriceq(filename); }
  else if (fileType == FASTAFILE){ rf = new ReadFileForWritingFasta(filename); }
  else if (fileType == ILLUMINAFILE){ rf = new ReadFileForWritingFastq(filename,fileType); }
  else { throw AssemblyException::ArgError("read file is of an inappropriate type."); }
  return rf;
}


ReadFileForWriting::~ReadFileForWriting(){}



/////////////////////////////////////////////
// ALL OF THE FORMATS ARE BELOW!!!!!!!!!!!
// 1. FASTA
// 2. FASTQ
// 3. PRICEQ
/////////////////////////////////////////////


/////////////////////////////////////////////
// 1. FASTA
/////////////////////////////////////////////



ReadFileForWritingFasta::ReadFileForWritingFasta(string filename) :
  _isOpen(false),
  _readIsWaiting(false)
{
  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());
}
ReadFileForWritingFasta::~ReadFileForWritingFasta(){
  if (_isOpen){ close(); }
  delete [] _filename;
}


string ReadFileForWritingFasta::getName(){ return string(_filename); }

void ReadFileForWritingFasta::open(){
  if (_isOpen){ throw AssemblyException::CallingError("don't open RFFWFa when it is already open."); }
  _isOpen = true;
  _getReadsFile.open(_filename,ifstream::in);
  _getReadsFile.close();
  if (_getReadsFile.fail()){ throw AssemblyException::FileError("RFFWFa doesn't exist"); }
  _getReadsFile.open(_filename,ifstream::in);
  tryAndGetRead();
}

bool ReadFileForWritingFasta::hasRead(){
  if (! _isOpen){ throw AssemblyException::CallingError("RFFWFa::hasRead file is not open"); }
  return _readIsWaiting;
}

WritableSeq* ReadFileForWritingFasta::getWritableRead(){
  WritableSeq* lastRead = _waitingRead;
  tryAndGetRead();
  return lastRead;
}

void ReadFileForWritingFasta::close(){
  if (_readIsWaiting){
    delete _waitingRead;
    _readIsWaiting = false;
  }
  _getReadsFile.close();
  _getReadsFile.clear();
  _isOpen = false;
}

FileType ReadFileForWritingFasta::getFileType(){ return FASTAFILE; }

void ReadFileForWritingFasta::tryAndGetRead(){
  // skip any initial empty lines
  while( (! _getReadsFile.eof()) and _getReadsFile.peek() != '>' ){
    string tempLine; // must be defined in this scope to initialize as null
    getline(_getReadsFile,tempLine); // the sequence lines
  }
  // now GO!!!
  if (_getReadsFile.eof()){ _readIsWaiting = false; }
  else {
    if (_getReadsFile.peek() != '>'){
      cerr << _filename << endl;
      throw AssemblyException::FileError("improperly formatted fasta file with leading dead space");
    }
    string nameLine;
    getline(_getReadsFile,nameLine);
    vector<string> seqLines;
    while( (! _getReadsFile.eof()) and _getReadsFile.peek() != '>' ){
      string tempLine; // must be defined in this scope to initialize as null
      getline(_getReadsFile,tempLine); // the sequence lines
      seqLines.push_back(tempLine);
    }
    _waitingRead = new WritableSeqFasta(nameLine, &seqLines);
    _readIsWaiting = true;
  }
}




/////////////////////////////////////////////
// 2. FASTQ
/////////////////////////////////////////////

ReadFileForWritingFastq::ReadFileForWritingFastq(string filename, FileType filetype) :
  _isOpen(false),
  _readIsWaiting(false),
  _filetype(filetype)
{
  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());
}
ReadFileForWritingFastq::~ReadFileForWritingFastq(){
  if (_isOpen){ close(); }
  delete [] _filename;
}


string ReadFileForWritingFastq::getName(){ return string(_filename); }

void ReadFileForWritingFastq::open(){
  if (_isOpen){ throw AssemblyException::CallingError("don't open RFFWFq when it is already open."); }
  _isOpen = true;
  _getReadsFile.open(_filename,ifstream::in);
  _getReadsFile.close();
  if (_getReadsFile.fail()){ throw AssemblyException::FileError("RFFWFq doesn't exist"); }
  _getReadsFile.open(_filename,ifstream::in);
  tryAndGetRead();
}

bool ReadFileForWritingFastq::hasRead(){
  if (! _isOpen){ throw AssemblyException::CallingError("RFFWFq::hasRead file is not open"); }
  return _readIsWaiting;
}

WritableSeq* ReadFileForWritingFastq::getWritableRead(){
  WritableSeq* lastRead = _waitingRead;
  tryAndGetRead();
  return lastRead;
}

void ReadFileForWritingFastq::close(){
  if (_readIsWaiting){
    delete _waitingRead;
    _readIsWaiting = false;
  }
  _getReadsFile.close();
  _getReadsFile.clear();
  _isOpen = false;
}

FileType ReadFileForWritingFastq::getFileType(){ return _filetype; }

void ReadFileForWritingFastq::tryAndGetRead(){
  // skip any initial empty lines
  while( (! _getReadsFile.eof()) and _getReadsFile.peek() != '@' ){
    string tempLine; // must be defined in this scope to initialize as null
    getline(_getReadsFile,tempLine); // the sequence lines
  }
  // now GO!!!
  if (_getReadsFile.eof()){ _readIsWaiting = false; }
  else {
    if (_getReadsFile.peek() != '@'){
      cerr << _filename << endl;
      throw AssemblyException::FileError("improperly formatted fastq file with leading dead space");
    }
    string nameLine;
    string seqLine;
    string tossLine;
    string scoreLine;
    getline(_getReadsFile,nameLine);
    getline(_getReadsFile,seqLine);
    getline(_getReadsFile,tossLine);
    getline(_getReadsFile,scoreLine);
    _waitingRead = new WritableSeqFastq(nameLine, seqLine, scoreLine, _filetype, string(_filename));
    _readIsWaiting = true;
  }
}





/////////////////////////////////////////////
// 3. PRICEQ
/////////////////////////////////////////////



ReadFileForWritingPriceq::ReadFileForWritingPriceq(string filename) :
  _isOpen(false),
  _readIsWaiting(false)
{
  // run checks on filename
  _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());
}
ReadFileForWritingPriceq::~ReadFileForWritingPriceq(){
  if (_isOpen){ close(); }
  delete [] _filename;
}


string ReadFileForWritingPriceq::getName(){ return string(_filename); }

void ReadFileForWritingPriceq::open(){
  if (_isOpen){ throw AssemblyException::CallingError("don't open RFFWPq when it is already open."); }
  _isOpen = true;
  _getReadsFile.open(_filename,ifstream::in);
  _getReadsFile.close();
  if (_getReadsFile.fail()){ throw AssemblyException::FileError("RFFWPq doesn't exist"); }
  _getReadsFile.open(_filename,ifstream::in);
  tryAndGetRead();
}

bool ReadFileForWritingPriceq::hasRead(){
  if (! _isOpen){ throw AssemblyException::CallingError("RFFWPq::hasRead file is not open"); }
  return _readIsWaiting;
}

WritableSeq* ReadFileForWritingPriceq::getWritableRead(){
  WritableSeq* lastRead = _waitingRead;
  tryAndGetRead();
  return lastRead;
}

void ReadFileForWritingPriceq::close(){
  if (_readIsWaiting){
    delete _waitingRead;
    _readIsWaiting = false;
  }
  _getReadsFile.close();
  _getReadsFile.clear();
  _isOpen = false;
}

FileType ReadFileForWritingPriceq::getFileType(){ return PRICEQFILE; }

void ReadFileForWritingPriceq::tryAndGetRead(){
  // skip any initial empty lines
  while( (! _getReadsFile.eof()) and _getReadsFile.peek() != '@' ){
    string tempLine; // must be defined in this scope to initialize as null
    getline(_getReadsFile,tempLine); // the sequence lines
  }
  // now GO!!!
  if (_getReadsFile.eof()){ _readIsWaiting = false; }
  else {
    if (_getReadsFile.peek() != '@'){
      cerr << _filename << endl;
      throw AssemblyException::FileError("improperly formatted fastq file with leading dead space");
    }
    string nameLine;
    string seqLine;
    string tossLine;
    string scoreLine;
    string tossLine2;
    string linkLine;
    getline(_getReadsFile,nameLine);
    getline(_getReadsFile,seqLine);
    getline(_getReadsFile,tossLine);
    getline(_getReadsFile,scoreLine);
    getline(_getReadsFile,tossLine2);
    getline(_getReadsFile,linkLine);
    _waitingRead = new WritableSeqPriceq(nameLine, seqLine, scoreLine, linkLine, string(_filename));
    _readIsWaiting = true;
  }
}








#endif

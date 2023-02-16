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


#ifndef READFILE_CPP
#define READFILE_CPP

#include "ReadFile.h"
#include "fileUtilities.h"
#include "AssemblyException.h"
#include "ReadFileFastaSingle.h"
#include "ReadFileFastqSingle.h"
#include "ReadFilePriceqSingle.h"
#include "ReadFileDual.h"
#include "ReadFilePaired.h"

using namespace::fileUtilities;

ReadFile* ReadFile::makeBasicReadFile(string filename, float countFactor){
  FileType fileType = getFileType(filename);
  ReadFile* rf;
  if (fileType == FASTQFILE){ rf = new ReadFileFastqSingle(filename,fileType,countFactor); }
  else if (fileType == PRICEQFILE){ rf = new ReadFilePriceqSingle(filename); }
  else if (fileType == FASTAFILE){ rf = new ReadFileFastaSingle(filename,countFactor); }
  else if (fileType == ILLUMINAFILE){ rf = new ReadFileFastqSingle(filename,fileType,countFactor); }
  else { throw AssemblyException::ArgError("read file is of an inappropriate type."); }
  return rf;
}
ReadFile* ReadFile::makeInvertedReadFile(string filename, float countFactor){
  FileType fileType = getFileType(filename);
  ReadFile* rf;
  if (fileType == FASTQFILE){ rf = new ReadFileFastqSingle(true,filename,fileType,countFactor); }
  else if (fileType == PRICEQFILE){ rf = new ReadFilePriceqSingle(true,filename); }
  else if (fileType == FASTAFILE){ rf = new ReadFileFastaSingle(filename,countFactor,true); }
  else if (fileType == ILLUMINAFILE){ rf = new ReadFileFastqSingle(true,filename,fileType,countFactor); }
  else { throw AssemblyException::ArgError("read file is of an inappropriate type."); }
  return rf;
}


ReadFile* ReadFile::makePairedReadFile(string filename, long ampliconSize, float countFactor){
  ReadFile* rfh = makeBasicReadFile(filename, countFactor);
  ReadFile* rf = new ReadFilePaired(rfh, ampliconSize);
  delete rfh;
  return rf;
}
ReadFile* ReadFile::makePairedReadFile(string filenameA, string filenameB, long ampliconSize, float countFactor){
  ReadFile* rfhA = makeBasicReadFile(filenameA, countFactor);
  ReadFile* rfhB = makeBasicReadFile(filenameB, countFactor);
  ReadFile* rf = new ReadFileDual(rfhA, rfhB, ampliconSize);
  delete rfhA;
  delete rfhB;
  return rf;
}


ReadFile* ReadFile::makeMatePairFile(string filename, long ampliconSize, float countFactor){
  ReadFile* rfh = makeInvertedReadFile(filename, countFactor);
  ReadFile* rf = new ReadFilePaired(rfh, ampliconSize);
  delete rfh;
  return rf;
}
ReadFile* ReadFile::makeMatePairFile(string filenameA, string filenameB, long ampliconSize, float countFactor){
  ReadFile* rfhA = makeInvertedReadFile(filenameA, countFactor);
  ReadFile* rfhB = makeInvertedReadFile(filenameB, countFactor);
  ReadFile* rf = new ReadFileDual(rfhA, rfhB, ampliconSize);
  delete rfhA;
  delete rfhB;
  return rf;
}


ReadFile* ReadFile::makeFalsePairFile(string filename, long readSize, long ampliconSize, float countFactor){
  FileType fileType = getFileType(filename);
  ReadFile* rfA;
  ReadFile* rfB;
  if (fileType == FASTQFILE){
    rfA = new ReadFileFastqSingle(false, readSize, filename,fileType, countFactor);
    rfB = new ReadFileFastqSingle(true, readSize, filename,fileType, countFactor);
  } else if (fileType == FASTAFILE){
    rfA = new ReadFileFastaSingle(filename, countFactor, readSize, false);
    rfB = new ReadFileFastaSingle(filename, countFactor, readSize, true);
  } else if (fileType == ILLUMINAFILE){
    rfA = new ReadFileFastqSingle(false, readSize, filename, fileType, countFactor);
    rfB = new ReadFileFastqSingle(true, readSize, filename, fileType, countFactor);
  } else { throw AssemblyException::ArgError("read file is of an inappropriate type for a false-pair file."); }
  ReadFile* rf = new ReadFileDual(rfA, rfB, ampliconSize);
  delete rfA;
  delete rfB;
  return rf;
}

/*
ReadFile::FileType ReadFile::getFileType(string filename){
  size_t lastDot = filename.rfind('.');
  if (lastDot == string::npos){ return NONEFILE; }
  else {
    string append = filename.substr(lastDot);
    if (append==".fa" or append==".fna" or append==".ffn" or append==".frn" or append==".fasta"){ return FASTAFILE; }
    else if (append == ".fq" or append == ".fastq"){ return FASTQFILE; }
    else if (append == ".pq" or append == ".priceq"){ return PRICEQFILE; }
    else if (append == ".txt"){
      size_t underscore = filename.rfind('_');
      if (underscore == string::npos){ return NONEFILE; }
      else {
	string uAppend = filename.substr(underscore);
	if (uAppend == "_sequence.txt"){ return ILLUMINAFILE; }
	else { return NONEFILE; }
      }
    }
    else { return NONEFILE; }
  }
}
*/

ReadFile::~ReadFile(){}




#endif

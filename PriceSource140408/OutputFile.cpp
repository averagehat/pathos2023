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


#ifndef OUTPUTFILE_CPP
#define OUTPUTFILE_CPP

#include "OutputFile.h"
#include "AssemblyException.h"
#include "fileUtilities.h"
#include "OutputFileNull.h"
#include "OutputFileFasta.h"
#include "OutputFilePriceq.h"
#include <algorithm>
#include <fstream>
#include <stdio.h>

using namespace::std;
using namespace::fileUtilities;

OutputFile::~OutputFile(){};

OutputFile* OutputFile::makeOutputFile(string filename){
  FileType fileType = fileUtilities::getFileType(filename);
  if (fileType == PRICEQFILE){ return new OutputFilePriceq( filename ); }
  else if (fileType == FASTAFILE){ return new OutputFileFasta( filename ); }
  else {
    throw AssemblyException::ArgError("output filenames must be fasta or priceq");
  }
}
OutputFile* OutputFile::makeNullOutputFile(){ return new OutputFileNull(); }

ScoredSeq** OutputFile::makeSortedContigArray(set<ScoredSeq*>* contigs){
  vector<ScoredSeq*> sortedSeqs;
  sortedSeqs.insert( sortedSeqs.end(), contigs->begin(), contigs->end() );
  SortSsByLength ssSorter; // sort by seq length, longest to shortest
  sort( sortedSeqs.begin(), sortedSeqs.end(), ssSorter );
  ScoredSeq** seqArray = new ScoredSeq*[ contigs->size() + 1 ];
  long index = 0;
  for (vector<ScoredSeq*>::iterator it = sortedSeqs.begin(); it != sortedSeqs.end(); ++it){
    seqArray[index] = *it;
    index++;
  }
  seqArray[index] = NULL;
  return seqArray;  
}

bool OutputFile::SortSsByLength::operator() (ScoredSeq* ssnA, ScoredSeq* ssnB){
  return ssnA->size() > ssnB->size();
}

void OutputFile::ensureWritability(const char* filename){
  // if the file already exists, then it will be over-written; that is OK
  ifstream inp;
  inp.open(filename, ifstream::in);
  inp.close();
  if(inp.fail()){
    // now make sure that the file can be created by creating a temp null version
    ofstream testOutfile;
    testOutfile.open( filename );
    testOutfile << "test\n";
    testOutfile.close();
    inp.open(filename, ifstream::in);
    inp.close();
    if(inp.fail()){
      throw AssemblyException::FileError("Output file could not be written.");
    }
    // delete the file; it will be re-created later
    remove( filename );
  }
}

#endif

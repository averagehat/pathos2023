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

/* these are static methods and typedefs that are useful when dealing with files
*/


#ifndef FILEUTILITIES_H
#define FILEUTILITIES_H

#include <stdexcept>
#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <set>
#include <map>
#include "AssemblyException.h"
using namespace::std;

namespace fileUtilities {

  // returns a string with the file type based on the append of the filename
  // NONEFILE: no described filetype
  // FASTAFILE: .fa or .fasta
  // FASTQFILE: .fq or .fastq or _sequence.txt
  enum FileType { NONEFILE, FASTAFILE, FASTQFILE, PRICEQFILE, ILLUMINAFILE };
  FileType getFileType(string filename);
  bool hasAsciiScoreEncoding(FileType fileType);


  //enum Encoding{ fastqSystem=0, illuminaSystem=1 };
  //static Encoding getEncoding(fileUtilities::FileType fileType);

  /////// HELP PARSING FILES /////////
  void getOkToSkipChars(set<char>* okToSkip, FileType fileType);
  void getWhitespaceChars(set<char>* okToSkip, FileType fileType);




  // basic string interpretation
  class ReadFileStringInterpreter{
  public:
    ReadFileStringInterpreter(char* rawSeq, long rawLen, long maxLen, set<char>* okToSkip, bool revComp);
    ~ReadFileStringInterpreter();
    bool okToAddSeq();
    void addSeq(char* rawSeq, long rawLen, set<char>* okToSkip);
    // does NOT get deleted, just created (values are publicly accessible)
    char* _seqString;
    long _seqLen;
    // sequence can only be added if _revComp is false
    bool _revComp;
  private:
    void modSeq(char* rawSeq, long rawLen, set<char>* okToSkip, bool revComp);

  };



  /* this is public so that other classes can use the same methods to interpret fastq files;
     implementing classes are private */
  class ScoreConverter {
  public:
    virtual ~ScoreConverter();
    // positions in the seqString with close-to-zero scores will be converted to 'N' chars
    // I will cut it off at less-than-one-bit-of-information
    virtual float* cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString) = 0;
    virtual FileType getEncoding() = 0;
    virtual float getCountFactor() = 0;
  };

  ScoreConverter* getScoreConverter(FileType encoding, float countFactor, string nameForErrMsg);

  float* cstringToScoresPriceq(char* scoreCstring, long readSize, bool invert);

  // these are public so that I can test them and also so that I
  // can have the two methods in the same place and still be accessible
  // to other classes that may use them (like the input file class).
  char scoreToCharPriceq(float score);
  float charToScorePriceq(char scoreChar);



  // IMPLEMENTATIONS
  // SCORE CONVERSION ENCODING
  // these will contain accelerated methods for converting score chars into floats; 
  // i will pre-compute the values and use switch statements
  // Phred+33 encoding
  class ScoreConverterFastq : public ScoreConverter {
  public:
    ScoreConverterFastq(float countFactor, string nameForErr);
    ~ScoreConverterFastq();
    float* cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString);
    FileType getEncoding();
    float getCountFactor();
  private:
    float charToScore(char nuc, char* seqString, long seqPos);
    float _scores[100];
    string _nameForErr;
    float _countFactor;
  };
  // Phred+64 encoding
  class ScoreConverterIllumina : public ScoreConverter {
  public:
    ScoreConverterIllumina(float countFactor, string nameForErr);
    ~ScoreConverterIllumina();
    float* cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString);
    FileType getEncoding();
    float getCountFactor();
  private:
    float charToScore(char nuc, char* seqString, long seqPos);
    float _scores[100];
    string _nameForErr;
    float _countFactor;
  };




}


#endif

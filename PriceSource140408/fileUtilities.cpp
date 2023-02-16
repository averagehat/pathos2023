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

#ifndef FILEUTILITIES_CPP
#define FILEUTILITIES_CPP

#include "fileUtilities.h"

#include <stdexcept>
#include <string>
#include <map>
#include <stdio.h>
#include <iostream>
#include "AssemblyException.h"
using namespace::std;



fileUtilities::FileType fileUtilities::getFileType(string filename){
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



bool fileUtilities::hasAsciiScoreEncoding(FileType fileType){
  if (fileType==FASTQFILE){ return true; }
  else if (fileType==ILLUMINAFILE){ return true; }
  else if (fileType==PRICEQFILE){ return true; }
  else { return false; }
}




void fileUtilities::getOkToSkipChars(set<char>* okToSkip, FileType fileType){
  getWhitespaceChars(okToSkip, fileType);
  if (fileType == FASTAFILE){
    okToSkip->insert('-');
  }
}
void fileUtilities::getWhitespaceChars(set<char>* okToSkip, FileType fileType){
  okToSkip->insert('\f');
  okToSkip->insert('\r');
  okToSkip->insert('\n');
  okToSkip->insert('\t');
  okToSkip->insert('\v');
  okToSkip->insert('\0');
  okToSkip->insert(' ');
}






fileUtilities::ReadFileStringInterpreter::ReadFileStringInterpreter(char* rawSeq, long rawLen, long maxLen, set<char>* okToSkip, bool revComp) :
  _revComp(revComp),
  _seqLen(0)
{
  _seqString = new char[maxLen+1];
  modSeq(rawSeq, rawLen, okToSkip, revComp);
}
fileUtilities::ReadFileStringInterpreter::~ReadFileStringInterpreter(){}


bool fileUtilities::ReadFileStringInterpreter::okToAddSeq(){ return (! _revComp); }
void fileUtilities::ReadFileStringInterpreter::addSeq(char* rawSeq, long rawLen, set<char>* okToSkip){
  if ( okToAddSeq() ){
    modSeq(rawSeq, rawLen, okToSkip, false);
  } else {
    throw AssemblyException::CallingError("RF::SI::addSeq can only be called if the seq is not RC");
  }
}

void fileUtilities::ReadFileStringInterpreter::modSeq(char* rawSeq, long rawLen, set<char>* okToSkip, bool revComp){
  // initially this is zero
  long n2 = _seqLen;
  if (revComp){
    // remember: this IS the RC case
    long n = rawLen;
    while (n > 0){
      --n;
      switch(rawSeq[n]) {
      // allowed nucleotide characters
      case 'A': _seqString[n2] = 'T'; break;
      case 'C': _seqString[n2] = 'G'; break;
      case 'G': _seqString[n2] = 'C'; break;
      case 'T': _seqString[n2] = 'A'; break;
      case 'U': _seqString[n2] = 'A'; break;
      case 'R': _seqString[n2] = 'N'; break;
      case 'Y': _seqString[n2] = 'N'; break;
      case 'S': _seqString[n2] = 'N'; break;
      case 'W': _seqString[n2] = 'N'; break;
      case 'K': _seqString[n2] = 'N'; break;
      case 'M': _seqString[n2] = 'N'; break;
      case 'B': _seqString[n2] = 'N'; break;
      case 'D': _seqString[n2] = 'N'; break;
      case 'H': _seqString[n2] = 'N'; break;
      case 'V': _seqString[n2] = 'N'; break;
      case 'N': _seqString[n2] = 'N'; break;
      case 'a': _seqString[n2] = 'T'; break;
      case 'c': _seqString[n2] = 'G'; break;
      case 'g': _seqString[n2] = 'C'; break;
      case 't': _seqString[n2] = 'A'; break;
      case 'u': _seqString[n2] = 'A'; break;
      case 'r': _seqString[n2] = 'N'; break;
      case 'y': _seqString[n2] = 'N'; break;
      case 's': _seqString[n2] = 'N'; break;
      case 'w': _seqString[n2] = 'N'; break;
      case 'k': _seqString[n2] = 'N'; break;
      case 'm': _seqString[n2] = 'N'; break;
      case 'b': _seqString[n2] = 'N'; break;
      case 'd': _seqString[n2] = 'N'; break;
      case 'h': _seqString[n2] = 'N'; break;
      case 'v': _seqString[n2] = 'N'; break;
      case 'n': _seqString[n2] = 'N'; break;
      case '.': _seqString[n2] = 'N'; break;
      default:
	if (okToSkip->count(rawSeq[n]) == 0){
	  throw AssemblyException::ArgError("nucleotide seq contains invalid char");
	} else { n2--; }
      // no default; any other letter is skipped
      }
      n2++;
    }
  } else {
    // remember: this is NOT the RC case
    long refN = 0;
    long refLen = rawLen;
    char* rawSeqRef = &rawSeq[refN];
    char* ssRef = &_seqString[_seqLen];
    long n = 0;
    while (n < refLen){
      switch(rawSeqRef[n]) {
      // allowed nucleotide characters
      case 'A': ssRef[n] = 'A'; break;
      case 'C': ssRef[n] = 'C'; break;
      case 'G': ssRef[n] = 'G'; break;
      case 'T': ssRef[n] = 'T'; break;
      case 'U': ssRef[n] = 'T'; break;
      case 'R': ssRef[n] = 'N'; break;
      case 'Y': ssRef[n] = 'N'; break;
      case 'S': ssRef[n] = 'N'; break;
      case 'W': ssRef[n] = 'N'; break;
      case 'K': ssRef[n] = 'N'; break;
      case 'M': ssRef[n] = 'N'; break;
      case 'B': ssRef[n] = 'N'; break;
      case 'D': ssRef[n] = 'N'; break;
      case 'H': ssRef[n] = 'N'; break;
      case 'V': ssRef[n] = 'N'; break;
      case 'N': ssRef[n] = 'N'; break;
      case 'a': ssRef[n] = 'A'; break;
      case 'c': ssRef[n] = 'C'; break;
      case 'g': ssRef[n] = 'G'; break;
      case 't': ssRef[n] = 'T'; break;
      case 'u': ssRef[n] = 'T'; break;
      case 'r': ssRef[n] = 'N'; break;
      case 'y': ssRef[n] = 'N'; break;
      case 's': ssRef[n] = 'N'; break;
      case 'w': ssRef[n] = 'N'; break;
      case 'k': ssRef[n] = 'N'; break;
      case 'm': ssRef[n] = 'N'; break;
      case 'b': ssRef[n] = 'N'; break;
      case 'd': ssRef[n] = 'N'; break;
      case 'h': ssRef[n] = 'N'; break;
      case 'v': ssRef[n] = 'N'; break;
      case 'n': ssRef[n] = 'N'; break;
      case '.': ssRef[n] = 'N'; break;
      default:
	if (okToSkip->count(rawSeqRef[n]) == 0){
	  throw AssemblyException::ArgError("nucleotide seq contains invalid char");
	} else {
	  ++refN;
	  rawSeqRef = &rawSeq[refN];
	  --refLen;
	  --n;
	}
      }
      n++;
    }
    n2 = refLen + _seqLen;
  }
  _seqString[n2] = '\0';
  _seqLen = n2;
}





fileUtilities::ScoreConverter::~ScoreConverter(){}
fileUtilities::ScoreConverter* fileUtilities::getScoreConverter(FileType encoding, float countFactor, string nameForErr){
  if (encoding == ILLUMINAFILE){ return new ScoreConverterIllumina(countFactor,nameForErr); }
  else if (encoding == FASTQFILE){ return new ScoreConverterFastq(countFactor,nameForErr); }
  else { throw AssemblyException::ArgError("RFFQS::getScoreConverter, encoding enum type is unknown."); }
}






fileUtilities::ScoreConverterFastq::ScoreConverterFastq(float countFactor, string nameForErr) :
  _nameForErr(nameForErr),
  _countFactor(countFactor)
{
  // calclulate the scores ahead of time
  for (int n = 0; n < 94; ++n){
    _scores[n] = countFactor * (float(1) - pow( float(10), float(0 - n) / float(10) ) );
  }
}
fileUtilities::ScoreConverterFastq::~ScoreConverterFastq(){}
fileUtilities::FileType fileUtilities::ScoreConverterFastq::getEncoding(){ return FASTQFILE; }
float fileUtilities::ScoreConverterFastq::getCountFactor(){ return _countFactor; }

float* fileUtilities::ScoreConverterFastq::cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString){
  float* scores = new float[ readSize + 1 ];
  if (invert){
    long readSizeM1 = readSize - 1;
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[readSizeM1 - n], seqString, n); }
  } else {
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[n], seqString, n); }
  }
  return scores;
}
inline float fileUtilities::ScoreConverterFastq::charToScore(char nuc, char* seqString, long seqPos){
  switch(nuc){
  case '!': seqString[seqPos] = 'N'; return 0;
  case '"': seqString[seqPos] = 'N'; return 0;
  case '#': seqString[seqPos] = 'N'; return 0;
  case '$': return _scores[3];
  case '%': return _scores[4];
  case '&': return _scores[5];
  case '\'': return _scores[6];
  case '(': return _scores[7];
  case ')': return _scores[8];
  case '*': return _scores[9];
  case '+': return _scores[10];
  case ',': return _scores[11];
  case '-': return _scores[12];
  case '.': return _scores[13];
  case '/': return _scores[14];
  case '0': return _scores[15];
  case '1': return _scores[16];
  case '2': return _scores[17];
  case '3': return _scores[18];
  case '4': return _scores[19];
  case '5': return _scores[20];
  case '6': return _scores[21];
  case '7': return _scores[22];
  case '8': return _scores[23];
  case '9': return _scores[24];
  case ':': return _scores[25];
  case ';': return _scores[26];
  case '<': return _scores[27];
  case '=': return _scores[28];
  case '>': return _scores[29];
  case '?': return _scores[30];
  case '@': return _scores[31];
  case 'A': return _scores[32];
  case 'B': return _scores[33];
  case 'C': return _scores[34];
  case 'D': return _scores[35];
  case 'E': return _scores[36];
  case 'F': return _scores[37];
  case 'G': return _scores[38];
  case 'H': return _scores[39];
  case 'I': return _scores[40];
  case 'J': return _scores[41];
  case 'K': return _scores[42];
  case 'L': return _scores[43];
  case 'M': return _scores[44];
  case 'N': return _scores[45];
  case 'O': return _scores[46];
  case 'P': return _scores[47];
  case 'Q': return _scores[48];
  case 'R': return _scores[49];
  case 'S': return _scores[50];
  case 'T': return _scores[51];
  case 'U': return _scores[52];
  case 'V': return _scores[53];
  case 'W': return _scores[54];
  case 'X': return _scores[55];
  case 'Y': return _scores[56];
  case 'Z': return _scores[57];
  case '[': return _scores[58];
  case '\\': return _scores[59];
  case ']': return _scores[60];
  case '^': return _scores[61];
  case '_': return _scores[62];
  case '`': return _scores[63];
  case 'a': return _scores[64];
  case 'b': return _scores[65];
  case 'c': return _scores[66];
  case 'd': return _scores[67];
  case 'e': return _scores[68];
  case 'f': return _scores[69];
  case 'g': return _scores[70];
  case 'h': return _scores[71];
  case 'i': return _scores[72];
  case 'j': return _scores[73];
  case 'k': return _scores[74];
  case 'l': return _scores[75];
  case 'm': return _scores[76];
  case 'n': return _scores[77];
  case 'o': return _scores[78];
  case 'p': return _scores[79];
  case 'q': return _scores[80];
  case 'r': return _scores[81];
  case 's': return _scores[82];
  case 't': return _scores[83];
  case 'u': return _scores[84];
  case 'v': return _scores[85];
  case 'w': return _scores[86];
  case 'x': return _scores[87];
  case 'y': return _scores[88];
  case 'z': return _scores[89];
  case '{': return _scores[90];
  case '|': return _scores[91];
  case '}': return _scores[92];
  case '~': return _scores[93];
  default:
    cerr << "Bad score ascii char (\"" << nuc << "\") in Fastq file " << _nameForErr << endl;
    throw AssemblyException::FileError("Fastq-encoded file contained an illegal score character.");
  }
}




fileUtilities::ScoreConverterIllumina::ScoreConverterIllumina(float countFactor, string nameForErr) :
  _nameForErr(nameForErr),
  _countFactor(countFactor)
{
  // calclulate the scores ahead of time
  for (int n = 0; n < 63; ++n){
    _scores[n] = countFactor * (float(1) - pow( float(10), float(0 - n) / float(10) ) );
  }
}
fileUtilities::ScoreConverterIllumina::~ScoreConverterIllumina(){}
fileUtilities::FileType fileUtilities::ScoreConverterIllumina::getEncoding(){ return ILLUMINAFILE; }
float fileUtilities::ScoreConverterIllumina::getCountFactor(){ return _countFactor; }

float* fileUtilities::ScoreConverterIllumina::cstringToScores(char* scoreCstring, long readSize, bool invert, char* seqString){
  float* scores = new float[ readSize + 1 ];
  if (invert){
    long readSizeM1 = readSize - 1;
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[readSizeM1 - n], seqString, n); }
  } else {
    for (long n = 0; n < readSize; n++){ scores[n] = charToScore(scoreCstring[n], seqString, n); }
  }
  return scores;
}
inline float fileUtilities::ScoreConverterIllumina::charToScore(char nuc, char* seqString, long seqPos){
  switch(nuc){
  // old illumina specifications (pre-CASAVA v1.3) ;<=>?
  case ';': seqString[seqPos] = 'N'; return 0;
  case '<': seqString[seqPos] = 'N'; return 0;
  case '=': seqString[seqPos] = 'N'; return 0;
  case '>': seqString[seqPos] = 'N'; return 0;
  case '?': seqString[seqPos] = 'N'; return 0;
  // this is where the really legit scores start
  case '@': seqString[seqPos] = 'N'; return 0;
  case 'A': seqString[seqPos] = 'N'; return 0;
  case 'B': seqString[seqPos] = 'N'; return 0;
  case 'C': return _scores[3];
  case 'D': return _scores[4];
  case 'E': return _scores[5];
  case 'F': return _scores[6];
  case 'G': return _scores[7];
  case 'H': return _scores[8];
  case 'I': return _scores[9];
  case 'J': return _scores[10];
  case 'K': return _scores[11];
  case 'L': return _scores[12];
  case 'M': return _scores[13];
  case 'N': return _scores[14];
  case 'O': return _scores[15];
  case 'P': return _scores[16];
  case 'Q': return _scores[17];
  case 'R': return _scores[18];
  case 'S': return _scores[19];
  case 'T': return _scores[20];
  case 'U': return _scores[21];
  case 'V': return _scores[22];
  case 'W': return _scores[23];
  case 'X': return _scores[24];
  case 'Y': return _scores[25];
  case 'Z': return _scores[26];
  case '[': return _scores[27];
  case '\\': return _scores[28];
  case ']': return _scores[29];
  case '^': return _scores[30];
  case '_': return _scores[31];
  case '`': return _scores[32];
  case 'a': return _scores[33];
  case 'b': return _scores[34];
  case 'c': return _scores[35];
  case 'd': return _scores[36];
  case 'e': return _scores[37];
  case 'f': return _scores[38];
  case 'g': return _scores[39];
  case 'h': return _scores[40];
  case 'i': return _scores[41];
  case 'j': return _scores[42];
  case 'k': return _scores[43];
  case 'l': return _scores[44];
  case 'm': return _scores[45];
  case 'n': return _scores[46];
  case 'o': return _scores[47];
  case 'p': return _scores[48];
  case 'q': return _scores[49];
  case 'r': return _scores[50];
  case 's': return _scores[51];
  case 't': return _scores[52];
  case 'u': return _scores[53];
  case 'v': return _scores[54];
  case 'w': return _scores[55];
  case 'x': return _scores[56];
  case 'y': return _scores[57];
  case 'z': return _scores[58];
  case '{': return _scores[59];
  case '|': return _scores[60];
  case '}': return _scores[61];
  case '~': return _scores[62];
  default:
    cerr << "Bad score ascii char (\"" << nuc << "\") in Solexa/Illumina file " << _nameForErr << endl;
    throw AssemblyException::FileError("Illumina-encoded (Phred+64) file contained an illegal score character.");
  }
}




float* fileUtilities::cstringToScoresPriceq(char* scoreCstring, long readSize, bool invert){
  float* scores = new float[ readSize + 1 ];
  if (invert){
    for (long n = 0; n < readSize; n++){ scores[readSize - n - 1] = charToScorePriceq( scoreCstring[n] ); }
  } else {
    for (long n = 0; n < readSize; n++){ scores[n] = charToScorePriceq( scoreCstring[n] ); }
  }
  return scores;
}

char fileUtilities::scoreToCharPriceq(float score){
  // 126 is the max value; 126 - 35 = 91
  int ascInt;
  if (score >= 1.0){
    float logged = log(score) / log(float(1.5));
    if (logged > float(91)){ ascInt = 126; }
    else { ascInt = int(logged) + 35; }
  }
  else if (score < 0.3){ ascInt = 33; }
  else if (score < 0.9){ ascInt = 34; }
  else { ascInt = 35; }
  return char(ascInt);
}
float fileUtilities::charToScorePriceq(char scoreChar){
  int ascInt = int(scoreChar);
  switch( ascInt ){
  case 33: return float(0.0);
  case 34: return float(0.3);
  default: return pow(float(1.5), float(ascInt - 35));
  }
}




#endif

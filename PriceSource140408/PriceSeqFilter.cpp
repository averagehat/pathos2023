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

/* This is the text-based interface for the PRICE assembler.
 *
 */


#include "SeqFilter.h"

#include <set>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <queue>
#include <string.h>
#include <bitset>
#include <cstddef>
#include <omp.h>
#include <stdlib.h>
using namespace::std;
using ::size_t;

string versionString = string(
"PRICE Sequence Filter v1.2\n"
);
string helpString = string(
"\n"
"These are the command options for the PRICE sequence filter. \n"
"For more details, see the accompanying README.txt file. \n"
"To run the PRICE assembler, use the executable: PriceTI. \n"
"Contact: price@derisilab.ucsf.edu \n"
"\n"
"Usage: ./PriceSeqFilter [args] \n"
"\n"
"INPUT/OUTPUT FILES: \n"
" accepted formats are fasta (appended .fa, .fasta, .fna, .ffn, .frn), fastq (.fq, .fastq, or _sequence.txt), or priceq (.pq or .priceq) \n"
"  NOTE ABOUT FASTQ ENCODING: multiple encodings are currently used for fastq quality scores.  The traditional encoding is Phred+33,\n"
"                             and PRICE will interpret scores from any .fq or .fastq file according to that encoding.  The Phred+64 \n"
"                             encoding has been used extensively by Illumina, and so it is applied to Illumina's commonly-used _sequence.txt\n"
"                             file append.  Please make sure that your encoding matches your file append.\n"
"INPUT FILES: \n"
" -f <a>: (a) input file of non-paired sequences. \n"
" -fp <a> <b> (a,b) two input files of sequences; the sequences in one file are the paired-ends of those in the other.\n"
"OUTPUT FILES: \n"
" -o <a>: (a) output file of non-paired sequences. \n"
" -op <a> <b> (a,b) two output files of sequences; the sequences in one file are the paired-ends of those in the other.\n"
"\n"
"OTHER PARAMS: \n"
" -r <a>: (a) alignment score reward for a nucleotide match; should be a positive integer (default=1)\n"
" -q <a>: (a) alignment score penalty for a nucleotide mismatch; should be a negative integer (default=-2)\n"
" -G <a>: (a) alignment score penalty for opening a gap; should be a negative integer (default=-5)\n"
" -E <a>: (a) alignment score penalty for extending a gap; should be a negative integer (default=-2)\n"
"\n"
"FILTERING RULES: \n"
"-pair <a>:  <a> is \"both\" or \"either\", describing whether both reads of a pair or either read of a pair must FAIL \n"
"            the filters in order for the entire pair to be removed.  If paired files are used as input, the read pairs \n"
"            will be retained or discarded together.  If you want the sequences to be individually evaluated, run twice \n"
"            using each file as an individual (non-paired) input file. (default=either)\n"
"\n"
"FILTERING SEQUENCES: \n"
" -rqf <a> <b>: filters out sequences with an unaccptably high number of low-quality nucleotides, as defined by the provided \n"
"               quality scores (only applies to files whose formats include quality score information). \n"
"               (a) the percentage of nucleotides in a read that must be high-quality; (b) the minimum allowed probability\n"
"               of a nucleotide being correct (must be between 0 and 1, and will usually be a decimal value close to 1);\n"
" -rnf <a>: filters pairs of reads if either has an unaccptably high number of uncalled nucleotides (Ns or other ambiguous\n"
"               IUPAC codes).  Like -rqf, but will also filter fasta-format data.  (a) the percentage of nucleotides in a read\n"
"               that must be called\n"
" -maxHp <a>: filters out a pair of reads if either read has a homo-polymer track >(a) nucleotides in length.\n"
" -maxDi <a>: filters out a pair of reads if either read has a repeating di-nucleotide track >(a) nucleotides in length.\n"
"           NOTE: this will also catch mono-nucleotide repeats of the specified length (a string of A's is also a string\n"
"           of AA's), so calling -maxHp in addition to -maxDi is superfluous unless -maxHp is given a smaller max value.\n"
" -badf <a> <b>: sequences fail if they match with at least (b)% identity to a sequence in file (a).\n"
" -goodf <a> <b>: sequences fail if they DON'T match with at least (b)% identity to a sequence in file (a).\n"
" -lenf <a>: filters out sequences shorter than (a) nt \n"
"\n"
"COMPUTATIONAL EFFICIENCY: \n"
" -a <x>: (x)num threads to use (default=1) \n"
"\n"
"USER INTERFACE: \n"
" -log <a>: determines the type of output:\n"
"         (a) = c: concise stdout (default)\n"
"         (a) = n: no stdout \n"
" -, -h, or --help: user interface info. \n"
);

// These are helpers for checking arg qualities that will need to be reviewed over and over
bool isFileMissing(string filename){
  // test that the file exists
  char* _filename = new char[ filename.size()+1 ];
  strcpy (_filename, filename.c_str());
  ifstream inp;
  inp.open(_filename, ifstream::in);
  inp.close();
  delete [] _filename;
  return inp.fail();
}

// don't need to check for zero-length error, but I will anyway
bool worksAsInt(string arg){
  bool argOk = true;
  if (arg.size() == 0){ argOk = false; }
  string::iterator it = arg.begin();
  bool negOk = arg.size() > 1 and (*it)=='-';
  if (negOk){ ++it; }
  while (argOk and it != arg.end()){
    if (! isdigit(*it)){ argOk = false; }
    ++it;
  }
  return argOk;
}
bool worksAsFloat(string arg){
  bool argOk = true;
  if (arg.size() == 0){ argOk = false; }
  string::iterator it = arg.begin();
  bool negOk = arg.size() > 1 and (*it)=='-';
  if (negOk){ ++it; }
  while (argOk and it != arg.end()){
    if (! (isdigit(*it) or (*it)=='.') ){ argOk = false; }
    ++it;
  }
  return argOk;
}
bool worksAsNumber(string arg){
  if (worksAsInt(arg)){ return true; }
  else if(worksAsFloat(arg)){ return true; }
  else { return false; }
}


// argn is the current flag being examined;
// argv and argc are constants from the system
// call being passed down.
// the flag itself is NOT included in the return count
// REQUIRES that the args themselves do not start with dashes
// unless they are numbers
int getNumFlagArgs(int argn, char* argv[], int argc){
  int argCount = 0;
  ++argn;
  while (argn < argc and (argv[argn][0] != '-' or worksAsNumber(argv[argn]))){ ++argn; ++argCount; }
  return argCount;
}
// last item in int array is zero - args are 1-indexed (flag is zero)
bool shouldBeIntArg(int argnStart, int plusN, char* argv[]){
  if ( worksAsInt(argv[argnStart+plusN]) ){ return true; }
  else {
    cerr << "Arg #" << plusN << " for the " << argv[argnStart] << " flag should be an integer: " << argv[argnStart+plusN] << endl;
    return false;
  }
}
bool shouldBeFloatArg(int argnStart, int plusN, char* argv[]){
  if ( worksAsFloat(argv[argnStart+plusN]) ){ return true; }
  else {
    cerr << "Arg #" << plusN << " for the " << argv[argnStart] << " flag should be a number: " << argv[argnStart+plusN] << endl;
    return false;
  }
}



int main(int argc, char *argv[]){

  cout << versionString << endl;

  int initialInput = 0;
  bool willDoFiltering = true; // allows full set of errors to be printed
  bool neitherFilteringNorError = false; // allows some calls to cause no assembly but also no fail msg
  bool listenerIsNull = false;

  // default is 1 core
  omp_set_num_threads( 1 );


  // same behavior as if --help had been called
  if (argc < 2){
    cout << helpString << endl;
    willDoFiltering = false;
    neitherFilteringNorError = true;
  }


  SeqFilter* seqFilter = new SeqFilter();
  // default behavior
  seqFilter->setEitherFailMode();



  // now go through and get the command-line args organized
  int argn = 1;
  while ( argn < argc ){
    string arg = string(argv[argn]);
    int numFlagArgs = getNumFlagArgs(argn,argv,argc);

    if (arg == "-h" or arg == "-" or arg == "--help"){
      cout << helpString << endl;
      willDoFiltering = false;
      neitherFilteringNorError = true;
    }

    else if (arg == "-fp"){
      if (numFlagArgs != 2){
	cerr << "-fp was given " << numFlagArgs << " args instead of the required 2 (filenames)." << endl;
	willDoFiltering = false;
      }
      string fileA = argv[argn+1];
      string fileB = argv[argn+2];
      if (isFileMissing(fileA)){
	cerr << "The following file input with the -fp flag could not be found: " << fileA << endl;
	willDoFiltering = false;
      }
      if (isFileMissing(fileB)){
	cerr << "The following file input with the -fp flag could not be found: " << fileB << endl;
	willDoFiltering = false;
      }
      if (willDoFiltering){ seqFilter->addInputFiles(fileA,fileB); }


    } else if (arg == "-f"){
      if (numFlagArgs != 1){
	cerr << "-f was given " << numFlagArgs << " args instead of the required 1 (filename)." << endl;
	willDoFiltering = false;
      }
      string file = argv[argn+1];
      if (isFileMissing(file)){
	cerr << "The following file input with the -f flag could not be found: " << file << endl;
	willDoFiltering = false;
      }
      if (willDoFiltering){ seqFilter->addInputFile(file); }


    } else if (arg == "-op"){
      if (numFlagArgs != 2){
	cerr << "-fp was given " << numFlagArgs << " args instead of the required 2 (filenames)." << endl;
	willDoFiltering = false;
      }
      string fileA = argv[argn+1];
      string fileB = argv[argn+2];
      if (willDoFiltering){ seqFilter->addOutputFiles(fileA,fileB); }


      } else if (arg == "-o"){
        if (numFlagArgs != 1){
	  cerr << "-o was given " << numFlagArgs << " args instead of the required 1 (filename)." << endl;
	  willDoFiltering = false;
	}
	string file = argv[argn+1];
	if (willDoFiltering){ seqFilter->addOutputFile(file); }


      } else if (arg == "-r" or arg == "-q" or arg == "-G" or arg == "-E"){
	// ALIGNMENT SCORING PARAMETERS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoFiltering = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoFiltering = false; }
	  if (willDoFiltering){
	    if (arg == "-r"){ seqFilter->setNucMatchScore( atol(argv[argn+1]) ); }
	    else if (arg == "-q"){ seqFilter->setNucMismatchPenalty( atol(argv[argn+1]) ); }
	    else if (arg == "-G"){ seqFilter->setOpenGapPenalty( atol(argv[argn+1]) ); }
	    else if (arg == "-E"){ seqFilter->setExtendGapPenalty( atol(argv[argn+1]) ); }
	  }
	}


      } else if (arg == "-pair"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoFiltering = false;
	} else if (willDoFiltering){
	  char answerChar = argv[argn+1][0];
	  if (answerChar=='b' or answerChar=='B' or answerChar=='e' or answerChar=='E'){
	    if (willDoFiltering){
	      if (answerChar=='b' or answerChar=='B'){ seqFilter->setBothFailMode(); }
	      else { seqFilter->setEitherFailMode(); }
	    }
	  } else {
	    cerr << "The following input was invalid for the -pair flag: " << argv[argn+1] << endl;
	    willDoFiltering = false;
	  }
	}


      } else if (arg == "-badf"){
        if (numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2." << endl;
	  willDoFiltering = false;
	} else if (willDoFiltering){
	  string filename = argv[argn+1];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -badf flag could not be found: " << filename << endl;
	    willDoFiltering = false;
	  }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoFiltering = false; }
	  float fractId = atof( argv[argn+2] ) / float(100.0);
	  if (willDoFiltering){ seqFilter->addBadSequenceFilter(filename, fractId); }
	}


      } else if (arg == "-goodf"){
        if (numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2." << endl;
	  willDoFiltering = false;
	} else if (willDoFiltering){
	  string filename = argv[argn+1];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -badf flag could not be found: " << filename << endl;
	    willDoFiltering = false;
	  }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoFiltering = false; }
	  float fractId = atof( argv[argn+2] ) / float(100.0);
	  if (willDoFiltering){ seqFilter->addGoodSequenceFilter(filename, fractId); }
	}


      } else if (arg=="-maxHp" or arg=="-maxDi"){
	// COMPLEXITY READ FILTERS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoFiltering = false;
	} else if (willDoFiltering){
	  if (arg == "-maxHp"){ seqFilter->addHomopolymerFilter( atol(argv[argn+1]) ); }
	  else if (arg == "-maxDi"){ seqFilter->addDinucRepeatFilter( atol(argv[argn+1]) ); }
	}


      } else if (arg == "-rqf"){
        if (numFlagArgs != 2 and numFlagArgs != 4){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2 or 4." << endl;
	  willDoFiltering = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoFiltering = false; }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoFiltering = false; }
	  float minFractGood = atof(argv[argn+1]) / float(100);
	  float minProbCorrect = atof(argv[argn+2]);
	  if (numFlagArgs==2){
	    if (willDoFiltering){ seqFilter->addReadQualityFilter(minFractGood, minProbCorrect); }
	  } else {
	    if (! shouldBeIntArg(argn,3,argv)){ willDoFiltering = false; }
	    if (! shouldBeIntArg(argn,4,argv)){ willDoFiltering = false; }
	    int numSkipCycles = atoi(argv[argn+3]);
	    int numRunCycles = atoi(argv[argn+4]);
	    if (willDoFiltering){ seqFilter->addReadQualityFilter(minFractGood, minProbCorrect, numSkipCycles, numRunCycles); }
	  }
	}


      } else if (arg == "-rnf"){
        if (numFlagArgs != 1 and numFlagArgs != 3){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1 or 3." << endl;
	  willDoFiltering = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoFiltering = false; }
	  float minFractGood = atof(argv[argn+1]) / float(100);
	  if (numFlagArgs==1){
	    if (willDoFiltering){ seqFilter->addReadCalledBasesFilter(minFractGood); }
	  } else {
	    if (! shouldBeIntArg(argn,2,argv)){ willDoFiltering = false; }
	    if (! shouldBeIntArg(argn,3,argv)){ willDoFiltering = false; }
	    int numSkipCycles = atoi(argv[argn+2]);
	    int numRunCycles = atoi(argv[argn+3]);
	    if (willDoFiltering){ seqFilter->addReadCalledBasesFilter(minFractGood, numSkipCycles, numRunCycles); }
	  }
	}


      } else if (arg == "-lenf"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoFiltering = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoFiltering = false; }
	  long minLength = atol( argv[argn+1] );
	  if (willDoFiltering){ seqFilter->addLengthFilter(minLength); }
	}


      } else if (arg=="-a"){
	// COMP EFFICIENCY ARGS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoFiltering = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoFiltering = false; }
          if (willDoFiltering){ omp_set_num_threads( atoi( argv[argn+1] ) ); }
	}


      } else if (arg == "-log"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoFiltering = false;
	} else {
	  string type = string(argv[argn+1]);
	  if (type == "n"){ listenerIsNull = true; }
	  else if (type == "c"){ }
	  else {
	    cerr << "'" << type << "' is not a valid stdout type ('n', or 'c')." << endl;
	    willDoFiltering = false;
	  }
	}


      } else {
	cerr << "'" << arg << "' is an invalid flag that was given " << numFlagArgs << " args." << endl;
	willDoFiltering = false;
      }
      argn += numFlagArgs + 1;

      // a double-check, but there should be an error already raised
      if (argn != argc and argv[argn][0] != '-'){
	cerr << arg << " was given an incorrect number of args." << endl;
	willDoFiltering = false;
      }
    }




    // check again!
    if ( willDoFiltering ){

      // TEMPORARY set-up for returning the command that launched the job
      ostringstream launchCommand;
      launchCommand << "Launch command:";
      for (int n = 0; n < argc; ++n){ launchCommand << " " << argv[n]; }
      cout << launchCommand.str().c_str() << endl;

      if (listenerIsNull){ seqFilter->makeLogNull(); }

      seqFilter->runFilter();

    }
    delete seqFilter;


  if (willDoFiltering){ return 0; }
  else if (neitherFilteringNorError){
    cerr << "No sequence filter was run; help message printed instead." << endl;
    return 0;
  } else {
    cerr << "SEQUENCE FILTER JOB EXITED WITH ERROR!!!" << endl;
    return 1; 
  }
}




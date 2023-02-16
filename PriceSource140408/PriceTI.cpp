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


#include "Assembler.h"

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
using namespace::std;
using ::size_t;

string versionString = string(
"PRICE Assembler v1.2\n"
);
string helpString = string(
"\n"
"These are the command options for the PRICE assembler. \n"
"For more details, see the accompanying README.txt file. \n"
"Contact: price@derisilab.ucsf.edu \n"
"\n"
"Usage: ./PriceTI [args] \n"
"\n"
"INPUT FILES: \n"
" accepted formats are fasta (appended .fa, .fasta, .fna, .ffn, .frn), fastq (.fq, .fastq, or _sequence.txt), or priceq (.pq or .priceq) \n"
"  NOTE ABOUT FASTQ ENCODING: multiple encodings are currently used for fastq quality scores.  The traditional encoding is Phred+33,\n"
"                             and PRICE will interpret scores from any .fq or .fastq file according to that encoding.  The Phred+64 \n"
"                             encoding has been used extensively by Illumina, and so it is applied to Illumina's commonly-used _sequence.txt\n"
"                             file append.  Please make sure that your encoding matches your file append.\n"
"INPUT READ FILES: \n"
"  NOTE: these flags can be used multiple times in the same command to include multiple read datasets. \n"
"  (default % ID = 90%) \n"
" PAIRED-END FILES (reads are 3p of one another on opposite strands, i.e. pointing towards one another)\n"
" -fp a b c [d e [f]]: (a,b)input file pair, (c)amplicon insert size (including read) \n"
"                  (d,e,f) are optional; (d)the num. cycles to be skipped before this file is used;\n"
"                  if (f) is provided, then the file will alternate between being used for (e) cycles and not used for (f) cycles;\n"
"                  otherwise, the file will be used for (e) cycles then will not be used again. \n"
" -fpp a b c d [e f [g]]: (a,b)input file pair, (c)amplicon insert size (including read), (d)required % identity for match (25-100 allowed) \n"
"                  (e,f,g) are optional; (e)the num. cycles to be skipped before this file is used;\n"
"                  if (g) is provided, then the file will alternate between being used for (f) cycles and not used for (g) cycles;\n"
"                  otherwise, the file will be used for (f) cycles then will not be used again. \n"
" -fs a b [c d [e]]: (a)input paired-end file (alternating entries are paired-end reads), (b)amplicon insert size (including read) \n"
"                  (c,d,e) are optional; (c)the num. cycles to be skipped before this file is used;\n"
"                  if (e) is provided, then the file will alternate between being used for (d) cycles and not used for (e) cycles;\n"
"                  otherwise, the file will be used for (d) cycles then will not be used again. \n"
" -fsp a b c [d e [f]]: (a)input paired-end file (alternating entries are paired-end reads), (b)amplicon insert size (including read),\n"
"                  (c)required % identity for match (25-100 allowed) \n"
"                  (d,e,f) are optional; (d)the num. cycles to be skipped before this file is used;\n"
"                  if (f) is provided, then the file will alternate between being used for (e) cycles and not used for (f) cycles;\n"
"                  otherwise, the file will be used for (e) cycles then will not be used again. \n"
" MATE-PAIR FILES (reads are 5p of one another on opposite strands, i.e. pointing away from one another)\n"
" -mp a b c [d e [f]]: like -fp above, but with reads in the opposite orientation.\n"
" -mpp a b c d [e f [g]]: like -fpp above, but with reads in the opposite orientation.\n"
" -ms a b [c d [e]]: like -fs above, but with reads in the opposite orientation.\n"
" -msp a b c [d e [f]]: like -fsp above, but with reads in the opposite orientation.\n"
" FALSE PAIRED-END FILES (unpaired reads are split into paired ends, with the scores of double-use nuceotides halved)\n"
" -spf a b c [d e [f]]: (a)input file, (b)the length of the 'reads' that will be taken from each side of the input reads,\n"
"                  (c)amplicon insert size (including read) \n"
"                  (d,e,f) are optional; (d)the num. cycles to be skipped before this file is used;\n"
"                  if (f) is provided, then the file will alternate between being used for (e) cycles and not used for (f) cycles;\n"
"                  otherwise, the file will be used for (e) cycles then will not be used again. \n"
" -spfp a b c d [e f [g]]: (a)input file, (b)the length of the 'reads' that will be taken from each side of the input reads,\n"
"                  (c)amplicon insert size (including read), (d)required % identity for match (25-100 allowed) \n"
"                  (e,f,g) are optional; (e)the num. cycles to be skipped before this file is used;\n"
"                  if (g) is provided, then the file will alternate between being used for (f) cycles and not used for (g) cycles;\n"
"                  otherwise, the file will be used for (f) cycles then will not be used again. \n"
"INPUT INITIAL CONTIG FILES: \n"
"  NOTE: these flags can be used multiple times in the same command to include multiple initial contig datasets. \n"
" -icf a b c d: (a)initial contig file, (b)number of addition steps, (c)number of cycles per step,\n"
"               (d)const by which to multiply quality scores \n"
" -picf a b c d e: (a)num of initial contigs from this file, (b)initial contig file, (c)num addition steps,\n"
"                  (d)num cycles per step (e)const by which to multiply quality scores \n"
" -icfNt/-picfNt: same as -icf/-picf, but if target mode is invoked, contigs with matches to these input sequences will not necessarily\n"
"                 be retained \n"
"OUTPUT FILES: \n"
" accepted formats are fasta (.fa or .fasta) or priceq (.pq or .priceq) \n"
" -o a: (a)output file name (.fasta or .priceq) \n"
" -nco a: (a)num. cycles that pass in between output files being written (default=1) \n"
"OTHER PARAMS: \n"
" -nc a: (a)num. of cycles \n"
" -link a: (a)max. number of contigs that are allowed to replace a read in a mini-assembly (default=2)\n"
" -mol a: (a)minimum overlap length for mini-assembly (default=35) \n"
" -tol a: (a)threshold seq num for scaling overlap for overhang assemblies (default=20)\n"
" NOTE: -mol and -tol do not affect the parameters for de-Bruijn-graph-based assembly.\n"
" -mpi a: (a)minimum % identity for mini-assembly (default=85) \n"
" -tpi a: (a)threshold seq num for scaling % ID for mini-assemblies (default=20)\n"
" -MPI, -TPI : same as above, but for meta-assembly (-MPI default=85, -TPI default=1000)\n"
" NOTE: there is no minimum overlap value for meta-assembly\n"
" -dbmax a: (a) the maximum length sequence that will be fed into de Bruijn assembly \n"
"           (default=100; recommended: max paired-end read length)\n"
" -dbk a: (a) the k-mer size for de Bruijn assembly (default=20; keep less than the read length)\n"
" -dbms a: (a) the minimum number of sequences to which de Bruijn assembly will be applied (default=3)\n"
" -r a: (a) alignment score reward for a nucleotide match; should be a positive integer (default=1)\n"
" -q a: (a) alignment score penalty for a nucleotide mismatch; should be a negative integer (default=-2)\n"
" -G a: (a) alignment score penalty for opening a gap; should be a negative integer (default=-5)\n"
" -E a: (a) alignment score penalty for extending a gap; should be a negative integer (default=-2)\n"
"FILTERING READS: \n"
" -rqf a b [c d]: filters pairs of reads if either has an unaccptably high number of low-quality nucleotides, as defined\n"
"                 by the provided quality scores (only applies to files whose formats include quality score information). \n"
"                 (a) the percentage of nucleotides in a read that must be high-quality; (b) the minimum allowed probability\n"
"                 of a nucleotide being correct (must be between 0 and 1, and will usually be a decimal value close to 1);\n"
"                 (c) and (d) optionally constrain this filter to use after (c) cycles have passed, to run for (d) cycles.\n"
"                 This flag may be called multiple times to generate variable behavior across a PRICE run.\n"
" -rnf a [b c]: filters pairs of reads if either has an unaccptably high number of uncalled nucleotides (Ns or other ambiguous\n"
"                 IUPAC codes).  Like -rqf, but will also filter fasta-format data.  (a) the percentage of nucleotides in a read\n"
"                 that must be called; (b) and (c) optionally constrain this filter to use after (b) cycles have passed, to run\n"
"                 for (c) cycles.  This flag may be called multiple times to generate variable behavior across a PRICE run.\n"
" -maxHp a: filters out a pair of reads if either read has a homo-polymer track >(a) nucleotides in length.\n"
" -maxDi a: filters out a pair of reads if either read has a repeating di-nucleotide track >(a) nucleotides in length.\n"
"           NOTE: this will also catch mono-nucleotide repeats of the specified length (a string of A's is also a string\n"
"           of AA's), so calling -maxHp in addition to -maxDi is superfluous unless -maxHp is given a smaller max value.\n"
" -badf a b: prevents reads with a match of at least (b)% identity to a sequence in file (a) from being mapped to contigs.\n"
" -repmask a b c d e f [g]: uses coverage levels of constructed and/or input contigs to find repetitive elements and mask them \n"
"                           as if they were sequences input using -badf.\n"
"                           (a) = cycle number (1-indexed) at which repeats will be detected.\n"
"                           (b) = 's' if repeats will be sought at the start of the cycle or 'f' if they will be sought at the finish.\n"
"                           (c) = the min. number of variance units above the median that will be counted as high-coverage.\n"
"                           (d) = the min. fold increase in coverage above the median that will be counted as high-coverage.\n"
"                           (e) = the min. size in nt for a detected repeat. \n"
"                           (f) = reads with a match of at least this % identity to a repeat will not be mapped to contigs.\n"
"                           (g) = an optional output file (.fasta or .priceq) to which the detected repeats will be written.\n"
" -reset a [b c d...]: re-introduces contigs that were previously not generating assembly jobs of their own\n"
"                      (a) is the one-indexed cycle where the contigs will be reset.  Same with b, c, d. \n"
"                      Any number of args may be added.\n"
"FILTERING INITIAL CONTIGS: \n"
" -icbf a b [c]: prevents input sequences with a match of at least (b)% identity to a sequence in file (a) from being used.\n"
"                This filter is optionally not applied to sequences of length greater than (c) nucleotides.\n"
" -icmHp a [b]: filters out an initial contig if it has a homo-polymer track >(a) nucleotides in length.\n"
"               This filter is optionally not applied to sequences of length greater than (b) nucleotides.\n"
" -icmDi a [b]: filters out an initial contig if it has a repeating di-nucleotide track >(a) nucleotides in length.\n"
"               NOTE: this will also catch mono-nucleotide repeats of the specified length (a string of A's is also a string\n"
"               of AA's), so calling -icmHp in addition to -icmDi is superfluous unless -icmHp is given a smaller max value.\n"
"               This filter is optionally not applied to sequences of length greater than (b) nucleotides.\n"
" -icqf a b [c]: filters out an initial contig if it has an unaccptably high number of low-quality nucleotides, as defined\n"
"                by the provided quality scores (only applies to files whose formats include quality score information). \n"
"                (a) the percentage of nucleotides in a read that must be high-quality; (b) the minimum allowed probability\n"
"                of a nucleotide being correct (must be between 0 and 1, and will usually be a decimal value close to 1);\n"
"                This filter is optionally not applied to sequences of length greater than (c) nucleotides.\n"
" -icnf a [b]: filters pairs of reads if either has an unaccptably high number of uncalled nucleotides (Ns or other ambiguous\n"
"              IUPAC codes).  Like -icqf, but will also filter fasta-format data.  (a) the percentage of nucleotides in a read\n"
"              that must be called. This filter is optionally not applied to sequences of length greater than (c) nucleotides.\n"
"FILTERING/PROCESSING ASSEMBLED CONTIGS: \n"
" -lenf a b: filters out contigs shorter than (a) nt at the end of every cycle, after skipping (b) cycles.\n"
"            NOTE: multiple -lenf commands can be entered; for any cycle, the most recently-initiated filter is used.\n"
"            Example: -lenf 50 2 -lenf 300 4 -lenf 200 6 => no filter for the first two cycles, then a 50nt filter for cycles\n"
"                     3 & 4, then a 300nt filter for cycles 5 & 6, then a 200nt filter for cycles 7 onwards.\n"
" -trim a b [c]: at the end of the (a)th cycle (indexed from 1), trim off the edges of conigs until reaching the minimum coverage\n"
"                level (b), optionally deleting contigs shorter than (c) after trimming; this flag may be used repeatedly.\n"
" -trimB a b [c]: basal trim; after skipping (a) cycles, trim off the edges of conigs until reaching the minimum coverage\n"
"                 level (b) at the end of EVERY cycle, optionally deleting contigs shorter than (c) after trimming.\n"
"                 -trimB may be called many times, and multiple calls will interact in the same way as multiple -lenf calls\n"
"                 (explained above). A call to -trim will override the basal trim values for that specified cycle only.\n"
" -trimI a [b]: initial trim; input initial contigs are trimmed before being used by PRICE to seed assemblies. Contigs are\n"
"               trimmed from their outside edges until reaching the minimum coverage level (a), optionally deleting contigs\n"
"               shorter than (b) after trimming. This flag is most appropriate for .priceq input, can be appropriate for\n"
"               .fastq input, and is inappropriate for .fasta input.  It will be equally applied to ALL input contigs.\n"
" -target a b [c d]: limit output contigs to those with matches to input initial contigs at the end of each cycle.\n"
"                    (a) % identity to an input initial contig to count as a match (ungapped); (b)num cycles to skip\n"
"                    before applying this filter.  [c and d are optional, but must both be provided if either is]\n"
"                    After target filtering has begin, target-filtered/-unfiltered cycles will alternate with (c)\n"
"                    filtered cycles followed by (d) unfiltered cycles.\n"
" -targetF a b [c d]: the same as -target, but now matches to all reads in the input set will be specified, not just\n"
"                    the ones that have been introduced up to that point (this is FullFile mode).\n"
"COMPUTATIONAL EFFICIENCY: \n"
" -a x: (x)num threads to use (default=1) \n"
" -mtpf a: (a)max threads per file (default=1) \n"
"USER INTERFACE: \n"
" -log a: determines the type of outputmakes the output verbose (lots of time stamp tags) \n"
"         (a) = c: concise stdout (default)\n"
"         (a) = n: no stdout \n"
"         (a) = v: verbose stdout \n"
" -logf a: (a)the name of an output file for verbose log info to be written (doesn't change stdout format) \n"
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
  bool willDoAssembly = true; // allows full set of errors to be printed
  bool neitherAssemblyNorError = false; // allows some calls to cause no assembly but also no fail msg
  bool distributeInitialContigs = false;
  int cyclesPerStep;
  int numSteps;

  // default is 1 core
  omp_set_num_threads( 1 );
  int maxThreadsPerFile = 1;

  bool listenerIsNull = false;
  bool listenerIsVerbose = false;
  bool includeListenerFile = false;
  string listenerFilename = "";

  // first, get the total number of cycles and create the cycle manager
  int numOfCyclesToDo;
  bool numCyclesSpecified = false;
  int argn = 1;
  // same behavior as if --help had been called
  if (argc < 2){
    cout << helpString << endl;
    willDoAssembly = false;
    neitherAssemblyNorError = true;
  }
  while ( argn < argc ){
    string arg = string(argv[argn]);
    int numFlagArgs = getNumFlagArgs(argn,argv,argc);

    if (arg == "-h" or arg == "-" or arg == "--help"){
      cout << helpString << endl;
      willDoAssembly = false;
      neitherAssemblyNorError = true;
    } else if (arg == "-nc"){
      if (numFlagArgs != 1){
	cerr << "-nc was given " << numFlagArgs << " args instead of the required 1." << endl;
	willDoAssembly = false;
      } else if (! shouldBeIntArg(argn,1,argv)){
	willDoAssembly = false;
      } else {
	numOfCyclesToDo = atoi( argv[argn+1] );
	numCyclesSpecified = true;
      }
    }
    argn += numFlagArgs + 1;
  }


  if (! numCyclesSpecified){
    cerr << "number of cycles MUST be specified with -nc flag" << endl;
    willDoAssembly = false;
  }


  if (willDoAssembly){

    Assembler* assembler = new Assembler(numOfCyclesToDo);
    assembler->setMetaFractIdContigThreshold(1000);

    // now go through again and get the command-line args organized
    argn = 1;
    while ( argn < argc ){
      string arg = string(argv[argn]);
      int numFlagArgs = getNumFlagArgs(argn,argv,argc);

      if (arg == "-fp"){
        if (numFlagArgs != 3 and numFlagArgs != 5 and numFlagArgs != 6){
	  cerr << "-fp was given " << numFlagArgs << " args instead of the required 3, 5, or 6." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 3){
	  string fileA = argv[argn+1];
	  string fileB = argv[argn+2];
	  if (isFileMissing(fileA)){
	    cerr << "The following file input with the -fp flag could not be found: " << fileA << endl;
	    willDoAssembly = false;
	  }
	  if (isFileMissing(fileB)){
	    cerr << "The following file input with the -fp flag could not be found: " << fileB << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+3]);
	  if (numFlagArgs == 3){
	    if (willDoAssembly){ assembler->addReadFile(fileA,fileB,ampSize); }
	  } else if (numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int numCycles = atoi(argv[argn+5]);
	    if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int onCycles = atoi(argv[argn+5]);
	    int skipCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize); }
	  }
	}


      } else if (arg == "-fpp"){
        if (numFlagArgs != 4 and numFlagArgs != 6 and numFlagArgs != 7){
	  cerr << "-fpp was given " << numFlagArgs << " args instead of the required 4, 6, or 7." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 4){
	  string fileA = argv[argn+1];
	  string fileB = argv[argn+2];
	  if (isFileMissing(fileA)){
	    cerr << "The following file input with the -fpp flag could not be found: " << fileA << endl;
	    willDoAssembly = false;
	  }
	  if (isFileMissing(fileB)){
	    cerr << "The following file input with the -fpp flag could not be found: " << fileB << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  if (! worksAsFloat(argv[argn+4])){
	    cerr << "the 4th arg for the " << arg << " flag should have been a number: " << argv[argn+4] << endl;
	    willDoAssembly = false;
	  }
	  long ampSize = atol(argv[argn+3]);
	  float fractId = float(atoi(argv[argn+4])) / float(100);
	  if (numFlagArgs == 4){
	    if (willDoAssembly){ assembler->addReadFile(fileA,fileB,ampSize,fractId); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+5]);
	    int numCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize); }
	  } else if (numFlagArgs == 7){
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,7,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+5]);
	    int onCycles = atoi(argv[argn+6]);
	    int skipCycles = atoi(argv[argn+7]);
	    if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize); }
	  }
	}


      } else if (arg == "-fs"){
        if (numFlagArgs != 2 and numFlagArgs != 4 and numFlagArgs != 5){
	  cerr << "-fs was given " << numFlagArgs << " args instead of the required 2, 4, or 5." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 2){
	  string file = argv[argn+1];
	  if (isFileMissing(file)){
	    cerr << "The following file input with the -fs flag could not be found: " << file << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+2]);
	  if (numFlagArgs == 2){
	    if (willDoAssembly){ assembler->addReadFile(file,ampSize); }
	  } else if (numFlagArgs == 4 or numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (numFlagArgs == 4){
	      int startCycle = atoi(argv[argn+3]);
	      int numCycles = atoi(argv[argn+4]);
	      if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,numCycles,file,ampSize); }
	    } else if (numFlagArgs == 5){
	      if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	      int startCycle = atoi(argv[argn+3]);
	      int onCycles = atoi(argv[argn+4]);
	      int skipCycles = atoi(argv[argn+5]);
	      if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize); }
	    }
	  }
	}


      } else if (arg == "-fsp"){
        if (numFlagArgs != 3 and numFlagArgs != 5 and numFlagArgs != 6){
	  cerr << "-fsp was given " << numFlagArgs << " args instead of the required 3, 5, or 6." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 3){
	  string file = argv[argn+1];
	  if (isFileMissing(file)){
	    cerr << "The following file input with the -fsp flag could not be found: " << file << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+2]);
	  float fractId = float(atoi(argv[argn+3])) / float(100);
	  if (numFlagArgs == 3){
	    if (willDoAssembly){ assembler->addReadFile(file,ampSize,fractId); }
	  } else if (numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int numCycles = atoi(argv[argn+5]);
	    if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,numCycles,file,ampSize); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int onCycles = atoi(argv[argn+5]);
	    int skipCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addReadFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize); }
	  }
	}


      } else if (arg == "-mp"){
        if (numFlagArgs != 3 and numFlagArgs != 5 and numFlagArgs != 6){
	  cerr << "-mp was given " << numFlagArgs << " args instead of the required 3, 5, or 6." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 3){
	  string fileA = argv[argn+1];
	  string fileB = argv[argn+2];
	  if (isFileMissing(fileA)){
	    cerr << "The following file input with the -mp flag could not be found: " << fileA << endl;
	    willDoAssembly = false;
	  }
	  if (isFileMissing(fileB)){
	    cerr << "The following file input with the -mp flag could not be found: " << fileB << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+3]);
	  if (numFlagArgs == 3){
	    if (willDoAssembly){ assembler->addMatePairFile(fileA,fileB,ampSize); }
	  } else if (numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int numCycles = atoi(argv[argn+5]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int onCycles = atoi(argv[argn+5]);
	    int skipCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize); }
	  }
	}


      } else if (arg == "-mpp"){
        if (numFlagArgs != 4 and numFlagArgs != 6 and numFlagArgs != 7){
	  cerr << "-mpp was given " << numFlagArgs << " args instead of the required 4, 6, or 7." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 4){
	  string fileA = argv[argn+1];
	  string fileB = argv[argn+2];
	  if (isFileMissing(fileA)){
	    cerr << "The following file input with the -mpp flag could not be found: " << fileA << endl;
	    willDoAssembly = false;
	  }
	  if (isFileMissing(fileB)){
	    cerr << "The following file input with the -mpp flag could not be found: " << fileB << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+3]);
	  float fractId = float(atoi(argv[argn+4])) / float(100);
	  if (numFlagArgs == 4){
	    if (willDoAssembly){ assembler->addMatePairFile(fileA,fileB,ampSize,fractId); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+5]);
	    int numCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,numCycles,fileA,fileB,ampSize); }
	  } else if (numFlagArgs == 7){
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,7,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+5]);
	    int onCycles = atoi(argv[argn+6]);
	    int skipCycles = atoi(argv[argn+7]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,fileA,fileB,ampSize); }
	  }
	}


      } else if (arg == "-ms"){
        if (numFlagArgs != 2 and numFlagArgs != 4 and numFlagArgs != 5){
	  cerr << "-ms was given " << numFlagArgs << " args instead of the required 2, 4, or 5." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 2){
	  string file = argv[argn+1];
	  if (isFileMissing(file)){
	    cerr << "The following file input with the -ms flag could not be found: " << file << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+2]);
	  if (numFlagArgs == 2){
	    if (willDoAssembly){ assembler->addMatePairFile(file,ampSize); }
	  } else if (numFlagArgs == 4){
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+3]);
	    int numCycles = atoi(argv[argn+4]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,numCycles,file,ampSize); }
	  } else if (numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+3]);
	    int onCycles = atoi(argv[argn+4]);
	    int skipCycles = atoi(argv[argn+5]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize); }
	  }
	}


      } else if (arg == "-msp"){
        if (numFlagArgs != 3 and numFlagArgs != 5 and numFlagArgs != 6){
	  cerr << "-msp was given " << numFlagArgs << " args instead of the required 3, 5, or 6." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 3){
	  string file = argv[argn+1];
	  if (isFileMissing(file)){
	    cerr << "The following file input with the -msp flag could not be found: " << file << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  long ampSize = atol(argv[argn+2]);
	  float fractId = float(atoi(argv[argn+3])) / float(100);
	  if (numFlagArgs == 3){
	    if (willDoAssembly){ assembler->addMatePairFile(file,ampSize,fractId); }
	  } else if (numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int numCycles = atoi(argv[argn+5]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,numCycles,file,ampSize); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int onCycles = atoi(argv[argn+5]);
	    int skipCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addMatePairFileLimitUse(startCycle,onCycles,skipCycles,file,ampSize); }
	  }
	}


      } else if (arg == "-spf"){
        if (numFlagArgs != 3 and numFlagArgs != 5 and numFlagArgs != 6){
	  cerr << "-spf was given " << numFlagArgs << " args instead of the required 3, 5, or 6." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 3){
	  string file = argv[argn+1];
	  if (isFileMissing(file)){
	    cerr << "The following file input with the -spf flag could not be found: " << file << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  long readSize = atol(argv[argn+2]);
	  long ampSize = atol(argv[argn+3]);
	  if (numFlagArgs == 3){
	    if (willDoAssembly){ assembler->addFalsePairFile(file,readSize,ampSize); }
	  } else if (numFlagArgs == 5){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int numCycles = atoi(argv[argn+5]);
	    if (willDoAssembly){ assembler->addFalsePairFileLimitUse(startCycle,numCycles,file,readSize,ampSize); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+4]);
	    int onCycles = atoi(argv[argn+5]);
	    int skipCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addFalsePairFileLimitUse(startCycle,onCycles,skipCycles,file,readSize,ampSize); }
	  }
	}


      } else if (arg == "-spfp"){
        if (numFlagArgs != 4 and numFlagArgs != 6 and numFlagArgs != 7){
	  cerr << "-spfp was given " << numFlagArgs << " args instead of the required 4, 6, or 7." << endl;
	  willDoAssembly = false;
	}
	if ( numFlagArgs >= 4){
	  string file = argv[argn+1];
	  if (isFileMissing(file)){
	    cerr << "The following file input with the -spfp flag could not be found: " << file << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  long readSize = atol(argv[argn+2]);
	  long ampSize = atol(argv[argn+3]);
	  float fractId = float(atoi(argv[argn+4])) / float(100);
	  if (numFlagArgs == 4){
	    if (willDoAssembly){ assembler->addFalsePairFile(file,readSize,ampSize,fractId); }
	  } else if (numFlagArgs == 6){
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+5]);
	    int numCycles = atoi(argv[argn+6]);
	    if (willDoAssembly){ assembler->addFalsePairFileLimitUse(startCycle,numCycles,file,readSize,ampSize); }
	  } else if (numFlagArgs == 7){
	    if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,6,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,7,argv)){ willDoAssembly = false; }
	    int startCycle = atoi(argv[argn+5]);
	    int onCycles = atoi(argv[argn+6]);
	    int skipCycles = atoi(argv[argn+7]);
	    if (willDoAssembly){ assembler->addFalsePairFileLimitUse(startCycle,onCycles,skipCycles,file,readSize,ampSize); }
	  }
	}


      } else if (arg == "-icf"){
        if (numFlagArgs != 4){
	  cerr << "-icf was given " << numFlagArgs << " args instead of the required 4." << endl;
	  willDoAssembly = false;
	} else {
	  string filename = argv[argn+1];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -icf flag could not be found: " << filename << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,4,argv)){ willDoAssembly = false; }
	  int numSteps = atoi( argv[argn+2] );
	  int cyclesPerStep = atoi( argv[argn+3] );
	  float countFactor = atof( argv[argn+4] );
	  if (willDoAssembly){ assembler->addInitialFile(filename,numSteps,cyclesPerStep,countFactor); }
	}


      } else if (arg == "-picf"){
        if (numFlagArgs != 5){
	  cerr << "-picf was given " << numFlagArgs << " args instead of the required 5." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  long initialInput = atol( argv[argn+1] );
	  string filename = argv[argn+2];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -picf flag could not be found: " << filename << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,5,argv)){ willDoAssembly = false; }
	  int numSteps = atoi( argv[argn+3] );
	  int cyclesPerStep = atoi( argv[argn+4] );
	  float countFactor = atof( argv[argn+5] );
	  if (willDoAssembly){ assembler->addInitialFile(initialInput,filename,numSteps,cyclesPerStep,countFactor); }
	}


      } else if (arg == "-icfNt"){
        if (numFlagArgs != 4){
	  cerr << "-icfNt was given " << numFlagArgs << " args instead of the required 4." << endl;
	  willDoAssembly = false;
	} else {
	  string filename = argv[argn+1];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -icfNt flag could not be found: " << filename << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,4,argv)){ willDoAssembly = false; }
	  int numSteps = atoi( argv[argn+2] );
	  int cyclesPerStep = atoi( argv[argn+3] );
	  float countFactor = atof( argv[argn+4] );
	  if (willDoAssembly){ assembler->addInitialFileNoTarget(filename,numSteps,cyclesPerStep,countFactor); }
	}


      } else if (arg == "-picfNt"){
        if (numFlagArgs != 5){
	  cerr << "-picfNt was given " << numFlagArgs << " args instead of the required 5." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  long initialInput = atol( argv[argn+1] );
	  string filename = argv[argn+2];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -picfNt flag could not be found: " << filename << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,5,argv)){ willDoAssembly = false; }
	  int numSteps = atoi( argv[argn+3] );
	  int cyclesPerStep = atoi( argv[argn+4] );
	  float countFactor = atof( argv[argn+5] );
	  if (willDoAssembly){ assembler->addInitialFileNoTarget(initialInput,filename,numSteps,cyclesPerStep,countFactor); }
	}


      } else if (arg == "-nc"){
	// nothing happens here because this info was already collected


      } else if (arg == "-repmask"){
        if (numFlagArgs != 6 and numFlagArgs != 7){
	  cerr << "-repmask was given " << numFlagArgs << " args instead of the required 6 or 7." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,3,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,4,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,5,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,6,argv)){ willDoAssembly = false; }
	  int cycleNum = atoi(argv[argn+1]);
	  string whenInCycle = argv[argn+2];
	  float minStdDev = atof(argv[argn+3]);
	  float minFoldUp = atof(argv[argn+4]);
	  long minSize = atol(argv[argn+5]);
	  float minFractId = atof(argv[argn+6]) / float(100);
	  if (minFractId < 0.25){
	    cerr << "-repmask min % ID must be at least 25%" << endl;
	    willDoAssembly = false;
	  }

	  // an output file may (7 args) or may not (6 args) be specified
          if (numFlagArgs == 7){
	    string outputFile = argv[argn+7];
	    if (whenInCycle[0] == 's'){
	      if (willDoAssembly){ assembler->findRepsBeforeCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId, outputFile); }
	    } else if (whenInCycle[0] == 'f'){
	      if (willDoAssembly){ assembler->findRepsAfterCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId, outputFile); }
	    } else {
	      cerr << "-repmask second arg must be 's' or 'f', not '" << whenInCycle << "'." << endl;
	      willDoAssembly = false;
	    }
	  } else {
	    if (whenInCycle[0] == 's'){
	      if (willDoAssembly){ assembler->findRepsBeforeCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId); }
	    } else if (whenInCycle[0] == 'f'){
	      if (willDoAssembly){ assembler->findRepsAfterCycle(cycleNum, minStdDev, minFoldUp, minSize, minFractId); }
	    } else {
	      cerr << "-repmask second arg must be 's' or 'f', not '" << whenInCycle << "'." << endl;
	      willDoAssembly = false;
	    }
	  }
	}
 

      } else if (arg == "-target" or arg == "-targetF"){
        if (numFlagArgs != 2 and numFlagArgs != 4){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2 or 4." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  bool fullFileMode = (arg == "-targetF");
	  float fractId = atof( argv[argn+1] ) / float(100.0);
	  int cyclesToSkip = atoi( argv[argn+2] );
	  if (numFlagArgs == 2){
	    if (willDoAssembly){ assembler->addTargetMode(fractId,cyclesToSkip,fullFileMode); }
	  } else {
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    int numFilterCycles = atoi( argv[argn+3] );
	    int numSkipCycles = atoi( argv[argn+4] );
	    if (willDoAssembly){ assembler->addTargetMode(fractId,cyclesToSkip,numFilterCycles,numSkipCycles,fullFileMode); }
	  }
	}


      } else if (arg == "-r" or arg == "-q" or arg == "-G" or arg == "-E"){
	// ALIGNMENT SCORING PARAMETERS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (willDoAssembly){
	    if (arg == "-r"){ assembler->setNucMatchScore( atol(argv[argn+1]) ); }
	    else if (arg == "-q"){ assembler->setNucMismatchPenalty( atol(argv[argn+1]) ); }
	    else if (arg == "-G"){ assembler->setOpenGapPenalty( atol(argv[argn+1]) ); }
	    else if (arg == "-E"){ assembler->setExtendGapPenalty( atol(argv[argn+1]) ); }
	  }
	}


      } else if (arg == "-badf"){
        if (numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2." << endl;
	  willDoAssembly = false;
	} else if (willDoAssembly){
	  string filename = argv[argn+1];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -badf flag could not be found: " << filename << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoAssembly = false; }
	  float fractId = atof( argv[argn+2] ) / float(100.0);
	  if (willDoAssembly){ assembler->addBadSequenceFilter(filename, fractId); }
	}


      } else if (arg=="-maxHp" or arg=="-maxDi"){
	// COMPLEXITY READ FILTERS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else if (willDoAssembly){
	  if (arg == "-maxHp"){ assembler->addHomopolymerFilter( atol(argv[argn+1]) ); }
	  else if (arg == "-maxDi"){ assembler->addDinucRepeatFilter( atol(argv[argn+1]) ); }
	}


      } else if (arg=="-icmHp" or arg=="-icmDi"){
	// COMPLEXITY INITIAL CONTIG FILTERS
        if (numFlagArgs != 1 and numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1 or 2." << endl;
	  willDoAssembly = false;
	} else if (numFlagArgs==1){
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  long maxRepLen = atol(argv[argn+1]);
	  if (willDoAssembly){
	    if (arg == "-icmHp"){ assembler->icHomopolymerFilter( maxRepLen ); }
	    else if (arg == "-icmDi"){ assembler->icDinucRepeatFilter( maxRepLen ); }
	  }
	} else if (numFlagArgs==2){
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  long maxRepLen = atol(argv[argn+1]);
	  long maxSeqLen = atol(argv[argn+2]);
	  if (willDoAssembly){
	    if (arg == "-icmHp"){ assembler->icHomopolymerFilter( maxRepLen, maxSeqLen ); }
	    else if (arg == "-icmDi"){ assembler->icDinucRepeatFilter( maxRepLen, maxSeqLen ); }
	  }
	}


      } else if (arg == "-rqf"){
        if (numFlagArgs != 2 and numFlagArgs != 4){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2 or 4." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoAssembly = false; }
	  float minFractGood = atof(argv[argn+1]) / float(100);
	  float minProbCorrect = atof(argv[argn+2]);
	  if (numFlagArgs==2){
	    if (willDoAssembly){ assembler->addReadQualityFilter(minFractGood, minProbCorrect); }
	  } else {
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,4,argv)){ willDoAssembly = false; }
	    int numSkipCycles = atoi(argv[argn+3]);
	    int numRunCycles = atoi(argv[argn+4]);
	    if (willDoAssembly){ assembler->addReadQualityFilter(minFractGood, minProbCorrect, numSkipCycles, numRunCycles); }
	  }
	}


      } else if (arg == "-rnf"){
        if (numFlagArgs != 1 and numFlagArgs != 3){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1 or 3." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  float minFractGood = atof(argv[argn+1]) / float(100);
	  if (numFlagArgs==1){
	    if (willDoAssembly){ assembler->addReadCalledBasesFilter(minFractGood); }
	  } else {
	    if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    int numSkipCycles = atoi(argv[argn+2]);
	    int numRunCycles = atoi(argv[argn+3]);
	    if (willDoAssembly){ assembler->addReadCalledBasesFilter(minFractGood, numSkipCycles, numRunCycles); }
	  }
	}


      } else if (arg == "-icbf"){
        if (numFlagArgs != 2 and numFlagArgs != 3){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2 or 3." << endl;
	  willDoAssembly = false;
	} else {
	  string filename = argv[argn+1];
	  if (isFileMissing(filename)){
	    cerr << "The following file input with the -icbf flag could not be found: " << filename << endl;
	    willDoAssembly = false;
	  }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoAssembly = false; }
	  float fractId = atof( argv[argn+2] ) / float(100.0);
	  if (numFlagArgs == 2){
	    if (willDoAssembly){ assembler->icBadSequenceFilter(filename, fractId); }
	  } else {
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    long maxSeqLen = atol( argv[argn+3] );
	    if (willDoAssembly){ assembler->icBadSequenceFilter(filename, fractId, maxSeqLen); }
	  }
	}


      } else if (arg == "-icqf"){
        if (numFlagArgs != 2 and numFlagArgs != 3){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2 or 3." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoAssembly = false; }
	  float minFractGood = atof(argv[argn+1]) / float(100);
	  float minProbCorrect = atof(argv[argn+2]);
	  if (numFlagArgs == 2){
	    if (willDoAssembly){ assembler->icReadQualityFilter(minFractGood, minProbCorrect); }
	  } else {
	    if (! shouldBeIntArg(argn,3,argv)){ willDoAssembly = false; }
	    long maxSeqLen = atol(argv[argn+3]);
	    if (willDoAssembly){ assembler->icReadQualityFilter(minFractGood, minProbCorrect, maxSeqLen); }
	  }
	}


      } else if (arg == "-icnf"){
        if (numFlagArgs != 1 and numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1 or 2." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  float minFractGood = atof(argv[argn+1]) / float(100);
	  if (numFlagArgs == 1){
	    if (willDoAssembly){ assembler->icReadCalledBasesFilter(minFractGood); }
	  } else {
	    if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	    long maxSeqLen = atol(argv[argn+2]);
	    if (willDoAssembly){ assembler->icReadCalledBasesFilter(minFractGood, maxSeqLen); }
	  }
	}


      } else if (arg == "-lenf"){
        if (numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	  long minLength = atol( argv[argn+1] );
	  int cyclesToSkip = atoi( argv[argn+2] );
	  if (willDoAssembly){ assembler->addLengthFilter(minLength, cyclesToSkip); }
	}


      } else if (arg == "-o"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  string filename = argv[argn+1];
	  if (willDoAssembly){ assembler->addOutputFile(filename); }
	}


      } else if (arg == "-nco"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (willDoAssembly){ assembler->setOutputInterval(atoi( argv[argn+1] )); }
	}


      } else if (arg=="-a" or arg=="mtpf"){
	// COMP EFFICIENCY ARGS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
          if (willDoAssembly){
	    if (arg == "-a"){ omp_set_num_threads( atoi( argv[argn+1] ) ); }
	    else if (arg == "-mtpf"){ maxThreadsPerFile = atoi( argv[argn+1] ); }
	  }
	}


      } else if (arg=="-dbk" or arg=="-dbmax" or arg=="-dbms"){
	// DE BRUIJN PARAMETERS
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (willDoAssembly){
	    if (arg == "-dbk"){ assembler->setDeBruijnKmerSize( atoi( argv[argn+1] ) ); }
	    else if (arg == "-dbmax"){ assembler->setDeBruijnMaxSeqLength( atol( argv[argn+1] ) ); }
	    else if (arg == "-dbms"){ assembler->setDeBruijnMinSeqNum( atol( argv[argn+1] ) ); }
	  }
	}


      } else if (arg=="-trim" or arg=="-trimB"){
	// CYCLIC TRIMMING
        if (numFlagArgs != 2 and numFlagArgs != 3){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 2 or 3." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (! shouldBeFloatArg(argn,2,argv)){ willDoAssembly = false; }
	  int cycleArg = atoi( argv[argn+1] );
	  float minCoverage = atof( argv[argn+2] );
	  if (numFlagArgs == 2){
	    if (willDoAssembly){
	      if (arg == "-trim"){ assembler->addTrim(cycleArg, minCoverage); }
	      else if (arg == "-trimB"){ assembler->addBasalTrim(cycleArg, minCoverage); }
	    }
	  } else if (numFlagArgs == 3){
	    int minLength = atol( argv[argn+3] );
	    if (willDoAssembly){
	      if (arg == "-trim"){ assembler->addTrim(cycleArg, minCoverage, minLength); }
	      else if (arg == "-trimB"){ assembler->addBasalTrim(cycleArg, minCoverage, minLength); }
	    }
	  }
	}


      } else if (arg == "-trimI"){
        if (numFlagArgs != 1 and numFlagArgs != 2){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1 or 2." << endl;
	  willDoAssembly = false;
	} else {
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  float minCoverage = atof( argv[argn+1] );
	  if (numFlagArgs == 1){
	    if (willDoAssembly){ assembler->addInitialTrim(minCoverage); }
	  } else {
	    if (! shouldBeIntArg(argn,2,argv)){ willDoAssembly = false; }
	    long minLength = atol( argv[argn+2] );
	    if (willDoAssembly){ assembler->addInitialTrim(minCoverage, minLength); }
	  }
	}


      } else if (arg == "-reset"){
	int numResets = numFlagArgs;
	if (numResets == 0 ){
	  cerr << "-reset should have at least one cycle number after it." << endl;
	  willDoAssembly = false;
	} else {
	  if (willDoAssembly){
	    for (int n = 1; n <=numResets; ++n){ 
	      if (! shouldBeIntArg(argn,n,argv)){ willDoAssembly = false; }
	      if (willDoAssembly){ assembler->addResetCycle( atoi( argv[argn+n] ) ); }
	    }
	  }
	}


      } else if (arg=="-link" or arg=="-mol" or arg=="-tol" or arg=="-mpi" or arg=="-tpi" or arg=="-MPI" or arg=="-TPI"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else if (arg=="-link" or arg=="-mol" or arg=="-tol" or arg=="-tpi" or arg=="-TPI"){
	  // int/long args
	  if (! shouldBeIntArg(argn,1,argv)){ willDoAssembly = false; }
	  if (arg == "-link"){ if (willDoAssembly){ assembler->setLinkMax( atol( argv[argn+1] ) ); } }
	  else if (arg == "-mol"){ if (willDoAssembly){ assembler->setMiniMinOverlap( atol( argv[argn+1] ) ); } }
	  else if (arg == "-tol"){ if (willDoAssembly){ assembler->setMiniMinOvlContigThreshold( atol( argv[argn+1] ) ); } }
	  else if (arg == "-tpi"){ if (willDoAssembly){ assembler->setMiniFractIdContigThreshold( atol( argv[argn+1] ) ); } }
	  else if (arg == "-TPI"){ if (willDoAssembly){ assembler->setMetaFractIdContigThreshold( atol( argv[argn+1] ) ); } }

	} else if (arg=="-mpi" or arg=="-MPI"){
	  // float args
	  if (! shouldBeFloatArg(argn,1,argv)){ willDoAssembly = false; }
	  else {
	    float fractId = atof(argv[argn+1]) / 100.0;
	    if (fractId < 0.25 or fractId > 1.0){
	      cerr << "% ID required for mini-assembly is out of range 25-100, was " << argv[argn+1] << endl;
	      willDoAssembly = false;
	    }
	    if (arg == "-mpi"){ // min percent ID for overhang assemblies
	      if (willDoAssembly){ assembler->setMiniMinFractId(fractId); }
	    } else if (arg == "-MPI"){ // min percent ID for meta assemblies
	      if (willDoAssembly){ assembler->setMetaMinFractId(fractId); }
	    }
	  }
	}


      } else if (arg == "-log"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  string type = string(argv[argn+1]);
	  if (type == "n"){ listenerIsNull = true; }
	  else if (type == "c"){ }
	  else if (type == "v"){ listenerIsVerbose = true; }
	  else {
	    cerr << "'" << type << "' is not a valid stdout type ('n', 'c', or 'v')." << endl;
	    willDoAssembly = false;
	  }
	}


      } else if (arg == "-logf"){
        if (numFlagArgs != 1){
	  cerr << arg << " was given " << numFlagArgs << " args instead of the required 1." << endl;
	  willDoAssembly = false;
	} else {
	  includeListenerFile = true;
	  listenerFilename = argv[argn+1];
	}


      } else {
	cerr << "'" << arg << "' is an invalid flag that was given " << numFlagArgs << " args." << endl;
	willDoAssembly = false;
      }
      argn += numFlagArgs + 1;

      // a double-check, but there should be an error already raised
      if (argn != argc and argv[argn][0] != '-'){
	cerr << arg << " was given an incorrect number of args." << endl;
	willDoAssembly = false;
      }
    }




    // check again!
    if ( willDoAssembly ){

      // TEMPORARY set-up for returning the command that launched the job
      ostringstream launchCommand;
      launchCommand << "Launch command:";
      for (int n = 0; n < argc; ++n){ launchCommand << " " << argv[n]; }
      cout << launchCommand.str().c_str() << endl;


      if (listenerIsNull){ assembler->makeLogNull(); }
      else if (listenerIsVerbose){ assembler->makeLogVerbose(); }
      if (includeListenerFile == true){ assembler->verboseLogFile(listenerFilename); }

      if (maxThreadsPerFile != 1){ assembler->setMaxThreadsPerFile(maxThreadsPerFile); }

      assembler->runAssembly();

    }
    delete assembler;
  }

  if (willDoAssembly){ return 0; }
  else if (neitherAssemblyNorError){
    cerr << "No assembly was run; help message printed instead." << endl;
    return 0;
  } else {
    cerr << "PRICE JOB EXITED WITH ERROR!!!" << endl;
    return 1; 
  }
}




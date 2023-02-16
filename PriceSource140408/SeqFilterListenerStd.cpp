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

/* implementations of this abstract class will monitor the
 * progress of Assembler and thereby handle things like
 * stdout, GUIs, or just take the place of the annoing cout
 * streams that I would otherwise have in place.
 */

#ifndef SEQFILTERLISTENERSTD_CPP
#define SEQFILTERLISTENERSTD_CPP
#include "SeqFilterListenerStd.h"
#include "AssemblyException.h"
#include <string>
#include <set>
#include <iostream>
using namespace::std;

  
SeqFilterListenerStd::SeqFilterListenerStd(){}
SeqFilterListenerStd::~SeqFilterListenerStd(){}

void SeqFilterListenerStd::startingFilter(bool isPaired){ 
  _isPaired = isPaired;
  if (_isPaired){ cout << "Sequence pairs processed/removed:\t0/0"; }
  else { cout << "Sequences processed/removed:\t0/0"; }
  cout.flush();
  _priorPrintSize = 3;
  _numProcessed = 0;
  _numRemoved = 0;
}
void SeqFilterListenerStd::updateProgress(long numProcessed, long numRemoved){
  _numProcessed += numProcessed;
  _numRemoved += numRemoved;
  for (int n = 0; n < _priorPrintSize; ++n){ cout << "\b"; }
  cout << _numProcessed << '/' << _numRemoved;
  cout.flush();
  _priorPrintSize = getNumChars(_numProcessed) + 1 + getNumChars(_numRemoved);
}
void SeqFilterListenerStd::finishingFilter(){
  cout << endl;
  if (_isPaired){ cout << "Num. sequence pairs kept (% of input):\t"; }
  else { cout << "Num. sequences kept (% of input):\t"; }
  long numKept = _numProcessed - _numRemoved;
  float keptPercent = float(100) * float(numKept) / float(_numProcessed);
  cout << numKept << " (" << keptPercent << "%)" << endl;
}

void SeqFilterListenerStd::message(string message){}
void SeqFilterListenerStd::verboseMessage(string message){}
void SeqFilterListenerStd::errorMessage(string message){}

void SeqFilterListenerStd::timeStamp(){}

int SeqFilterListenerStd::getNumChars(long num){
  if (num < 0){ throw AssemblyException::ArgError("SSLS::getNumChars can't use a negative number"); }
  int numChars = 1;
  long currentMin = 10;
  while (num >= currentMin){
    numChars += 1;
    currentMin *= 10;
  }
  return numChars;
}

#endif

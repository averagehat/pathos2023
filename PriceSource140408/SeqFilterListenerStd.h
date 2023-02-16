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

#ifndef SEQFILTERLISTENERSTD_H
#define SEQFILTERLISTENERSTD_H
#include "SeqFilterListener.h"
#include <string>
#include <set>
using namespace::std;


class SeqFilterListenerStd : public SeqFilterListener {

 public:
  
  SeqFilterListenerStd();
  ~SeqFilterListenerStd();

  void startingFilter(bool isPaired);
  void updateProgress(long numProcessed, long numRemoved);
  void finishingFilter();

  void message(string message);
  void verboseMessage(string message);
  void errorMessage(string message);

  void timeStamp();

 private:
  static int getNumChars(long num);

  bool _isPaired;
  long _numProcessed;
  long _numRemoved;
  int _priorPrintSize;
};

#endif

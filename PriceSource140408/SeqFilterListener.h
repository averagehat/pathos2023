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

#ifndef SEQFILTERLISTENER_H
#define SEQFILTERLISTENER_H
#include "ScoredSeq.h"
#include <string>
#include <set>
using namespace::std;


class SeqFilterListener {

 public:
  static SeqFilterListener* makeNullListener();
  static SeqFilterListener* makeStdListener();
  static SeqFilterListener* makeVerboseListener();

  virtual ~SeqFilterListener();

  virtual void startingFilter(bool isPaired) = 0;
  virtual void updateProgress(long numProcessed, long numRemoved) = 0;
  virtual void finishingFilter() = 0;


  virtual void message(string message) = 0;
  virtual void verboseMessage(string message) = 0;
  virtual void errorMessage(string message) = 0;

  virtual void timeStamp() = 0;

  static long getTotalSize(set<ScoredSeq*>* contigs);
  static long getMaxSize(set<ScoredSeq*>* contigs);
  static long getN50(set<ScoredSeq*>* contigs);

};

#endif

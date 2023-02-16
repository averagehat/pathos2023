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

#ifndef SEQFILTERLISTENERNULL_CPP
#define SEQFILTERLISTENERNULL_CPP
#include "SeqFilterListenerNull.h"
#include <string>
#include <set>
using namespace::std;

  
SeqFilterListenerNull::SeqFilterListenerNull(){}
SeqFilterListenerNull::~SeqFilterListenerNull(){}

void SeqFilterListenerNull::startingFilter(bool isPaired){}
void SeqFilterListenerNull::updateProgress(long numProcessed, long numRemoved){}
void SeqFilterListenerNull::finishingFilter(){}

void SeqFilterListenerNull::message(string message){}
void SeqFilterListenerNull::verboseMessage(string message){}
void SeqFilterListenerNull::errorMessage(string message){}

void SeqFilterListenerNull::timeStamp(){}

#endif

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

/*
This is an interface for implementations of ScoredSeq that is designed for efficient file parsing and
writing.  The added methods deal with storage of the file-type-specific information about the sequence
(like the full name line) and the file-type-specific output-writing options.  These need not be visible
when writing to a file.

Currently, these are immutable, but I may expand them in the future to allow for things like sequence
truncation.  This class is the appropriate context for this as the requirements for shortening a sequence
are format-specific.
 */

#ifndef WRITABLESEQ_H
#define WRITABLESEQ_H
# include <vector>
# include <string>
#include "ScoredSeq.h"
#include "fileUtilities.h"
using std::string;
using std::vector;


class WritableSeq : public ScoredSeq {

 public:

  virtual ~WritableSeq();

  virtual char* getFileString() = 0;
  virtual fileUtilities::FileType fileType() = 0;


};

#endif

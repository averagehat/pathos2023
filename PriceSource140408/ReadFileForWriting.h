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
This is an abstract class for a file of reads.  I will create an implementation
for each format of input file, paired-end or not, etc.

This subtype of a ReadFile guarantees output of a special type of ScoredSeq that
makes it easy to write that seq to a file, retaining its original name/etc.
*/


#ifndef READFILEFORWRITING_H
#define READFILEFORWRITING_H

# include <vector>
# include <string>
# include <set>
#include <fstream>
# include "fileUtilities.h"
# include "WritableSeq.h"

using namespace::std;
using namespace::fileUtilities;


// I will maybe develop this into a full ReadFile at some point, but for now,
// it is its own class that just implements some of the ReadFile methods
class ReadFileForWriting {
 public:

  // factory methods
  static ReadFileForWriting* makeBasicReadFile(string filename);

  // these define an interface for use and re-use of read files

  virtual ~ReadFileForWriting();

  virtual string getName() = 0;
  virtual void open() = 0;
  virtual bool hasRead() = 0;
  // and alternative to getRead that returns the writable type.
  // NOTE: does not respect any of the wrapper-type open args;
  // one reason why this is not currently a ReadFile subclass
  virtual WritableSeq* getWritableRead() = 0;
  virtual void close() = 0;

  virtual FileType getFileType() = 0;   
};




class ReadFileForWritingFasta : public ReadFileForWriting {
 public:
  ReadFileForWritingFasta(string filename);
  ~ReadFileForWritingFasta();
  string getName();
  void open();
  bool hasRead();
  WritableSeq* getWritableRead();
  void close();
  FileType getFileType();

 private:
  void tryAndGetRead();
  char* _filename;
  bool _isOpen;
  WritableSeq* _waitingRead;
  bool _readIsWaiting;
  ifstream _getReadsFile; // the file object for obtaining reads
};




class ReadFileForWritingFastq : public ReadFileForWriting {
 public:
  ReadFileForWritingFastq(string filename, FileType fileType);
  ~ReadFileForWritingFastq();
  string getName();
  void open();
  bool hasRead();
  WritableSeq* getWritableRead();
  void close();
  FileType getFileType();

 private:
  void tryAndGetRead();
  char* _filename;
  bool _isOpen;
  WritableSeq* _waitingRead;
  bool _readIsWaiting;
  ifstream _getReadsFile; // the file object for obtaining reads
  FileType _filetype;
};



class ReadFileForWritingPriceq : public ReadFileForWriting {
 public:
  ReadFileForWritingPriceq(string filename);
  ~ReadFileForWritingPriceq();
  string getName();
  void open();
  bool hasRead();
  WritableSeq* getWritableRead();
  void close();
  FileType getFileType();

 private:
  void tryAndGetRead();
  char* _filename;
  bool _isOpen;
  WritableSeq* _waitingRead;
  bool _readIsWaiting;
  ifstream _getReadsFile; // the file object for obtaining reads
};



#endif

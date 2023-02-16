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

/* A collection of ScoredSeq objects indexed for quick retrieval

This version is optimized for high speed when high fractId is required.
The set of ScoredSeqs that it contains cannot be modified.

The contents of this collection are typed.
 */


#ifndef SCOREDSEQCOLLECTIONBWT_H
#define SCOREDSEQCOLLECTIONBWT_H

#include <omp.h>
#include <set>
#include <map>
#include <vector>
#include "ScoredSeq.h"
#include "Alignment.h"
#include "DynamicProgrammingAligner.h"
#include "AlignmentScoreMatrix.h"
#include "burrowsWheelerUtilities.h"
#include "MatchSeqTest.h"
#include "MatchOffsetTest.h"
using namespace::std;


class ScoredSeqCollectionBwt{

 public:
  // for constructors, note that there would be no advantage to passing in an array

  // these two will be used if ungapped alignments are to be returned
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, bool penalizeEdgeGaps=true);
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap);
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, float fractId, long minOverlap, long maxOverlap);
  // these two will be used if gapped alignments are to be returned
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DyProAlignerFactory * asif);
  ScoredSeqCollectionBwt(set<ScoredSeq*>* inputSeqs, DyProAlignerFactory * asif, long globalMaxOverlap);
  ~ScoredSeqCollectionBwt();

 
  enum GapMode { gapped=0, ungapped=1 };
  enum MinOvlStringency { hardMinOvl=0, softMinOvl=1 };
  enum AlignmentType { semiGlobal=0, fullAlignment=1 };
  enum AlignmentThreadedness{ notThreaded=0, threaded=1 };
  enum MatchesSought{ allMatches=0, bestMatches=1, firstMatch=2 };

  // DEFAULT: edge scaling is enabled
  // that means that matches to smaller words from the edges of query contigs
  // will be sought in order to find matches with small overlaps.
  void enableEdgeScaling();
  void disableEdgeScaling();

  // returns a deep copy of the collection that uses only the
  // shallow versions of the contained objects so that there are
  // no pointers that could be used in some other part of the
  // software (excepting singletons).
  ScoredSeqCollectionBwt * copy();

  // methods for dealing with contents
  bool contains(ScoredSeq* s);
  void getSeqs( set<ScoredSeq*>* );

  // learn about the collection
  float getFractId();
  long getMinOverlap();

  // In the methods below, a pre-created set of alignments is passed in as a reference.
  // It will be modified by having the new matches added to it.

  // these methods return just the best matches
  // (according to score; all equal matches are returned)
  //void getBestFullMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense = '.', MatchSeqTest* seqTest = NULL);

  // these allow a new minOverlap to be set
  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'

  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, MatchSeqTest* seqTest,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getBestMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* seqTest,
		      long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'

  // THREADED
  /*
  void getBestFullMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
			  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);
  void getBestFullMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
			  AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);
  */
  // these allow a new minOverlap to be set
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, MatchSeqTest** seqTestArray,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  void getBestMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray,
		      long* minOvlArray, AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);


  // these just determine if there is a full match (either full query or full target)
  bool hasFullMatch(ScoredSeq* seq, char sense = '.', MatchSeqTest* seqTest = NULL);
  // minOverlap MUST be set
  bool hasMatch(ScoredSeq* seq, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  bool hasMatch(ScoredSeq* seq, char sense, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  bool hasMatch(ScoredSeq* seq, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  bool hasMatch(ScoredSeq* seq, char sense, MatchSeqTest* seqTest, long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);
  // THREADED

  bool* hasFullMatch(long numQueries, ScoredSeq** seqArray,
		     AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);
  bool* hasFullMatch(long numQueries, ScoredSeq** seqArray, char sense,
		     AlignmentThreadedness threadedness, MatchSeqTest** seqTestArray = NULL);

  bool* hasMatch(long numQueries, ScoredSeq** seqArray, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  bool* hasMatch(long numQueries, ScoredSeq** seqArray, char sense, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  bool* hasMatch(long numQueries, ScoredSeq** seqArray, MatchSeqTest** seqTestArray, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);
  bool* hasMatch(long numQueries, ScoredSeq** seqArray, char sense, MatchSeqTest** seqTestArray, long* minOvlArray,
		 AlignmentThreadedness threadedness, MinOvlStringency ovlStringency = hardMinOvl);




  // GETMATCHES

  // and i have equivalents that return all matches
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq, MatchSeqTest* matchTest,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(vector<Alignment*>* matches, ScoredSeq* seq, char sense, MatchSeqTest* matchTest,
		  long minOverlap, MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'
  // these will be threaded
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, MatchSeqTest** matchTestArray,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // both senses are returned
  void getMatches(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, MatchSeqTest** matchTestArray,
		  long* minOvlArray, AlignmentThreadedness threadedness,
		  MinOvlStringency ovlStringency = hardMinOvl);  // sense can be '+' or '-'


  // these two methods deal with the buffering for all the member seqs.  they 
  // don't change the abstract state; they are provided as tools for saving memory
  // without destroying the collection; their converse operations are private and 
  // will be called based on demand (i.e. when matches are requested).
  void unbufferSeqs();
  void unindexSeqs();

  long size(); // number of member ScoredSeqs
  void OK();

  void bufferSeqs();

  // this is only public so that I can test it
  //Alignment* reverseAlignmentSeqA(Alignment* oldAl, ScoredSeq* newSeqA);



 private:
  void constructorHelper(set<ScoredSeq*>* inputSeqs);

  bool _penalizeEdgeGaps;
  static long _defaultMinScore;

  // this should be true during assembly and false during read mapping
  bool _edgeScalingEnabled;

  // this is used in several places
  static long getOverlapFromOffset(long offset, long seqSizeA, long seqSizeB);

  // private versions of these methods; the above methods will forward to them
  void runSearchesConverter(vector<Alignment*>* matches, ScoredSeq* seq, char sense,
			    MatchSeqTest* seqTest, long minOverlap,
			    AlignmentType alType, MatchesSought matchesSought,
			    MinOvlStringency ovlStringency, AlignmentThreadedness threadedness);
  void runSearchesPrivate(long numQueries, vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense,
			  MatchSeqTest** seqTestArray, long* minOvlArray,
			  AlignmentType alType, MatchesSought matchesSought,
			  MinOvlStringency ovlStringency, AlignmentThreadedness threadedness);

  // these methods create an internal representation with a bigger memory
  // footprint and so are called based on demand (i.e. when matches are 
  // requested).
  void indexSeqs();
  void indexSeq(ScoredSeq * seq); // some functional abstraction

  float _fractId;
  float _maxFractMis;
  long _maxOverlap;
  long _minOverlap;
  bool _usingMaxOverlap;

  GapMode _gapMode;
  // a possible alternative to _minOverlap for full alignments; will be
  // set to the length of the shortest sequence or _minOverlap, whichever
  // is longer
  long _minFragmentSize;

  // these are for the gapped alignment mode
  DyProAlignerFactory * _asif;


  // BWT stuff
  long* _sortedSuffixes; // also an array of text length
  long* _bwTransform; // also an array of text length
  long* _occCounts; // also an array of text length
  long* _tableC; // an array of alphabet length
  long _textSize;

  // I NEED A CLASS OF SCOREDSEQ WRAPPERS THAT DOES THE OVERLAP FUNCTION AND KEEPS TRACK
  // OF COORD POSITION IN THE META-TEXT USED FOR BWT MAPPING.
  // used to store the ScoredSeqs so they can be queried
  class ScoredSeqLocus {
  public:
    ScoredSeqLocus();
    ScoredSeqLocus(ScoredSeq* seq, long start, long end, long index);
    bool overlaps(long otherStart, long otherEnd);
    bool contains(long otherStart, long otherEnd);
    ScoredSeq* _seq;
    long _start;
    long _end;
    long _index;
  private:
  };



  // helps enforce the limits of offsets that may be generated from
  // a given query's word match
  class WordQueryWithIndex;
  class WordQuery{
  public:
    virtual ~WordQuery();
    virtual WordQuery* copy() = 0;
    virtual long size() = 0;
    virtual long start() = 0;
    // end is pos+1, good for iterating
    virtual long end() = 0;
    virtual bool acceptableOverlap(long targetStart, ScoredSeq* target) = 0;
    virtual bool notContinuous(WordQuery* otherWq) = 0;
    //virtual WordQueryWithIndex* indexedCopy(long index) = 0;
    virtual long getIndex() = 0;
    virtual void setIndex(long index) = 0;
    // adjusts the start/end/size aspects of the query
    // REQUIRES continuity, MODIFIES self
    virtual void merge(long start, long end) = 0;
  };


  // VERSION 2: WITH INDEX
  // for indexedCopy methods, new indexes replace the old

  class WordQuerySimple : public WordQuery {
  public:
    WordQuerySimple(long start, long size);
    ~WordQuerySimple();
    WordQuerySimple* copy();
    void merge(long start, long end);
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
    long getIndex();
    void setIndex(long index);
  private:
    WordQuerySimple(long start, long end, long size, long index);
    // define the souce seq and position on it of the query word
    long _start;
    long _end;
    long _size;
    long _index;
  };

  class WordQueryNoEdgeLimit : public WordQuery {
  public:
    WordQueryNoEdgeLimit(long start, long size, char sense);
    ~WordQueryNoEdgeLimit();
    WordQueryNoEdgeLimit* copy();
    void merge(long start, long end);
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
    long getIndex();
    void setIndex(long index);
  private:
    WordQueryNoEdgeLimit(long start, long end, long size, char sense, long index);
    // define the souce seq and position on it of the query word
    long _start;
    long _end;
    long _size;
    char _sense;
    long _index;
  };

  class WordQueryFiveEdgeLimit : public WordQuery {
  public:
    WordQueryFiveEdgeLimit(long start, long size, char sense, long edgeLimit);
    ~WordQueryFiveEdgeLimit();
    WordQueryFiveEdgeLimit* copy();
    void merge(long start, long end);
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
    long getIndex();
    void setIndex(long index);
  private:
    WordQueryFiveEdgeLimit(long start, long end, long size, char sense, long edgeLimitX2, long index);
    // define the souce seq and position on it of the query word
     long _start;
    long _end;
    long _size;
    char _sense;
    long _index;
    // define the distance to the edge of a contig within
    // which the word match must be found (or its relevance);
    // stored as X2 to speed up operations (calculated by external
    // constructor, copied by private constructor
    long _edgeLimitX2;
  };
  class WordQueryThreeEdgeLimit : public WordQuery {
  public:
    WordQueryThreeEdgeLimit(long start, long size, char sense, long edgeLimit);
    ~WordQueryThreeEdgeLimit();
    WordQueryThreeEdgeLimit* copy();
    void merge(long start, long end);
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
    long getIndex();
    void setIndex(long index);
  private:
    WordQueryThreeEdgeLimit(long start, long end, long size, char sense, long edgeLimitX2, long index);
    // define the souce seq and position on it of the query word
    long _start;
    long _end;
    long _size;
    char _sense;
    long _index;
    // define the distance to the edge of a contig within
    // which the word match must be found (or its relevance)
    // stored as X2 to speed up operations (calculated by external
    // constructor, copied by private constructor
    long _edgeLimitX2;
  };

  class WordQueryFullTarget : public WordQuery {
  public:
    WordQueryFullTarget(long sourceSize, long start, long size, char sense, float fractId, long maxTargetSize);
    ~WordQueryFullTarget();
    WordQueryFullTarget* copy();
    void merge(long start, long end);
    bool acceptableOverlap(long targetStart, ScoredSeq* target);
    long size();
    long start();
    long end();
    bool notContinuous(WordQuery* otherWq);
    long getIndex();
    void setIndex(long index);
  private:
    WordQueryFullTarget(long sourceSize, long start, long end, long size, char sense, float fractId, long maxTargetSize, long index);
    // define the souce seq and position on it of the query word
    long _sourceSize;
    long _start;
    long _end;
    long _size;
    char _sense;
    long _index;
    float _fractId;
    long _maxTargetSize;
  };




  // used to deal with redundant word queries
  class WordQueryBinOrganizer {
  public:
    WordQueryBinOrganizer(vector<WordQuery*>* wordQueries);
    ~WordQueryBinOrganizer();

    void deleteWithWordQueries();
    // allows the object to be re-used without being re-created
    void reset();

    long numIndexes();
    long numValidQueries();
    bool queryStillValid(long qIndex);
    WordQuery* getQuery(long qIndex);
    WordQuery* getIndexedQueryCopy(long qIndex);

    void removeQuery(long qIndex);
    void removeQuery(long qIndex, long currentBlock);
    void removeSupersetQueries(long qIndex);

    // last valid entry is followed by a NULL entry
    WordQuery** getQueries();
    long numBlocks();
    bool blockHasQueries(long block);

    long blockStart(long block);
    long blockEnd(long block);
    bool blockHasQuery(long block, long qIndex);

    //void getBlockQueryIndexes(long block, set<long>* indexes);
    // these two methods replace the one above
    long getBlockQueryIndexCount(long block);
    long* getBlockQueryIndexes(long block);

  private:
    long* viewBlockQueryIndexes(long block);
    WordQuery** _queries;
    bool* _queryIsValid;
    long _numIndexes;
    long _numValidQueries;
    long _numBlocks;
    // this next one is for optimizing the reset function; loads from the back!!
    long* _invalidQueries;
    // for these, the N value in End equals the N+1 value in start; pre-computing
    // is faster than adding one every time the boundaries are checked
    long* _blockToStart;
    long* _blockToEnd;
    long* _startToBlock;

    // sparse: for all of the qIndexes in the bin, there is a value showing
    // where they appear in the _blockToQisInBin array
    //long** _blockToQiBinIndex;

    // organizes the data by block, each of which is defined by a segment of one big array
    long* _blockToBinStart;
    long* _blockToBinEnd;
    // all the blocks in one big array
    // useful elements bounded by [_blockToBinStart,_blockToBinEnd)
    // each sub-array section has the qIndexes of the queries in that bin in an arbitrary order;
    long* _blockToAllQisInBin;
    long* _blockToAllQiBinIndex;
  };


  class OffsetAndQueryCarrier {
  public:
    OffsetAndQueryCarrier(long offset, long numQueries);
    ~OffsetAndQueryCarrier();
    void deleteWithQueries();
    long _offset;
    long _numQueries;
    WordQuery** _queries;
    OffsetAndQueryCarrier* copy();
  };

  // for keeping track of offsets and passing them around
  class OffsetTrackerImmutable;
  class OffsetTracker {
  public:
    virtual ~OffsetTracker();
    virtual long numOffsets() = 0;
    virtual long maxOffsetCount() = 0;
    virtual long maxCountOffset() = 0;
    virtual ScoredSeq* targetContig() = 0;
    virtual OffsetTrackerImmutable* immutableCopy(bool includeQueries) = 0;
    // getting offsets with queries included is only allowed if this is TRUE
    virtual bool queriesIncluded() = 0;
    virtual OffsetAndQueryCarrier** getOffsets(bool includeQueries) = 0;
  };

  // two implemeting classes; the second saves memory by deleting information
  // that is no longer needed by the first if it will no longer be mutated
  class OffsetTrackerMutable : public OffsetTracker {
  public:
    OffsetTrackerMutable(ScoredSeq* target);
    ~OffsetTrackerMutable();
    void addOffset(long offset, WordQuery* seedingQuery, long normFactor=1);
    long maxOffsetCount();
    long maxCountOffset();
    long numOffsets();
    OffsetAndQueryCarrier** getOffsets(bool includeQueries);
    ScoredSeq* targetContig();
    OffsetTrackerImmutable* immutableCopy(bool includeQueries);
    bool queriesIncluded();
  private:
    class OffsetFields {
    public:
      OffsetFields(float count, WordQuery* maxWord);
      ~OffsetFields();
      void addOffset(float localAddition, WordQuery* seedingQuery);
      float _count;
      WordQuery* _maxWord;
      vector<WordQuery*> _otherWords;
    };
    map<long,OffsetFields*> _offsetToFields;
    long _maxCountOffset;
    float _maxOffsetCount;
    ScoredSeq* _targetContig;
    // for sorting word queries by start coord
    struct SortWqByStart {
      bool operator() (WordQuery* wqA, WordQuery* wqB);
    };
  };

  class OffsetTrackerImmutable : public OffsetTracker {
  public:
    // NOTE: the array is rep-exposed, so it and its component offsets should already be copies;
    // they will be deleted along with this object!!!
    OffsetTrackerImmutable(ScoredSeq* targetContig, long numOffsets, OffsetAndQueryCarrier** offsets,
			   long maxCountOffset, long maxOffsetCount, bool qIncluded);
    ~OffsetTrackerImmutable();
    long maxOffsetCount();
    long maxCountOffset();
    long numOffsets();
    OffsetAndQueryCarrier** getOffsets(bool includeQueries);
    ScoredSeq* targetContig();
    bool queriesIncluded();
    OffsetTrackerImmutable* immutableCopy(bool includeQueries);
  private:
    bool _queriesIncluded;
    long _maxCountOffset;
    float _maxOffsetCount;
    ScoredSeq* _targetContig;
    long _numOffsets;
    OffsetAndQueryCarrier** _offsets;
  };




  // pure abstract class that will replace the one above; I will change
  // the name when I am done with it
  class AlignmentMaker {
  public:
    virtual ~AlignmentMaker();
    virtual Alignment* makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense, OffsetAndQueryCarrier* offsetCarrier) = 0;
    virtual long getScore(Alignment* al) = 0;
    virtual long getMinScore() = 0;
    virtual void setMinScore(long minScore) = 0;
    virtual AlignmentScoreMatrix* getScoreMatrix() = 0;
    // returns an ASM that has the scores (0,1,1,1) - prevents need for "new" statment + deletion in static methods
    // REQUIRES that the accessing object is one of the static methods below.  "accessed" object is rep-exposed!
    //virtual AlignmentScoreMatrix* accessMisCountAsm() = 0;
    //static bool makeAlignmentHelperPlus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis, long alLength);
    static bool makeAlignmentHelperPlus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis,
					long exBlockCount, long* examineStarts, long* examineEnds);
    //static bool makeAlignmentHelperMinus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis, long alLength);
    static bool makeAlignmentHelperMinus(ScoredSeq* seqA, long aStart, ScoredSeq* seqB, long bStart, long maxNumMis,
					 long exBlockCount, long* examineStarts, long* examineEnds);
    static inline int scoreNucMisMatch(char nA, char nB);
  };
  // implementing classes
  class AlMakerNoEdge : public AlignmentMaker {
  public:
    //AlMakerNoEdge();
    AlMakerNoEdge(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis);
    ~AlMakerNoEdge();
    Alignment* makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense, OffsetAndQueryCarrier* offsetCarrier);
    long getScore(Alignment* al);
    long getMinScore();
    void setMinScore(long minScore);
    AlignmentScoreMatrix* getScoreMatrix();
    //AlignmentScoreMatrix* accessMisCountAsm();
  private:
    AlignmentScoreMatrix* _scoreMatrix;
    //AlignmentScoreMatrix* _misCountAsm;
    long _minScore;
    float _maxFractMis;
  };
  class AlMakerPenalizeEdgeA : public AlignmentMaker {
  public:
    //AlMakerPenalizeEdgeA();
    AlMakerPenalizeEdgeA(AlignmentScoreMatrix* scoreMatrix, long minScore, float maxFractMis);
    ~AlMakerPenalizeEdgeA();
    Alignment* makeAlignment(ScoredSeq* seqA, ScoredSeq* seqB, char sense, OffsetAndQueryCarrier* offsetCarrier);
    long getScore(Alignment* al);
    long getMinScore();
    void setMinScore(long minScore);
    AlignmentScoreMatrix* getScoreMatrix();
    //AlignmentScoreMatrix* accessMisCountAsm();
  private:
    AlignmentScoreMatrix* _scoreMatrix;
    //AlignmentScoreMatrix* _misCountAsm;
    long _minScore;
    float _maxFractMis;
  };



  // I will use this to find the ScoredSeqs to which queries were mapped quickly
  // indexes are 0 -> total text length / _binSize (continuous)
  // this will be an array of ScoredSeqLocus* sets.
  vector<ScoredSeqLocus*>* _binToSeqs;

  // an array of loci
  ScoredSeqLocus** _orderedContigLoci;
  long _numContigLoci;
  // for speeding up use of contig loci; should be one element per thread (up to max count)
  bool** _contigIndexWasHit;
  long** _localContigIndex;
  bool* _contigIndexWasHitMemory;
  long* _localContigIndexMemory;

  long _binSize;
  long _maxBin; // for safety checks
  long getOverlapContigIndex(long start, long end);
  AlignmentScoreMatrix* _matchCountAsm;
  AlignmentScoreMatrix* _misCountAsm; // gaps count as mismatches here (terminal gaps will exist)
  AlignmentScoreMatrix* _scoreAsm;

  // helpers for functions above that use private datatypes
  void runAlignmentsHelper(long numQueries, OffsetTracker*** contigToOffsetsArray,
			   vector<Alignment*>** matchesArray, ScoredSeq** seqArray, char sense, GapMode gapMode,
			   MatchSeqTest** seqTestArray, long* minOvlArray, long* bestScoreArray,
			   MatchesSought matchesSought, bool* seekMatch, AlignmentThreadedness threadedness);
  AlignmentMaker* makeAlMaker(GapMode gapMode, long bestScoreSoFar);
  void dealWithMatches(vector<Alignment*>* tempMatches, vector<Alignment*>* realMatches,
		       bool getAllMatches, bool justFirstMatch,
		       long index, long* bestScoreArray, long newBestScore, bool* matchSoughtArray);


  // this is used in many contexts; its calling is controlled by the hasMatch, getBestMatches, and getMatches
  // wrapper methods.  note that the default min score so far is zero, but i am going to force that to be
  // input so that i don't fuck up and forget to include it somewhere
  // NOTE: if onlyBestAlignment==false, then the input bestScoreSoFar will just be spit back without being updated
  long makeAlignmentHelper(OffsetTracker* tracker, ScoredSeq* seq, vector<Alignment*>* matches, long minOverlap,
			   AlignmentMaker* alMaker, ScoredSeq* seqFlip, char sense, long bestScoreSoFar,
			   bool onlyBestAlignment, GapMode gapMode, MatchesSought matchesSought);


  // used to define the set of subtext strings for which perfect matches will be sought
  static int _specialCaseLenDenom[6];
  static int _specialCaseStepDenom[6];
  static int _specialCaseWindowNum[170];
  // this next one is derivative of the first two, but I'll compute manually ahead of time
  static int _specialCaseStartDenom[6];

  WordQueryBinOrganizer* getSubtextWordQueries(ScoredSeq* querySeq, char* seqSeq, char sense, long minOverlap);
  WordQueryBinOrganizer* getSubtextWordQueriesFullAlignment(ScoredSeq* querySeq, char* seqSeq, char sense, long minTargetLength);
  void getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
				   long fragmentLength);
  void getSubtextWordQueriesHelper(ScoredSeq* querySeq, vector<WordQuery*>* wordQuerySet, char sense,
				   long blockStartPos, long blockEndPos, bool useLimits);
  void subtextWordQueriesFilterNs(ScoredSeq* querySeq, char* seqSeq, WordQueryBinOrganizer* wqbo);


  long getNumWindowsHelper(long numAllowedMis);

  // a helper for getMatches/getBestMatches
  // the startQn makes threading easier
  void getOffsets(long numQ, long startQn, ScoredSeq** seqs, char sense, long* minOverlaps, MatchSeqTest** matchTests,
		  bool* lookForMatches, AlignmentType alType, OffsetTracker*** resultArray);

  // these are just because "getOffsets" is too massive on its own to keep track of what's going on easily
  void getOffsetsSub1(long currentIndex, OffsetTrackerMutable** hitArrayMutable, OffsetTracker** sortedHitWithOffsets, 
		      bool queriesNeeded, long* hitIndexes, bool* localIndexWasHit);
  long getOffsetsSub2(long currentIndex, WordQuery* wordQuery, ScoredSeq* seq, long* matchArray, OffsetTrackerMutable** hitArrayMutable, 
		      char* hitArrayMutableMemory, long* hitIndexes, long* localIndex, bool* localIndexWasHit, MatchSeqTest* matchTest);

  // for sorting word queries by length, shortest to longest
  struct SortWqByLength {
    bool operator() (WordQuery* wqA, WordQuery* wqB);
  };

  // these keep track of the total content, binned by whether or not the
  // work that needs to be done to get the collection ready to use has
  // already been performed.
  set<ScoredSeq*> _allSeqs;

};

#endif

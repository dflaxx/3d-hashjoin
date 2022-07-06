#pragma once

#include "util/standard_includes.hh"
#include "util/aggregate.hh"

#include <ostream>

// hash table bucket statistics
struct HtBucketStatistics {
  size_t _bucketIndex;
  size_t _numEntries;     // #stored key/data pairs
  size_t _chainLen;       // collision chain length

  std::string toCsvString() const;
};

// hash table statistics
struct HtStatistics {
  size_t _numBuckets;
  size_t _numEmptyBuckets;
  size_t _numEntries;           // #stored key/data pairs
  size_t _numDistinctKeys;      // <= _numEntries
  Aggregate<size_t> _collisionChainLen;
  Aggregate<size_t> _collisionChainLenNonempty;
  Aggregate<size_t> _numDistinctKeysPerBucket;
  Aggregate<size_t> _numDistinctKeysPerNonemptyBucket;

  std::vector<HtBucketStatistics> _bucketStats;

  inline double numEntriesPerKey() const { return (_numEntries + 0.0) / _numDistinctKeys; }

  inline double fracEmptyBuckets() const { return (_numEmptyBuckets + 0.0) / _numBuckets; }

  //inline HtStatistics()
  //  : _numBuckets(), _numEmptyBuckets(), _numEntries(), _numDistinctKeys(),
  //    _numEntriesPerBucket(), _numEntriesPerNonemptyBucket(),
  //    _numDistinctKeysPerBucket(), _numDistinctKeysPerNonemptyBucket() {}
  
  inline HtStatistics() { reset(); }

  void reset();

  void print(std::ostream& os = std::cout) const;

  std::string toCsvString() const;
  static std::string toCsvStringHeader();

  std::string bucketStatToCsvString() const;
  static std::string bucketStatToCsvStringHeader();

  void printCsv(std::ostream& os = std::cout) const;
  static void printCsvHeader(std::ostream& os = std::cout);

};





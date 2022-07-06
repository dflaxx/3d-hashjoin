#include "ht_statistics.hh"
#include <sstream>

void
HtStatistics::reset() {
  _numBuckets = 0;
  _numEmptyBuckets = 0;
  _numEntries = 0;
  _numDistinctKeys = 0;
  _collisionChainLen.init();
  _collisionChainLenNonempty.init();
  _numDistinctKeysPerBucket.init();
  _numDistinctKeysPerNonemptyBucket.init();
}

void
HtStatistics::print(std::ostream& os) const {
  os << "#buckets                  = " << _numBuckets << "\n";
  os << "#empty buckets            = " << _numEmptyBuckets << "\n";
  os << "#entries                  = " << _numEntries << "\n";
  os << "#distinct keys            = " << _numDistinctKeys << "\n";
  os << "cc length:                  "
    << _collisionChainLen.avg() << " | "
    << _collisionChainLen.min() << " | " 
    << _collisionChainLen.max() << "\n";
  os << "cc length nonempty:         "
    << _collisionChainLenNonempty.avg() << " | "
    << _collisionChainLenNonempty.min() << " | " 
    << _collisionChainLenNonempty.max() << "\n";
  //os << "#keys/bucket              : "
  //  << _numDistinctKeysPerBucket.avg() << " | "
  //  << _numDistinctKeysPerBucket.min() << " | " 
  //  << _numDistinctKeysPerBucket.max() << "\n";
  //os << "#keys/non-empty bucket    : "
  //  << _numDistinctKeysPerNonemptyBucket.avg() << " | "
  //  << _numDistinctKeysPerNonemptyBucket.min() << " | " 
  //  << _numDistinctKeysPerNonemptyBucket.max() << "\n";
}

std::string
HtStatistics::toCsvString() const {
  std::stringstream res;
  res
    << _numBuckets << ";"
    << _numEmptyBuckets << ";"
    << _numEntries << ";"
    << _numDistinctKeys << ";";
  res
    << _collisionChainLen.avg() << ";"
    << _collisionChainLen.min() << ";" 
    << _collisionChainLen.max() << ";";
  res
    << _collisionChainLenNonempty.avg() << ";"
    << _collisionChainLenNonempty.min() << ";" 
    << _collisionChainLenNonempty.max() << ";";
  //res
  //  << _numDistinctKeysPerBucket.avg() << ";"
  //  << _numDistinctKeysPerBucket.min() << ";" 
  //  << _numDistinctKeysPerBucket.max() << ";";
  //res
  //  << _numDistinctKeysPerNonemptyBucket.avg() << ";"
  //  << _numDistinctKeysPerNonemptyBucket.min() << ";" 
  //  << _numDistinctKeysPerNonemptyBucket.max() << ";";
  return res.str();
}

std::string
HtStatistics::toCsvStringHeader() {
  std::stringstream res;
  res
    << "#buckets;#empty_buckets;#entries;#distinct_keys;"
    << "#e/b_avg;#e/b_min;#e/b_max;"
    << "#e/neb_avg;#e/neb_min;#e/neb_max;";
    //<< "#k/b_avg;#k/b_min;#k/b_max;"
    //<< "#k/neb_avg;#k/neb_min;#k/neb_max;";
  return res.str();
}

std::string
HtStatistics::bucketStatToCsvString() const {
  std::stringstream res;
  res << "bucket_idx;num_entries;chain_len\n";
  for (const auto& bs : _bucketStats) {
    res << bs._bucketIndex << ";" << bs._numEntries << ";" << bs._chainLen << "\n";
  }
  return res.str();
}

std::string
HtStatistics::bucketStatToCsvStringHeader() {
  return "bucket_idx;num_entries;chain_len";

}

void
HtStatistics::printCsv(std::ostream& os) const {
  os << toCsvString();
}

void
HtStatistics::printCsvHeader(std::ostream& os) {
  os << toCsvStringHeader();
}

std::string
HtBucketStatistics::toCsvString() const {
  std::stringstream res;
  res << _bucketIndex << ";" << _numEntries << ";" << _chainLen;
  return res.str();
}



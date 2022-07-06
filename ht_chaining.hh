#pragma once

#include "util/standard_includes.hh"

#include <cmath>
#include <functional>
#include <unordered_set>

#include "util/math.hh"
#include "util/output_helpers.hh"
#include "util/reservoir.hh"

#include "concepts.hh"
#include "ht_iterators.hh"
#include "ht_statistics.hh"

#define PADDING_BYTESIZE 0   // this is just for testing

/*
 * A templated hash table that (conceptually) stores (key, data) pairs
 * and resolves collisions by chaining.
 * Data can, e.g., be 1) a pointer to a tuple (tuple_t*) or 2) a row id of a tuple (uint).
 *
 * The hash function and the key comparison function are given as template paramters.
 *
 * Template parameters:
 * - Tdata: the data stored in the hash table, e.g., a tuple pointer or a tuple identifier
 * - Thashfun: a class/type implementing the hash function that maps keys to hash table buckets.
 *   Needs types input_t (the key type) and output_t (an unsigned int type for the bucket index).
 *   Also needs a static member function eval(input_t*) -> output_t.
 * - Tcontenteqfun: a class/type implementing the key comparison function for inserting.
 *   Needs types left_t and right_t and a static member function eval(left_t*, right_t*) -> bool.
 *   TODO: not needed, remove
 *
 * This implementation uses C++20 concepts to enforce compile-time predicates on the template parameters.
 * Concepts are specified in concepts.hh.
 */
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
class HtChaining1 {
  // some type asserts
  static_assert(std::same_as<Tdata, typename Thashfun::input_t>, "Thashfun::input_t does not match Tdata");
  static_assert( std::same_as<Tdata, typename Tcontenteqfun::left_t> &&
                 std::same_as<Tdata, typename Tcontenteqfun::right_t>,
                 "Tcontenteqfun::left_t/right_t do not match Tdata" );

  public:
    struct Node;  // fwd decl of inner struct
    friend NodeIterator<Node, false>;  // iterator class may access this class's private members
    friend NodeIterator<Node, true>;

  public:
    using data_t      = Tdata;
    using hashfun_t   = Thashfun;
    using hashvalue_t = typename hashfun_t::output_t;
    using hashdir_t   = std::vector<Node>;
    using eqfun_t     = Tcontenteqfun;

    using hashdir_iterator = typename hashdir_t::iterator;
    using const_hashdir_iterator = typename hashdir_t::const_iterator;
    using node_iterator       = NodeIterator<Node, false>;
    using const_node_iterator = NodeIterator<Node, true>;

    using stats_t = HtStatistics;

    inline static Node* EMPTY_ENTRY = reinterpret_cast<Node*>(0x1);

  public:
    /* Hash directory and collision chain node */
    struct Node {
      Node*       _next;       // next node in the collision chain
      data_t*     _data;       // pointer to tuple
      hashvalue_t _hashvalue;  // = hashfun_t::eval(t), for fast comparisons, e.g. strings
#if PADDING_BYTESIZE > 0
#warning Compiling HtChaining1::Node with memory padding
      std::byte   __padding[PADDING_BYTESIZE];   // force node to be 32 B
#endif

      inline Node(data_t* aData)
        : _next(nullptr), _data(aData), _hashvalue() {}
      inline Node() : _next(EMPTY_ENTRY), _data(), _hashvalue() {}

      inline void init(data_t* aData, const hashvalue_t aHashvalue, Node* aNext) {
        _data = aData;
        _next = aNext;
        _hashvalue = aHashvalue;
      }
      inline void init(data_t* aData, const hashvalue_t aHashvalue) {
        init(aData, aHashvalue, nullptr);
      }

      inline data_t*     data()      const { return _data; }
      inline Node*       next()      const { return isEmpty() ? nullptr : _next; }
      inline hashvalue_t hashvalue() const { return _hashvalue; }
      inline bool        isEmpty()   const { return _next == EMPTY_ENTRY; }
      inline bool        hasNext()   const { return !isEmpty() && (_next != nullptr); }

      inline void        clear()           { _next = EMPTY_ENTRY; }  // !! danger of dangling ptrs

      inline node_iterator begin() { return node_iterator(this); }
      inline node_iterator end()   { return node_iterator(nullptr); }
      inline const_node_iterator cbegin() const { return const_node_iterator(this); }
      inline const_node_iterator cend()   const { return const_node_iterator(nullptr); }
    };

  public:
    inline HtChaining1(const size_t aNumBuckets, const uint32_t aReservoirLog2ChunkSize)
      : _hashdir(aNumBuckets, Node()), _reservoir(aReservoirLog2ChunkSize), _size() {}

  public:
    inline size_t      numBuckets()              const { return _hashdir.size(); }
    inline hashvalue_t hash(const data_t* aData) const { return hashfun_t::eval(aData); }
    inline size_t      size()                    const { return _size; }
    inline size_t      getRsvSize()              const { return _reservoir.cardinality(); }

    size_t             memoryConsupmtion()       const;
    size_t             memoryConsupmtionDir()    const;
    size_t             memoryConsupmtionChains() const;

  public:
    void insert(data_t* aData);
    void print(uint64_t aIndent, uint64_t aTab, const bool aSkipEmpty, std::ostream& os = std::cout) const;

    template <typename Tprobedata, alg_hashfun_c Tprobehashfun>
    const_node_iterator findDirEntryByOther(const Tprobedata* aProbeTuple) const;

    void clear();

  public:
    // statistics
    stats_t makeStatistics() const;

    // get iterator of _hashdir_t aka std::vector<Node>
    inline typename hashdir_t::iterator begin()              { return _hashdir.begin(); }
    inline typename hashdir_t::iterator end()                { return _hashdir.end(); }
    inline typename hashdir_t::const_iterator cbegin() const { return _hashdir.cbegin(); }
    inline typename hashdir_t::const_iterator cend()   const { return _hashdir.cend(); }

  private:
    inline       size_t getDirIndex(const data_t* aData) const { return hash(aData) % numBuckets(); }
    inline       size_t getDirIndex(const hashvalue_t h) const { return h % numBuckets(); }
    inline const Node*  getDirEntry(const data_t* aData) const { return &_hashdir.at(getDirIndex(aData)); }
    inline       Node*  getDirEntry(const data_t* aData)       { return &_hashdir.at(getDirIndex(aData)); }
    inline const Node*  getDirEntry(const hashvalue_t h) const { return &_hashdir.at(getDirIndex(h)); }

    // given data_t*, return its directory index, hash value and ptr to directory entry (as a tuple)
    inline std::tuple<size_t, hashvalue_t, Node*>
    getDirIdxHashEntry(const data_t* aData) {
      return { getDirIndex(aData), hash(aData), getDirEntry(aData) };
    }

    inline bool isEqual(const data_t* x, const data_t* y) { return eqfun_t::eval(x, y); }  // TODO: not needed, remove

  private:
    hashdir_t       _hashdir;
    Reservoir<Node> _reservoir;  // reservoir stores objects (here: Nodes) in chunks of consecutive
                                 // memory with constant memory addresses (unlike vectors)
    size_t          _size;       // number of data items stored
};

// member functions returning the memory consumption in bytes
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtion() const {
  return memoryConsupmtionDir() + memoryConsupmtionChains();
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtionDir() const {
  return numBuckets() * sizeof(Node);
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtionChains() const {
  return getRsvSize() * sizeof(Node);
}


// insert a new data item
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::insert(data_t* aData) {
  auto [lDirIdx, lHashvalue, lDirEntry] = getDirIdxHashEntry(aData);
  if (lDirEntry->isEmpty()) {
    // bucket is empty -> insert into directory
    assert(nullptr == lDirEntry->next());
    lDirEntry->init(aData, lHashvalue);
  } else {
    // bucket has at least one entry -> insert into collision chain (head)
    Node* lNewNode = _reservoir.get_new_entry();
    lNewNode->init(aData, lHashvalue, lDirEntry->next());
    lDirEntry->_next = lNewNode;
  }
  ++_size;
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::print(uint64_t aIndent, uint64_t aTab,
                                                   const bool aSkipEmpty,
                                                   std::ostream& os) const {
  using df::infra::indent;
  using df::infra::number_of_digits;
  const size_t lDirDigits = std::log10(numBuckets());
  const size_t lHashvalDigits = number_of_digits(std::numeric_limits<hashvalue_t>::max());
  for (size_t i = 0; i < numBuckets(); ++i) {
    const Node& lDirEntry = _hashdir.at(i);
    if (lDirEntry.isEmpty() && aSkipEmpty) { continue; }
    os << indent().lvl(aTab);
    os << std::setw(lDirDigits) << i << ": ";

    if (lDirEntry.isEmpty()) {
      os << "(empty)\n";
      continue;
    }
    os << "\n";

    // lDirEntry's k/d pair
    os
      << indent().marg(aIndent).lvl(aTab)
      << std::setw(lHashvalDigits) << lDirEntry.hashvalue() <<  ": "
      << lDirEntry.data() << "\n";

    // go through collision chain entries
    for (const Node* it = lDirEntry.next(); it != nullptr; it = it->next()) {
      os <<
        indent().marg(aIndent).lvl(aTab)
        << std::setw(lHashvalDigits) << it->hashvalue() <<  ": "
        << it->data() << "\n";
    }

  }
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
template <typename Tprobedata, alg_hashfun_c Tprobehashfun>
auto
HtChaining1<Tdata, Thashfun, Tcontenteqfun>
::findDirEntryByOther(const Tprobedata* aProbeTuple) const
-> const_node_iterator {
  using other_hashfun_t = Tprobehashfun;
  static_assert(std::is_same_v<hashvalue_t, typename other_hashfun_t::output_t>);

  hashvalue_t lProbeHashValue = other_hashfun_t::eval(aProbeTuple);
  const_node_iterator lIt = getDirEntry(lProbeHashValue);
  return lIt;
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::clear() {
  _reservoir.erase();  // empty reservoir and release memory
  // set hash dir buckets to empty
  for (size_t i = 0; i < numBuckets(); ++i) {
    _hashdir.at(i).clear();
  }
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
auto
HtChaining1<Tdata, Thashfun, Tcontenteqfun>::makeStatistics() const
-> stats_t {
  stats_t s;
  s._numBuckets = numBuckets();
  s._numEntries = size();
  std::unordered_set<key_t> lDistinctKeys;

  // iterate buckets
  for (const_hashdir_iterator dir_it = cbegin(); dir_it != cend(); ++dir_it) {
    if (dir_it->isEmpty()) {
      ++s._numEmptyBuckets;
      s._collisionChainLen.step(0);
      continue;
    }
    // iterate collision chain
    uint32_t lCurrChainLen = 0;  // len of current bucket's collision chain (incl dir entry)
    // collision chain length also counts the dir entry, i.e. non-empty buckets have len >= 1
    for (const_node_iterator node_it = dir_it->cbegin(); node_it.valid(); ++node_it) {
      assert(!node_it->isEmpty());
      ++lCurrChainLen;
      lDistinctKeys.insert(node_it->_hashvalue);  // !! cannot access key, must use hash!
                                                  //    therefore, unaccuracy possible!
    }
    s._collisionChainLen.step(lCurrChainLen);
    s._collisionChainLenNonempty.step(lCurrChainLen);
  }

  s._numDistinctKeys = lDistinctKeys.size();

  return s;
}



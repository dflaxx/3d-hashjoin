#pragma once

#include "util/standard_includes.hh"

#include <cmath>
#include <concepts>
#include <functional>
#include <ios>
#include <type_traits>
#include <utility>
#include <variant>

#include "util/reservoir.hh"

#include "concepts.hh"
#include "ht_iterators.hh"
#include "ht_statistics.hh"


// ## HtNested1

/*
 * A templated nested hash table that (conceptually) stores (key, data) pairs.
 * data can, e.g., be 1) a pointer to a tuple (tuple_t*) or 2) a row id of a tuple (uint).
 *
 * Collision resolution by chaining. Collision chains are "partitioned" by the key:
 * Two distict keys k1, k2 that hash to the same bucket form a collision chain (linked list).
 * Each distinct key k results in its own "collision sub-chain":
 * 
 *     dir
 *     +--+    +--+    +--+
 *     |  + -> |k1+ -> |k2+ -> ...
 *     +--+    ++-+    ++-+
 *     |  |     |       |
 *     +--+    +v-+    +v-+
 *     |  |    |d1|    |d3|
 *     +--+    ++-+    +--+
 *     |  |     |
 *     +--+    +v-+ 
 *             |d2| 
 *             +--+ 
 *
 * k1, k2 are MainNodes; d1, d2 and d3 are SubNodes.
 * The figure is not entirely correct:
 * Each hash directory entry also contains a key/data pair and has its own SubNode chain.
 *
 * There are two types of nodes: MainNodes and SubNodes.
 * - MainNode: hash directory and main collision chain entry,
 *             represents a distinct key,
 *             contains one key/data pair,
 *             is the "head" of a chain of SubNodes,
 *             points to the next MainNode in the same hash bucket (linked list),
 *             points to its own list of SubNodes (linked list).
 * - SubNode:  chain entries belonging to one MainNode, i.e. one distinct key,
 *             contains one data element,
 *             points to the next SubNode of the same chain (linked list)
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
 *
 * This implementation uses C++20 concepts to enforce compile-time predicates on the template parameters.
 * Concepts are specified in concepts.hh.
 */
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
class HtNested1 {
  // some type asserts
  static_assert(std::same_as<Tdata, typename Thashfun::input_t>, "Thashfun::input_t does not match Tdata");
  static_assert(std::same_as<Tdata, typename Tcontenteqfun::left_t> &&
                std::same_as<Tdata, typename Tcontenteqfun::right_t>,
                "Tcontenteqfun::left_t/right_t do not match Tdata" );
  
  public:
    // forward decl of inner structs
    struct MainNode;
    struct SubNode;
    // iterator class may access this class's private members
    friend NodeIterator<MainNode, false>;  // false = non-const
    friend NodeIterator<MainNode, true>;   // true  = const
    friend NodeIterator<SubNode, false>;
    friend NodeIterator<SubNode, true>;

  public:
    using data_t      = Tdata;
    using hashfun_t   = Thashfun;
    using hashvalue_t = typename hashfun_t::output_t;
    using hashdir_t   = std::vector<MainNode>;
    using eqfun_t     = Tcontenteqfun;

    using hashdir_iterator = typename hashdir_t::iterator;
    using const_hashdir_iterator = typename hashdir_t::const_iterator;
    using main_node_iterator = NodeIterator<MainNode, false>;
    using const_main_node_iterator = NodeIterator<MainNode, true>;
    using sub_node_iterator = NodeIterator<SubNode, false>;
    using const_sub_node_iterator = NodeIterator<SubNode, true>;
    
    using stats_t = HtStatistics;

    // using the least-significant bit of main node and sub node pointers to indicate their state.
    inline static MainNode* EMPTY_ENTRY           = reinterpret_cast<MainNode*>(0x1);
    inline static SubNode*  UNINITIALIZED_SUBNODE = reinterpret_cast<SubNode*>(0x1);

  public:
    /* Hash directory + main collision chain entry, represents one distinct value */
    struct MainNode {
      /*
       * MainNode States
       *                   (1)          (2)              (3)              (4)
       * _next             0x1          nullptr          nullptr          valid pointer
       * _subchain_head    nullptr      nullptr          valid pointer    valid pointer
       * _data             nullptr      valid pointer    valid pointer    valid pointer
       * _hashvalue        NA/0         valid value      valid value      valid value
       *
       * (1) empty directory entry
       * (2) directory entry or last collision chain entry with a single stored value
       * (3) directory entry or last collision chain entry with a stored value and a subchain
       * (4) directory or collision chain entry (not the last one) with a stored value and a subchain
       *
       * Note: 0x1 == EMPTY_ENTRY
       */
      MainNode*   _next;           // next main node in the collision chain
      SubNode*    _subchain_head;  // sub-chain head
      data_t*     _data;           // pointer to first element (tuple)  // TODO const?
      hashvalue_t _hashvalue;      // = hashfun_t::eval(t), for fast comparisons, e.g. strings

      inline MainNode(data_t* aDataPtr)
        : _next(nullptr), _subchain_head(nullptr), _data(aDataPtr), _hashvalue(0) {}

      inline MainNode()
        : _next(EMPTY_ENTRY), _subchain_head(nullptr), _data(nullptr), _hashvalue(0) {}

      inline const data_t*     data()  const { return _data; }
      inline       hashvalue_t hashvalue() const { return _hashvalue; }
      inline       MainNode*   next()  const { return isEmpty() ? nullptr : _next; }
      inline       SubNode*    child() const { return _subchain_head; }

      inline bool isEmpty()  const { return _next == EMPTY_ENTRY; }
      inline bool hasNext()  const { return !isEmpty() && (next() != nullptr); }
      inline bool hasChild() const { return child() != nullptr; }

      inline void init(data_t* aData, const hashvalue_t aHashvalue) {  // required because of Reservoir
        _next = nullptr;
        _subchain_head = nullptr;
        _data = aData;
        _hashvalue = aHashvalue;
      }

      inline void        clear()           { _next = EMPTY_ENTRY; }  // !! danger of dangling ptrs

      inline main_node_iterator       begin()        { return main_node_iterator(this); }
      inline main_node_iterator       end()          { return main_node_iterator(nullptr); }
      inline const_main_node_iterator cbegin() const { return const_main_node_iterator(this); }
      inline const_main_node_iterator cend()   const { return const_main_node_iterator(nullptr); }
    };

    /* Chain below one MainNode entry, belongs to one distinct value */
    struct SubNode {
      SubNode*   _next;  // next node in the linked list of SubNodes
      data_t*    _data;

      inline SubNode() : _next(UNINITIALIZED_SUBNODE), _data(nullptr) {}
      inline SubNode(data_t* aData, SubNode* aNext) : _next(aNext), _data(aData) {}
      inline SubNode(data_t* aData) : SubNode(aData, nullptr) {}

      inline const data_t*  data()    const { return _data; }

      inline       bool     hasNext() const { return _next != nullptr; }
      inline       SubNode* next()    const { return _next; }

             void init(data_t* aData, SubNode* aNext);  // required because of Reservoir
      inline void init(data_t* aData) { init(aData, nullptr); }

      inline sub_node_iterator begin()              { return sub_node_iterator(this); }
      inline sub_node_iterator end()                { return sub_node_iterator(nullptr); }
      inline const_sub_node_iterator cbegin() const { return const_sub_node_iterator(this); }
      inline const_sub_node_iterator cend()   const { return const_sub_node_iterator(nullptr); }
    };

  public:
    HtNested1(const size_t aNumBuckets, const uint32_t aMainRsvLog2ChunkSize, const uint32_t aSubRsvLog2ChunkSize);

  public:
    inline size_t      numBuckets()              const { return _hashdir.size(); }
    inline hashvalue_t hash(const data_t* aData) const { return hashfun_t::eval(aData); }
    inline size_t      size()                    const { return _size; }
    inline size_t      getRsvMainSize()          const { return _mainReservoir.cardinality(); }
    inline size_t      getRsvSubSize()           const { return _subReservoir.cardinality(); }

    size_t             memoryConsupmtion()           const;
    size_t             memoryConsupmtionDir()        const;
    size_t             memoryConsupmtionMainChains() const;
    size_t             memoryConsupmtionSubChains()  const;

  public:
    void insert(data_t* aData);
    void print(uint64_t aIndent, uint64_t aTab, const bool aSkipEmpty, std::ostream& os = std::cout) const;
    //inline bool     containsKey(const key_t aKey) { return findKey(aKey).first; }  // TODO
    
    template <typename Tprobedata, alg_hashfun_c Tprobehashfun, alg_binary_predicate_c Tjoinpred>
    const std::tuple<const MainNode*, const uint64_t> findMainNodeByOther(const Tprobedata* aProbeTuple) const;

    void clear();  // empty the hash table (clear reservoirs and reset directory nodes)

  public:
    // statistics
    stats_t makeStatistics() const;

    // get iterator of std::vector<MainNode>
    inline hashdir_iterator       begin()        { return _hashdir.begin(); }
    inline hashdir_iterator       end()          { return _hashdir.end(); }
    inline const_hashdir_iterator cbegin() const { return _hashdir.cbegin(); }
    inline const_hashdir_iterator cend()   const { return _hashdir.cend(); }

  private:
    inline       size_t     getDirIndex(const data_t* aData) const { return hash(aData) % numBuckets(); }
    inline       size_t     getDirIndex(const hashvalue_t h) const { return h % numBuckets(); }
    inline const MainNode*  getDirEntry(const data_t* aData) const { return &_hashdir.at(getDirIndex(aData)); }
    inline       MainNode*  getDirEntry(const data_t* aData)       { return &_hashdir.at(getDirIndex(aData)); }
    inline const MainNode*  getDirEntry(const hashvalue_t h) const { return &_hashdir.at(getDirIndex(h)); }

    // given data_t*, return its directory index, hash value and ptr to directory entry (as a tuple)
    inline std::tuple<size_t, hashvalue_t, MainNode*>
    getDirIdxHashEntry(const data_t* aData) {
      return { getDirIndex(aData), hash(aData), getDirEntry(aData) };
    }

  private:
    void                       insertAtMainNode(data_t* aData, const hashvalue_t aHashvalue, MainNode* aMainNode);
    SubNode*                   insertIntoSubchain(data_t* aData, MainNode* aMainNode);

    std::pair<bool, MainNode*> findMainNode(const data_t* aData);
    std::pair<bool, MainNode*> findMainNode(const data_t* aData, MainNode* aDirEntry);

    inline bool                isEqual(const data_t* x, const data_t* y) { return eqfun_t::eval(x, y); }
    inline bool                isMainNodeMatch(const data_t* aData, const MainNode* aNode) {
      return (aNode->hashvalue() == hash(aData)) && isEqual(aData, aNode->data());
    }

  private:
    hashdir_t           _hashdir;
    Reservoir<MainNode> _mainReservoir;  // reservoir stores objects (here: MainNodes/SubNodes) in chunks of
    Reservoir<SubNode>  _subReservoir;   // consecutive memory with constant memory addresses (unlike vectors)
    size_t              _size = 0;       // number of data items stored

};

// c'tor - parameters: number of hash table buckets, size of reservoirs in log2
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
HtNested1<Tdata, Thashfun, Tcontenteqfun> ::HtNested1(const uint64_t aNumBuckets,
                                                      const uint32_t aMainRsvLog2ChunkSize,
                                                      const uint32_t aSubRsvLog2ChunkSize)
  : _hashdir(aNumBuckets, MainNode()),
    _mainReservoir(aMainRsvLog2ChunkSize), _subReservoir(aSubRsvLog2ChunkSize) {}

// member functions returning the memory consumption in bytes
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtNested1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtion() const {
  return memoryConsupmtionDir() + memoryConsupmtionMainChains() + memoryConsupmtionSubChains();
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtNested1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtionDir() const {
  return numBuckets() * sizeof(MainNode);
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtNested1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtionMainChains() const {
  return getRsvMainSize() * sizeof(MainNode);
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
size_t
HtNested1<Tdata, Thashfun, Tcontenteqfun>::memoryConsupmtionSubChains() const {
  return getRsvSubSize() * sizeof(SubNode);
}

// insert a new data item
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtNested1<Tdata, Thashfun, Tcontenteqfun>::insert(data_t* aData) {
  const hashvalue_t lHashvalue = hash(aData);
  MainNode* lDirEntry = getDirEntry(aData);

  // check the directory entry of the bucket that aData hashes to
  if (lDirEntry->isEmpty() || isMainNodeMatch(aData, lDirEntry)) {
    insertAtMainNode(aData, lHashvalue, lDirEntry);

  } else {
    // must find the right MainNode to insert our data
    auto [lMatch, lCandidate] = findMainNode(aData);                    // TODO: findMainNode(aData, lDirEntry)
    if (lMatch) {
      // a matching MainNode for the data already exists
      insertAtMainNode(aData, lHashvalue, lCandidate);
    } else {
      // no matching MainNode exists, create one
      assert(!lCandidate->hasNext());
      lCandidate->_next = _mainReservoir.get_new_entry();
      insertAtMainNode(aData, lHashvalue, lCandidate->_next);
    }
  }
  ++_size;
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtNested1<Tdata, Thashfun, Tcontenteqfun>::print(uint64_t aIndent, uint64_t aTab, const bool aSkipEmpty, std::ostream& os) const {
  const std::string lIndentStr(aIndent, ' ');
  const std::string lTabStr(aTab, ' ');
  const size_t lDirDigits = std::log10(numBuckets());
  for (size_t i = 0; i < numBuckets(); ++i) {
    const MainNode& lDirEntry = _hashdir.at(i);
    if (lDirEntry.isEmpty() && aSkipEmpty) { continue; }
    os << lIndentStr;
    os << std::setw(lDirDigits) << i << ": ";
    if (lDirEntry.isEmpty()) {
      os << "(empty)" << std::endl;
      continue;
    }
    os << "\n";

    // go through dir entry's sub chain (if exists)
    os << lIndentStr << lTabStr << lDirEntry._hashvalue << " -> {" << lDirEntry._data;
    for (const SubNode* it = lDirEntry._subchain_head; it != nullptr; it = it->_next) {
      os << ", " << it->_data;
    }
    os << "}\n";

    // go through collision chain entries and their sub chains
    for (const MainNode* itMain = lDirEntry._next; itMain != nullptr; itMain = itMain->_next) {
      os << lIndentStr << lTabStr << itMain->_hashvalue << " -> {" << itMain->_data;
      for (const SubNode* itSub = itMain->_subchain_head; itSub != nullptr; itSub = itSub->_next) {
        os << ", " << itSub->_data;
      }
      os << "}\n";
    }

    //os << lIndentStr << "\n";
  }
}

// find a matching main node for a probe tuple
// requires a separate probe hash function that can compute the hash value of a data_t tuple and a Tprobedata tuple;
// also requires a probe equality function that can compare a data_t and a Tprobedata tuple for equality on some 
// join attribute (i.e., this function is the join predicate).
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
template <typename Tprobedata, alg_hashfun_c Tprobehashfun, alg_binary_predicate_c Tjoinpred>
auto
HtNested1<Tdata, Thashfun, Tcontenteqfun>::findMainNodeByOther(const Tprobedata* aProbeTuple) const
-> const std::tuple<const MainNode*, const uint64_t> {
  using other_hashfun_t = Tprobehashfun;
  using joinpred_t = Tjoinpred;
  static_assert(std::is_same_v<hashvalue_t, typename other_hashfun_t::output_t>);
  static_assert(std::is_same_v<typename joinpred_t::left_t, Tprobedata>);
  static_assert(std::is_same_v<typename joinpred_t::right_t, data_t>);

  hashvalue_t lProbeHashValue = other_hashfun_t::eval(aProbeTuple);
  const MainNode* lResultNode = getDirEntry(lProbeHashValue);

  uint64_t lNumCmps = 0;  // number of comparisons along the collision chain

  // walk the collision chain
  do {
    if (lResultNode->isEmpty()) { return {nullptr, lNumCmps}; }
    ++lNumCmps;
    if ((lResultNode->hashvalue() == lProbeHashValue) &&
        (joinpred_t::eval(aProbeTuple, lResultNode->data()))) {
      return {lResultNode, lNumCmps};
    }
    lResultNode = lResultNode->next();
  } while (nullptr != lResultNode);

  return {nullptr, lNumCmps};
}


// insert aData into a given MainNode (if its still empty) or into this MainNode's sub chain.
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtNested1<Tdata, Thashfun, Tcontenteqfun>::insertAtMainNode(data_t* aData, const hashvalue_t aHashvalue, MainNode* aMainNode) {
  if (aMainNode->isEmpty()) {  // MainNode doesn't contain any data yet
    aMainNode->init(aData, aHashvalue);
  } else {
    assert(isMainNodeMatch(aData, aMainNode));
    [[maybe_unused]] const SubNode* lNew = insertIntoSubchain(aData, aMainNode);
    assert(lNew->_data == aData);
  }
}

// insert aData into the sub chain of a given MainNode
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
typename HtNested1<Tdata, Thashfun, Tcontenteqfun>::SubNode*
HtNested1<Tdata, Thashfun, Tcontenteqfun>::insertIntoSubchain(data_t* aData, MainNode* aParent) {
  assert(isMainNodeMatch(aData, aParent));
  if (aParent->hasChild()) {  // at least one node in the subchain
    SubNode* lNewFirstChild = _subReservoir.get_new_entry();
    lNewFirstChild->init(aData, aParent->_subchain_head);
    aParent->_subchain_head = lNewFirstChild;
  } else {  // no node in the subchain -> create one
    aParent->_subchain_head = _subReservoir.get_new_entry();
    aParent->_subchain_head->init(aData);
  }
  return aParent->_subchain_head;
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
std::pair<bool, typename HtNested1<Tdata, Thashfun, Tcontenteqfun>::MainNode*>
HtNested1<Tdata, Thashfun, Tcontenteqfun>::findMainNode(const data_t* aData) {
  MainNode* lDirEntry = getDirEntry(aData);
  return findMainNode(aData, lDirEntry);
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
std::pair<bool, typename HtNested1<Tdata, Thashfun, Tcontenteqfun>::MainNode*>
HtNested1<Tdata, Thashfun, Tcontenteqfun>::findMainNode(const data_t* aData, MainNode* aDirEntry) {
  if (isMainNodeMatch(aData, aDirEntry)) {
    return {true, aDirEntry};
  } else {
    // iterator provides hasNext(), MainNode* provides next()
    main_node_iterator it = aDirEntry;
    for (; it.hasNext(); ++it) {
      if (isMainNodeMatch(aData, it->next())) {
        return {true, it->next()};
      }
    }
    return { false, &*it };
  }
}

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtNested1<Tdata, Thashfun, Tcontenteqfun>::clear() {
  _subReservoir.erase();  // empty reservoir and release memory
  _mainReservoir.erase();  // empty reservoir and release memory
  // set hash dir buckets to empty
  for (size_t i = 0; i < numBuckets(); ++i) {
    _hashdir.at(i).clear();
  }
}

// create hash table statistics
template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
auto
HtNested1<Tdata, Thashfun, Tcontenteqfun>::makeStatistics() const
-> stats_t {
  stats_t s;
  s._numBuckets = numBuckets();
  s._numEntries = size();

  // iterate buckets
  for (const_hashdir_iterator dir_it = cbegin(); dir_it != cend(); ++dir_it) {
    if (dir_it->isEmpty()) {
      ++s._numEmptyBuckets;
      s._collisionChainLen.step(0);
      continue;
    }
    // iterate main collision chain
    uint32_t lCurrChainLen = 0;  // len of current bucket's collision chain (incl dir entry)
    // collision chain length also counts the dir entry, i.e. non-empty buckets have len >= 1
    for (const_main_node_iterator main_it = dir_it->cbegin(); main_it.valid(); ++main_it) {
      assert(!main_it->isEmpty());
      ++s._numDistinctKeys;
      ++lCurrChainLen;
      // iterate sub chain
      for (const_sub_node_iterator sub_it{main_it->child()}; sub_it.valid(); ++sub_it) {
        // nop
      }
    }
    s._collisionChainLen.step(lCurrChainLen);
    s._collisionChainLenNonempty.step(lCurrChainLen);
  }

  return s;
}

// nested classes

template <typename Tdata, alg_hashfun_c Thashfun, alg_binary_predicate_c Tcontenteqfun>
void
HtNested1<Tdata, Thashfun, Tcontenteqfun>::SubNode::init(data_t* aData, SubNode* aNext) {
  _data = aData;
  _next = aNext;
}


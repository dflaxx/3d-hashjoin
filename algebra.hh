#pragma once

#include "util/standard_includes.hh"

#include <chrono>
#include <functional>
#include <thread>

#include "util/debugging_helpers.hh"
#include "util/chrono_helpers.hh"
#include "ht_nested.hh"
#include "ht_chaining.hh"

/*
 * A simple algebra to demonstrate hash joins with a nested hashtable.
 * 
 * Properties: templated, push-based (producers and consumers), tuple-wise processing
 * Underlying storage: row store (see class RelationRS) 
 *
 * Operators:
 * - Table Scan
 * - Selection
 * - Hashjoin Build for regular chaining hashtable
 * - Hashjoin Probe with regular chaining hashtable
 * - Hashjoin Build for nested hashtable
 * - Hashjoin Probe with nested hashtable
 * - Unnest
 * - Top
 *
 * All operators inherit from a common base class AlgBase which provides some
 * debugging functionalities. Otherwise, there is no inheritance and also no
 * overriding of virtual functions.
 * Instead, the algebra is templated, i.e., all operators are class templates.
 * To ensure that only compatible types are combined, C++20 concepts are used.
 *
 * TODO
 * - Rethink const-ness of members and member function parameters.
 *   is overload a solution?
 *     void step(const input_t*, globstat_t*)
 *     void step(      input_t*, globstat_t*)
 *
 */

class AlgBase;  // fwd decl for definition of concepts

// Concepts (C++20)
#if __cplusplus >= 202002L
#include <concepts>
#include "concepts.hh"  // own concepts from separate file

// concept for every algebra operator
template <typename T>
concept alg_operator_c =
  std::derived_from<T, AlgBase> &&
  requires(T c) {
    typename T::globstat_t;
    typename T::input_t;
    typename T::output_t;
  };

// concept for consumer operators like HashJoinBuild
template <typename T>
concept alg_consumer_c =
  alg_operator_c<T> &&
  requires(T c) {
    { c.init(static_cast<typename T::globstat_t*>(nullptr)) }
      -> std::same_as<void>;
    { c.step(static_cast<typename T::input_t*>(nullptr),
             static_cast<typename T::globstat_t*>(nullptr)) }
      -> std::same_as<void>;
    { c.fin(static_cast<typename T::globstat_t*>(nullptr)) }
      -> std::same_as<void>;
  };

// concept for producers operators like TableScan
template <typename T>
concept alg_producer_c =
  alg_operator_c<T> &&
  requires(T c) {
    { c.run(static_cast<typename T::globstat*>(nullptr)) } -> std::same_as<void>;
  };

// hash table build operators must have a hashtable() getter member function
// to return a reference to the hashtable to the respective probe operator.
template <typename T>
concept alg_buildop_c =
  alg_consumer_c<T> &&
  requires(T t) {
    typename T::hashtable_t;
    { t.hashtable() } -> std::same_as<const typename T::hashtable_t&>;
  };
#else
#error Compiling with C++ < 20, concepts not supported.
#endif  // __cplusplus >= 202002L


// a row store relation
template <typename Ttuple>
struct RelationRS {
  using tuple_t  = Ttuple;
  using tuple_vt = std::vector<tuple_t>;

  tuple_vt _tuples;

  inline size_t card() const { return _tuples.size(); }
};

template <typename Ttuple>
std::ostream& operator<<(std::ostream& os, const RelationRS<Ttuple>& aRel) {
  for (const auto& t : aRel._tuples) {
    os << t << "\n";
  }
  return os;
}


// global state struct -- currently not used
struct GlobStat0 {
  size_t _ht_num_buckets;
  size_t _ht_rsv_log2_chunksize_main;
  size_t _ht_rsv_log2_chunksize_sub;
  size_t _ht_rsv_log2_chunksize;
};


// get runtime of operator without parent (consumer) operator runtimes
// XXX Does NOT return the exclusive runtime as one might expect due to pipelining of tuples.
// XXX Don't use this.
template <alg_operator_c Toperator>
auto
get_runtime_excl(const Toperator* aOp) -> typename Toperator::duration_t {
  assert (aOp != nullptr);
  if constexpr (requires (Toperator t) { t.consumer(); }) {
    return (aOp->getRuntime() - aOp->consumer()->getRuntime());
  } else {
    return aOp->getRuntime();
  }
}

/*
 * Print a plan strand from the top-most the bottom-most operator:
 * operator name|count|exclusive runtime|runs.
 * e.g. given strand Scan -> Selection -> HashJoinProbe -> Top,
 * print_strand(&Scan) prints: Top \n HashJoinProbe \n Selection \n Scan.
 *
 * XXX See warning for function get_runtime_excl.
 */
template <alg_operator_c Toperator>
void print_strand(const Toperator* aOp, const size_t aIndentLvl = 0, std::ostream& os = std::cout)
{
  assert (aOp != nullptr);
  if constexpr (requires (Toperator t) { t.consumer(); }) {
    print_strand(aOp->consumer(), aIndentLvl, os);
  }
  os
    << df::infra::indent(aIndentLvl)
    << aOp->name() << "|"
    << aOp->count() << "|"
    << df::infra::to_string(get_runtime_excl(aOp), true) << "|"
    << aOp->runs()
    << "\n";
}


// common base class for all operators, provides some statistics functions
class AlgBase {
  public:
    using clock_t = std::chrono::steady_clock;
    using time_point_t = std::chrono::time_point<clock_t>;
    using duration_t = std::chrono::nanoseconds;
  public:
    explicit AlgBase(const std::string& aName)
      : _count(0), _ok(true), _startTime(), _stopTime(), _name(aName), _runs(0) {}
    AlgBase() : AlgBase("") {}
  public:
    inline       void         reset() { _count = 0; _ok = true; startTimer(); ++_runs; }
    inline       void         inc() { ++_count; }
    inline       uint64_t     count() const { return _count; }
    inline       bool         ok() const { return _ok; }
    inline       bool         ok(const bool b) { return (_ok = b); }
    inline       void         startTimer() { _startTime = clock_t::now(); }
    inline       void         stopTimer()  { _stopTime = clock_t::now(); }
    inline const std::string& name() const { return _name; }
    inline       uint64_t     runs() const { return _runs; }

    // note: runtime is inclusive, i.e. 'earlier' operator includes parent op runtimes
    inline duration_t getRuntime() const {
      return std::chrono::duration_cast<duration_t>(_stopTime - _startTime);
    }
    inline std::string getRuntimeStr() const {
      return df::infra::to_string(getRuntime(), true);
    }
  protected:
    uint64_t     _count;     // e.g. number of step() calls for consumers
    bool         _ok;        // e.g. to indicate error states
    time_point_t _startTime; // beginning and end of processing
    time_point_t _stopTime;
    std::string  _name;      // operator name like "AlgTop"
    uint64_t     _runs;      // number of times the operator was re-run (e.g. for repeated experiments)
};


// Top operator
template <typename Tinput, typename Tglobstat>
class AlgTop : public AlgBase {
  public:
    using globstat_t  = Tglobstat;
    using input_t     = Tinput;
    using output_t    = void;  // top operator is the root and has no consumer
    using print_fun_t = std::function<void (const input_t*, std::ostream& os)>;
  public:
    inline AlgTop()
      : AlgTop(std::cout, true) {}
    inline AlgTop(std::ostream& aOs, const bool aPrintResult)
      : AlgBase("AlgTop"), _os(aOs), _print_result(aPrintResult) {}
    inline AlgTop(std::ostream& aOs, const bool aPrintResult, print_fun_t aPrintFunction)
      : AlgBase("AlgTop"), _os(aOs), _print_result(aPrintResult), _print_fun(aPrintFunction) {}
  public:
    inline void init([[maybe_unused]] globstat_t* aGlobstat) {
      //std::cout << "AlgTop::init()\n";
      reset();
    }
    inline void step(input_t* aInput, [[maybe_unused]] globstat_t* aGlobstat) {
      inc();
      if (_print_result && runs() == 1) {  // only print result the first time the operator is run
        _print_fun(aInput, _os);
        _os << "\n";
      }
    }
    inline void fin([[maybe_unused]] globstat_t* aGlobstat) {
      //std::cout << "AlgTop::fin()\n";
      stopTimer();
    }
  public:
    bool printResult() const { return _print_result; }
    void printResult(const bool aPrint) { _print_result = aPrint; }
  private:
    std::ostream& _os;
    bool          _print_result;
    print_fun_t   _print_fun = [] (const input_t* aInput, std::ostream& aOs) {
      aOs << aInput;
    };
};


// Table scan operator
template <alg_consumer_c Tconsumer> 
class AlgScan : public AlgBase {
  public:
    using consumer_t = Tconsumer;
    using globstat_t = typename consumer_t::globstat_t;
    using input_t    = typename consumer_t::input_t;  // push input as is
    using output_t   = typename consumer_t::input_t;
    using input_rel_t = RelationRS<input_t>;
  public:
    inline AlgScan(consumer_t* aConsumer, input_rel_t* aRelation)
      : AlgBase("AlgScan"), _consumer(aConsumer), _relation(aRelation) {}
  public:
    inline void run(globstat_t *aGlobstat) {
      //std::cout << "AlgScan::run()\n";
      reset();
      _consumer->init(aGlobstat);
      for (auto &t : _relation->_tuples) {
        inc();
        _consumer->step(&t, aGlobstat);
      }
      _consumer->fin(aGlobstat);
      stopTimer();
    }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
  private:
    consumer_t*        _consumer;
    input_rel_t* _relation;
};


// Selection operator
template <alg_consumer_c Tconsumer, alg_predicate_c Tpredicate>
class AlgSelection : public AlgBase {
  public:
    using consumer_t  = Tconsumer;
    using globstat_t  = typename consumer_t::globstat_t;
    using input_t     = typename consumer_t::input_t;
    using output_t    = typename consumer_t::input_t;
    using predicate_t = Tpredicate;
  public:
    inline AlgSelection(Tconsumer* aConsumer)
      : AlgBase("AlgSelection"), _consumer(aConsumer) {}
  public:
    inline void init(globstat_t* aGlobstat) {
      reset();
      _consumer->init(aGlobstat);
    }
    inline void step(input_t* aInput, globstat_t* aGlobstat) {
      if (predicate_t::eval(aInput)) {
        inc();
        _consumer->step(aInput, aGlobstat);
      }
    }
    inline void step(const input_t* aInput, globstat_t* aGlobstat) {
      if (predicate_t::eval(aInput)) {
        inc();
        _consumer->step(aInput, aGlobstat);
      }
    }
    inline void fin(globstat_t* aGlobstat) {
      _consumer->fin(aGlobstat);
      stopTimer();
    }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
  private:
    consumer_t* _consumer;
};


// Dynamic selection operator (runtime predicate instead of compile-time)
template <alg_consumer_c Tconsumer, alg_dyn_predicate_c Tpredicate>
class AlgDynSelection : public AlgBase {
  public:
    using consumer_t  = Tconsumer;
    using globstat_t  = typename consumer_t::globstat_t;
    using input_t     = typename consumer_t::input_t;
    using output_t    = typename consumer_t::input_t;
    using predicate_t = Tpredicate;
  public:
    inline AlgDynSelection(consumer_t* aConsumer, predicate_t aPredicate)
      : AlgBase("AlgDynSelection"), _consumer(aConsumer), _pred(aPredicate) {}
    inline AlgDynSelection(consumer_t* aConsumer)
      : AlgDynSelection(aConsumer, predicate_t()) {}
  public:
    inline void init(globstat_t* aGlobstat) {
      reset();
      _consumer->init(aGlobstat);
    }
    inline void step(input_t* aInput, globstat_t* aGlobstat) {
      if (_pred(aInput)) {
        inc();
        _consumer->step(aInput, aGlobstat);
      }
    }
    inline void step(const input_t* aInput, globstat_t* aGlobstat) {
      if (_pred(aInput)) {
        inc();
        _consumer->step(aInput, aGlobstat);
      }
    }
    inline void fin(globstat_t* aGlobstat) {
      _consumer->fin(aGlobstat);
      stopTimer();
    }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
  private:
    consumer_t* _consumer;
    predicate_t _pred;
};


// 3D hash join build operator (nested hash table)
template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, typename Tglobstat>
class AlgNestJoinBuild : public AlgBase {
  public:
    using globstat_t  = Tglobstat;
    using hashfun_t   = Thashfun;
    using input_t     = typename hashfun_t::input_t;
    using output_t    = void;
    using eqfun_t     = Tequalfun;
    using hashtable_t = HtNested1<input_t, hashfun_t, eqfun_t>;
  public:
    AlgNestJoinBuild(const size_t aHashDirSize,
                     const uint32_t aHtLog2ChunkSizeMain,
                     const uint32_t aHtLog2ChunkSizeSub)
      : AlgBase("AlgNestJoinBuild"),
        _hashtable(aHashDirSize, aHtLog2ChunkSizeMain, aHtLog2ChunkSizeSub) {}
    AlgNestJoinBuild(const globstat_t* aGlobstat)
      : AlgNestJoinBuild(aGlobstat->_ht_num_buckets,
                                    aGlobstat->_ht_rsv_log2_chunksize_main,
                                    aGlobstat->_ht_rsv_log2_chunksize_sub) {}
  public:
    inline void init([[maybe_unused]] globstat_t* aGlobstat) {
      //std::cout << "AlgNestJoinBuild::init()\n";
      reset();
    }
    inline void step(input_t* aInput, [[maybe_unused]] globstat_t* aGlobstat) {
      inc();
      _hashtable.insert(aInput);
    }
    inline void fin([[maybe_unused]] globstat_t* aGlobstat) {
      //std::cout << "AlgNestJoinBuild::fin()\n";
      stopTimer();
      //std::cout << "  Build Hash Table:\n";
      //_hashtable.print(4, 2, true);
    }
  public:
    inline const hashtable_t& hashtable() const { return _hashtable; }
    inline void               clear_ht()        { _hashtable.clear(); }
  private:
    hashtable_t _hashtable;
};


/*
 * 3D hash join probe operator
 *
 * Note:
 * Probe into the nested hashtable provided by the respective AlgNestJoinBuild operator
 * AlgBase::_count counts the number of matches when probing
 */
template <alg_consumer_c Tconsumer, alg_buildop_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg_concatfun_c Tconcatfun>
class AlgNestJoinProbe : public AlgBase {
  public:
    using consumer_t     = Tconsumer;
    using build_t        = Tbuild;
    using hashfun_t      = Thashfun;
    using globstat_t     = typename consumer_t::globstat_t;
    using input_t        = typename hashfun_t::input_t;
    using output_t       = typename consumer_t::input_t;
    using joinpred_t     = Tjoinpred;
    using concatfun_t    = Tconcatfun;
  public:
    AlgNestJoinProbe(consumer_t* aConsumer, build_t* aBuildOperator)
      : AlgBase("AlgNestJoinProbe"),
        _consumer(aConsumer), _buildOperator(aBuildOperator), _outputTuple(), _numCmps() {}
  public:
    inline void init(globstat_t* aGlobstat) {
      //std::cout << "AlgNestJoinProbe::init()\n";
      reset();
      _numCmps = 0;
      _consumer->init(aGlobstat);
    }
    inline void step(input_t* aProbeTuple, [[maybe_unused]] globstat_t* aGlobstat) {
      // probe into _buildOperator->hashtable() with the probetuple.
      constexpr bool lTrace = false;
      if (lTrace) { std::cout << "  AlgNestJoinProbe::step() -> probe tuple: " << *aProbeTuple << "\n"; }

      //using MainNode = typename build_t::hashtable_t::MainNode;

      // Signature of HtNested1::findMainNodeByOther:
      // template <typename Tprobedata, typename Tprobehashfun, typename Tjoinpred>
      // const std::tuple<const MainNode*, const uint64_t>
      // findMainNodeByOther(const Tprobedata* aProbeTuple) const;
      const auto [lJoinMatch, lNumComparisons] = _buildOperator->hashtable()
        .template findMainNodeByOther<input_t, hashfun_t, joinpred_t>(aProbeTuple);

      _numCmps += lNumComparisons;

      if (nullptr == lJoinMatch) {
        if (lTrace) { std::cout << "    NO join partner found.\n"; }
      } else {
        if (lTrace) { std::cout << "    join partner(s) found: " << lJoinMatch->data() << "\n"; }
        _outputTuple = concatfun_t::eval(aProbeTuple, lJoinMatch);
        inc();
        _consumer->step(&_outputTuple, aGlobstat);
      }
    }
    inline void fin(globstat_t* aGlobstat) {
      //std::cout << "AlgNestJoinProbe::fin()\n";
      _consumer->fin(aGlobstat);
      stopTimer();
    }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline       uint64_t    numCmps()  const { return _numCmps; }
  private:
    consumer_t* _consumer;
    build_t*    _buildOperator;  // AlgNestJoinBuild, has getter method hashtable()
    output_t    _outputTuple;    // nested output tuple: <tuple_left_t*, HtNested::MainNode*>
    uint64_t    _numCmps;        // total number of collision chain comparisons
};


/*
 * 3D Hash Join Unnest Operator
 *
 * Expands/unnests nested tuples to concatenated tuples from left and right relation
 * input:  1x     (probe_side_tuple*, hashtable_main_node*)
 * output: [0,n]x (probe_side_tuple*, build_side_tuple*)
 *
 * Note: Tunnestfun is responsible for knowing the layout of the nested tuple.
 *       AlgUnnestHt does not have to know it.
 *
 * AlgBase::_count counts the output size,
 * i.e., the number of *unnested* tuples leaving the operator
 */
template <alg_consumer_c Tconsumer, alg_unnestfun_c Tunnestfun, typename Thtnested>
class AlgUnnestHt : public AlgBase {
  public:
    using consumer_t  = Tconsumer;
    using globstat_t  = typename consumer_t::globstat_t;
    using unnestfun_t = Tunnestfun;
    using output_t    = typename consumer_t::input_t;
    using input_t     = typename unnestfun_t::input_t;  // a nested tuple
    using ht_nested_t = Thtnested;
    static_assert(std::is_same_v<output_t, typename unnestfun_t::output_t>,
        "AlgUnnestHt::output_t (aka consumer_t::input_t) does not match unnestfun_t::output_t");
  public:
    inline AlgUnnestHt(consumer_t* aConsumer)
      : AlgBase("AlgUnnest"),
        _consumer(aConsumer), _outputTuple() {}
  public:
    inline void init([[maybe_unused]] globstat_t* aGlobstat) { 
      //std::cout << "AlgUnnestHt::init()\n";
      reset();
      _consumer->init(aGlobstat);
    }
    inline void step(input_t* aNestedTuple, [[maybe_unused]] globstat_t* aGlobstat) {
      const bool lTrace = false;
      if (lTrace) { std::cout << "AlgUnnestHt::step() -> nested tuple: " << aNestedTuple << "\n"; }
      // aNestedTuple is a nested tuple with _left (probe side tuple ptr)
      // and _right (ptr to hashtable main node)
      // --> walk the subchain of the main node and produce one unnested output tuple
      //     per sub node
      using MainNode = typename ht_nested_t::MainNode;
      using SubNode  = typename ht_nested_t::SubNode;

      const MainNode* lNestedRight = unnestfun_t::getMainNode(aNestedTuple);
      assert(!lNestedRight->isEmpty()); // a valid nested input tuple is present,
                                        // so the nested part must at least contain
                                        // one join match

      // first data from _right is in the MainNode
      unnestfun_t::eval_left(&_outputTuple, aNestedTuple);
      unnestfun_t::eval_right(&_outputTuple, aNestedTuple, lNestedRight->data());
      if (lTrace) { std::cout << "  unnested tuple: " << _outputTuple << "\n"; }
      _consumer->step(&_outputTuple, aGlobstat);
      inc();

      SubNode* lSubNode = lNestedRight->child();
      while (nullptr != lSubNode) {
        unnestfun_t::eval_right(&_outputTuple, aNestedTuple, lSubNode->data());
        if (lTrace) { std::cout << "  unnested tuple: " << _outputTuple << "\n"; }
        _consumer->step(&_outputTuple, aGlobstat);
        inc();
        lSubNode = lSubNode->next();
      }

    }
    inline void fin([[maybe_unused]] globstat_t* aGlobstat) {
      //std::cout << "AlgUnnestHt::fin()\n";
      _consumer->fin(aGlobstat);
      stopTimer();
    }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
  private:
    consumer_t* _consumer;
    output_t    _outputTuple;   // unnested output tuple for pushing to the consumer
};


// Regular hash join build operator (chaining hash table)
template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, typename Tglobstat>
class AlgHashJoinBuild : public AlgBase {
  public:
    using globstat_t  = Tglobstat;
    using hashfun_t   = Thashfun;
    using input_t     = typename hashfun_t::input_t;
    using output_t    = void;
    using eqfun_t     = Tequalfun;
    using hashtable_t = HtChaining1<input_t, hashfun_t, eqfun_t>;
  public:
    AlgHashJoinBuild(const size_t aHashDirSize, const uint32_t aHtLog2ChunkSize)
      : AlgBase("AlgHashJoinBuild"), _hashtable(aHashDirSize, aHtLog2ChunkSize) {}
    AlgHashJoinBuild(const globstat_t* aGlobstat)
      : AlgHashJoinBuild(aGlobstat->_ht_num_buckets, aGlobstat->_ht_rsv_log2_chunksize) {}
  public:
    inline void init([[maybe_unused]] globstat_t* aGlobstat) {
      reset();
    }
    inline void step(input_t* aTuple, [[maybe_unused]] globstat_t* aGlobstat) {
      inc();
      _hashtable.insert(aTuple);
    }
    inline void fin([[maybe_unused]] globstat_t* aGlobstat) {
      stopTimer();
    }
  public:
    inline const hashtable_t& hashtable() const { return _hashtable; }
    inline void               clear_ht()        { _hashtable.clear(); }
  private:
    hashtable_t _hashtable;
};


/*
 * Regular hash join probe operator
 *
 * Note
 * Probe into the chaining hashtable provided by the respective AlgHashJoinBuild operator
 * AlgBase::_count counts the number of matches when probing
 *
 * Template parameters:
 * - bool IsBuildKeyUnique: If the hash table was built on a unique (key) attribute,
 *                          probe can terminate after the first match.
 */
template <alg_consumer_c Tconsumer, alg_buildop_c Tbuild,
          alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
          alg_concatfun_c Tconcatfun,
          bool IsBuildKeyUnique = false>
class AlgHashJoinProbe : public AlgBase {
  public:
    using consumer_t  = Tconsumer;
    using build_t     = Tbuild;
    using hashfun_t   = Thashfun;
    using globstat_t  = typename consumer_t::globstat_t;
    using input_t     = typename hashfun_t::input_t;
    using output_t    = typename consumer_t::input_t;
    using hashvalue_t = typename hashfun_t::output_t;
    using joinpred_t  = Tjoinpred;
    using concatfun_t = Tconcatfun;
  public:
      inline AlgHashJoinProbe(consumer_t* aConsumer, build_t* aBuildOperator)
        : AlgBase("AlgHashJoinProbe"),
          _consumer(aConsumer), _buildOperator(aBuildOperator), _outputTuple(), _numCmps() {}
  public:
    inline void init([[maybe_unused]] globstat_t* aGlobstat) {
      reset();
      _numCmps = 0;
      _consumer->init(aGlobstat);
    }
    inline void step(input_t* aTuple, [[maybe_unused]] globstat_t* aGlobstat) {
      constexpr bool lTrace = false;
      if (lTrace) { std::cout << "  AlgHashJoinProbe::step() -> probe tuple: " << *aTuple << "\n"; }

      using const_node_iterator = typename build_t::hashtable_t::const_node_iterator;

      hashvalue_t lProbeHashValue = hashfun_t::eval(aTuple);
      uint64_t lNumComparisons = 0;

      // Signature of HtChaining1::findDirEntryByOther:
      // template <typename Tprobedata, alg_hashfun_c Tprobehashfun>
      // const_node_iterator
      // findDirEntryByOther(const Tprobedata* aProbeTuple) const;
      const_node_iterator it = _buildOperator->hashtable()
        .template findDirEntryByOther<input_t, hashfun_t>(aTuple);
      if (it->isEmpty()) {
        if (lTrace) { std::cout << "    NO join partner found.\n"; }
        return;
      }
      for (; it != nullptr; ++it) {
        assert(!it->isEmpty());  // none of the collision chain nodes should be empty
        ++lNumComparisons;
        if ((it->hashvalue() == lProbeHashValue) &&
            (joinpred_t::eval(aTuple, it->data()))) {
          if (lTrace) { std::cout << "    join partner found: " << it->data() << "\n"; }
          _outputTuple = concatfun_t::eval(aTuple, it->data());
          inc();
          _consumer->step(&_outputTuple, aGlobstat);
          if constexpr (IsBuildKeyUnique) {    // this is the one and only match -> can terminate
            break;
          }
        }
      }
      _numCmps += lNumComparisons;
    }
    inline void fin([[maybe_unused]] globstat_t* aGlobstat) {
      _consumer->fin(aGlobstat);
      stopTimer();
    }
  public:
    inline const consumer_t* consumer() const { return _consumer; }
    inline       uint64_t    numCmps()  const { return _numCmps; }
  private:
    consumer_t* _consumer;
    build_t*    _buildOperator; // AlgHashJoinBuild, has getter method hashtable()
    output_t    _outputTuple;   // concat'd output tuple: <tuple_left_t*, tuple_right_t*>
    uint64_t    _numCmps;       // total number of collision chain comparisons
};


#include "util/standard_includes.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <limits>
#include <ostream>
#include <random>
#include <unordered_set>

#include "lib/CLI11.hpp"

#include "util/csv_writer.hh"
#include "util/chrono_helpers.hh"
#include "util/measure_helpers.hh"
#include "util/output_helpers.hh"
#include "util/string_helpers.hh"
#include "util/GenRandIntVec.hh"
#include "util/hasht.hh"
#include "util/zipf_distribution.hh"

#include "algebra.hh"
#include "ht_chaining.hh"
#include "ht_nested.hh"

/*
 * Measure runtime of different query execution plans involving a single 
 * key/foreign key hash join between two base relations R and S on integer attributes.
 * There are two implementations for the join:
 * one using a chaining hash table and one using a nested hash table.
 *
 * Relations:
 * - R := {[k int, a int, b int]}   -> key relation
 *   R.k \in [0, |R| - 1]
 *   R.a and b unused
 * - S := {[k int, a int, b int]}   -> foreign key relation, S.a references R.k
 *   S.a \in A \subseteq [0, |R| - 1]
 *   S.k \in [0, |S| - 1] used as unique tuple id; b unused
 *
 * Hash function: murmur32 (from hasht.hh)
 *
 * Dimensions:
 * - cardinality of R and S \in {2^10, 2^11, ..., 2^25}
 *   => 16 possibilities per relation => 16^2 = 256 combinations
 *
 * - skew: no = uniform distribution, yes = zipf distribution with param = 1.0
 *   => 2 combinations
 *
 * - values for S.a
 *   S.a \in [0, a_max(t)]
 *   with a_max(t) = (|R| / 2^t) - 1, for t \in {0, ..., 9}
 *   => 10 combinations
 *
 * - number of hash table buckets
 *   In the default case, choose number of HT buckets equal to
 *   the number of distinct values (#dv) to store.
 *   Parameter b scales down this value linearly:
 *   #bucket := #dv / b, but at least 1.
 *   b \in {1, 2, 3, 4}
 * 
 * Total: 256 * 2 * 10 * 4 = 20480 parameter combinations
 *
 * Want to evaluate four plans:
 * - Csr: chaining hash table with build on R and probe on S
 * - Crs: chaining hash table with build on S and probe on R
 * - Nrs: nested hash table with build on S and probe on R
 * - Nsr: nested hash table with build on R and probe on S
 *
 * We're also interested in the runtime of the build, probe and unnest operators
 * alone. Therefore, the following plans are also evaluated to be able to 
 * compute the relevant runtimes by subtraction:
 * - scr and scs: scan R and scan S (scan op -> top op).
 *   Note: independent from the hash table implementation
 * - NrsNU: like Nrs, but without the unnest operator,
 *   i.e., the probe result (nested tuples) goes straight to the top operator.
 *
 * TODO
 * - output hash table size and memory consumption
 */


// a generic tuple of 3 unsigned 32-bit integers
struct tuple_uint32_3_t { uint32_t k, a, b; };

std::ostream& operator<<(std::ostream& os, const tuple_uint32_3_t& t) {
  os << "[" << t.k << "," << t.a << "]";
  return os;
}

// one instance = run of one parameter combination
class Experiment1 {
  public:
    enum class plans_e {
      NONE     = 0,
      scr      = 1,
      scs      = 2,
      Csr      = 4,
      Crs      = 8,
      Nrs      = 16,
      Nsr      = 32,
      NrsNU    = 64,
      //CsrSel   = 128,
      //CrsSel   = 256,
      //NrsSel   = 512,
      //CsrSelLf = 1024,  // unused/meaning unclear
      //CrsSelLf = 2048,
      //NrsSelLf = 4096,
      CsrUU    = 8192,
      ALL    = (CsrUU << 1) - 1
    };
  public:
    using base_tuple_t = tuple_uint32_3_t;
    using value_t      = uint32_t;
    using value_vt     = std::vector<value_t>;
    using clock_t      = std::chrono::steady_clock;
    using time_point_t = std::chrono::time_point<clock_t>;
    using plans_map_t  = std::unordered_map<std::string, plans_e>;
  public:
    struct GlobStat;
    struct nested_tuple_RS_t;  // [r, {s | joinpred(r,s)}]
    struct nested_tuple_SR_t;  // [s, {r | joinpred(r,s)}]
    struct result_tuple_t;
    struct HashfunOnR;
    struct HashfunOnS;
    struct EqfunBuildR;
    struct EqfunBuildS;
    struct EqfunJoinpredRS;
    struct EqfunJoinpredSR;
    struct ConcatfunChaining;
    struct ConcatfunNestedRS;
    struct ConcatfunNestedSR;
    struct UnnestFunRS;
    struct UnnestFunSR;
    struct SelectionRkMostFreq;
    struct SelectionRkLeastFreq;
  public:
    inline Experiment1(const uint32_t aLog2CardR, const uint32_t aLog2CardS,
                       const bool aSkew, const uint32_t aParakT, const uint32_t aParamB,
                       const std::filesystem::path aMeasureCsvFile, const std::vector<std::string>& aPlans)
      : _log2CardR(aLog2CardR), _log2CardS(aLog2CardS), _skew(aSkew), _t(aParakT), _b(aParamB),
        _plans(plans_e::NONE), _measureCsv(aMeasureCsvFile) {
      assert(64 > aLog2CardR);
      assert(64 > aLog2CardS);
      plansFromVec(aPlans);
    }
    inline Experiment1(const uint32_t aLog2CardR, const uint32_t aLog2CardS,
                       const bool aSkew, const uint32_t aParakT, const uint32_t aParamB,
                       const std::filesystem::path aMeasureCsvFile)
      : Experiment1(aLog2CardR, aLog2CardS, aSkew, aParakT, aParamB,
                    aMeasureCsvFile, std::vector<std::string>({"ALL"})) {
      assert(64 > aLog2CardR);
      assert(64 > aLog2CardS);
    }
  public:
    void init();
    void run();
  public:
    inline size_t   cardR() const { return (1U << _log2CardR); }
    inline size_t   cardS() const { return (1U << _log2CardS); }
    inline uint32_t t()     const { return _t; }
    inline uint32_t b()     const { return _b; }
    
    inline std::chrono::milliseconds minRuntime() const { return _minRuntime; }
    inline size_t                    minRepeat()  const { return _minRepeat; }

    inline static const plans_map_t& getSupportedPlans() { return _supportedPlansMap; }

  public:
    void printRelations() const;
    void printTimers() const;
    void printTypes() const;


  private:
    void generateR(const value_vt& aKeys);
    void generateS(const value_vt& aKeys, const value_vt& aForeignKeys);

    void runCsr();
    void runCrs();
    void runNrs();
    void runNsr();
    void runScanR();  // scan -> top
    void runScanS();
    void runNrsNoUnnest();  // scan(S) -> build && scan(R) -> probe -> top
    void runCsrUnknownUnique();

    inline uint32_t getFkMax() const { return (1 << (_log2CardR - _t)); }  // |R| / 2^t

    //inline uint32_t getSelectionLimitLeastFreq() const {
    //  return (getFkMax() > _withSelectBottomM) ? (getFkMax() - _withSelectBottomM) : 0;
    //}

    void writeCsvHeader();
    void writeExpParamsToCsv();  // common experiment parameters to CSV

    inline void timerStart(const std::string& aDesc) { timerCtl(aDesc); }
    inline void timerStop(const std::string& aDesc)  { timerCtl(aDesc, true); }
    void timerCtl(const std::string& aDesc, const bool aIsStop = false);

    void plansFromVec(const std::vector<std::string>&);

  private:
    // experiment parameters
    uint32_t _log2CardR;
    uint32_t _log2CardS;
    bool     _skew;
    uint32_t _t;  // scale parameter for S.a (foreign key)
    uint32_t _b;  // scale parameter for #buckets in hash table

    std::chrono::milliseconds _minRuntime{300};  // repeat running plan until this bound is reached
    size_t                    _minRepeat{8};     // min number of repeated measurements per plan

    plans_e  _plans = plans_e::ALL;  // plans to execute in run()

    const uint32_t _log2RsvChunkSize = 10;  // reservoir chunk size for all hashtables

    RelationRS<base_tuple_t> _R;  // key relation
    RelationRS<base_tuple_t> _S;  // foreign key relation

    size_t _numDvSa;  // number of distinct values in S.a

    df::infra::CSVWriter _measureCsv;

    std::map<std::string, std::pair<time_point_t, time_point_t>> _timePoints{};  // experiment timers
    bool _trace = true;

    // the base hash function used by all Hashfun classes
    inline static std::function<uint32_t(value_t)> _hashfunc = ht::murmur_hash<value_t>;

    static std::function<void(const result_tuple_t*, std::ostream&)> _result_tuple_printer;

    static const plans_map_t _supportedPlansMap;

  public:
    // return type deduction: http://mochan.info/c++/2019/06/12/returntype-deduction.html
    using hashvalue_t = decltype(_hashfunc(std::declval<value_t>()));
};

// enum operators
inline Experiment1::plans_e
operator|(const Experiment1::plans_e lhs, const Experiment1::plans_e rhs) {
  using type_t = std::underlying_type_t<Experiment1::plans_e>;
  return static_cast<Experiment1::plans_e>(static_cast<type_t>(lhs) | static_cast<type_t>(rhs));
}
inline Experiment1::plans_e&
operator|=(Experiment1::plans_e& lhs, const Experiment1::plans_e rhs) {
  lhs = lhs | rhs;
  return lhs;
}
inline Experiment1::plans_e
operator&(const Experiment1::plans_e lhs, const Experiment1::plans_e rhs) {
  using type_t = std::underlying_type_t<Experiment1::plans_e>;
  return static_cast<Experiment1::plans_e>(static_cast<type_t>(lhs) & static_cast<type_t>(rhs));
}
inline Experiment1::plans_e&
operator&=(Experiment1::plans_e& lhs, const Experiment1::plans_e rhs) {
  lhs = lhs & rhs;
  return lhs;
}

const Experiment1::plans_map_t Experiment1::_supportedPlansMap = {
  {"none",     plans_e::NONE},
  {"NONE",     plans_e::NONE},
  {"scr",      plans_e::scr},
  {"scs",      plans_e::scs},
  {"Csr",      plans_e::Csr},
  {"CsrUU",    plans_e::CsrUU},
  {"Crs",      plans_e::Crs},
  {"Nrs",      plans_e::Nrs},
  {"Nsr",      plans_e::Nsr},
  {"NrsNU",    plans_e::NrsNU},
  //{"CsrSel",   plans_e::CsrSel},
  //{"CrsSel",   plans_e::CrsSel},
  //{"NrsSel",   plans_e::NrsSel},
  //{"CsrSelLf", plans_e::CsrSelLf},
  //{"CrsSelLf", plans_e::CrsSelLf},
  //{"NrsSelLf", plans_e::NrsSelLf},
  {"all",      plans_e::ALL},
  {"ALL",      plans_e::ALL}
};

// nested classes

struct Experiment1::GlobStat {};
struct Experiment1::HashfunOnR {
  using input_t = base_tuple_t;
  using output_t = hashvalue_t;
  inline static output_t eval(const input_t* aBaseTuple) {
    return _hashfunc(aBaseTuple->k);
  }
};
struct Experiment1::HashfunOnS {
  using input_t = base_tuple_t;
  using output_t = hashvalue_t;
  inline static output_t eval(const input_t* aBaseTuple) {
    return _hashfunc(aBaseTuple->a);
  }
};
struct Experiment1::EqfunBuildR {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->k == r->k);
  }
};
struct Experiment1::EqfunBuildS {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->a == r->a);
  }
};
struct Experiment1::EqfunJoinpredRS {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->k == r->a);
  }
};
struct Experiment1::EqfunJoinpredSR {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->a == r->k);
  }
};
// for Nrs
struct Experiment1::nested_tuple_RS_t {
  base_tuple_t* _left;
  const HtNested1<base_tuple_t, HashfunOnS, EqfunBuildS>::MainNode* _right;
};
// for Nsr
struct Experiment1::nested_tuple_SR_t {
  base_tuple_t* _left;
  const HtNested1<base_tuple_t, HashfunOnR, EqfunBuildR>::MainNode* _right;
};

struct Experiment1::result_tuple_t {
  const base_tuple_t* _left;
  const base_tuple_t* _right;
};
std::ostream& operator<<(std::ostream& os, const Experiment1::result_tuple_t& t) {
  os << "[" << *t._left << "," << *t._right << "]";
  return os;
}
std::function<void(const Experiment1::result_tuple_t*, std::ostream&)>
Experiment1::_result_tuple_printer = [](const result_tuple_t* t, std::ostream& os) {
        os << "[" << *t->_left << "," << *t->_right << "]";
};
struct Experiment1::ConcatfunChaining {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  using output_t = result_tuple_t;
  inline static output_t eval(left_t* l, const right_t* r) {
    return {l, r}; 
  }
};
struct Experiment1::ConcatfunNestedRS {
  using left_t = base_tuple_t;
  using right_t = HtNested1<base_tuple_t, HashfunOnS, EqfunBuildS>::MainNode;
  using output_t = nested_tuple_RS_t;
  inline static output_t eval(left_t* l, const right_t* r) {
    return {l, r}; 
  }
};
struct Experiment1::ConcatfunNestedSR {
  using left_t = base_tuple_t;
  using right_t = HtNested1<base_tuple_t, HashfunOnR, EqfunBuildR>::MainNode;
  using output_t = nested_tuple_SR_t;
  inline static output_t eval(left_t* l, const right_t* r) {
    return {l, r}; 
  }
};
struct Experiment1::UnnestFunRS {
  public:
    using input_t  = nested_tuple_RS_t;
    using output_t = result_tuple_t;
    using MainNode = HtNested1<base_tuple_t, HashfunOnS, EqfunBuildS>::MainNode;
    using data_t   = HtNested1<base_tuple_t, HashfunOnS, EqfunBuildS>::data_t;
  public:
    inline static const MainNode* getMainNode(input_t* aNestedTuple) {
      return aNestedTuple->_right;
    }
    inline static void eval_left(output_t* out, input_t* in) {
      out->_left = in->_left;
    }
    inline static void eval_right(output_t* out, [[maybe_unused]] input_t* in, const data_t* data) {
      out->_right = data;
    }
};
struct Experiment1::UnnestFunSR {
  public:
    using input_t  = nested_tuple_SR_t;
    using output_t = result_tuple_t;
    using MainNode = HtNested1<base_tuple_t, HashfunOnR, EqfunBuildR>::MainNode;
    using data_t   = HtNested1<base_tuple_t, HashfunOnR, EqfunBuildR>::data_t;
  public:
    inline static const MainNode* getMainNode(input_t* aNestedTuple) {
      return aNestedTuple->_right;
    }
    inline static void eval_left(output_t* out, input_t* in) {
      out->_left = in->_left;
    }
    inline static void eval_right(output_t* out, [[maybe_unused]] input_t* in, const data_t* data) {
      out->_right = data;
    }
};


// public

void
Experiment1::init() {
  timerStart(__FUNCTION__);

  std::mt19937 lRng;

  // unique keys (R.k, S.k)
  timerStart(std::string(__FUNCTION__) + "::prepare key vectors");
  value_vt lKeysR;
  value_vt lKeysS;
  lKeysR.resize(cardR());
  lKeysS.resize(cardS());
  for (value_t i = 0; i < cardR(); ++i) { lKeysR.at(i) = i; }
  std::shuffle(lKeysR.begin(), lKeysR.end(), lRng);
  for (value_t i = 0; i < cardS(); ++i) { lKeysS.at(i) = i; }
  timerStop(std::string(__FUNCTION__) + "::prepare key vectors");

  // foreign key S.a
  timerStart(std::string(__FUNCTION__) + "::prepare foreign key vec");
  value_vt lForeignKeys;
  uint32_t lFkMax = getFkMax();  // max value for foreign key S.a
  GenRandIntVec lGriv;
  GenRandIntVec::param_t lGrivParam;
  // param_t(const dist_t aDist, const int aMax, const int aShift, const double aParam,
  //         const int aFlags, const int aOrder)
  // max is exclusive; order is -1 for permute and +1 for sort
  if (!_skew) {
    lGrivParam = GenRandIntVec::param_t(GenRandIntVec::dist_t::kUni, lFkMax, 0, 0.0, 0, -1);
  } else {
    lGrivParam = GenRandIntVec::param_t(GenRandIntVec::dist_t::kZipf, lFkMax, 0, 1.0, 0, -1);
  }
  lGriv.generate(lForeignKeys, cardS(), lGrivParam, lRng);  // resizes vec + generates values
  timerStop(std::string(__FUNCTION__) + "::prepare foreign key vec");

  // copy key and foreign key vector into relations/tuples
  generateR(lKeysR);
  generateS(lKeysS, lForeignKeys);

  // count number of generated distinct values in S.a
  _numDvSa = std::unordered_set(lForeignKeys.cbegin(), lForeignKeys.cend()).size();

  timerStop(__FUNCTION__);
}

void
Experiment1::printRelations() const {
  std::cout << "-- R --\n";
  for (const auto& t : _R._tuples) {
    std::cout << t.k << "|" << t.a << "|" << t.b << "\n";
  }
  std::cout << "-- S --\n";
  for (const auto& t : _S._tuples) {
    std::cout << t.k << "|" << t.a << "|" << t.b << "\n";
  }
}

void
Experiment1::printTimers() const {
  for (const auto& [desc, tp_pair] : _timePoints) {
    std::cout
      << df::infra::indent(1)
      << desc << "|"
      << df::infra::to_string(
           std::chrono::duration_cast<df::infra::milliseconds>(tp_pair.second - tp_pair.first)
         )
      << std::endl;
  }
}

void
Experiment1::printTypes() const {
  std::cout
    << "  value_t          = " << type_name<value_t>() << " (" << sizeof(value_t) << " B)\n"
    << "  hashvalue_t      = " << type_name<hashvalue_t>() << " (" << sizeof(hashvalue_t) << " B)\n"
    << "  tuple_uint32_3_t = " << type_name<tuple_uint32_3_t>() << " (" << sizeof(tuple_uint32_3_t) << " B)\n";
}

// private

void
Experiment1::generateR(const value_vt& aKeys) {
  timerStart(__FUNCTION__);
  assert(aKeys.size() >= cardR());
  _R._tuples.resize(cardR());
  for (size_t i = 0; i < cardR(); ++i) {
    _R._tuples.at(i).k = aKeys.at(i);
  }
  timerStop(__FUNCTION__);
}

void
Experiment1::generateS(const value_vt& aKeys, const value_vt& aForeignKeys) {
  timerStart(__FUNCTION__);
  assert(aKeys.size() >= cardS());
  _S._tuples.resize(cardS());
  for (size_t i = 0; i < cardS(); ++i) {
    _S._tuples.at(i).k = aKeys.at(i);
    _S._tuples.at(i).a = aForeignKeys.at(i);
  }
  timerStop(__FUNCTION__);
}

void
Experiment1::run() {
  timerStart(__FUNCTION__);
  writeCsvHeader();
  if ((plans_e::scr & _plans) != plans_e::NONE)    runScanR();
  if ((plans_e::scs & _plans) != plans_e::NONE)    runScanS();
  if ((plans_e::Csr & _plans) != plans_e::NONE)    runCsr();
  if ((plans_e::CsrUU & _plans) != plans_e::NONE)  runCsrUnknownUnique();
  if ((plans_e::Crs & _plans) != plans_e::NONE)    runCrs();
  if ((plans_e::Nsr & _plans) != plans_e::NONE)    runNsr();
  if ((plans_e::Nrs & _plans) != plans_e::NONE)    runNrs();
  if ((plans_e::NrsNU & _plans) != plans_e::NONE)  runNrsNoUnnest();
  timerStop(__FUNCTION__);
}

void
Experiment1::runScanR() {
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using top_t = AlgTop<base_tuple_t, GlobStat>;
  using scan_t = AlgScan<top_t>;
  top_t lOpTop(std::cout, false, [](const base_tuple_t* t, std::ostream& os) { os << t; } );
  scan_t lOpScanR(&lOpTop, &_R);

  time_point_t lTotalRunStart{clock_t::now()};
  lOpScanR.run(&lGs);
  time_point_t lTotalRunStop{clock_t::now()};

  // write results to csv
  writeExpParamsToCsv();
  _measureCsv
    .writeField("scr") // plan
    .writeField("NA")  // impl
    .writeField("NA")  // build
    .writeField("NA")  // probe
    .writeField("NA")  // ht_buckets
    .writeField("NA")  // ht_fracEmpty
    .writeField("NA")  // ht_cc0_avg
    .writeField("NA")  // ht_cc0_min
    .writeField("NA")  // ht_cc0_max
    .writeField("NA")  // ht_cc1_avg
    .writeField("NA")  // ht_cc1_min
    .writeField("NA")  // ht_cc1_max
    .writeField((lTotalRunStop - lTotalRunStart).count())  // t_total
    .writeField("NA")                                      // t_buildStr
    .writeField("NA")                                      // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanR.count())
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField(lOpTop.count())
    .newline();

  timerStop(__FUNCTION__);
}

void
Experiment1::runScanS() {
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using top_t = AlgTop<base_tuple_t, GlobStat>;
  using scan_t = AlgScan<top_t>;
  top_t lOpTop(std::cout, false, [](const base_tuple_t* t, std::ostream& os) { os << t; } );
  scan_t lOpScanS(&lOpTop, &_S);
  time_point_t lTotalRunStart{clock_t::now()};
  lOpScanS.run(&lGs);
  time_point_t lTotalRunStop{clock_t::now()};

  // write results to csv
  writeExpParamsToCsv();
  _measureCsv
    .writeField("scs") // plan
    .writeField("NA")  // impl
    .writeField("NA")  // build
    .writeField("NA")  // probe
    .writeField("NA")  // ht_buckets
    .writeField("NA")  // ht_fracEmpty
    .writeField("NA")  // ht_cc0_avg
    .writeField("NA")  // ht_cc0_min
    .writeField("NA")  // ht_cc0_max
    .writeField("NA")  // ht_cc1_avg
    .writeField("NA")  // ht_cc1_min
    .writeField("NA")  // ht_cc1_max
    .writeField((lTotalRunStop - lTotalRunStart).count())  // t_total
    .writeField("NA")                                      // t_buildStr
    .writeField("NA")                                      // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanS.count())
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField("NA")
    .writeField(lOpTop.count())
    .newline();

  timerStop(__FUNCTION__);
}

void
Experiment1::runCsr() {
  /*
   * Chaining hashtable w/ build on R, probe on S
   * Probe terminates early as key property of R.k from build side is known.
   *
   *          Top
   *           |
   *   [result_tuple_t]
   *           |
   *         Probe.....(S.a = R.k).....Build
   *           |                         |
   *         Scan                      Scan
   *           |                         |
   *           S                         R
   *    [base_tuple_t]            [base_tuple_t]
   */
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using build_t = AlgHashJoinBuild<HashfunOnR, EqfunBuildR, GlobStat>;
  using scan_R_t = AlgScan<build_t>;
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using probe_t = AlgHashJoinProbe
    <top_t, build_t, HashfunOnS, EqfunJoinpredSR, ConcatfunChaining, true>;
  using scan_S_t = AlgScan<probe_t>;

  // build on R.k => #dv(R.k) = |R| (key)
  // choose #buckets for HT as #dv(R.k) / b
  uint32_t lHtNumBuckets = std::max((cardR() / _b), 1UL);
  build_t lOpBuild(lHtNumBuckets, _log2RsvChunkSize);
  scan_R_t lOpScanR(&lOpBuild, &_R);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  probe_t lOpProbe(&lOpTop, &lOpBuild);
  scan_S_t lOpScanS(&lOpProbe, &_S);

  std::chrono::nanoseconds lDurationBuild{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(), [&lDurationTotal, &lDurationBuild, &lDurationProbe, &lOpScanR, &lOpScanS, &lGs] () {
        time_point_t lTotalRunStart{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTotalRunMid{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTotalRunStop{clock_t::now()};

        lDurationBuild += (lTotalRunMid - lTotalRunStart);
        lDurationProbe += (lTotalRunStop - lTotalRunMid);
        lDurationTotal += (lTotalRunStop - lTotalRunStart);
      },
      [&lOpBuild] () { lOpBuild.clear_ht(); }, false, minRepeat()
    );

  //size_t lNumRepeat = 1;
  //for (uint32_t i = 0; i < lNumRepeat; ++i) {
  //  time_point_t lTotalRunStart{clock_t::now()};
  //  lOpScanR.run(&lGs);
  //  time_point_t lTotalRunMid{clock_t::now()};
  //  lOpScanS.run(&lGs);
  //  time_point_t lTotalRunStop{clock_t::now()};

  //  lDurationBuild += (lTotalRunMid - lTotalRunStart);
  //  lDurationProbe += (lTotalRunStop - lTotalRunMid);
  //  lDurationTotal += (lTotalRunStop - lTotalRunStart);

  //  // empty the hash table for the next iteration, otherwise adding tuples twice!
  //  if (i != repeat() - 1) {
  //    // in every but the last iteration -> retain last iter's HT for statistics etc.
  //    lOpBuild.clear_ht();
  //  }
  //}

  lDurationBuild /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan Csr\n";
    std::cout << "  Build Strand\n";
    print_strand(&lOpScanR, 2);
    std::cout << "  Probe Strand\n";
    print_strand(&lOpScanS, 2);
    std::cout << "  sizeof(Node): " << sizeof(build_t::hashtable_t::Node) << "\n";
  }

  const HtStatistics lHtStat = lOpBuild.hashtable().makeStatistics();

  // write results to csv
  writeExpParamsToCsv();
  _measureCsv
    .writeField("Csr")
    .writeField("chaining")
    .writeField("R")
    .writeField("S")
    .writeField(lOpBuild.hashtable().numBuckets())
    .writeField(lHtStat.fracEmptyBuckets())                // ht_fracEmpty
    .writeField(lHtStat._collisionChainLen.avg())          // ht_cc0_avg
    .writeField(lHtStat._collisionChainLen.min())          // ht_cc0_min
    .writeField(lHtStat._collisionChainLen.max())          // ht_cc0_max
    .writeField(lHtStat._collisionChainLenNonempty.avg())  // ht_cc1_avg
    .writeField(lHtStat._collisionChainLenNonempty.min())  // ht_cc1_min
    .writeField(lHtStat._collisionChainLenNonempty.max())  // ht_cc1_max
    .writeField(lIterations)                               // reps
    .writeField(lDurationTotal.count())                    // t_total
    .writeField(lDurationBuild.count())                    // t_buildStr
    .writeField(lDurationProbe.count())                    // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanR.count())
    .writeField("NA")
    .writeField(lOpBuild.count())
    .writeField(lOpScanS.count())
    .writeField("NA")
    .writeField(lOpProbe.count())
    .writeField(lOpProbe.numCmps())
    .writeField("NA")
    .writeField(lOpTop.count())
    .newline();

  timerStop(__FUNCTION__);
}

void
Experiment1::runCsrUnknownUnique() {
  /*
   * Chaining hashtable w/ build on R, probe on S
   * Probe does not terminate early (R.k key property from build side unknown)
   *
   *          Top
   *           |
   *   [result_tuple_t]
   *           |
   *         Probe.....(S.a = R.k).....Build
   *           |                         |
   *         Scan                      Scan
   *           |                         |
   *           S                         R
   *    [base_tuple_t]            [base_tuple_t]
   */
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using build_t = AlgHashJoinBuild<HashfunOnR, EqfunBuildR, GlobStat>;
  using scan_R_t = AlgScan<build_t>;
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using probe_t = AlgHashJoinProbe
    <top_t, build_t, HashfunOnS, EqfunJoinpredSR, ConcatfunChaining>;
  using scan_S_t = AlgScan<probe_t>;

  // build on R.k => #dv(R.k) = |R| (key)
  // choose #buckets for HT as #dv(R.k) / b
  uint32_t lHtNumBuckets = std::max((cardR() / _b), 1UL);
  build_t lOpBuild(lHtNumBuckets, _log2RsvChunkSize);
  scan_R_t lOpScanR(&lOpBuild, &_R);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  probe_t lOpProbe(&lOpTop, &lOpBuild);
  scan_S_t lOpScanS(&lOpProbe, &_S);

  std::chrono::nanoseconds lDurationBuild{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(), [&lDurationTotal, &lDurationBuild, &lDurationProbe, &lOpScanR, &lOpScanS, &lGs] () {
        time_point_t lTotalRunStart{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTotalRunMid{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTotalRunStop{clock_t::now()};

        lDurationBuild += (lTotalRunMid - lTotalRunStart);
        lDurationProbe += (lTotalRunStop - lTotalRunMid);
        lDurationTotal += (lTotalRunStop - lTotalRunStart);
      },
      [&lOpBuild]() { lOpBuild.clear_ht(); }, false, minRepeat()
    );

  lDurationBuild /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan CsrUU\n";
    std::cout << "  Build Strand\n";
    print_strand(&lOpScanR, 2);
    std::cout << "  Probe Strand\n";
    print_strand(&lOpScanS, 2);
    std::cout << "  sizeof(Node): " << sizeof(build_t::hashtable_t::Node) << "\n";
  }

  const HtStatistics lHtStat = lOpBuild.hashtable().makeStatistics();

  // write results to csv
  writeExpParamsToCsv();
  _measureCsv
    .writeField("CsrUU")
    .writeField("chaining")
    .writeField("R")
    .writeField("S")
    .writeField(lOpBuild.hashtable().numBuckets())
    .writeField(lHtStat.fracEmptyBuckets())                // ht_fracEmpty
    .writeField(lHtStat._collisionChainLen.avg())          // ht_cc0_avg
    .writeField(lHtStat._collisionChainLen.min())          // ht_cc0_min
    .writeField(lHtStat._collisionChainLen.max())          // ht_cc0_max
    .writeField(lHtStat._collisionChainLenNonempty.avg())  // ht_cc1_avg
    .writeField(lHtStat._collisionChainLenNonempty.min())  // ht_cc1_min
    .writeField(lHtStat._collisionChainLenNonempty.max())  // ht_cc1_max
    .writeField(lIterations)                               // reps
    .writeField(lDurationTotal.count())                    // t_total
    .writeField(lDurationBuild.count())                    // t_buildStr
    .writeField(lDurationProbe.count())                    // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanR.count())
    .writeField("NA")
    .writeField(lOpBuild.count())
    .writeField(lOpScanS.count())
    .writeField("NA")
    .writeField(lOpProbe.count())
    .writeField(lOpProbe.numCmps())
    .writeField("NA")
    .writeField(lOpTop.count())
    .newline();

  timerStop(__FUNCTION__);
}

void
Experiment1::runCrs() {
  /*
   * Chaining hashtable w/ build on S, probe on R
   *
   *          Top
   *           |
   *   [result_tuple_t]
   *           |
   *         Probe.....(R.k = S.a).....Build
   *           |                         |
   *         Scan                      Scan
   *           |                         |
   *           R                         S
   *    [base_tuple_t]            [base_tuple_t]
   */
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using build_t = AlgHashJoinBuild<HashfunOnS, EqfunBuildS, GlobStat>;
  using scan_S_t = AlgScan<build_t>;
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using probe_t = AlgHashJoinProbe
    <top_t, build_t, HashfunOnR, EqfunJoinpredRS, ConcatfunChaining>;
  using scan_R_t = AlgScan<probe_t>;

  uint32_t lHtNumBuckets = std::max((_numDvSa / _b), 1UL);
  build_t lOpBuild(lHtNumBuckets, _log2RsvChunkSize);  // build on S.a, choose #buckets := #dv(S.a) / b
  scan_S_t lOpScanS(&lOpBuild, &_S);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  probe_t lOpProbe(&lOpTop, &lOpBuild);
  scan_R_t lOpScanR(&lOpProbe, &_R);

  std::chrono::nanoseconds lDurationBuild{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(), [&lDurationTotal, &lDurationBuild, &lDurationProbe, &lOpScanR, &lOpScanS, &lGs] () {
        time_point_t lTotalRunStart{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTotalRunMid{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTotalRunStop{clock_t::now()};

        lDurationBuild += (lTotalRunMid - lTotalRunStart);
        lDurationProbe += (lTotalRunStop - lTotalRunMid);
        lDurationTotal += (lTotalRunStop - lTotalRunStart);
      },
      [&lOpBuild] () { lOpBuild.clear_ht(); }, false, minRepeat()
    );

  //for (uint32_t i = 0; i < repeat(); ++i) {
  //  time_point_t lTotalRunStart{clock_t::now()};
  //  lOpScanS.run(&lGs);
  //  time_point_t lTotalRunMid{clock_t::now()};
  //  lOpScanR.run(&lGs);
  //  time_point_t lTotalRunStop{clock_t::now()};

  //  lDurationBuild += (lTotalRunMid - lTotalRunStart);
  //  lDurationProbe += (lTotalRunStop - lTotalRunMid);
  //  lDurationTotal += (lTotalRunStop - lTotalRunStart);

  //  // empty the hash table for the next iteration, otherwise adding tuples twice!
  //  if (i != repeat() - 1) {
  //    // in every but the last iteration -> retain last iter's HT for statistics etc.
  //    lOpBuild.clear_ht();
  //  }
  //}

  lDurationBuild /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan Crs\n";
    std::cout << "  Build Strand\n";
    print_strand(&lOpScanS, 2);
    std::cout << "  Probe Strand\n";
    print_strand(&lOpScanR, 2);
    std::cout << "  sizeof(Node): " << sizeof(build_t::hashtable_t::Node) << "\n";
  }

  const HtStatistics lHtStat = lOpBuild.hashtable().makeStatistics();

  // write results to csv
  writeExpParamsToCsv();
  _measureCsv
    .writeField("Crs")
    .writeField("chaining")
    .writeField("S")
    .writeField("R")
    .writeField(lOpBuild.hashtable().numBuckets())
    .writeField(lHtStat.fracEmptyBuckets())                // ht_fracEmpty
    .writeField(lHtStat._collisionChainLen.avg())          // ht_cc0_avg
    .writeField(lHtStat._collisionChainLen.min())          // ht_cc0_min
    .writeField(lHtStat._collisionChainLen.max())          // ht_cc0_max
    .writeField(lHtStat._collisionChainLenNonempty.avg())  // ht_cc1_avg
    .writeField(lHtStat._collisionChainLenNonempty.min())  // ht_cc1_min
    .writeField(lHtStat._collisionChainLenNonempty.max())  // ht_cc1_max
    .writeField(lIterations)                               // reps
    .writeField(lDurationTotal.count())                    // t_total
    .writeField(lDurationBuild.count())                    // t_buildStr
    .writeField(lDurationProbe.count())                    // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanS.count())
    .writeField("NA")
    .writeField(lOpBuild.count())
    .writeField(lOpScanR.count())
    .writeField("NA")
    .writeField(lOpProbe.count())
    .writeField(lOpProbe.numCmps())
    .writeField("NA")
    .writeField(lOpTop.count())
    .newline();

  timerStop(__FUNCTION__);
}

void
Experiment1::runNrs() {
  /*
   * Nested hashtable w/ build on S, probe on R
   *
   *          Top
   *           |
   *   [result_tuple_t]
   *           |
   *        Unnest
   *           |
   *  [nested_tuple_RS_t]
   *           |
   *        Nested                    Nested
   *         Probe.....(R.k = S.a).....Build
   *           |                         |
   *         Scan                      Scan
   *           |                         |
   *           R                         S
   *    [base_tuple_t]            [base_tuple_t]
   */
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using build_t = AlgNestJoinBuild<HashfunOnS, EqfunBuildS, GlobStat>;
  using scan_S_t = AlgScan<build_t>;
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using unnest_t = AlgUnnestHt<top_t, UnnestFunRS, build_t::hashtable_t>;
  using probe_t = AlgNestJoinProbe
    <unnest_t, build_t, HashfunOnR, EqfunJoinpredRS, ConcatfunNestedRS>;
  using scan_R_t = AlgScan<probe_t>;

  // build on S.a, choose #buckets := #dv(S.a) / b
  uint32_t lHtNumBuckets = std::max((_numDvSa / _b), 1UL);
  build_t lOpBuild(lHtNumBuckets, _log2RsvChunkSize, _log2RsvChunkSize);
  scan_S_t lOpScanS(&lOpBuild, &_S);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  unnest_t lOpUnnest(&lOpTop);
  probe_t lOpProbe(&lOpUnnest, &lOpBuild);
  scan_R_t lOpScanR(&lOpProbe, &_R);

  std::chrono::nanoseconds lDurationBuild{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(), [&lDurationTotal, &lDurationBuild, &lDurationProbe, &lOpScanR, &lOpScanS, &lGs] () {
        time_point_t lTotalRunStart{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTotalRunMid{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTotalRunStop{clock_t::now()};

        lDurationBuild += (lTotalRunMid - lTotalRunStart);
        lDurationProbe += (lTotalRunStop - lTotalRunMid);
        lDurationTotal += (lTotalRunStop - lTotalRunStart);
      },
      [&lOpBuild] () { lOpBuild.clear_ht(); }, false, minRepeat()
    );

  lDurationBuild /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;


  if (_trace) {
    std::cout << "Plan Nrs\n";
    std::cout << "  Build Strand\n";
    print_strand(&lOpScanS, 2);
    std::cout << "  Probe Strand\n";
    print_strand(&lOpScanR, 2);
    std::cout << "  sizeof(MainNode): " << sizeof(build_t::hashtable_t::MainNode) << "\n";
    std::cout << "  sizeof(SubNode):  " << sizeof(build_t::hashtable_t::SubNode) << "\n";
  }

  const HtStatistics lHtStat = lOpBuild.hashtable().makeStatistics();

  writeExpParamsToCsv();
  _measureCsv
    .writeField("Nrs")
    .writeField("nested")
    .writeField("S")
    .writeField("R")
    .writeField(lOpBuild.hashtable().numBuckets())
    .writeField(lHtStat.fracEmptyBuckets())                // ht_fracEmpty
    .writeField(lHtStat._collisionChainLen.avg())          // ht_cc0_avg
    .writeField(lHtStat._collisionChainLen.min())          // ht_cc0_min
    .writeField(lHtStat._collisionChainLen.max())          // ht_cc0_max
    .writeField(lHtStat._collisionChainLenNonempty.avg())  // ht_cc1_avg
    .writeField(lHtStat._collisionChainLenNonempty.min())  // ht_cc1_min
    .writeField(lHtStat._collisionChainLenNonempty.max())  // ht_cc1_max
    .writeField(lIterations)                               // reps
    .writeField(lDurationTotal.count())                    // t_total
    .writeField(lDurationBuild.count())                    // t_buildStr
    .writeField(lDurationProbe.count())                    // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanS.count())
    .writeField("NA")
    .writeField(lOpBuild.count())
    .writeField(lOpScanR.count())
    .writeField("NA")
    .writeField(lOpProbe.count())
    .writeField(lOpProbe.numCmps())
    .writeField(lOpUnnest.count())
    .writeField(lOpTop.count())
    .newline();
  timerStop(__FUNCTION__);
}

void
Experiment1::runNsr() {
  /*
   * Nested hashtable w/ build on R, probe on S
   *
   *          Top
   *           |
   *   [result_tuple_t]
   *           |
   *        Unnest
   *           |
   *  [nested_tuple_SR_t]
   *           |
   *        Nested                    Nested
   *         Probe.....(S.a = R.k).....Build
   *           |                         |
   *         Scan                      Scan
   *           |                         |
   *           S                         R
   *    [base_tuple_t]            [base_tuple_t]
   */
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using build_t = AlgNestJoinBuild<HashfunOnR, EqfunBuildR, GlobStat>;
  using scan_R_t = AlgScan<build_t>;
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using unnest_t = AlgUnnestHt<top_t, UnnestFunSR, build_t::hashtable_t>;
  using probe_t = AlgNestJoinProbe
    <unnest_t, build_t, HashfunOnS, EqfunJoinpredSR, ConcatfunNestedSR>;
  using scan_S_t = AlgScan<probe_t>;

  // build on R.k => #dv(R.k) = |R| (key)
  // choose #buckets for HT as #dv(R.k) / b
  uint32_t lHtNumBuckets = std::max((cardR() / _b), 1UL);
  build_t lOpBuild(lHtNumBuckets, _log2RsvChunkSize, _log2RsvChunkSize);
  scan_R_t lOpScanR(&lOpBuild, &_R);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  unnest_t lOpUnnest(&lOpTop);
  probe_t lOpProbe(&lOpUnnest, &lOpBuild);
  scan_S_t lOpScanS(&lOpProbe, &_S);

  std::chrono::nanoseconds lDurationBuild{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(), [&lDurationTotal, &lDurationBuild, &lDurationProbe, &lOpScanR, &lOpScanS, &lGs] () {
        time_point_t lTotalRunStart{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTotalRunMid{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTotalRunStop{clock_t::now()};

        lDurationBuild += (lTotalRunMid - lTotalRunStart);
        lDurationProbe += (lTotalRunStop - lTotalRunMid);
        lDurationTotal += (lTotalRunStop - lTotalRunStart);
      },
      [&lOpBuild] () { lOpBuild.clear_ht(); }, false, minRepeat()
    );

  lDurationBuild /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan Nsr\n";
    std::cout << "  Build Strand\n";
    print_strand(&lOpScanR, 2);
    std::cout << "  Probe Strand\n";
    print_strand(&lOpScanS, 2);
    std::cout << "  sizeof(MainNode): " << sizeof(build_t::hashtable_t::MainNode) << "\n";
    std::cout << "  sizeof(SubNode):  " << sizeof(build_t::hashtable_t::SubNode) << "\n";
  }

  const HtStatistics lHtStat = lOpBuild.hashtable().makeStatistics();

  writeExpParamsToCsv();
  _measureCsv
    .writeField("Nsr")                                     // plan
    .writeField("nested")                                  // ht_impl
    .writeField("R")                                       // build
    .writeField("S")                                       // probe
    .writeField(lOpBuild.hashtable().numBuckets())
    .writeField(lHtStat.fracEmptyBuckets())                // ht_fracEmpty
    .writeField(lHtStat._collisionChainLen.avg())          // ht_cc0_avg
    .writeField(lHtStat._collisionChainLen.min())          // ht_cc0_min
    .writeField(lHtStat._collisionChainLen.max())          // ht_cc0_max
    .writeField(lHtStat._collisionChainLenNonempty.avg())  // ht_cc1_avg
    .writeField(lHtStat._collisionChainLenNonempty.min())  // ht_cc1_min
    .writeField(lHtStat._collisionChainLenNonempty.max())  // ht_cc1_max
    .writeField(lIterations)                               // reps
    .writeField(lDurationTotal.count())                    // t_total
    .writeField(lDurationBuild.count())                    // t_buildStr
    .writeField(lDurationProbe.count())                    // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanR.count())                          // c_scanBuild
    .writeField("NA")
    .writeField(lOpBuild.count())                          // c_htBuild
    .writeField(lOpScanS.count())                          // c_scanProbe
    .writeField("NA")
    .writeField(lOpProbe.count())                          // c_htProbe
    .writeField(lOpProbe.numCmps())                        // c_htProbeCmp
    .writeField(lOpUnnest.count())                         // c_unnest
    .writeField(lOpTop.count())                            // c_top
    .newline();
  timerStop(__FUNCTION__);
}

void
Experiment1::runNrsNoUnnest() {
  /*
   * Nested hashtable w/ build on S, probe on R, no unnesting
   *
   *          Top
   *           |
   *  [nested_tuple_RS_t]
   *           |
   *        Nested                    Nested
   *         Probe.....(R.k = S.a).....Build
   *           |                         |
   *         Scan                      Scan
   *           |                         |
   *           R                         S
   *    [base_tuple_t]            [base_tuple_t]
   */
  timerStart(__FUNCTION__);
  GlobStat lGs;
  using build_t = AlgNestJoinBuild<HashfunOnS, EqfunBuildS, GlobStat>;
  using scan_S_t = AlgScan<build_t>;
  using top_t = AlgTop<nested_tuple_RS_t, GlobStat>;
  using probe_t = AlgNestJoinProbe
    <top_t, build_t, HashfunOnR, EqfunJoinpredRS, ConcatfunNestedRS>;
  using scan_R_t = AlgScan<probe_t>;

  // build on S.a, choose #buckets := #dv(S.a) / b
  uint32_t lHtNumBuckets = std::max((_numDvSa / _b), 1UL);
  build_t lOpBuild(lHtNumBuckets, _log2RsvChunkSize, _log2RsvChunkSize);
  scan_S_t lOpScanS(&lOpBuild, &_S);

  top_t lOpTop(std::cout, false);
  probe_t lOpProbe(&lOpTop, &lOpBuild);
  scan_R_t lOpScanR(&lOpProbe, &_R);

  std::chrono::nanoseconds lDurationBuild{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(), [&lDurationTotal, &lDurationBuild, &lDurationProbe, &lOpScanR, &lOpScanS, &lGs] () {
        time_point_t lTotalRunStart{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTotalRunMid{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTotalRunStop{clock_t::now()};

        lDurationBuild += (lTotalRunMid - lTotalRunStart);
        lDurationProbe += (lTotalRunStop - lTotalRunMid);
        lDurationTotal += (lTotalRunStop - lTotalRunStart);
      },
      [&lOpBuild] () { lOpBuild.clear_ht(); }, false, minRepeat()
    );

  lDurationBuild /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan Nrs w/o Unnesting\n";
    std::cout << "  Build Strand\n";
    print_strand(&lOpScanS, 2);
    std::cout << "  Probe Strand\n";
    print_strand(&lOpScanR, 2);
  }

  const HtStatistics lHtStat = lOpBuild.hashtable().makeStatistics();

  writeExpParamsToCsv();
  _measureCsv
    .writeField("NrsNU")
    .writeField("nested")
    .writeField("S")
    .writeField("R")
    .writeField(lOpBuild.hashtable().numBuckets())
    .writeField(lHtStat.fracEmptyBuckets())                // ht_fracEmpty
    .writeField(lHtStat._collisionChainLen.avg())          // ht_cc0_avg
    .writeField(lHtStat._collisionChainLen.min())          // ht_cc0_min
    .writeField(lHtStat._collisionChainLen.max())          // ht_cc0_max
    .writeField(lHtStat._collisionChainLenNonempty.avg())  // ht_cc1_avg
    .writeField(lHtStat._collisionChainLenNonempty.min())  // ht_cc1_min
    .writeField(lHtStat._collisionChainLenNonempty.max())  // ht_cc1_max
    .writeField(lIterations)                               // reps
    .writeField(lDurationTotal.count())                    // t_total
    .writeField(lDurationBuild.count())                    // t_buildStr
    .writeField(lDurationProbe.count())                    // t_probeStr
    .writeField(get_runtime_excl(&lOpTop).count())         // t_top
    .writeField(lOpScanS.count())
    .writeField("NA")
    .writeField(lOpBuild.count())
    .writeField(lOpScanR.count())
    .writeField("NA")
    .writeField(lOpProbe.count())
    .writeField(lOpProbe.numCmps())
    .writeField("NA")
    .writeField(lOpTop.count())
    .newline();
  timerStop(__FUNCTION__);
}


void
Experiment1::writeCsvHeader() {
  _measureCsv
    .writeField("mintime")
    .writeField("minreps")
    .writeField("log2CardR")
    .writeField("log2CardS")
    .writeField("skew")
    .writeField("t")
    .writeField("fkMax")
    .writeField("numDvSa")
    .writeField("b")
    // --
    .writeField("plan")
    .writeField("ht_impl")
    .writeField("build")
    .writeField("probe")
    .writeField("ht_buckets")
    .writeField("ht_fracEmpty")  // fraction of empty buckets in the hashtable
    .writeField("cc0_avg")  // hash table collision chain length avg/min/max
    .writeField("cc0_min")  // cc0: all buckets
    .writeField("cc0_max")
    .writeField("cc1_avg")  // cc1: nonempty buckets
    .writeField("cc1_min")
    .writeField("cc1_max")
    .writeField("reps")     // number of repetitions of plan execution
    .writeField("t_total")
    .writeField("t_buildStr")  // build strand
    .writeField("t_probeStr")  // probe strand
    //.writeField("t_scanBuild")
    //.writeField("t_htBuild")
    //.writeField("t_scanProbe")
    //.writeField("t_htProbe")
    //.writeField("t_unnest")
    .writeField("t_top")
    .writeField("c_scanBuild")
    .writeField("c_selBuild")
    .writeField("c_htBuild")
    .writeField("c_scanProbe")
    .writeField("c_selProbe")
    .writeField("c_htProbe")
    .writeField("c_htProbeCmp")
    .writeField("c_unnest")
    .writeField("c_top")
    .newline();
}

void
Experiment1::writeExpParamsToCsv() {
  _measureCsv
    .writeField(df::infra::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(minRuntime())))
    .writeField(minRepeat())
    .writeField(_log2CardR)
    .writeField(_log2CardS)
    .writeField(_skew)
    .writeField(_t)
    .writeField(getFkMax())
    .writeField(_numDvSa)
    .writeField(_b);
}

void
Experiment1::timerCtl(const std::string& aDesc, const bool aIsStop) {
  time_point_t lNow = clock_t::now();
  if (!aIsStop) {
    _timePoints.insert({aDesc, {lNow, time_point_t{}}});
  } else {
    assert(_timePoints.contains(aDesc));
    auto& [desc, tp_pair] = *(_timePoints.find(aDesc));  // find returns pair<key_t, value_t>
    tp_pair.second = lNow;
  }
}

void
Experiment1::plansFromVec(const std::vector<std::string>& aPlanVec) {
  using df::infra::to_lower;
  _plans = plans_e::NONE;
  for (auto p : aPlanVec) {
    if (_supportedPlansMap.contains(p)) {
      _plans |= _supportedPlansMap.at(p);
    }
  }
}

// ############################################################################

int main(int argc, char** argv) {
  CLI::App lApp{"Hash Table Experiment 1"};


  uint32_t                 lLog2CardR;
  uint32_t                 lLog2CardS;
  bool                     lSkew;
  uint32_t                 lParamT;  // scaledown of dist values in S.a
  uint32_t                 lParamB = 1;  // scaledown of buckets of hashtable
  std::filesystem::path    lMeasureFile;
  std::vector<std::string> lPlansToRun = {"all"};
  bool                     lPrintTimers = false;
  bool                     lPrintRelations = false;


  lApp.add_option("-R,--card-R", lLog2CardR, "Cardinality of key relation R as log2")
    ->required()
    ->check(CLI::Range(0, 30));
  lApp.add_option("-S,--card-S", lLog2CardS, "Cardinality of foreign key relation S as log2")
    ->required()
    ->check(CLI::Range(0, 30));
  lApp.add_flag("--skew,!--no-skew", lSkew, "Skewed distribution of foreign keys (Zipf) or not (uniform)")
    ->required();
  lApp.add_option("-t,--param-t", lParamT, "Parameter t, scales the max distinct value of the foreign key column based on --card-R")
    ->required()
    ->check(CLI::Range(0, 9));
  lApp.add_option("-b,--param-b", lParamB,
      "Parameter b, scales down the number of buckets of the build hash table \
      wrt the number of distict values of the build side")
    ->check(CLI::Range(1, 4));
  lApp.add_option("--measure-file", lMeasureFile, "CSV file where the measurements are written to")
    ->required();
  lApp.add_option("-p,--plans",
      [&lPlansToRun](std::vector<std::string> val) {
	if (val.size() > 0 && val.at(0).size() > 0) { lPlansToRun.clear(); }  // remove default "all"
        std::string item;
        for (const auto &v : val) {
          std::stringstream ss(v);
          while (std::getline(ss, item, ',')) {
            if (!item.empty()) { lPlansToRun.push_back(item); }
          }
        }
        return true;
      },
      "Comma-separated list of the plans to run");
  //->check(CLI::IsMember(Experiment1::getSupportedPlans()));  // XXX doesn't work yet: checks whole string (with commas) instead of parsed one
  //
  lApp.add_flag("--print-timers,!--no-print-timers", lPrintTimers, "Print timers after experiment run");
  lApp.add_flag("--print-relations,!--no-print-relations", lPrintRelations, "Print generated relations R and S");
  
  CLI11_PARSE(lApp, argc, argv);
  if (lParamT > lLog2CardR) {
    std::cerr << "--param-t must not be greater than --card-R\n";
    return EXIT_FAILURE;
  }

  std::cout << "Running Experiment 1 with the following config:\n"
    << df::infra::indent(1) << "--card-R " << lLog2CardR << "\n"
    << df::infra::indent(1) << "--card-S " << lLog2CardS << "\n"
    << std::boolalpha
    << df::infra::indent(1) << "--skew " << lSkew << "\n"
    << df::infra::indent(1) << "--param-t " << lParamT << "\n"
    << df::infra::indent(1) << "--param-b " << lParamB << "\n"
    << df::infra::indent(1) << "--measure-file " << lMeasureFile << "\n"
    << df::infra::indent(1) << "--plans ";
  std::for_each(lPlansToRun.cbegin(), lPlansToRun.cend(), [](const std::string& x) { std::cout << x << ","; });
  std::cout << "\n";


  Experiment1 lEx1(lLog2CardR, lLog2CardS, lSkew, lParamT, lParamB, lMeasureFile, lPlansToRun);
  lEx1.init();
  if (lPrintRelations) {
    lEx1.printRelations();
  }

  lEx1.run();

  if (lPrintTimers) {
    std::cout << "Timers:\n";
    lEx1.printTimers();
  }

  std::cout << "Types:\n";
  lEx1.printTypes();

  std::cout << "----" << std::endl;
  return EXIT_SUCCESS;
}

#include "util/math.hh"
#include "util/standard_includes.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <limits>
#include <ostream>
#include <random>
#include <unordered_set>

#include "lib/CLI11.hpp"

#include "util/GenRandIntVec.hh"
#include "util/chrono_helpers.hh"
#include "util/csv_writer.hh"
#include "util/hasht.hh"
#include "util/measure_helpers.hh"
#include "util/output_helpers.hh"
#include "util/string_helpers.hh"
#include "util/zipf_distribution.hh"

#include "algebra.hh"
#include "ht_chaining.hh"
#include "ht_nested.hh"

/*
 * Measure the runtime of different query execution plans involving two key/foreign key joins
 * between three relations on integer attributes.
 * Compare a plan with standard hash joins to plans with 3D hash joins and deferred unnesting.
 *
 * Goal: Show that using the 3D hash join with DEFERRED UNNESTING is beneficial.
 *
 * Design/Idea:
 * - R "large"; S, T "small(er)"
 * - join(R, S) and join(R, T) "large", join(R, S, T) "small"
 *   -> each single join produces a certain amount of output tuples,
 *      many of which are dropped by the second join.
 * - Use a hash join: build hash tables on S, T; probe with R.
 *
 * Relations/Schema: 'inverted star'
 * - R := {[k int]}           => central key relation
 *   R.k \in [0, |R| - 1], key
 * - S := {[k int, a int]}   => foreign key relation
 *   S.k \in [0, |S| - 1], key
 *   S.a foreign key references R.k
 * - T := {[k int, b int]}   => foreign key relation
 *   T.k \in [0, |S| - 1], key
 *   T.b foreign key references R.k
 *
 * Hash function: murmur32 (from hasht.hh)
 *
 * Generation of foreign keys S.a, T.b
 *   (1)     Fraction |R| * alpha is referenced by both S and T (S.a = T.b).
 *   (2)+(3) Fraction |R| * beta is referenced only by either S or T (disjoint: S.a \cap T.b = \emptyset)
 *   (4)     Fraction |R| * (1 - (alpha + 2 * beta) is not referenced at all.
 * 
 * (2)+(3) have a multiplicity of c,
 * (1) has a multiplicity of 1.
 *
 * 0                                                   |R|-1
 * |-------------------------------------------------------|
 * |     |              .              |                   |
 * | (1) |      (2)     .      (3)     |        (4)        |
 * |     |              .              |                   |
 * |-------------------------------------------------------|
 * |.....| <- alpha
 *       |..............:..............| <- 2*beta
 *                                     |...................| <- 1 - (alpha + 2*beta)
 *
 * Plans:
 * - Ndu (nested hash join w/ deferred unnesting), Nnu (no unnesting at all, i.e., skip (*))
 *
 *            Top
 *             |
 *     [result_tuple_t]
 *             |
 *          Unnest*
 *             |
 *          Unnest*
 *             |
 *   [nested_tuple_RST_t]
 *             |
 *          Nested                                              Nested
 *           Probe.................(R.k = T.a)..................Build
 *             |                                                   |
 *    [nested_tuple_RS_t]                                          |
 *             |                                                   |
 *          Nested                    Nested                       |
 *           Probe.....(R.k = S.a).....Build                       |
 *             |                         |                         |
 *           Scan                      Scan                      Scan
 *             |                         |                         |
 *             R                         S                         T
 *      [base_tuple_t]            [base_tuple_t]            [base_tuple_t]
 *
 *
 * - Niu (nested hash join w/ immediate unnesting)
 *
 *            Top
 *             |
 *     [result_tuple_t]
 *             |
 *          Unnest
 *             |
 *   [partial_nested_tuple_RST_t]
 *             |
 *          Nested                                              Nested
 *           Probe.................(R.k = T.a)..................Build
 *             |                                                   |
 *    [result_tuple_RS_t]                                          |
 *             |                                                   |
 *          Unnest                                                 |
 *             |                                                   |
 *    [nested_tuple_RS_t]                                          |
 *             |                                                   |
 *          Nested                    Nested                       |
 *           Probe.....(R.k = S.a).....Build                       |
 *             |                         |                         |
 *           Scan                      Scan                      Scan
 *             |                         |                         |
 *             R                         S                         T
 *      [base_tuple_t]            [base_tuple_t]            [base_tuple_t]
 *
 *
 * - C (normal 'chaining' hash join)
 *
 *            Top
 *             |
 *     [result_tuple_t]
 *             |
 *           Probe.................(R.k = T.a)..................Build
 *             |                                                   |
 *    [result_tuple_RS_t]                                          |
 *             |                                                   |
 *          Nested                    Nested                       |
 *           Probe.....(R.k = S.a).....Build                       |
 *             |                         |                         |
 *           Scan                      Scan                      Scan
 *             |                         |                         |
 *             R                         S                         T
 *      [base_tuple_t]            [base_tuple_t]            [base_tuple_t]
 */

// all base tuples are just pairs of uint32_t's, i.e., this is how the tuples in R, S and T look like
struct tuple_uint32_2_t { uint32_t k, a; };
std::ostream& operator<<(std::ostream& os, const tuple_uint32_2_t& t) {
  os << "[" << t.k << "," << t.a << "]";
  return os;
}

class Experiment4 {
  public:
    enum class plans_e {
      NONE     = 0,
      Ndu      = 1,
      Nnu      = 2,
      Chj      = 4,
      ALL      = (Chj << 1) - 1
    };
  public:
    using base_tuple_t = tuple_uint32_2_t;
    using relation_t   = RelationRS<base_tuple_t>;
    using attrvalue_t  = uint32_t;
    using attrvalue_vt = std::vector<uint32_t>;
    using clock_t      = std::chrono::steady_clock;
    using time_point_t = std::chrono::time_point<clock_t>;
    using plans_map_t  = std::unordered_map<std::string, plans_e>;
  public:
    struct result_tuple_t;        // [r, s0, s1]
    struct result_tuple_RS_t;     // C: [r, s0]
    struct nested_tuple_RS_t;     // N: [r, {s0 | joinpred(r, s0)}]
    struct nested_tuple_RST_t;    // N: [r, {s0 | joinpred(r, s0)}, {s1 | joinpred(r, s1)}]
    struct tuple_R_nS_xT_t;       // N: [r, {s0 | joinpred(r, s0)}, s1] 
  public:
    struct GlobStat;
    struct HashfunR;
    struct HashfunRS;
    struct HashfunNestedRS;
    struct HashfunFkRel;
    struct EqfunBuildFkRel;
    struct JoinpredRS;
    struct Joinpred_RS_T;
    struct JoinpredRTnested;
    struct ConcatfunChaining_RS;
    struct ConcatfunChaining_RS_T;
    struct ConcatfunNested_RS;
    struct ConcatfunNested_RST;
    struct Unnestfun_R_nS_xT;
    struct Unnestfun_R_xS_xT;
  public:
    inline Experiment4(const uint32_t aLog2CardR,
                       const uint32_t aAlpha, const uint32_t aMultAlpha,
                       const uint32_t aBeta, const uint32_t aMultBeta,
                       const std::filesystem::path aMeasureCsvFile, const std::vector<std::string>& aPlans)
      : _log2CardR(aLog2CardR), _alpha(aAlpha), _beta(aBeta), _multiplicityAlpha(aMultAlpha), _multiplicityBeta(aMultBeta), 
        _plans(plans_e::NONE), _measureCsv(aMeasureCsvFile) {
      plansFromVec(aPlans);
    }
  public:
    void init(const bool aShuffle = true);
    void run();
  public:
    inline size_t   cardR()             const { return (1U << _log2CardR); }
    inline size_t   cardS()             const { return cardFkRelations(); }
    inline size_t   cardT()             const { return cardFkRelations(); }
    inline uint32_t alpha()             const { return _alpha; }
    inline uint32_t beta()              const { return _beta; }
    inline uint32_t multiplicityAlpha() const { return _multiplicityAlpha; }
    inline uint32_t multiplicityBeta()  const { return _multiplicityBeta; }

    inline size_t numFkCommon()     const { return cardR() / (1U << alpha()); }
    inline size_t numFkExclusive()  const { return cardR() / (1U << beta()); }
    inline size_t cardFkCommon()    const { return numFkCommon() * multiplicityAlpha(); }
    inline size_t cardFkExclusive() const { return numFkExclusive() * multiplicityBeta(); }
    inline size_t cardFkRelations() const { return cardFkCommon() + cardFkExclusive(); }

    size_t calcJoinCard1() const;  // |join(R, S)| = |join(R, T)|
    size_t calcJoinCard2() const;  // |join(R, S, T)|

    inline size_t   cardRGenerated() const { return _R.card(); }
    inline size_t   cardSGenerated() const { return _S.card(); }
    inline size_t   cardTGenerated() const { return _T.card(); }

    inline std::chrono::milliseconds minRuntime() const { return _minRuntime; }
    inline size_t                    minRepeat()  const { return _minRepeat; }

  public:
    void printRelations(const bool aWideFmt = false, const size_t aIndent = 0) const;
    void printTimers() const;
    void printTypes() const;
    static void printParamTable();

  private:
    void generateR(const attrvalue_vt& aKeys);
    void generateFkRel(relation_t& aTargetRelation, const attrvalue_vt& aKeys,
        const attrvalue_vt& aForeignKeysCommon, const attrvalue_vt aForeignKeysExclusive);

    void plansFromVec(const std::vector<std::string>&);

    void runNdu();
    void runChj();

    void writeCsvHeader();       // header for measurement file
    void writeExpParamsToCsv();  // common experiment parameters to CSV

  private:
    // generation parameters
    uint32_t _log2CardR;    // base relation cardinalities in log2
    // Note: cardS := (cardR() / 2^alpha) * multiplicityAlpha + (cardR / 2^beta) * multiplicityBeta

    uint32_t _alpha;        // fraction 1/2^alpha of keys in R that are referenced by both S and T
    uint32_t _beta;         // fraction 1/2^beta  of keys in R that are referenced by either S or T

    uint32_t _multiplicityAlpha;   // how often is each 'alpha key' referenced by a FK?
    uint32_t _multiplicityBeta;    // how often is each 'beta key' referenced by a FK?

    // the base hash function used by all Hashfun classes
    inline static std::function<uint32_t(attrvalue_t)> _hashfunc = ht::murmur_hash<attrvalue_t>;

    // run parameters
    std::chrono::milliseconds _minRuntime{300};  // repeat running plan until this bound is reached
    size_t                    _minRepeat{8};     // min number of repeated measurements per plan

    plans_e  _plans = plans_e::ALL;  // plans to execute in run()

    const uint32_t _log2RsvChunkSize = 10;  // reservoir chunk size for all hashtables

    // base relations
    RelationRS<base_tuple_t> _R;  // key relation
    RelationRS<base_tuple_t> _S;  // foreign key S relation
    RelationRS<base_tuple_t> _T;  // foreign key T relation

    // bookkeeping & measurement
    df::infra::CSVWriter _measureCsv;
    std::map<std::string, std::pair<time_point_t, time_point_t>> _timePoints{};  // experiment timers
    bool _trace = true;

    static std::function<void(const result_tuple_t*, std::ostream&)> _result_tuple_printer;

    static const plans_map_t _supportedPlansMap;

  public:
    // return type deduction: http://mochan.info/c++/2019/06/12/returntype-deduction.html
    using hashvalue_t = decltype(_hashfunc(std::declval<attrvalue_t>()));

};


/* enum operators */
inline Experiment4::plans_e
operator|(const Experiment4::plans_e lhs, const Experiment4::plans_e rhs) {
  using type_t = std::underlying_type_t<Experiment4::plans_e>;
  return static_cast<Experiment4::plans_e>(static_cast<type_t>(lhs) | static_cast<type_t>(rhs));
}
inline Experiment4::plans_e&
operator|=(Experiment4::plans_e& lhs, const Experiment4::plans_e rhs) {
  lhs = lhs | rhs;
  return lhs;
}
inline Experiment4::plans_e
operator&(const Experiment4::plans_e lhs, const Experiment4::plans_e rhs) {
  using type_t = std::underlying_type_t<Experiment4::plans_e>;
  return static_cast<Experiment4::plans_e>(static_cast<type_t>(lhs) & static_cast<type_t>(rhs));
}
inline Experiment4::plans_e&
operator&=(Experiment4::plans_e& lhs, const Experiment4::plans_e rhs) {
  lhs = lhs & rhs;
  return lhs;
}

// map of all supported plans (poor man's reflection for enums)
const Experiment4::plans_map_t Experiment4::_supportedPlansMap = {
  {"none", plans_e::NONE},
  {"NONE", plans_e::NONE},
  {"Ndu",  plans_e::Ndu},
  {"Nnu",  plans_e::Nnu},
  {"Chj",  plans_e::Chj},
  {"all",  plans_e::ALL},
  {"ALL",  plans_e::ALL}
};


/* nested classes/structs */

struct Experiment4::result_tuple_t {
  // [r, s, t]
  const base_tuple_t* _r;
  const base_tuple_t* _s;
  const base_tuple_t* _t;
};
std::function<void(const Experiment4::result_tuple_t*, std::ostream&)>
Experiment4::_result_tuple_printer = [](const result_tuple_t* t, std::ostream& os) {
        os << "[" << *(t->_r) << "," << *(t->_s) << "," << *(t->_t) << "]";
};
struct Experiment4::result_tuple_RS_t {
  // C: [r, s]
  const base_tuple_t* _r;
  const base_tuple_t* _s;
};

struct Experiment4::GlobStat {};

struct Experiment4::HashfunR {
  using input_t = base_tuple_t;
  using output_t = hashvalue_t;
  inline static output_t eval(const input_t* aBaseTuple) {
    return _hashfunc(aBaseTuple->k);
  }
};
struct Experiment4::HashfunRS {
  using input_t = result_tuple_RS_t;
  using output_t = hashvalue_t;
  inline static output_t eval(const input_t* aBaseTuple) {
    return _hashfunc(aBaseTuple->_r->k);
  }
};
struct Experiment4::HashfunFkRel {
  using input_t = base_tuple_t;
  using output_t = hashvalue_t;
  inline static output_t eval(const input_t* aBaseTuple) {
    return _hashfunc(aBaseTuple->a);
  }
};
struct Experiment4::EqfunBuildFkRel {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->a == r->a);
  }
};
// for first probe
struct Experiment4::JoinpredRS {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->k == r->a);
  }
};
// result of first nested probe
struct Experiment4::nested_tuple_RS_t {
  // N: [r, {s0 | joinpred(r, s0)}]
  base_tuple_t* _r;
  const HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode* _s;
};
struct Experiment4::ConcatfunNested_RS {
  using left_t = base_tuple_t;
  using right_t = const HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode;
  using output_t = nested_tuple_RS_t;
  inline static output_t eval(left_t* l, const right_t* r) {
    return {l, r}; 
  }
};
// for second nested probe
struct Experiment4::JoinpredRTnested {
  using left_t = nested_tuple_RS_t;  // [r, list of s]
  using right_t = base_tuple_t;      // t
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->_r->k == r->a);
  }
};
// result of second nested probe
struct Experiment4::nested_tuple_RST_t {
  // N: [r, {s0 | joinpred(r, s0)}, {s1 | joinpred(r, s1)}]
  base_tuple_t* _r;
  const HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode* _s;
  const HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode* _t;
};
struct Experiment4::HashfunNestedRS {
  using input_t = nested_tuple_RS_t;
  using output_t = hashvalue_t;
  inline static output_t eval(const input_t* aNestedTuple) {
    return _hashfunc(aNestedTuple->_r->k);
  }
};
struct Experiment4::ConcatfunNested_RST {
  using left_t = nested_tuple_RS_t;
  using right_t = const HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode;
  using output_t = nested_tuple_RST_t;
  inline static output_t eval(left_t* left, const right_t* right) {
    return {left->_r, left->_s, right}; 
  }
};
// for unnest
struct Experiment4::tuple_R_nS_xT_t {
  // N: [r, {s | joinpred(r, s)}, t] 
  base_tuple_t* _r;
  const HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode* _s;
  const base_tuple_t* _t;
};
struct Experiment4::Unnestfun_R_nS_xT {
  using input_t  = nested_tuple_RST_t;
  using output_t = tuple_R_nS_xT_t;
  using MainNode = HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode;
  using data_t   = HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::data_t;
  inline static const MainNode* getMainNode(input_t* aNestedTuple) {
    return aNestedTuple->_t;
  }
  inline static void eval_left(output_t* out, input_t* in) {
    out->_r = in->_r;
    out->_s = in->_s;
  }
  inline static void eval_right(output_t* out, [[maybe_unused]] input_t* in, const data_t* data) {
    out->_t = data;
  }
};
struct Experiment4::Unnestfun_R_xS_xT {
  using input_t  = tuple_R_nS_xT_t;
  using output_t = result_tuple_t;
  using MainNode = HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::MainNode;
  using data_t   = HtNested1<base_tuple_t, HashfunFkRel, EqfunBuildFkRel>::data_t;
  inline static const MainNode* getMainNode(input_t* aNestedTuple) {
    return aNestedTuple->_s;
  }
  inline static void eval_left(output_t* out, input_t* in) {
    out->_r = in->_r;
    out->_t = in->_t;
  }
  inline static void eval_right(output_t* out, [[maybe_unused]] input_t* in, const data_t* data) {
    out->_s = data;
  }
};

// for chaining ht
struct Experiment4::Joinpred_RS_T {
  using left_t = result_tuple_RS_t;  // [r, s]
  using right_t = base_tuple_t;      // t
  inline static bool eval(const left_t* l, const right_t* r) {
    return (l->_r->k == r->a);
  }
};
struct Experiment4::ConcatfunChaining_RS {
  using left_t = base_tuple_t;
  using right_t = base_tuple_t;
  using output_t = result_tuple_RS_t;
  inline static output_t eval(left_t* l, const right_t* r) {
    return {l, r}; 
  }
};
struct Experiment4::ConcatfunChaining_RS_T {
  using left_t = result_tuple_RS_t;
  using right_t = base_tuple_t;
  using output_t = result_tuple_t;
  inline static output_t eval(left_t* left, const right_t* right) {
    return {left->_r, left->_s, right}; 
  }
};

// output operators
std::ostream& operator<<(std::ostream& os, const Experiment4::result_tuple_t& t) {
  os << "[" << *(t._r) << "," << *(t._s) << "," << *(t._t) << "]";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Experiment4::result_tuple_RS_t& t) {
  os << "[" << *(t._r) << "," << *(t._s) << "]";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Experiment4::nested_tuple_RS_t& t) {
  os << "[" << *(t._r) << "," << (t._s) << "]";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Experiment4::nested_tuple_RST_t& t) {
  os << "[" << *(t._r) << "," << (t._s) << (t._t) << "]";
  return os;
}
std::ostream& operator<<(std::ostream& os, const Experiment4::tuple_R_nS_xT_t& t) {
  os << "[" << *(t._r) << "," << (t._s) << "," << *(t._t) << "]";
  return os;
}

/* public member functions */

void
Experiment4::init(const bool aShuffle) {
  assert(cardR() >= numFkCommon() + 2 * numFkExclusive());
  std::mt19937 lRng;

  // unique keys: R.k, S.k, T.k
  attrvalue_vt lKeys;
  lKeys.resize(std::max(cardR(), cardFkRelations()));
  std::iota(lKeys.begin(), lKeys.end(), 0);

  // foreign keys S.a, T.b
  attrvalue_vt lFkCommon, lFkExclusiveS, lFkExclusiveT;
  lFkCommon.resize(cardFkCommon());
  lFkExclusiveS.resize(cardFkExclusive());
  lFkExclusiveT.resize(cardFkExclusive());

  attrvalue_t lVal = 0;
  // fill lFkCommon
  size_t index = 0;
  for (; lVal < numFkCommon(); ++lVal) {
    for (size_t i = 0; i < multiplicityAlpha(); ++i) {
      // std::cout << "common: [" << std::setw(2) << index << "] = " << std::setw(3) << lVal << std::endl;
      lFkCommon.at(index) = lVal;
      ++index;
    }
  }
  // fill lFkExclusiveS
  index = 0;
  for (; lVal < numFkCommon() + numFkExclusive(); ++lVal) {
    for (size_t i = 0; i < multiplicityBeta(); ++i) {
      // std::cout << "excl S: [" << std::setw(2) << index << "] = " << std::setw(3) << lVal << std::endl;
      lFkExclusiveS.at(index) = lVal;
      ++index;
    }
  }
  // fill lFkExclusiveT
  index = 0;
  for (; lVal < numFkCommon() + 2 * numFkExclusive(); ++lVal) {
    for (size_t i = 0; i < multiplicityBeta(); ++i) {
      // std::cout << "excl T: [" << std::setw(2) << index << "] = " << std::setw(3) << lVal << std::endl;
      lFkExclusiveT.at(index) = lVal;
      ++index;
    }
  }

  generateR(lKeys);

  // shuffle the FKs and insert into relations S and T
  if (aShuffle) {
    std::shuffle(lFkExclusiveS.begin(), lFkExclusiveS.end(), lRng);
    std::shuffle(lFkExclusiveT.begin(), lFkExclusiveT.end(), lRng);
    std::shuffle(lFkCommon.begin(), lFkCommon.end(), lRng);
  }
  generateFkRel(_S, lKeys, lFkCommon, lFkExclusiveS);
  if (aShuffle) {
    std::shuffle(lFkCommon.begin(), lFkCommon.end(), lRng);
  }
  generateFkRel(_T, lKeys, lFkCommon, lFkExclusiveT);
}

void
Experiment4::run() {
  writeCsvHeader();
  if ((plans_e::Ndu & _plans) != plans_e::NONE) runNdu();
  if ((plans_e::Chj & _plans) != plans_e::NONE) runChj();
}

size_t
Experiment4::calcJoinCard1() const {
  // cardinality after first join: join(R, S), join(R, T)
  // Let a := 1/2^alpha, b := 1/2^beta; mA, mB := multiplicity alpha, beta
  // |join(R,S)| |join(R,T) = |R| * a * mA + |R| * b * mB
  // = |S| = |T| (by design)
  return cardFkRelations();
}

size_t
Experiment4::calcJoinCard2() const {
  // cardinality after second join: join(R, S, T)
  return numFkCommon() * multiplicityAlpha() * multiplicityAlpha();
}

void
Experiment4::printRelations(const bool aWideFmt, const size_t aIndent) const {
  if (!aWideFmt) {
    std::cout << df::infra::indent(aIndent) << "-- R --\n";
    for (const auto& t : _R._tuples) {
      std::cout << df::infra::indent(aIndent) << t.k << "|" << t.a << "\n";
    }
    std::cout << df::infra::indent(aIndent) << "-- S --\n";
    for (const auto& t : _S._tuples) {
      std::cout << df::infra::indent(aIndent) << t.k << "|" << t.a << "\n";
    }
    std::cout << df::infra::indent(aIndent) << "-- T --\n";
    for (const auto& t : _T._tuples) {
      std::cout << df::infra::indent(aIndent) << t.k << "|" << t.a << "\n";
    }
  } else {
    size_t lMaxCard = std::max({cardR(), cardS(), cardT()});
    size_t lColW = std::max(6UL, df::infra::number_of_digits(lMaxCard)+2);
    std::cout
      << df::infra::indent(aIndent) << "|"
      << std::setw(2*lColW + 1) << std::left << "R" << "|"
      << std::setw(2*lColW + 1) << std::left << "S" << "|"
      << std::setw(2*lColW + 1) << std::left << "T" << "|\n";
    std::cout
      << df::infra::indent(aIndent) << "|"
      << std::setw(lColW) << std::left << "R.k" << "|"
      << std::setw(lColW) << std::left << "R.a" << "|"
      << std::setw(lColW) << std::left << "S.k" << "|"
      << std::setw(lColW) << std::left << "S.a" << "|"
      << std::setw(lColW) << std::left << "T.k" << "|"
      << std::setw(lColW) << std::left << "T.a" << "|\n";
    std::cout
      << df::infra::indent(aIndent) << "|"
      << std::setfill('-')
      << std::setw(lColW) << std::left << "" << "|"
      << std::setw(lColW) << std::left << "" << "|"
      << std::setw(lColW) << std::left << "" << "|"
      << std::setw(lColW) << std::left << "" << "|"
      << std::setw(lColW) << std::left << "" << "|"
      << std::setw(lColW) << std::left << "" << "|\n";
    std::cout << std::setfill(' ') << std::right;
    for (size_t i = 0; i < lMaxCard; ++i) {
      std::cout << df::infra::indent(aIndent) << "|";
      for (const auto& lRel : {_R, _S, _T}) {
        if (i >= lRel.card()) {
          std::cout
            << std::setw(lColW) << "" << "|"
            << std::setw(lColW) << "" << "|";
        } else {
          std::cout
            << std::setw(lColW) << std::to_string(lRel._tuples.at(i).k) << "|"
            << std::setw(lColW) << std::to_string(lRel._tuples.at(i).a) << "|";
        }
      }
      if (i == cardFkCommon()) { std::cout << " <-"; }
      std::cout << "\n";
    }
  }
}

void
Experiment4::printParamTable() {
  //lExp.cardR()
  //lExp.cardS()
  //lExp.cardT()
  //lExp.numFkCommon()
  //lExp.multiplicityAlpha() 
  //lExp.cardFkCommon()
  //lExp.numFkExclusive()
  //lExp.multiplicityBeta() <
  //lExp.cardFkExclusive()
  //lExp.calcJoinCard1()
  //lExp.calcJoinCard2()

  using std::setw;

  std::cout
    << setw(2) << "r" << " "
    << setw(2) << "a" << " "
    << setw(2) << "am" << " "
    << setw(2) << "b" << " "
    << setw(2) << "bm" << " "
    << setw(10) << "cardR" << " "
    << setw(10) << "cardS" << " "
    << setw(10) << "cardT" << " "
    << setw(10) << "FkC_dv" << " "
    << setw(10) << "FkC_mul" << " "
    << setw(10) << "FkC_card" << " "
    << setw(10) << "FkE_dv" << " "
    << setw(10) << "FkE_mul" << " "
    << setw(10) << "FkE_card" << " "
    << setw(10) << "c(RS)" << " "
    << setw(10) << "c(RST)" << " "
    << "\n";
  size_t lCtr = 0;
  for (uint32_t lLog2CardR = 10; lLog2CardR <= 25; ++lLog2CardR) {
    for (uint32_t lAlpha = 0; lAlpha <= lLog2CardR; ++lAlpha) {
      for (uint32_t lAlphaMult = 1; lAlphaMult < 10; ++lAlphaMult) {
        for (uint32_t lBeta = 0; lBeta <= lLog2CardR; ++lBeta) {
          for (uint32_t lBetaMult = 1; lBetaMult < 10; ++lBetaMult) {
            Experiment4 lExp(lLog2CardR, lAlpha, lAlphaMult, lBeta, lBetaMult, "/tmp/exp4_dummy.csv", {"NONE"});
            //if (lExp.cardS() > lExp.cardR() || lExp.cardT() > lExp.cardR()) { continue; }
            std::cout
              << setw(2) << lLog2CardR << " "
              << setw(2) << lAlpha << " "
              << setw(2) << lAlphaMult << " "
              << setw(2) << lBeta << " "
              << setw(2) << lBetaMult << " "
              << setw(10) << lExp.cardR() << " "
              << setw(10) << lExp.cardS() << " "
              << setw(10) << lExp.cardT() << " "
              << setw(10) << lExp.numFkCommon() << " "
              << setw(10) << lExp.multiplicityAlpha()  << " "
              << setw(10) << lExp.cardFkCommon() << " "
              << setw(10) << lExp.numFkExclusive() << " "
              << setw(10) << lExp.multiplicityBeta() << " "
              << setw(10) << lExp.cardFkExclusive() << " "
              << setw(10) << lExp.calcJoinCard1() << " "
              << setw(10) << lExp.calcJoinCard2() << " "
              << "\n";
            ++lCtr;
          }
        }
      }
    }
  }
  std::cout << "Counter: " << lCtr << "\n";
}

/* private member functions */

void
Experiment4::generateR(const attrvalue_vt& aKeys) {
  assert(cardR() <= aKeys.size());
  _R._tuples.resize(cardR());
  for (size_t i = 0; i < cardR(); ++i) {
    _R._tuples.at(i).k = aKeys.at(i);
    _R._tuples.at(i).a = 0;
  }
}
void
Experiment4::generateFkRel(relation_t& aTargetRelation,
                           const attrvalue_vt& aKeys, const attrvalue_vt& aForeignKeysCommon,
                           const attrvalue_vt aForeignKeysExclusive) {
  assert(cardFkRelations() <= aKeys.size());
  assert(cardFkRelations() <= aForeignKeysCommon.size() + aForeignKeysExclusive.size());
  aTargetRelation._tuples.resize(cardFkRelations());
  for (size_t i = 0; i < cardFkRelations(); ++i) {
    aTargetRelation._tuples.at(i).k = aKeys.at(i);
  }
  size_t i = 0;
  for (; i < aForeignKeysCommon.size(); ++i) {
    aTargetRelation._tuples.at(i).a = aForeignKeysCommon.at(i);
  }
  for (; i < aForeignKeysCommon.size() + aForeignKeysExclusive.size(); ++i) {
    aTargetRelation._tuples.at(i).a = aForeignKeysExclusive.at(i - aForeignKeysCommon.size());
  }
}

void
Experiment4::plansFromVec(const std::vector<std::string>& aPlanVec) {
  using df::infra::to_lower;
  _plans = plans_e::NONE;
  for (auto p : aPlanVec) {
    if (_supportedPlansMap.contains(p)) {
      _plans |= _supportedPlansMap.at(p);
    }
  }
}

void
Experiment4::writeCsvHeader() {
  _measureCsv
    .writeField("mintime")
    .writeField("minreps")
    // --
    .writeField("log2CardR")
    .writeField("a")      // alpha
    .writeField("aM")     // multiplicity for alpha FKs
    .writeField("b")      // beta
    .writeField("bM")     // multiplicity for beta FKs
    .writeField("cardR")
    .writeField("cardS")
    .writeField("cardT")
    // --
    .writeField("plan")
    .writeField("ht_impl")
    //.writeField("cc0_avg")  // hash table collision chain length avg/min/max
    //.writeField("cc0_min")  // cc0: all buckets
    //.writeField("cc0_max")
    //.writeField("cc1_avg")  // cc1: nonempty buckets
    //.writeField("cc1_min")
    //.writeField("cc1_max")
    .writeField("reps")     // number of repetitions of plan execution
    // --
    .writeField("t_total")
    .writeField("t_build_S")  // build on S (strand)
    .writeField("t_build_T")  // build on T (strand)
    .writeField("t_probe_R")  // probe with R
    // --
    .writeField("c_sc_R")
    .writeField("c_sc_S")
    .writeField("c_sc_T")
    .writeField("c_build_S")
    .writeField("c_build_T")
    .writeField("c_probe_RS")
    .writeField("c_probe_RS_cmp")
    .writeField("c_probe_RT")
    .writeField("c_probe_RT_cmp")
    .writeField("c_unnest_S")
    .writeField("c_unnest_T")
    .writeField("c_top")
    .newline();
}

void
Experiment4::writeExpParamsToCsv() {
  _measureCsv
    .writeField(df::infra::to_string(std::chrono::duration_cast<std::chrono::milliseconds>(minRuntime())))
    .writeField(minRepeat())
    .writeField(_log2CardR)
    .writeField(alpha())
    .writeField(multiplicityAlpha())
    .writeField(beta())
    .writeField(multiplicityBeta())
    .writeField(cardR())
    .writeField(cardS())
    .writeField(cardT());
}

// run functions

void
Experiment4::runNdu() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;

  // TYPEDEFS
  // Build Strands
  // -> both build strands are the same; name them differently anyway for better distinction/readaibility
  using build_S_t = AlgNestJoinBuild<HashfunFkRel, EqfunBuildFkRel, GlobStat>;
  using build_T_t = AlgNestJoinBuild<HashfunFkRel, EqfunBuildFkRel, GlobStat>;
  using scan_S_t = AlgScan<build_S_t>;
  using scan_T_t = AlgScan<build_T_t>;

  // Probe Strand
  // Unnest #2 is the "upper"/second unnest and unpacks S
  // Unnest #1 is the "lower"/first unnest and unpacks T
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using unnest_2_t = AlgUnnestHt<top_t, Unnestfun_R_xS_xT, build_S_t::hashtable_t>;
  using unnest_1_t = AlgUnnestHt<unnest_2_t, Unnestfun_R_nS_xT, build_T_t::hashtable_t>;
  using probe_RT_t = AlgNestJoinProbe<unnest_1_t, build_T_t, HashfunNestedRS, JoinpredRTnested, ConcatfunNested_RST>;
  using probe_RS_t = AlgNestJoinProbe<probe_RT_t, build_S_t, HashfunR, JoinpredRS, ConcatfunNested_RS>;
  using scan_R_t = AlgScan<probe_RS_t>;

  // PLAN
  GlobStat lGs;
  size_t lNumDvFk = numFkCommon() + numFkExclusive();
  build_S_t lOpBuildS(lNumDvFk, _log2RsvChunkSize, _log2RsvChunkSize);
  scan_S_t  lOpScanS(&lOpBuildS, &_S);
  build_T_t lOpBuildT(lNumDvFk, _log2RsvChunkSize, _log2RsvChunkSize);
  scan_T_t  lOpScanT(&lOpBuildT, &_T);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  unnest_2_t lOpUnnest2(&lOpTop);
  unnest_1_t lOpUnnest1(&lOpUnnest2);
  probe_RT_t lOpProbeRT(&lOpUnnest1, &lOpBuildT);
  //probe_RT_t lOpProbeRT(&lOpTop, &lOpBuildT);
  probe_RS_t lOpProbeRS(&lOpProbeRT, &lOpBuildS);
  scan_R_t   lOpScanR(&lOpProbeRS, &_R);

  // RUN
  std::chrono::nanoseconds lDurationBuildS{0};
  std::chrono::nanoseconds lDurationBuildT{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(),
      [&lOpTop, &lDurationTotal, &lDurationBuildS, &lDurationBuildT, &lDurationProbe, &lOpScanR, &lOpScanS, &lOpScanT, &lGs] () {
        time_point_t lTpStart{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTpMid1{clock_t::now()};
        lOpScanT.run(&lGs);
        time_point_t lTpMid2{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTpStop{clock_t::now()};

        lDurationBuildS += (lTpMid1 - lTpStart);
        lDurationBuildT += (lTpMid2 - lTpMid1);
        lDurationProbe += (lTpStop - lTpMid2);
        lDurationTotal += (lTpStop - lTpStart);

        lOpTop.printResult(false);  // disable printing after the first run
      },
      [&lOpBuildS, &lOpBuildT] () {
        lOpBuildS.clear_ht();
        lOpBuildT.clear_ht();
      },
      false, minRepeat()
    );

  lDurationBuildS /= lIterations;
  lDurationBuildT /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan Ndu\n";
    std::cout << "  S: sizeof(MainNode): " << sizeof(build_S_t::hashtable_t::MainNode) << "\n";
    std::cout << "  S: sizeof(SubNode):  " << sizeof(build_S_t::hashtable_t::SubNode) << "\n";
    std::cout << "  T: sizeof(MainNode): " << sizeof(build_T_t::hashtable_t::MainNode) << "\n";
    std::cout << "  T: sizeof(SubNode):  " << sizeof(build_T_t::hashtable_t::SubNode) << "\n";
  }

  const HtStatistics lHtStatS = lOpBuildS.hashtable().makeStatistics();
  const HtStatistics lHtStatT = lOpBuildT.hashtable().makeStatistics();

  // write result to csv
  writeExpParamsToCsv();
  _measureCsv
    .writeField("Ndu")                    // plan
    .writeField("nested")                 // ht_impl
    .writeField(lIterations)              // reps
    //
    .writeField(lDurationTotal.count())   // t_total
    .writeField(lDurationBuildS.count())  // t_build_S
    .writeField(lDurationBuildT.count())  // t_build_T
    .writeField(lDurationProbe.count())   // t_probe_R
    //
    .writeField(lOpScanR.count())         // c_sc_R
    .writeField(lOpScanS.count())         // c_sc_S
    .writeField(lOpScanT.count())         // c_sc_T
    .writeField(lOpBuildS.count())        // c_build_S
    .writeField(lOpBuildT.count())        // c_build_T
    .writeField(lOpProbeRS.count())       // c_probe_R
    .writeField(lOpProbeRS.numCmps())     // c_probe_R
    .writeField(lOpProbeRT.count())       // c_probe_R
    .writeField(lOpProbeRT.numCmps())     // c_probe_R
    .writeField(lOpUnnest1.count())       // c_unnest_S
    .writeField(lOpUnnest2.count())       // c_unnest_T
    .writeField(lOpTop.count())           // c_top
    .newline();
}

void
Experiment4::runChj() {
  std::cout << __PRETTY_FUNCTION__ << std::endl;

  // TYPEDEFS
  // Build Strands
  // -> both build strands are the same; name them differently anyway for better distinction/readaibility
  using build_S_t = AlgHashJoinBuild<HashfunFkRel, EqfunBuildFkRel, GlobStat>;
  using build_T_t = AlgHashJoinBuild<HashfunFkRel, EqfunBuildFkRel, GlobStat>;
  using scan_S_t = AlgScan<build_S_t>;
  using scan_T_t = AlgScan<build_T_t>;

  // Probe Strand
  // Unnest #2 is the "upper"/second unnest and unpacks S
  // Unnest #1 is the "lower"/first unnest and unpacks T
  using top_t = AlgTop<result_tuple_t, GlobStat>;
  using probe_RT_t = AlgHashJoinProbe<top_t, build_T_t, HashfunRS, Joinpred_RS_T, ConcatfunChaining_RS_T>;
  using probe_RS_t = AlgHashJoinProbe<probe_RT_t, build_S_t, HashfunR, JoinpredRS, ConcatfunChaining_RS>;
  using scan_R_t = AlgScan<probe_RS_t>;

  // PLAN
  GlobStat lGs;
  size_t lNumDvFk = numFkCommon() + numFkExclusive();
  build_S_t lOpBuildS(lNumDvFk, _log2RsvChunkSize);
  scan_S_t  lOpScanS(&lOpBuildS, &_S);
  build_T_t lOpBuildT(lNumDvFk, _log2RsvChunkSize);
  scan_T_t  lOpScanT(&lOpBuildT, &_T);

  top_t lOpTop(std::cout, false, _result_tuple_printer);
  probe_RT_t lOpProbeRT(&lOpTop, &lOpBuildT);
  //probe_RT_t lOpProbeRT(&lOpTop, &lOpBuildT);
  probe_RS_t lOpProbeRS(&lOpProbeRT, &lOpBuildS);
  scan_R_t   lOpScanR(&lOpProbeRS, &_R);

  // RUN
  std::chrono::nanoseconds lDurationBuildS{0};
  std::chrono::nanoseconds lDurationBuildT{0};
  std::chrono::nanoseconds lDurationProbe{0};
  std::chrono::nanoseconds lDurationTotal{0};

  auto [_, lIterations] = df::infra::repeat_mintime(
      minRuntime(),
      [&lOpTop, &lDurationTotal, &lDurationBuildS, &lDurationBuildT, &lDurationProbe, &lOpScanR, &lOpScanS, &lOpScanT, &lGs] () {
        time_point_t lTpStart{clock_t::now()};
        lOpScanS.run(&lGs);
        time_point_t lTpMid1{clock_t::now()};
        lOpScanT.run(&lGs);
        time_point_t lTpMid2{clock_t::now()};
        lOpScanR.run(&lGs);
        time_point_t lTpStop{clock_t::now()};

        lDurationBuildS += (lTpMid1 - lTpStart);
        lDurationBuildT += (lTpMid2 - lTpMid1);
        lDurationProbe += (lTpStop - lTpMid2);
        lDurationTotal += (lTpStop - lTpStart);

        lOpTop.printResult(false);  // disable printing after the first run
      },
      [&lOpBuildS, &lOpBuildT] () {
        lOpBuildS.clear_ht();
        lOpBuildT.clear_ht();
      },
      false, minRepeat()
    );

  lDurationBuildS /= lIterations;
  lDurationBuildT /= lIterations;
  lDurationProbe /= lIterations;
  lDurationTotal /= lIterations;

  if (_trace) {
    std::cout << "Plan Chj\n";
    std::cout << "  S: sizeof(Node): " << sizeof(build_S_t::hashtable_t::Node) << "\n";
    std::cout << "  T: sizeof(Node): " << sizeof(build_T_t::hashtable_t::Node) << "\n";
  }

  writeExpParamsToCsv();
  _measureCsv
    .writeField("Chj")                    // plan
    .writeField("chaining")               // ht_impl
    .writeField(lIterations)              // reps
    //
    .writeField(lDurationTotal.count())   // t_total
    .writeField(lDurationBuildS.count())  // t_build_S
    .writeField(lDurationBuildT.count())  // t_build_T
    .writeField(lDurationProbe.count())   // t_probe_R
    //
    .writeField(lOpScanR.count())         // c_sc_R
    .writeField(lOpScanS.count())         // c_sc_S
    .writeField(lOpScanT.count())         // c_sc_T
    .writeField(lOpBuildS.count())        // c_build_S
    .writeField(lOpBuildT.count())        // c_build_T
    .writeField(lOpProbeRS.count())       // c_probe_R
    .writeField(lOpProbeRS.numCmps())     // c_probe_R
    .writeField(lOpProbeRT.count())       // c_probe_R
    .writeField(lOpProbeRT.numCmps())     // c_probe_R
    .writeField("NA")                     // c_unnest_S
    .writeField("NA")                     // c_unnest_T
    .writeField(lOpTop.count())           // c_top
    .newline();
}


// ############################################################################

int main(int argc, char** argv) {
  CLI::App lApp{"Hash Table Experiment 4: Two Joins w/ Deferred Unnesting"};

  bool lRunExp = true;
  bool lPrintTimers = true;
  bool lPrintRelations = false;
  bool lPrintParamTable = false;

  uint32_t lLog2CardR;
  uint32_t lAlpha;
  uint32_t lBeta;
  uint32_t lMultiplicityAlpha;
  uint32_t lMultiplicityBeta;
  std::filesystem::path lMeasureFile;
  std::vector<std::string> lPlansToRun = {"all"};

  // specify CLI arguments
  lApp.add_flag("--run,!--no-run", lRunExp, "Run the experiment");
  lApp.add_flag("--print-timers,!--no-print-timers", lPrintTimers, "Print timers after experiment run");
  lApp.add_flag("--print-relations,!--no-print-relations", lPrintRelations, "Print generated relations R and S");
  lApp.add_flag("--print-paramtable,!--no-print-paramtable", lPrintParamTable, "Print parameter table and values");

  lApp.add_option("-R,--card-R", lLog2CardR, "Cardinality of key relation R as log2")
    ->required()
    ->check(CLI::Range(0, 30));
  lApp.add_option("-a,--alpha", lAlpha, "Fraction of distinct foreign keys that survive both joins, as log2")
    ->required()
    ->check(CLI::Range(0, 30));
  lApp.add_option("-b,--beta", lBeta, "Fraction of distinct foreign keys that survive one join but not the other, as log2")
    ->required()
    ->check(CLI::Range(0, 30));
  lApp.add_option("-A,--alpha-mult", lMultiplicityAlpha, "Multiplicity of foreign keys from alpha")
    ->required()
    ->check(CLI::PositiveNumber);
  lApp.add_option("-B,--beta-mult", lMultiplicityBeta, "Multiplicity of foreign keys from beta")
    ->required()
    ->check(CLI::PositiveNumber);
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

  // parse args
  CLI11_PARSE(lApp, argc, argv);

  // print what we've read
  std::cout << "Running Experiment 4 with the following config:\n"
    << df::infra::indent(1) << "--card-R " << lLog2CardR << "\n"
    << df::infra::indent(1) << "--alpha " << lAlpha << "\n"
    << df::infra::indent(1) << "--beta " << lBeta << "\n"
    << df::infra::indent(1) << "--alpha-mult " << lMultiplicityAlpha << "\n"
    << df::infra::indent(1) << "--beta-mult " << lMultiplicityBeta << "\n"
    << df::infra::indent(1) << "--measure-file " << lMeasureFile << "\n"
    << df::infra::indent(1) << "--plans ";
  std::for_each(lPlansToRun.cbegin(), lPlansToRun.cend(), [](const std::string& x) { std::cout << x << ","; });
  std::cout << "\n";
  std::cout << "\n";

  // create experiment instance and generate data
  Experiment4 lExp(lLog2CardR, lAlpha, lMultiplicityAlpha, lBeta, lMultiplicityBeta, lMeasureFile, lPlansToRun);
  std::cout << "Data Generation Config\n"
    << df::infra::indent(1) << "Base Relations" << "\n"
    << df::infra::indent(1) << std::setw(14) << "|R|"            << ": " << std::setw(10) << lExp.cardR() << "\n"
    << df::infra::indent(1) << std::setw(14) << "|S|"            << ": " << std::setw(10) << lExp.cardS() << "\n"
    << df::infra::indent(1) << std::setw(14) << "|T|"            << ": " << std::setw(10) << lExp.cardT() << "\n"
    << df::infra::indent(1) << "FK_Common" << "\n"
    << df::infra::indent(1) << std::setw(14) << "dv"             << ": " << std::setw(10) << lExp.numFkCommon() << "\n"
    << df::infra::indent(1) << std::setw(14) << "mult"           << ": " << std::setw(10) << lExp.multiplicityAlpha() << "\n"
    << df::infra::indent(1) << std::setw(14) << "card"           << ": " << std::setw(10) << lExp.cardFkCommon() << "\n"
    << df::infra::indent(1) << "FK_Excl" << "\n"
    << df::infra::indent(1) << std::setw(14) << "dv"             << ": " << std::setw(10) << lExp.numFkExclusive() << "\n"
    << df::infra::indent(1) << std::setw(14) << "mult"           << ": " << std::setw(10) << lExp.multiplicityBeta() << "\n"
    << df::infra::indent(1) << std::setw(14) << "card"           << ": " << std::setw(10) << lExp.cardFkExclusive() << "\n"
    << df::infra::indent(1) << "Joins" << "\n"
    << df::infra::indent(1) << std::setw(14) << "card(j(R,Fk))"  << ": " << std::setw(10) << lExp.calcJoinCard1() << "\n"
    << df::infra::indent(1) << std::setw(14) << "card(j(R,S,T))" << ": " << std::setw(10) << lExp.calcJoinCard2() << "\n"
    << "\n";

  lExp.init();
  
  std::cout << "Data Generation Result\n"
    << df::infra::indent(1) << std::setw(14) << "|R|" << ": " << std::setw(10) << lExp.cardRGenerated() << "\n"
    << df::infra::indent(1) << std::setw(14) << "|S|" << ": " << std::setw(10) << lExp.cardSGenerated() << "\n"
    << df::infra::indent(1) << std::setw(14) << "|T|" << ": " << std::setw(10) << lExp.cardTGenerated() << "\n\n";

  if (lPrintRelations) {
    std::cout << "# Relations\n";
    lExp.printRelations(true, 1);
    std::cout << "\n";
  }

  if (lRunExp) {
    lExp.run();
  }

  if (lPrintParamTable) {
    std::cout << "# Param Table\n";
    Experiment4::printParamTable();
  }


  return EXIT_SUCCESS;
}

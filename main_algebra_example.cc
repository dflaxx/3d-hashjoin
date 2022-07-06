#include "util/standard_includes.hh"

#include "algebra.hh"
#include "ht_chaining.hh"
#include "ht_nested.hh"
#include "util/hasht.hh"


/*
 * for algebra from algebra.hh
 */

using attrval_t = int;

struct tuple_L_t { attrval_t a, b; };
struct tuple_R_t { attrval_t c, d; };

std::ostream& operator<< (std::ostream& os, const tuple_L_t t) {
  os << "(" << t.a << "," << t.b << ")"; return os;
}
std::ostream& operator<< (std::ostream& os, const tuple_R_t t) {
  os << "(" << t.c << "," << t.d << ")"; return os;
}
std::ostream& operator<< (std::ostream& os, const tuple_L_t* t) {
  os << "(" << t->a << "," << t->b << ")"; return os;
}
std::ostream& operator<< (std::ostream& os, const tuple_R_t* t) {
  os << "(" << t->c << "," << t->d << ")"; return os;
}

class SelectionL {
  public:
    using input_t = tuple_L_t;
  public:
    inline static bool eval(const input_t* aInput) {
      return (aInput->b < 40);
    }
};
class DynSelectionL {
  public:
    using input_t = tuple_L_t;
  public:
    inline bool operator()(const input_t* aInput) {
      return (aInput->b < 40);
    }
};

class HashfunBuild {
  public:
    using input_t = tuple_R_t;
    using output_t = uint64_t;
  public:
    inline static output_t eval(const input_t* aInput) {
      return ht::murmur_hash<uint64_t>(aInput->c);
    }
};

class HashfunProbe {
  public:
    using input_t = tuple_L_t;
    using output_t = uint64_t;
  public:
    inline static output_t eval(const input_t* aInput) {
      return ht::murmur_hash<uint64_t>(aInput->a);
    }
};

class EqFunBuild {
  public:
    using left_t = tuple_R_t;
    using right_t = tuple_R_t;
  public:
    inline static bool eval(const left_t* l, const right_t* r) {
      return (l->c == r->c);
    }
};

// the join predicate: L.a = R.c
class EqFunProbe {
  public:
    using left_t = tuple_L_t;
    using right_t = tuple_R_t;
  public:
    inline static bool eval(const left_t* l, const right_t* r) {
      return (l->a == r->c);
    }
};

struct tuple_nested_t {
  tuple_L_t* _left;
  const HtNested1<tuple_R_t, HashfunBuild, EqFunBuild>::MainNode* _right;
};

struct tuple_LR_t {
  const tuple_L_t* _left;
  const tuple_R_t* _right;
};

class ConcatFunNested {
  public:
    using left_t = tuple_L_t;
    using right_t = HtNested1<tuple_R_t, HashfunBuild, EqFunBuild>::MainNode;
    using output_t = tuple_nested_t;
  public:
    inline static output_t eval(left_t* l, const right_t* r) {
      return tuple_nested_t{l, r}; 
    }
};

class ConcatFunChaining {
  public:
    using left_t = tuple_L_t;
    using right_t = tuple_R_t;
    using output_t = tuple_LR_t;
  public:
    inline static output_t eval(left_t* l, const right_t* r) {
      return {l, r}; 
    }
};

// print output tuple
std::ostream& operator<< (std::ostream& os, const tuple_LR_t t) {
  os << "(" << t._left->a << "," << t._left->b << "," << t._right->c << "," << t._right->d << ")"; return os;
}

class UnnestFun {
  public:
    using input_t  = tuple_nested_t;
    using output_t = tuple_LR_t;
    using MainNode = HtNested1<tuple_R_t, HashfunBuild, EqFunBuild>::MainNode;
    using data_t   = HtNested1<tuple_R_t, HashfunBuild, EqFunBuild>::data_t;
  public:
    inline static const MainNode* getMainNode(input_t* aNestedTuple) {
      return aNestedTuple->_right;
    }
    inline static void eval_left(output_t* out, input_t* in) {
      out->_left = in->_left;  // just copy the pointer to tuple_L_t
    }
    inline static void eval_right(output_t* out, input_t* in, const data_t* data) {
      out->_right = data;
    }
};

/* some query evaluation plans */

// scan --> selection --> top
void algebra_test0() {
  std::cout << "### " << __PRETTY_FUNCTION__ << " ###" << std::endl;

  GlobStat0 lGs0;

  RelationRS<tuple_L_t> lRelL = { ._tuples{
  // a  b
    {1, 11},
    {2, 21},
    {3, 31},
    {4, 41}}
  };

  std::cout
    << "-- Relation L --\n"
    << lRelL
    << "\n";

  using top_t = AlgTop<tuple_L_t, GlobStat0>;
  //using sel_t = AlgSelection<top_t, SelectionL>;
  using sel_t = AlgDynSelection<top_t, DynSelectionL>;
  using scan_t = AlgScan<sel_t>;

  top_t lTop(std::cout, true,
      [](const auto* t, std::ostream& os) {
        os << "(" << t->a << "," << t->b << ") @ " << t;
      });
  //sel_t lSelL(&lTop);
  sel_t lSelL(&lTop, DynSelectionL());
  scan_t lScanL(&lSelL, &lRelL);

  std::cout << "Output tuples\n";
  lScanL.run(&lGs0);

  std::cout << std::endl;
  std::cout << "count Top:  " << lTop.count() << "\n";
  std::cout << "count Sel:  " << lSelL.count() << "\n";
  std::cout << "count Scan: " << lScanL.count() << "\n";
}

// scan, select, build+probe, top (output nested tuples)
void algebra_test1() {
  std::cout << "### " << __PRETTY_FUNCTION__ << " ###" << std::endl;

  /*
   *      Top
   *       |
   *     Probe.....(L.a = R.c).....Build
   *       |                         |
   *    Select                       |
   *       |                         |
   *     Scan                      Scan
   *       |                         |
   *       L                         R
   */

  GlobStat0 lGs0;

  RelationRS<tuple_L_t> lRelL = { ._tuples{
  // a  b
    {1, 11},
    {2, 21},
    {3, 31},
    {4, 41}}
  };
  RelationRS<tuple_R_t> lRelR = { ._tuples{
  // c  d
    {1, -1},
    {1, -2},
    {1, -3},
    {2, -1},
    {2, -2},
    {3, -1}}
  };

  std::cout << "-- Relation L --\n" << lRelL << "\n";
  std::cout << "-- Relation R --\n" << lRelR << "\n";

  //template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, typename Tglobstat>
  using build_t  = AlgNestJoinBuild<HashfunBuild, EqFunBuild, GlobStat0>;
  using scan_R_t = AlgScan<build_t>;

  using top_t = AlgTop<tuple_nested_t, GlobStat0>;
  // template <alg_consumer_c Tconsumer, alg_buildop_c Tbuild,
  //           alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred, alg_concatfun_c Tconcatfun>
  using probe_t = AlgNestJoinProbe<top_t, build_t, HashfunProbe, EqFunProbe, ConcatFunNested>;
  using sel_t = AlgSelection<probe_t, SelectionL>;
  using scan_L_t = AlgScan<sel_t>;

  build_t lBuildOperator(5, 4, 4);
  scan_R_t lScanR(&lBuildOperator, &lRelR);

  top_t lTop(std::cout, true,
      [](const auto* t, std::ostream& os) {
        os << "(" << t->_left->a << "," << t->_left->b << "," << t->_right->data()->c << "," << t->_right->data()->d << ") @ " << t;
      }
      );
  probe_t lProbeOperator(&lTop, &lBuildOperator);
  sel_t lSelL(&lProbeOperator);
  scan_L_t lScanL(&lSelL, &lRelL);

  lScanR.run(&lGs0);

  std::cout << "Output tuples\n";
  lScanL.run(&lGs0);

  std::cout << std::endl;
  std::cout << "Build Strand:\n";
  std::cout << "  count Build: " << lBuildOperator.count() << "\n";
  std::cout << "  count Scan:  " << lScanR.count() << "\n";
  std::cout << "Probe Strand:\n";
  std::cout << "  count Top:   " << lTop.count() << "\n";
  std::cout << "  count Probe: " << lProbeOperator.count() << "\n";
  std::cout << "  count Sel:   " << lSelL.count() << "\n";
  std::cout << "  count Scan:  " << lScanL.count() << "\n";
}

// nested join of two relations, subsequent unnest + output
void algebra_test2() {
  std::cout << "### " << __PRETTY_FUNCTION__ << " ###" << std::endl;

  /*
   *      Top
   *       |
   *    Unnest
   *       |
   *    NProbe.....(L.a = R.c)....NBuild  (probe & build w/ nested ht)
   *       |                         |
   *    Select                       |
   *       |                         |
   *     Scan                      Scan
   *       |                         |
   *       L                         R
   */

  GlobStat0 lGs0;

  RelationRS<tuple_L_t> lRelL = { ._tuples{
  // a  b
    {1, 11},
    {2, 21},
    {3, 31},
    {4, 41}}
  };
  RelationRS<tuple_R_t> lRelR = { ._tuples{
  // c  d
    {1, -1},
    {1, -2},
    {1, -3},
    {2, -1},
    {2, -2},
    {3, -1}}
  };

  std::cout << "-- Relation L --\n" << lRelL << "\n";
  std::cout << "-- Relation R --\n" << lRelR << "\n";

  //template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, typename Tglobstat>
  using build_t  = AlgNestJoinBuild<HashfunBuild, EqFunBuild, GlobStat0>;
  using scan_R_t = AlgScan<build_t>;

  using top_t = AlgTop<tuple_LR_t, GlobStat0>;
  // AlgUnnestHt: template <alg_consumer_c Tconsumer, alg_unnestfun_c Tunnestfun, typename Thtnested>
  using unnest_t = AlgUnnestHt<top_t, UnnestFun, typename build_t::hashtable_t>;
  // AlgNestJoinProbe:
  // template <alg_consumer_c Tconsumer, alg_buildop_c Tbuild,
  //           alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred, alg_concatfun_c Tconcatfun>
  using probe_t = AlgNestJoinProbe<unnest_t, build_t, HashfunProbe, EqFunProbe, ConcatFunNested>;
  using sel_t = AlgSelection<probe_t, SelectionL>;
  using scan_L_t = AlgScan<sel_t>;

  build_t lBuildOperator(5, 4, 4);
  scan_R_t lScanR(&lBuildOperator, &lRelR);

  top_t lTop(std::cout, true,
      [](const auto* t, std::ostream& os) {
        os << "(" << t->_left->a << "," << t->_left->b << "," << t->_right->c << "," << t->_right->d << ") @ " << t;
      }
      );
  unnest_t lUnnestOp(&lTop);
  probe_t lProbeOperator(&lUnnestOp, &lBuildOperator);
  sel_t lSelL(&lProbeOperator);
  scan_L_t lScanL(&lSelL, &lRelL);

  lScanR.run(&lGs0);

  std::cout << "Output tuples\n";
  lScanL.run(&lGs0);

  std::cout << std::endl;
  std::cout
    << "Build Strand:\n"
    << "  Build:  " << std::setw(6) << lBuildOperator.count() << "  " << std::setw(10) << lBuildOperator.getRuntimeStr() << "\n"
    << "  Scan:   " << std::setw(6) << lScanR.count()         << "  " << std::setw(10) << lScanR.getRuntimeStr() << "\n"
    << "Probe Strand:\n"
    << "  Top:    " << std::setw(6) << lTop.count()           << "  " << std::setw(10) << lTop.getRuntimeStr() << "\n"
    << "  Unnest: " << std::setw(6) << lUnnestOp.count()      << "  " << std::setw(10) << lUnnestOp.getRuntimeStr() << "\n"
    << "  Probe:  " << std::setw(6) << lProbeOperator.count() << "  " << std::setw(10) << lProbeOperator.getRuntimeStr() << "\n"
    << "  Sel:    " << std::setw(6) << lSelL.count()          << "  " << std::setw(10) << lSelL.getRuntimeStr() << "\n"
    << "  Scan:   " << std::setw(6) << lScanL.count()         << "  " << std::setw(10) << lScanL.getRuntimeStr() << "\n";
}

/*
 * Conventional hash join (not nested)
 */
void algebra_test3() {
  std::cout << "### " << __PRETTY_FUNCTION__ << " ###" << std::endl;

  /*
   *      Top
   *       |
   *     Probe.....(L.a = R.c).....Build
   *       |                         |
   *    Select                       |
   *       |                         |
   *     Scan                      Scan
   *       |                         |
   *       L                         R
   */

  GlobStat0 lGs0;

  RelationRS<tuple_L_t> lRelL = { ._tuples{
  // a  b
    {1, 11},
    {2, 21},
    {3, 31},
    {4, 41}}
  };
  RelationRS<tuple_R_t> lRelR = { ._tuples{
  // c  d
    {1, -1},
    {1, -2},
    {1, -3},
    {2, -1},
    {2, -2},
    {3, -1}}
  };

  std::cout << "-- Relation L --\n" << lRelL << "\n";
  std::cout << "-- Relation R --\n" << lRelR << "\n";

  //template <alg_hashfun_c Thashfun, alg_binary_predicate_c Tequalfun, typename Tglobstat>
  using build_t  = AlgHashJoinBuild<HashfunBuild, EqFunBuild, GlobStat0>;
  using scan_R_t = AlgScan<build_t>;

  using top_t = AlgTop<tuple_LR_t, GlobStat0>;
  // AlgHashJoinProbe:
  // template <alg_consumer_c Tconsumer, alg_buildop_c Tbuild,
  //           alg_hashfun_c Thashfun, alg_binary_predicate_c Tjoinpred,
  //           alg_concatfun_c Tconcatfun>
  using probe_t = AlgHashJoinProbe<top_t, build_t, HashfunProbe, EqFunProbe, ConcatFunChaining>;
  using sel_t = AlgSelection<probe_t, SelectionL>;
  using scan_L_t = AlgScan<sel_t>;

  build_t lOpBuild(5, 4);
  scan_R_t lOpScanR(&lOpBuild, &lRelR);

  top_t lOpTop(std::cout, true,
      [](const auto* t, std::ostream& os) {
        os << "(" << t->_left->a << "," << t->_left->b << "," << t->_right->c << "," << t->_right->d << ") @ " << t;
      }
      );
  probe_t lOpProbe(&lOpTop, &lOpBuild);
  sel_t lOpSelL(&lOpProbe);
  scan_L_t lOpScanL(&lOpSelL, &lRelL);

  lOpScanR.run(&lGs0);

  std::cout << "Output tuples\n";
  lOpScanL.run(&lGs0);

  std::cout << std::endl;
  std::cout
    << "Build Strand:\n"
    << "  Build:  " << std::setw(6) << lOpBuild.count() << "  " << std::setw(10) << lOpBuild.getRuntimeStr() << "\n"
    << "  Scan:   " << std::setw(6) << lOpScanR.count() << "  " << std::setw(10) << lOpScanR.getRuntimeStr() << "\n"
    << "Probe Strand:\n"
    << "  Top:    " << std::setw(6) << lOpTop.count()   << "  " << std::setw(10) << lOpTop.getRuntimeStr() << "\n"
    << "  Probe:  " << std::setw(6) << lOpProbe.count() << "  " << std::setw(10) << lOpProbe.getRuntimeStr() << "\n"
    << "  Sel:    " << std::setw(6) << lOpSelL.count()  << "  " << std::setw(10) << lOpSelL.getRuntimeStr() << "\n"
    << "  Scan:   " << std::setw(6) << lOpScanL.count() << "  " << std::setw(10) << lOpScanL.getRuntimeStr() << "\n";

  std::cout << std::endl;
  std::cout << "Build Strand:\n";
  print_strand(&lOpScanR, 1);
  std::cout << "Probe Strand:\n";
  print_strand(&lOpScanL, 1);
}


int
main() {
  std::cout << "__cplusplus == " << __cplusplus << std::endl;

  std::cout << "- - - - - - - - - - - - - - - - - - -" << std::endl;

  algebra_test0();
  std::cout << "- - - - - - - - - - - - - - - - - - -" << std::endl;
  algebra_test1();
  std::cout << "- - - - - - - - - - - - - - - - - - -" << std::endl;
  algebra_test2();
  std::cout << "- - - - - - - - - - - - - - - - - - -" << std::endl;
  algebra_test3();


  return 0;
}

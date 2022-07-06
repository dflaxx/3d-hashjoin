#ifndef INFRA_GEN_RANDOM_UINT_VECTOR_HH
#define INFRA_GEN_RANDOM_UINT_VECTOR_HH

#include "standard_includes.hh"
#include <random>
#include "zipf_distribution.hh"

/*
 * class GenRandIntVec
 * generates integer vector of length aSize
 * with values in [0,max[
 * according to some distribution dist_t;
 *
 * Author: Guido Moerkotte (moerkotte@uni-mannheim.de)
 * Modifications by Daniel Flachs (flachs@uni-mannheim.de)
 * 
 */

typedef std::mt19937 rng32_t;

class GenRandIntVec {
  public:
    // supported distributions
    enum dist_t  {
      kKey    = 0, // key column
      kDiv    = 1,
      kUni    = 2,
      kExp    = 3,
      kNorm   = 4,
      kZipf   = 5,
      kSelf   = 6,
      kPois   = 7,
      kNoDist = 8
    };
    enum flag_t {
      kFill    = 1, // all values will occur at least once
      kShuffle = 2, // shuffle frequency vector after generating it
    };
    struct param_t {
      dist_t _dist;  // the distribution
      uint   _max;   // maximum value, or for kDiv: i/_max
      uint   _shift; // shift, not for kUni, kDiv, kUni
      double _param; // lambda for kExp, sigma for kNorm, theta for kZipf, h for kSelf, etc
      uint   _flags; // the flags
      int    _order; // -1: permute, 0: leave as is, +1: sort
      param_t() : _dist(kNoDist), _max(0), _shift(0), _param(0), _flags(0), _order(0) {}
      param_t(const dist_t aDist, 
              const int aMax, const int aShift, const double aParam,
              const int aFlags, const int aOrder)
             : _dist(aDist), _max(aMax), _shift(aShift), _param(aParam),
               _flags(aFlags), _order(aOrder) {}
      param_t(const std::string& aDistName, 
              const int aMax, const int aShift, const double aParam,
              const int aFlags, const int aOrder)
             : _dist(GenRandIntVec::str2dist(aDistName)), 
                _max(aMax), _shift(aShift), _param(aParam),
               _flags(aFlags), _order(aOrder)  {}
      inline uint   max()     const { return _max; }
      inline uint   mod()     const { return _max; }
      inline uint   div()     const { return _max; }
      inline dist_t dist()    const { return _dist; }
      inline uint   shift()   const { return _shift; }
      inline double param()   const { return _param; }
      inline double lambda()  const { return _param; } // for exponential distribution
      inline uint   flags()   const { return _flags; }
      inline bool   fill()    const { return (flags() & kFill); }
      inline bool   shuffle() const { return (flags() & kShuffle); }
      inline int    order()   const { return _order; }
      inline bool   permute() const { return (-1 == order()); }
      inline bool   sort()    const { return (+1 == order()); }
      inline const std::string& name() const { return GenRandIntVec::dist2str(_dist); }
    };
    typedef std::vector<std::string> string_vt;
    typedef std::exponential_distribution<double> dist_exp_t;
    typedef std::normal_distribution<double>      dist_norm_t;
    typedef std::poisson_distribution<uint>       dist_pois_t;
    typedef zipf_distribution<uint,double>        dist_zipf_t;
    // could add gamma distribution
  private:
    GenRandIntVec(const GenRandIntVec&) = delete;
    GenRandIntVec operator=(const GenRandIntVec&) = delete;
  public:
    GenRandIntVec();
  public:
    void generate(uint_vt& aVecOut, const uint aCard, const param_t& aParam, rng32_t& aRng);
  public:
    static dist_t             str2dist(const std::string& aDist);
    static const std::string& dist2str(const dist_t aDist);
  public:
    void generate_key(uint_vt& aVecOut,  const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_div(uint_vt& aVecOut,  const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_uni(uint_vt& aVecOut,  const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_exp(uint_vt& aVecOut,  const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_norm(uint_vt& aVecOut, const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_zipf(uint_vt& aVecOut, const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_self(uint_vt& aVecOut, const uint aCard, const param_t& aParam, rng32_t& aRng);
    void generate_pois(uint_vt& aVecOut, const uint aCard, const param_t& aParam, rng32_t& aRng);
  public:
    uint  genval_exp (const param_t& aParam, rng32_t& aRng);
    uint  genval_norm(const param_t& aParam, rng32_t& aRng);
    uint  genval_zipf(const param_t& aParam, rng32_t& aRng);
    uint  genval_self(const param_t& aParam, rng32_t& aRng);
    uint  genval_pois(const param_t& aParam, rng32_t& aRng);
  public:
    static void vec_permute(uint_vt& aVec, rng32_t& aRng);
    static void vec_sort(uint_vt& aVec);
    static void vec_revsort(uint_vt& aVec);
    static void freq_expand(uint_vt& aVecOut, const uint_vt& aFreqVec);
  private:
    uint_vt       _freq; // frequency vector, used if flag kFill is set then max < aSize must hold
    dist_exp_t   _dist_exp;
    dist_norm_t  _dist_norm;
    dist_pois_t  _dist_pois;
    dist_zipf_t* _dist_zipf;
  private:
    static const string_vt _distname;
};


#endif

#ifndef AGGREGATE_HH
#define AGGREGATE_HH

#include <iostream>
#include <iomanip>
#include <limits>
#include <math.h>

/*
 * Author: Guido Moerkotte (moerkotte@uni-mannheim.de)
 */
template<class Num>
class Aggregate {
  public:
    enum kind_t {
      count_e = 1,
      sum_e   = 2,
      avg_e   = 4,
      min_e   = 8,
      max_e   = 16,
      dist_e  = 32,
      qmid_e  = 64,
      qerr_e  = 128, // qerror of qmiddle
      middle_e = 256
    };
  public:
    Aggregate(int i = 0) : _dist(0), _min(std::numeric_limits<Num>::max()), _max(0), _sum(0), _sumsq(0), _count(0) { init(i); }
    Aggregate(const Num aMin, const Num aMax) : _dist(0), _min(aMin), _max(aMax), _sum(0), _sumsq(0), _count(0) {}
    ~Aggregate() {}
  public:
    inline
    void init(int i = 0) { _dist = i; 
                       _min =  std::numeric_limits<Num>::max();
                       _max =  std::numeric_limits<Num>::min();
                       _sum = _count = (Num) 0; 
                       _sumsq = (Num) 0;
                     };
    inline
    void init(Num aMin, Num aMax) { _dist = 0;
                                    _min = aMin;
                                    _max = aMax;
                                    _sum = _count = (Num) 0;
                                    _sumsq = (Num) 0;
                                  };
    inline
    void step(Num x) {
      if(x < _min) _min = x;
      if(x > _max) _max = x;
      _sum += x;
      _sumsq += (x*x);
      _count += (Num) 1;
    }
    inline void fin() {}
  public:
     inline int    dist() const { return _dist; }
     inline Num    count() const { return _count; }
     inline Num    min() const { return _min; }
     inline double avg() const { return ((double) sum() / (double) count()); }
     inline Num    max() const { return _max; }
     inline Num    sum() const { return _sum; }
     inline Num    sumsq() const { return _sumsq; }
     inline Num    span() const { return (max() - min()); }
     inline double middle() const { return (((double) (min() + max())) / 2.0); }
     inline double qmiddle() const { return (0.0 == (double) min()) ? sqrt((double) max()) : 
                                                                      sqrt((double) min() * (double) max()); }
     inline double qErrorOfQmiddle() const { return (max() / qmiddle()); }
     inline double qErrorOfAvg() const { return std::max<double>(avg()/min(), max()/avg()); }
     inline double qSpread() const { return (((double) max()) / ((double) min())); }
     inline Aggregate& operator=(const Aggregate& x) {
                          _dist   = x._dist;
                          _min    = x._min;
                          _max    = x._max;
                          _sum    = x._sum;
                          _sumsq  = x._sumsq;
                          _count  = x._count;
                          return (*this);
                       }
  public:
    std::ostream& print(std::ostream& os, int aKindSet, int aFieldWidth = 8) const {
      if(0 != (aKindSet & dist_e))  os << std::setw(aFieldWidth) << dist() << ' ';
      if(0 != (aKindSet & count_e)) os << std::setw(aFieldWidth) << count() << ' ';
      if(0 != (aKindSet & sum_e))   os << std::setw(aFieldWidth) << sum() << ' ';
      if(0 != (aKindSet & min_e))   os << std::setw(aFieldWidth) << min() << ' ';
      if(0 != (aKindSet & avg_e))   os << std::setw(aFieldWidth) << avg() << ' ';
      if(0 != (aKindSet & max_e))   os << std::setw(aFieldWidth) << max() << ' ';
      if(0 != (aKindSet & qmid_e))   os << std::setw(aFieldWidth) << qmiddle() << ' ';
      if(0 != (aKindSet & qerr_e))   os << std::setw(aFieldWidth) << qErrorOfQmiddle() << ' ';
      if(0 != (aKindSet & middle_e))   os << std::setw(aFieldWidth) << middle() << ' ';
      return os;
    }
  private:
    int _dist;
    Num _min;
    Num _max;
    Num _sum;
    Num _sumsq; // sum squared
    Num _count;
  private:
    static int _kindset;
  public:
    static int  setToPrint() { return _kindset; }
    static void setToPrint(int aKindset) { _kindset = aKindset; }
};

template<class Num>
std::ostream& 
operator<<(std::ostream& os, const Aggregate<Num>& aAggr) {
  if(0 == Aggregate<Num>::setToPrint()) {
    os << "min: " << aAggr.min() << ", ";
    os << "max: " << aAggr.max() << ", ";
    os << "avg: " << aAggr.avg() << ", ";
    os << "qmid: " << aAggr.qmiddle();
  } else {
    aAggr.print(os, Aggregate<Num>::setToPrint());
  }
  return os;
}


template<class Num>
int Aggregate<Num>::_kindset = 0;

#endif


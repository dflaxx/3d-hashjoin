#include "GenRandIntVec.hh"
#include <algorithm>

GenRandIntVec::GenRandIntVec() : _freq(), _dist_exp(), _dist_norm(), _dist_pois(), _dist_zipf(0) {}


GenRandIntVec::dist_t
GenRandIntVec::str2dist(const std::string& aDist) {
  for(int i = 0; i < kNoDist; ++i) {
    if(_distname[i] == aDist) {
      return ((dist_t) i);
    }
  }
  return kNoDist;
}

const std::string&
GenRandIntVec::dist2str(const dist_t aDist) {
  return _distname[aDist < kNoDist ? aDist : kNoDist];
}

void 
GenRandIntVec::generate(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  if(aVec.size() != aCard) {
    aVec.resize(aCard);
  }
  if(aParam.fill()) {
    assert((uint) aParam.max() <= aCard);
  }

  switch(aParam.dist()) {
    case kKey:  generate_key(aVec, aCard, aParam, aRng); break;
    case kDiv:  generate_div(aVec, aCard, aParam, aRng); break;
    case kUni:  generate_uni(aVec, aCard, aParam, aRng); break;
    case kExp:  generate_exp(aVec, aCard, aParam, aRng); break;
    case kNorm: generate_norm(aVec, aCard, aParam, aRng); break;
    case kZipf: generate_zipf(aVec, aCard, aParam, aRng); break;
    case kSelf: generate_self(aVec, aCard, aParam, aRng); break;
    case kPois: generate_pois(aVec, aCard, aParam, aRng); break;
    default: assert(0 == 1);
  }
}

void 
GenRandIntVec::generate_key(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  for(uint i = 0; i < aCard; ++i) {
    aVec[i] = i;
  }
  if(aParam.permute()) {
    vec_permute(aVec, aRng);
  } else
  if(aParam.sort()) {
    // is sorted already
  }
}

void 
GenRandIntVec::generate_div(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  const uint d = aParam.div();
  for(uint i = 0; i < aCard; ++i) {
    aVec[i] = i / d;
  }
  if(aParam.permute()) {
    vec_permute(aVec, aRng);
  } else
  if(aParam.sort()) {
    // is sorted already
  }
}

void 
GenRandIntVec::generate_uni(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  std::uniform_int_distribution<int> lDist(0, aParam.max() - 1);
  if(aParam.fill() || aParam.shuffle()) {
    _freq.resize(aParam.max());
    uint n = aCard;
    if(aParam.fill()) {
      for(size_t i = 0; i < _freq.size(); ++i) { _freq[i] = 1; } // fill
      n = aCard - _freq.size();
    }
    for(uint i = 0; i < n; ++i) {
      ++_freq.at(lDist(aRng));
    }
    freq_expand(aVec, _freq);
  } else {
    for(uint i = 0; i < aCard; ++i) {
      aVec[i] = lDist(aRng);
      // aVec[i] = aRng() % aParam.max();
    }
  }

  if(aParam.permute()) {
    vec_permute(aVec, aRng); // this is not really necessary
  } else
  if(aParam.sort()) {
    vec_sort(aVec);
  }
}

void 
GenRandIntVec::generate_exp(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  _dist_exp.param(dist_exp_t::param_type(aParam.lambda()));
  if(aParam.fill() || aParam.shuffle()) {
    _freq.resize(aParam.max());
    uint n = aCard;
    if(aParam.fill()) {
      for(size_t i = 0; i < _freq.size(); ++i) { _freq[i] = 1; } // fill
      n = aCard - _freq.size();
    }
    uint v = 0;
    for(uint i = 0; i < n; ++i) {
      v = genval_exp(aParam, aRng);
      assert((0 <= v) && (v < aParam.max()));
      ++_freq[v];
    }
    if(aParam.shuffle()) {
      vec_permute(_freq, aRng);
    }
    freq_expand(aVec, _freq);
  } else {
    for(uint i = 0; i < aCard; ++i) {
      aVec[i] = genval_exp(aParam, aRng);
    }
  }
  if(aParam.permute()) {
    vec_permute(aVec, aRng);
  } else
  if(aParam.sort()) {
    vec_sort(aVec);
  }
}

void 
GenRandIntVec::generate_norm(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
   const double lMu = ((double) aParam.max()) / 2;
  _dist_norm.param(dist_norm_t::param_type(lMu, aParam.param()));
  if(aParam.fill() || aParam.shuffle()) {
    _freq.resize(aParam.max());
    uint n = aCard;
    if(aParam.fill()) {
      for(size_t i = 0; i < _freq.size(); ++i) { _freq[i] = 1; } // fill
      n = aCard - _freq.size();
    }
    uint v = 0;
    for(uint i = 0; i < n; ++i) {
      v = genval_norm(aParam, aRng);
      assert((0 <= v) && (v < aParam.max()));
      ++_freq[v];
    }
    if(aParam.shuffle()) {
      vec_permute(_freq, aRng);
    }
    freq_expand(aVec, _freq);
  } else {
    for(uint i = 0; i < aCard; ++i) {
      aVec[i] = genval_norm(aParam, aRng);
    }
  }
  if(aParam.permute()) {
    vec_permute(aVec, aRng);
  } else
  if(aParam.sort()) {
    vec_sort(aVec);
  }
}

void 
GenRandIntVec::generate_zipf(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  assert(0 <= aParam.param());
  _dist_zipf = new dist_zipf_t(aParam.max(), aParam.param());
  if(aParam.fill() || aParam.shuffle()) {
    _freq.resize(aParam.max());
    uint n = aCard;
    if(aParam.fill()) {
      for(size_t i = 0; i < _freq.size(); ++i) { _freq[i] = 1; } // fill
      n = aCard - _freq.size();
    }
    uint v = 0;
    for(uint i = 0; i < n; ++i) {
      v = genval_zipf(aParam, aRng);
      assert((0 <= v) && (v < aParam.max()));
      ++_freq[v];
    }
    if(aParam.shuffle()) {
      vec_permute(_freq, aRng);
    }
    freq_expand(aVec, _freq);
  } else {
    for(uint i = 0; i < aCard; ++i) {
      aVec[i] = genval_zipf(aParam, aRng);
    }
  }
  if (aParam.permute()) {
    vec_permute(aVec, aRng);
  } else if (aParam.sort()) {
    vec_sort(aVec);
  }
  delete _dist_zipf;
  _dist_zipf = 0;
}

void 
GenRandIntVec::generate_self(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  assert((0 < aParam.param()) && (aParam.param() < 1));
  if(aParam.fill() || aParam.shuffle()) {
    _freq.resize(aParam.max());
    uint n = aCard;
    if(aParam.fill()) {
      for(size_t i = 0; i < _freq.size(); ++i) { _freq[i] = 1; } // fill
      n = aCard - _freq.size();
    }
    uint v = 0;
    for(uint i = 0; i < n; ++i) {
      v = genval_self(aParam, aRng);
      assert((0 <= v) && (v < aParam.max()));
      ++_freq[v];
    }
    if(aParam.shuffle()) {
      vec_permute(_freq, aRng);
    }
    freq_expand(aVec, _freq);
  } else {
    for(uint i = 0; i < aCard; ++i) {
      aVec[i] = genval_self(aParam, aRng);
    }
  }
  if(aParam.permute()) {
    vec_permute(aVec, aRng);
  } else
  if(aParam.sort()) {
    vec_sort(aVec);
  }
}

void 
GenRandIntVec::generate_pois(uint_vt& aVec, const uint aCard, const param_t& aParam, rng32_t& aRng) {
  assert(1 <=  aParam.param());
  _dist_pois.param(dist_pois_t::param_type(aParam.param()));

  if(aParam.fill() || aParam.shuffle()) {
    _freq.resize(aParam.max());
    uint n = aCard;
    if(aParam.fill()) {
      for(size_t i = 0; i < _freq.size(); ++i) { _freq[i] = 1; } // fill
      n = aCard - _freq.size();
    }
    uint v = 0;
    for(uint i = 0; i < n; ++i) {
      v = genval_pois(aParam, aRng);
      assert((0 <= v) && (v < aParam.max()));
      ++_freq[v];
    }
    if(aParam.shuffle()) {
      vec_permute(_freq, aRng);
    }
    freq_expand(aVec, _freq);
  } else {
    for(uint i = 0; i < aCard; ++i) {
      aVec[i] = genval_pois(aParam, aRng);
    }
  }
  if(aParam.permute()) {
    vec_permute(aVec, aRng);
  } else
  if(aParam.sort()) {
    vec_sort(aVec);
  }
}


uint  
GenRandIntVec::genval_exp(const param_t& aParam, rng32_t& aRng) {
  return (((int) std::floor((_dist_exp(aRng) * aParam.max()) + aParam.shift())) % aParam.max());
}

uint  
GenRandIntVec::genval_norm(const param_t& aParam, rng32_t& aRng) {
  int    lRes = 0;
  double v = 0;
  while(true) {
    v = _dist_norm(aRng);
    lRes = v;
    if((0 <= v) && (v < aParam.max())) {
      break;
    }
  }
  return ((lRes + aParam.shift()) % aParam.max());
}

uint
GenRandIntVec::genval_zipf(const param_t& aParam, rng32_t& aRng) {
  return (((*_dist_zipf)(aRng) - 1 + aParam.shift()) % aParam.max()); // dist_zipf returns v in [1,n]
}

/*
 *  implementation of SelfSimilar (80-20-rule) distribution by 
 *  Gray et al.: Quickly Generating Billion-Record Synthetic Databases.
 */

uint  
GenRandIntVec::genval_self(const param_t& aParam, rng32_t& aRng) {
  const double n = aParam.max();
  const double h = aParam.param();
  const double u = ((double) aRng()) / (((double) std::numeric_limits<uint32_t>::max())); // randf
  return ((int) (n * std::pow(u, std::log(h) / std::log(1.0 - h))));
}

/*
 *  implementation of poisson distribution by 
 *  Gray et al.: Quickly Generating Billion-Record Synthetic Databases.
 */
/*
int  
GenRandIntVec::genval_pois(const param_t& aParam, rng32_t& aRng) {
  const double lLambda = aParam.lambda();
  const double c = std::exp(-lLambda);
  double u = 0;
  double p = 1.0;
  int    n = 0;
  while(p > c) {
    u = ((double) aRng()) / (((double) std::numeric_limits<uint32_t>::max()));
    p *= u;
    ++n;
  }
  return (((n-1) + aParam.shift()) % aParam.max()); 
}
*/

uint
GenRandIntVec::genval_pois(const param_t& aParam, rng32_t& aRng) {
  return ((_dist_pois(aRng) + aParam.shift()) % aParam.max());
}


void
GenRandIntVec::vec_permute(uint_vt& aVec, rng32_t& aRng) {
  for(size_t i = aVec.size() - 1; i > 0; --i) {
    std::swap(aVec[i], aVec[aRng() % i]);
  }
}

void
GenRandIntVec::vec_sort(uint_vt& aVec) {
  std::sort(aVec.begin(), aVec.end());
}

void
GenRandIntVec::vec_revsort(uint_vt& aVec) {
  std::sort(aVec.begin(), aVec.end(), std::greater<int>());
}
void
GenRandIntVec::freq_expand(uint_vt& aVecOut, const uint_vt& aFreqVec) {
  size_t k = 0;
  for(size_t i = 0; i < aFreqVec.size(); ++i) {
    for(uint j = 0; j < aFreqVec[i]; ++j) {
      aVecOut[k++] = i;
    }
  }
}


const GenRandIntVec::string_vt GenRandIntVec::_distname = {
  "key",
  "div",
  "uni",
  "exp",
  "norm",
  "zipf",
  "self",
  "pois",
  "<invalid_distribution>"
};

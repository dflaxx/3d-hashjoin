#ifndef INFRA_HASH_TEMPLATES_HH
#define INFRA_HASH_TEMPLATES_HH

/*
 * Author: Guido Moerkotte (moerkotte@uni-mannheim.de)
 * Modifications by Daniel Flachs (flachs@uni-mannheim.de)
 */

#include <inttypes.h>
#include <cmath>
#include <limits>

#ifdef __SSE4_2__
  #include <nmmintrin.h>
#endif

namespace ht {

template<typename Tuint>
inline Tuint fibhash(const Tuint x, const Tuint n);

template<>
inline uint16_t
fibhash<uint16_t>(const uint16_t x, const uint16_t n) {
  const double b = 40503;
  const double a = b * (((double) 1.0)  / ((double) ((uint32_t) 1 << 16)));
  const double z = a * x;
  return (uint16_t) std::floor(n * (z - std::floor(z)));
}

template<>
inline uint32_t
fibhash<uint32_t>(const uint32_t x, const uint32_t n) {
  const double b = 2654435769;
  const double a = b * (((double) 1.0)  /   ( ((double) ((uint32_t) 1 << 16)) * ((double) ((uint32_t) 1 << 16)) ));
  const double z = a * x;
  return (uint32_t) std::floor(n * (z - std::floor(z)));
}

template<>
inline uint64_t
fibhash<uint64_t>(const uint64_t x, const uint64_t n) {
  const double b = 11400714819323198485UL;  // problematic, check!
  const double a = b * (((double) 1.0)  /   ( ((double) ((uint64_t) 1 << 32)) * ((double) ((uint64_t) 1 << 32)) ));
  const double z = a * x;
  return (uint64_t) std::floor(n * (z - std::floor(z)));
}

template<typename Tuint>
inline Tuint murmur_hash(Tuint x);

template<>
inline uint32_t
murmur_hash(uint32_t x) {
  x ^= x >> 16;
  x *= 0x85ebca6b;
  x ^= x >> 13;
  x *= 0xc2b2ae35;
  x ^= x >> 16;
  return x;
}

template<>
inline uint64_t 
murmur_hash(uint64_t x) {
  x ^= (x >> 33);
  x *= 0xFF51AFD7ED558CCD;
  x ^= (x >> 33);
  x *= 0xC4CEB9FE1A95EC63;
  x ^= (x >> 33);
  return x;
}

// Larson string hash
template<typename Tuint>
inline Tuint larson_hash(const char* x, const Tuint salt) {
  Tuint h = salt;
  while (*x) {
    h = h * 101 + static_cast<uint8_t>(*x++);
  }
  return h;
}
template<typename Tuint>
inline Tuint larson_hash(const char* x) {
  return larson_hash<Tuint>(x, 0);
}


template<typename Tuint>
class MultiplicativeHashing {
  public:
    MultiplicativeHashing(Tuint a, Tuint b) : _a(a), _b(b) {}
  public:
    inline Tuint hash(const Tuint x) const { return (_a * x + _b); }
  private:
    Tuint _a;
    Tuint _b;
};

template<typename Tuint>
inline Tuint boncz_hash(const Tuint);

template<>
inline uint32_t
boncz_hash(const uint32_t x) {
  return ((x >> 21) ^ (x >> 13) ^ (x >> 7) ^ (x));
}

template<>
inline uint64_t
boncz_hash(const uint64_t x) {
  return ((x >> 7) ^ (x >> 13) ^ (x >> 17) ^ (x >> 23) ^ (x >> 43));
}

#ifdef __SSE4_2__

template<typename Tuint>
inline Tuint hash_crc32(const Tuint aCurr, const Tuint aVal);

template<>
inline uint32_t
hash_crc32(const uint32_t aCurr, const uint32_t aVal) {
  return _mm_crc32_u32(aCurr, aVal);
}


// achtung: produziert nur 32 bit bei 64 bit inputs
template<>
inline uint64_t
hash_crc32(const uint64_t aCurr, const uint64_t aVal) {
  return _mm_crc32_u64(aCurr, aVal);
}

#endif

// same as classes


// identity
template<typename Tuint>
class HashId {
  public:
    typedef Tuint item_t;
    inline Tuint operator()(const Tuint x) { return x; }
    inline Tuint operator()(const Tuint x, const Tuint aMod) { return (x % aMod); }
};


template<typename Tuint>
class HashMurmur {
  public:
    typedef Tuint item_t;
    inline Tuint operator() (const Tuint x) const { return murmur_hash<Tuint>(x); }
    inline Tuint operator() (const Tuint x, const Tuint aMod) const { return (murmur_hash<Tuint>(x) % aMod); }
    // DF: changed operator() to be const member functions
};

template<typename Tuint>
class HashBoncz {
  public:
    typedef Tuint item_t;
    inline Tuint operator()(const Tuint x) { return boncz_hash<Tuint>(x); }
    inline Tuint operator()(const Tuint x, const Tuint aMod) { return (boncz_hash<Tuint>(x) % aMod); }
};

template<typename Tuint>
class HashFib {
  public:
    typedef Tuint item_t;
    inline Tuint operator()(const Tuint x) { return fibhash<Tuint>(x, std::numeric_limits<Tuint>::max()); }
    inline Tuint operator()(const Tuint x, const Tuint aMax) { return fibhash<Tuint>(x, aMax); }
};

#ifdef __SSE4_2__
template<typename Tuint>
class HashCrc32 {
  public:
    HashCrc32() : _curr((Tuint) _seed_initial) {}
    HashCrc32(const Tuint aSeed) : _curr(aSeed) {}
  public:
    inline void  seed(const Tuint aSeed) { _curr = aSeed; }
    inline Tuint hash(const Tuint aValu) { 
                   return (_curr = hash_crc32<Tuint>(_curr, aValu)); 
                 }
    inline Tuint hash() const { return _curr; }
  public:
    inline Tuint operator()(const Tuint aValu) const { 
                   return hash_crc32<Tuint>(_curr, aValu); 
                 }
    inline Tuint operator()(const Tuint aValu, const Tuint aMod) const { 
                   return (hash_crc32<Tuint>(_curr, aValu) % aMod); 
                 }
  private:
    Tuint _curr;
  public:
    static constexpr uint64_t _seed_initial = 0x6ca55437ae08fe14LL;
};

class HashCrc64 {
  public:
    HashCrc64() : _curr1(_seed_initial_1), _curr2(_seed_initial_2) {}
    HashCrc64(const uint64_t aSeed1, const uint64_t aSeed2) : _curr1(aSeed1), _curr2(aSeed2) {}
  public:
    inline uint64_t hash(const uint64_t aValu) const {
                      uint64_t lRes = hash_crc32<uint64_t>(_curr1, aValu) ^
                                      (hash_crc32<uint64_t>(_curr2, aValu) << 32);
                      return lRes;
                    }
    inline uint64_t operator()(const uint64_t aValu) const {
                      return hash(aValu);
                    }
  private:
    uint64_t _curr1;
    uint64_t _curr2;
  public:
    static constexpr uint64_t _seed_initial_1 = 0x832ca348c6f5dae9LL;
    static constexpr uint64_t _seed_initial_2 = 0x53c5a37a580308e3LL;
};

template<typename Tuint>
class HashCrc32x2 {
  public:
    HashCrc32x2() : _curr_1((Tuint) _seed_initial_1), 
                    _curr_2((Tuint) _seed_initial_2) {}
    HashCrc32x2(const Tuint aSeed1, const Tuint aSeed2) : _curr_1(aSeed1), _curr_2(aSeed2) {}
  public:
    inline void  seed(const Tuint aSeed1, const Tuint aSeed2) { _curr_1 = aSeed1; _curr_2 = aSeed2; }
    inline Tuint hash(const Tuint aValu) { 
                   _curr_1 = hash_crc32<Tuint>(_curr_1, aValu); 
                   _curr_2 = hash_crc32<Tuint>(_curr_2, aValu); 
                   return (_curr_1 ^ _curr_2);
                 }
    inline Tuint hash() const { return (_curr_1 ^ _curr_2); }
  public:
    inline Tuint operator()(const Tuint aValu) const { 
                   // return (hash_crc32<Tuint>(_curr_1, aValu) ^ hash_crc32<Tuint>(_curr_2, aValu));
                   // return hash_crc32<Tuint>(_curr_1, aValu);
                   // uint32_t h1 = hash_crc32<Tuint>(_curr_1, aValu);
                   // uint32_t h2 = hash_crc32<Tuint>(_curr_2, aValu);
                   // std::cout << "val: " << aValu << std::endl;
                   // std::cout << std::hex;
                   // std::cout << " h1: " << h1 << std::endl;
                   // std::cout << " h2: " << h2 << std::endl;
                   // std::cout << " hx: " << (h1 ^ h2) << std::endl;
                   // std::cout << std::dec;
                   return hash_crc32<Tuint>(_curr_1, aValu);
                 }
    inline Tuint operator()(const Tuint aValu, const Tuint aMod) const { 
                   // return ((hash_crc32<Tuint>(_curr_1, aValu) ^ hash_crc32<Tuint>(_curr_2, aValu)) % aMod); 
                   const uint64_t h1 = hash_crc32<Tuint>(_curr_1, aValu);
                   const uint64_t h2 = hash_crc32<Tuint>(_curr_2, aValu);
                   const uint64_t h  = (h1 << 32 | h2);
                   return (Tuint) (h % aMod);
                 }
  private:
    Tuint _curr_1;
    Tuint _curr_2;
  public:
    static constexpr uint64_t _seed_initial_1 = 0x832ca348c6f5dae9LL;
    static constexpr uint64_t _seed_initial_2 = 0x53c5a37a580308e3LL;
};
#endif

} // end namespace


#endif

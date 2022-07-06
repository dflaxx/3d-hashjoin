#ifndef INFRA_RESERVOIR_HH
#define INFRA_RESERVOIR_HH

#include <string.h>
#include <vector>
#include <queue>

/*
 * Author: Guido Moerkotte (moerkotte@uni-mannheim.de)
 */
template<class T>
class Reservoir {
  public:
    typedef unsigned int uint;
    typedef std::vector<T*> chunk_vt;
    typedef std::queue<T*> chunk_qt;
  private:
    Reservoir(const Reservoir&);
    Reservoir& operator=(const Reservoir&);
  public:
    class iterator {
      public:
        iterator() : _reservoir(0), _curr(0), _chunkend(0), _currchunkno(0), _ok(false) {}
        iterator(Reservoir<T>* aReservoir, 
                 T* aCurr, 
                 T* aChunkEnd,
                 uint aCurrChunkNo,
                 bool aOk) : _reservoir(aReservoir), _curr(aCurr), _chunkend(aChunkEnd),
                             _currchunkno(aCurrChunkNo), _ok(aOk) {}
      public:
        inline T* operator*() { return _curr; }
        inline bool ok() const { return _ok; }
        inline iterator& operator++();
      private:
        Reservoir<T>* _reservoir;
        T*   _curr;
        T*   _chunkend;
        uint _currchunkno;
        bool _ok;
    };
  public:
    Reservoir();
    Reservoir(const uint aLog2ChunkSize);
    ~Reservoir();
  public:
    void init(uint aLog2ChunkSize); // if called with Reservoir(), call init() before doing anything else
  public:
    T* get_new_entry(); // without initial T element
    T* get_new_entries(const uint k); // get k entries
    T* push_back(T*);
    inline const T* operator[](const uint i) const { return &(_chunks[idx_chunk(i)][idx_in_chunk(i)]); }
    inline       T* operator[](const uint i)       { return &(_chunks[idx_chunk(i)][idx_in_chunk(i)]); }
    inline uint cardinality() const { return _cardinality; }
    inline uint noChunks() const { return _chunks.size(); }
  public:
    void clear(); // empties container but keeps memory
    void erase(); // empties container and releases memory
  public:
    iterator begin() {
               return iterator(this,
                               first_of_chunk(0), 
                               end_of_chunk(0),
                               0,
                               (0 < cardinality()));
             }
  public:
    inline T* first_of_chunk(uint aChunkNo) { return (&(_chunks[aChunkNo][0])); }
    inline T* end_of_chunk(uint aChunkNo);
  private:
    T* allocChunk();
    inline uint idx_chunk(uint i) const { return (i >> _log2chunksize); }
    inline uint idx_in_chunk(uint i) const { return (i & _mask); }
  private:
    uint _log2chunksize;
    uint _mask;
    T*   _firstemptyinchunk; // first empty entry in last chunk
    T*   _endofchunk;        // end of last chunk
    uint _cardinality;
    chunk_vt _chunks;
    chunk_qt _freechunks;
};


template<class T>
Reservoir<T>::Reservoir() 
             : _log2chunksize(0), _mask(0),
               _firstemptyinchunk(0), _endofchunk(0),
               _cardinality(0), _chunks(), _freechunks() {
}

template<class T>
Reservoir<T>::Reservoir(const uint aLog2ChunkSize) 
             : _log2chunksize(aLog2ChunkSize), _mask((1 << aLog2ChunkSize) - 1), 
               _firstemptyinchunk(0), _endofchunk(0),
               _cardinality(0), _chunks(), _freechunks() {
  // allocChunk(); // XXX
}

template<class T>
Reservoir<T>::~Reservoir() {
  for(typename chunk_vt::iterator iter = _chunks.begin(); iter != _chunks.end(); ++iter) {
    delete[] (*iter);
  }
  for(;!_freechunks.empty(); _freechunks.pop()) {
    delete[] (_freechunks.front());
  }
}

template<class T>
void
Reservoir<T>::init(uint aLog2ChunkSize) {
  _log2chunksize = aLog2ChunkSize;
  _mask = ((1 << aLog2ChunkSize) - 1);
  _firstemptyinchunk = 0;
  _endofchunk = 0;
  _cardinality = 0;
  // allocChunk(); // XXX
}

template<class T>
T*
Reservoir<T>::get_new_entry() {
  if(_endofchunk <= _firstemptyinchunk) { // XXX was at the end XXX
    allocChunk();
  }
  T* lRes = _firstemptyinchunk;
  ++_firstemptyinchunk;
  ++_cardinality;
  return lRes;
}

template<class T>
T*
Reservoir<T>::get_new_entries(const uint k) {
  if(_endofchunk <= (_firstemptyinchunk + k)) { // XXX was at the end XXX
    allocChunk();
    if(_endofchunk <= (_firstemptyinchunk + k)) {
      return 0; // requested size larger than chunk size
    }
  }
  T* lRes = _firstemptyinchunk;
  _firstemptyinchunk += k;
  _cardinality += k;
  return lRes;
}

template<class T>
T*
Reservoir<T>::push_back(T* aT) {
  T* lRes = get_new_entry();
  (*lRes) = (*aT);
  return lRes;
}

template<class T>
T*
Reservoir<T>::allocChunk() {
  T* lChunk = 0;
  if(_freechunks.empty()) {
    lChunk = new T[1 << _log2chunksize];
    // lChunk = reinterpret_cast<T*>(new char[(1 << _log2chunksize) * sizeof(T)]);
    // memset(lChunk, 0, sizeof(T) * (1 << _log2chunksize));
  } else {
    lChunk = _freechunks.front();
    _freechunks.pop();
    memset((void*) lChunk, 0, sizeof(T) * (1 << _log2chunksize));
  }
  _chunks.push_back(lChunk);
  _firstemptyinchunk = lChunk;
  _endofchunk = (&(_chunks[_chunks.size() - 1][_mask + 1]));
  return _firstemptyinchunk;
}

template<class T>
void
Reservoir<T>::clear() {
  for(typename chunk_vt::iterator iter = _chunks.begin(); iter != _chunks.end(); ++iter) {
    _freechunks.push((*iter));
  }
  _chunks.clear();
  _firstemptyinchunk = 0;
  _endofchunk = 0;
  _cardinality = 0;
  // allocChunk(); // XXX
}

template<class T>
void
Reservoir<T>::erase() {
  for(;!_freechunks.empty(); _freechunks.pop()) {
    delete[] (_freechunks.front());
  }

  for(typename chunk_vt::iterator iter = _chunks.begin(); iter != _chunks.end(); ++iter) {
    delete[] (*iter);
  }
  _chunks.clear();
  // _freechunks.clear(); // done above
  _endofchunk = 0;
  _cardinality = 0;
  // allocChunk(); // XXX
}


template<class T>
T*
Reservoir<T>::end_of_chunk(uint aChunkNo) { 
  if((aChunkNo  + 1) == noChunks()) {
    return _firstemptyinchunk;
  } else {
    return (&(_chunks[aChunkNo][_mask + 1]));
  };
}

template<class T>
typename Reservoir<T>::iterator &
Reservoir<T>::iterator::operator++() {
   ++_curr;
   if(_chunkend <= _curr) {
     ++_currchunkno;
     if(_currchunkno < _reservoir->noChunks()) {
       _curr = _reservoir->first_of_chunk(_currchunkno);
       _chunkend = _reservoir->end_of_chunk(_currchunkno);
       _ok = (_curr < _chunkend);
     } else {
       _ok = false;
     }
   }
   return (*this);
}


#endif

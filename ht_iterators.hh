#pragma once

#include <functional>
#include <type_traits>

/*
 * Tnode must provide a method next() that returns a pointer_t
 * TODO: define concept for this property
 */
template <typename Tnode, bool IsConst>
class NodeIterator {
  public:
    using self_t = NodeIterator<Tnode, IsConst>;
    using node_t = Tnode;
    /*
     * avoid code duplication for const and non-const iterators using compile-time conditional                                
     * Sources:
     *   https://stackoverflow.com/questions/15182963/avoid-code-duplication-when-writing-iterator-and-const-iterator
     *   https://stackoverflow.com/questions/2150192/how-to-avoid-code-duplication-implementing-const-and-non-const-iterators
     *   https://www.drdobbs.com/the-standard-librarian-defining-iterato/184401331
     *   http://www.sj-vs.net/c-implementing-const_iterator-and-non-const-iterator-without-code-duplication/
     */
    using pointer_t = typename std::conditional<IsConst, const node_t*, node_t*>::type;
    using reference_t = typename std::conditional<IsConst, const node_t&, node_t&>::type;

  public:
    inline NodeIterator(pointer_t aNodePtr)
      : _cur(aNodePtr) {}

    // copy constructor: allow for implicit conversion from a regular iterator to a const_iterator 
    inline NodeIterator(const NodeIterator<Tnode, false>& aOther)
      : _cur(aOther._cur) {}

    inline self_t& operator++() {
      _cur = _cur->next();
      return *this;
    }

    inline self_t operator++(int) {
      self_t tmp = *this;
      ++(*this);
      return tmp;
    }

    inline reference_t operator*() {
      return *_cur;
    }

    inline pointer_t operator->() {
      return _cur;
    }

    inline bool operator==(const self_t& other) const {
      return _cur == other._cur;
    }

    inline bool operator!=(const self_t& other) const {
      return !(*this == other);
    }

    inline bool valid() const {
      return _cur != nullptr;
    }

    inline bool hasNext() const {
      return valid() && ((_cur->next()) != nullptr);
    }

    // allow const NodeIterator to access private members of
    // non-const NodeIterator (needed for copy constructor)
    friend NodeIterator<Tnode, true>;

  private:
    pointer_t _cur;
};

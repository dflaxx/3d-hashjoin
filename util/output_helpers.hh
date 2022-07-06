#pragma once

#include "standard_includes.hh"

#include <iomanip>
#include <iostream>
#include <ostream>


namespace df::infra {

/*
 * Modifier to be piped ('<<') into a std::ostream& os
 * to indent the next output to os by _margin + _level * _tabWidth characters.
 * Not sticky, i.e. does not change state of os and doesn't influence subsequent outputs.
 *
 * Usage examples:
 *   std::cout << indent(2) << "this is indented by two levels\n"
 *   std::cout << indent().marg(2).lvl(1) << "this has a margin of 2 chars, indented by 1 lvl\n"
 *
 * Source/inspiration: https://stackoverflow.com/a/29337924
 */
class indent {
  public:
    using self_t = indent;
  public:
             indent() : indent(0) {}
    explicit indent(const size_t aLvl, const size_t aMargin = 0, const size_t aTabWidth = 2,
                    const char aFill = ' ')
      : _level(aLvl), _tabWidth(aTabWidth), _margin(aMargin), _fill(aFill) {}
    friend inline std::ostream& operator<<(std::ostream& os, const indent& obj) {
      os << std::setw(obj._margin) << "";
      const char lOldFillChar = os.fill();
      os.fill(obj._fill);
      const size_t lTotalIndent = obj._level * obj._tabWidth;
      os << std::setw(lTotalIndent) << "";
      os.fill(lOldFillChar);
      return os;
    }
    inline indent& lvl (const size_t aLvl)  { _level = aLvl; return *this; }
    inline indent& tabw(const size_t aTabW) { _tabWidth = aTabW; return *this; }
    inline indent& marg(const size_t aMarg) { _margin = aMarg; return *this; }
    inline indent& fill(const char aFill)   { _fill = aFill; return *this; }

  private:
    size_t _level;
    size_t _tabWidth;
    size_t _margin;
    char   _fill;
};



};

#pragma once

#include "standard_includes.hh"

#include <filesystem>
#include <fstream>
#include <memory>
#include <ostream>


namespace df::infra {

class CSVWriter {
  public:
    enum class FlushPolicy {
      kDefault,  // the ostream's default
      kNewline,  // flush after each newline
      kField     // flush after each field
    };
  public:
    CSVWriter();
    CSVWriter(const std::string& aFile,
                       const bool aAppend = false,
                       const std::string& aSep = ";");
    CSVWriter(std::ostream& aStream,
                       const std::string& aSep = ";");
    CSVWriter(const CSVWriter&) = delete;
    CSVWriter(CSVWriter&&) = delete;

    inline ~CSVWriter() = default;

  public:
    //CSVWriter& writeLine(const std::string& aLine);

    template <typename T>
    inline
    CSVWriter& writeField(const T& aField, const bool aLastInLine = false) {
      if (_colIndex > 0) { *_os << _sep; }
      *_os << aField;
      if (_flushPolicy == FlushPolicy::kField) { _os->flush() ; }
      ++_colIndex;
      if (aLastInLine) { newline(); }
      return *this;
    }

    inline
    CSVWriter& newline() {
      *_os << NEWLINE;
      if (_flushPolicy == FlushPolicy::kField || _flushPolicy == FlushPolicy::kNewline) {
        _os->flush();
      }
      _colIndex = 0;
      return *this;
    }

    inline void setFlushPolicy(const FlushPolicy aPolicy) {
      _flushPolicy = aPolicy;
    }

  private:
    std::shared_ptr<std::ostream> _os;
    std::string _sep = ";";
    std::string _commentPrefix = "#";  // unused

    const static char NEWLINE = '\n';
    
    uint32_t _colIndex = 0;

    FlushPolicy _flushPolicy = FlushPolicy::kDefault;


};  // CSVWriter

}  // df::infra

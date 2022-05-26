#ifndef JULY_SYMBOLNAME_H
#define JULY_SYMBOLNAME_H

#include <string>

namespace model::symbols {

// TODO: Does this class make sense, if public constructor is available?
struct SymbolName {
  public:
    explicit SymbolName(std::string);
    SymbolName() = default;
    const std::string& get_name() const;
    bool operator<(const SymbolName& rhs) const;
    bool operator>(const SymbolName& rhs) const;
    bool operator<=(const SymbolName& rhs) const;
    bool operator>=(const SymbolName& rhs) const;
    bool operator==(const SymbolName& rhs) const;
    bool operator!=(const SymbolName& rhs) const;

  private:
    std::string name_;
};

}  // namespace symbols

#endif  //JULY_SYMBOLNAME_H

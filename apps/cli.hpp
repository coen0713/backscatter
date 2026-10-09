#pragma once

// Tiny "--key value" / "--flag" parser shared by the command-line tools.

#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace bsar::cli {

class Args {
 public:
  Args(int argc, char** argv, const std::set<std::string>& flags) {
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      if (a.rfind("--", 0) == 0) {
        const std::string key = a.substr(2);
        if (flags.count(key) != 0) {
          values_[key] = "1";
        } else if (i + 1 < argc) {
          values_[key] = argv[++i];
        } else {
          throw std::invalid_argument("missing value for --" + key);
        }
      } else {
        positional_.push_back(a);
      }
    }
  }

  [[nodiscard]] bool has(const std::string& key) const {
    used_.insert(key);
    return values_.count(key) != 0;
  }
  [[nodiscard]] std::string str(const std::string& key, const std::string& fallback) const {
    used_.insert(key);
    const auto it = values_.find(key);
    return it == values_.end() ? fallback : it->second;
  }
  [[nodiscard]] double num(const std::string& key, double fallback) const {
    used_.insert(key);
    const auto it = values_.find(key);
    if (it == values_.end()) {
      return fallback;
    }
    try {
      return std::stod(it->second);
    } catch (const std::exception&) {
      throw std::invalid_argument("--" + key + " expects a number, got '" + it->second + "'");
    }
  }
  [[nodiscard]] std::vector<double> nums(const std::string& key) const {
    std::vector<double> out;
    std::stringstream ss(str(key, ""));
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      out.push_back(std::stod(tok));
    }
    return out;
  }
  [[nodiscard]] const std::vector<std::string>& positional() const { return positional_; }

  /// Throw on options that no code path asked about (typos).
  void reject_unknown() const {
    for (const auto& [key, value] : values_) {
      if (used_.count(key) == 0) {
        throw std::invalid_argument("unknown option --" + key);
      }
    }
  }

 private:
  std::map<std::string, std::string> values_;
  std::vector<std::string> positional_;
  mutable std::set<std::string> used_;
};

}  // namespace bsar::cli

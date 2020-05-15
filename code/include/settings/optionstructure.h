/*!
 * \file OptionStructure.h
 * \brief Classes for different Options in rtsn
 * \author S. Schotthoefer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef OPTION_STRUCTURE_H
#define OPTION_STRUCTURE_H

#include <vector>
#include <string>
#include <map>
#include <sstream>

// Base Class for all kinds of options
class OptionBase {

private:
  std::vector<std::string> _value;
public:
  OptionBase() {};
  virtual  ~OptionBase() = 0;

  virtual std::string SetValue(std::vector<std::string> value){this->_value = value; return "";}
  std::vector<std::string> GetValue() {return _value;}
  virtual void SetDefault() = 0;

  std::string optionCheckMultipleValues(std::vector<std::string> & option_value, std::string type_id, std::string option_name) {
    if (option_value.size() != 1) {
      std::string newString;
      newString.append(option_name);
      newString.append(": multiple values for type ");
      newString.append(type_id);
      return newString;
    }
    return "";
  }

  std::string badValue(std::vector<std::string> & option_value, std::string type_id, std::string option_name) {
    std::string newString;
    newString.append(option_name);
    newString.append(": improper option value for type ");
    newString.append(type_id);
    return newString;
  }
};

inline OptionBase::~OptionBase() {}

template <class Tenum>
class OptionEnum : public OptionBase {

  std::map<std::string, Tenum> _map;
  Tenum & _field; // Reference to the fieldname
  Tenum _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionEnum(std::string option_field_name, const std::map<std::string, Tenum> m, Tenum & option_field, Tenum default_value) : _field(option_field) {
    this->_map = m;
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionEnum() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    // Check if there is more than one string
    std::string out = optionCheckMultipleValues(option_value, "enum", this->_name);
    if (out.compare("") != 0) {
      return out;
    }

    // Check to see if the enum value is in the map
    if (this->_map.find(option_value[0]) == _map.end()) {
      std::string str;
      str.append(this->_name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current RTSN options in config_template.cfg.");
      return str;
    }
    // If it is there, set the option value
    Tenum val = this->_map[option_value[0]];
    this->_field = val;
    return "";
  }

  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionDouble : public OptionBase {
  double & _field; // Reference to the fieldname
  double _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionDouble(std::string option_field_name, double & option_field, double default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionDouble() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    // check if there is more than one value
    std::string out = optionCheckMultipleValues(option_value, "double", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    double val;
    if (is >> val) {
      this->_field = val;
      return "";
    }
    return badValue(option_value, "double", this->_name);
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionString : public OptionBase {
  std::string & _field; // Reference to the fieldname
  std::string _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionString(std::string option_field_name, std::string & option_field, std::string default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionString() override {};

  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    // check if there is more than one value
    std::string out = optionCheckMultipleValues(option_value, "double", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    this->_field.assign(option_value[0]);
    return "";
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionInt : public OptionBase {
  int & _field; // Reference to the feildname
  int _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionInt(std::string option_field_name, int & option_field, int default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionInt() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "int", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    int val;
    if (is >> val) {
      this->_field = val;
      return "";
    }
    return badValue(option_value, "int", this->_name);
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionULong : public OptionBase {
  unsigned long & _field; // Reference to the feildname
  unsigned long _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionULong(std::string option_field_name, unsigned long & option_field, unsigned long default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionULong() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "unsigned long", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    unsigned long val;
    if (is >> val) {
      this->_field = val;
      return "";
    }
    return badValue(option_value, "unsigned long", this->_name);
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class COptionUShort : public OptionBase {
  unsigned short & _field; // Reference to the feildname
  unsigned short _def; // Default value
  std::string _name; // identifier for the option

public:
  COptionUShort(std::string option_field_name, unsigned short & option_field, unsigned short default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~COptionUShort() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "unsigned short", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    unsigned short val;
    if (is >> val) {
      this->_field = val;
      return "";
    }
    return badValue(option_value, "unsigned short", this->_name);
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionLong : public OptionBase {
  long & _field; // Reference to the feildname
  long _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionLong(std::string option_field_name, long & option_field, long default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionLong() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "long", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    long val;
    if (is >> val) {
      this->_field = val;
      return "";
    }
    return badValue(option_value, "long", this->_name);
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionBool : public OptionBase {
  bool & _field; // Reference to the feildname
  bool _def; // Default value
  std::string _name; // identifier for the option

public:
  OptionBool(std::string option_field_name, bool & option_field, bool default_value) : _field(option_field) {
    this->_def = default_value;
    this->_name = option_field_name;
  }

  ~OptionBool() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    // check if there is more than one value
    std::string out = optionCheckMultipleValues(option_value, "bool", this->_name);
    if (out.compare("") != 0) {
      return out;
    }
    if (option_value[0].compare("YES") == 0) {
      this->_field = true;
      return "";
    }
    if (option_value[0].compare("NO") == 0) {
      this->_field = false;
      return "";
    }
    return badValue(option_value, "bool", this->_name);
  }
  void SetDefault() override {
    this->_field = this->_def;
  }
};

class OptionStringList : public OptionBase {
  std::vector<std::string> & _field; // Reference to the fieldname
  std::string _name; // identifier for the option
  unsigned short & _size;

public:
  OptionStringList(std::string option_field_name, unsigned short & list_size, std::vector<std::string> & option_field) : _field(option_field), _size(list_size) {
    this->_name = option_field_name;
  }

  ~OptionStringList() override {};

  std::string SetValue(std::vector<std::string> option_value) override {
    OptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      this->_size = 0;
      return "";
    }
    this->_size = option_size;

    // Parse all of the options
    this->_field.resize(this->_size);
    for (unsigned long i  = 0; i < option_size; i++) {
      this->_field.insert(this->_field.begin()+i,option_value[i]);
    }
    return "";
  }

  void SetDefault() override {
    this->_size = 0; // There is no default value for list
  }
};

#endif // OPTION_STRUCTURE_H

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
class COptionBase {

private:
  std::vector<std::string> value;
public:
  COptionBase() {};
  virtual  ~COptionBase() = 0;

  virtual std::string SetValue(std::vector<std::string> value){this->value = value; return "";}
  std::vector<std::string> GetValue() {return value;}
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

inline COptionBase::~COptionBase() {}


template <class Tenum>
class COptionEnum : public COptionBase {

  std::map<std::string, Tenum> m;
  Tenum & field; // Reference to the feildname
  Tenum def; // Default value
  std::string name; // identifier for the option

public:
  COptionEnum(std::string option_field_name, const std::map<std::string, Tenum> m, Tenum & option_field, Tenum default_value) : field(option_field) {
    this->m = m;
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionEnum() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    // Check if there is more than one string
    std::string out = optionCheckMultipleValues(option_value, "enum", this->name);
    if (out.compare("") != 0) {
      return out;
    }

    // Check to see if the enum value is in the map
    if (this->m.find(option_value[0]) == m.end()) {
      std::string str;
      str.append(this->name);
      str.append(": invalid option value ");
      str.append(option_value[0]);
      str.append(". Check current RTSN options in config_template.cfg.");
      return str;
    }
    // If it is there, set the option value
    Tenum val = this->m[option_value[0]];
    this->field = val;
    return "";
  }

  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionDouble : public COptionBase {
  double & field; // Reference to the fieldname
  double def; // Default value
  std::string name; // identifier for the option

public:
  COptionDouble(std::string option_field_name, double & option_field, double default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionDouble() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    std::string out = optionCheckMultipleValues(option_value, "double", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    double val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "double", this->name);
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionString : public COptionBase {
  std::string & field; // Reference to the fieldname
  std::string def; // Default value
  std::string name; // identifier for the option

public:
  COptionString(std::string option_field_name, std::string & option_field, std::string default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionString() override {};

  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    std::string out = optionCheckMultipleValues(option_value, "double", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    this->field.assign(option_value[0]);
    return "";
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionInt : public COptionBase {
  int & field; // Reference to the feildname
  int def; // Default value
  std::string name; // identifier for the option

public:
  COptionInt(std::string option_field_name, int & option_field, int default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionInt() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "int", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    int val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "int", this->name);
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionULong : public COptionBase {
  unsigned long & field; // Reference to the feildname
  unsigned long def; // Default value
  std::string name; // identifier for the option

public:
  COptionULong(std::string option_field_name, unsigned long & option_field, unsigned long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionULong() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "unsigned long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    unsigned long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned long", this->name);
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionUShort : public COptionBase {
  unsigned short & field; // Reference to the feildname
  unsigned short def; // Default value
  std::string name; // identifier for the option

public:
  COptionUShort(std::string option_field_name, unsigned short & option_field, unsigned short default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionUShort() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "unsigned short", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    unsigned short val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "unsigned short", this->name);
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionLong : public COptionBase {
  long & field; // Reference to the feildname
  long def; // Default value
  std::string name; // identifier for the option

public:
  COptionLong(std::string option_field_name, long & option_field, long default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionLong() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    std::string out = optionCheckMultipleValues(option_value, "long", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    std::istringstream is(option_value[0]);
    long val;
    if (is >> val) {
      this->field = val;
      return "";
    }
    return badValue(option_value, "long", this->name);
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionBool : public COptionBase {
  bool & field; // Reference to the feildname
  bool def; // Default value
  std::string name; // identifier for the option

public:
  COptionBool(std::string option_field_name, bool & option_field, bool default_value) : field(option_field) {
    this->def = default_value;
    this->name = option_field_name;
  }

  ~COptionBool() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    // check if there is more than one value
    std::string out = optionCheckMultipleValues(option_value, "bool", this->name);
    if (out.compare("") != 0) {
      return out;
    }
    if (option_value[0].compare("YES") == 0) {
      this->field = true;
      return "";
    }
    if (option_value[0].compare("NO") == 0) {
      this->field = false;
      return "";
    }
    return badValue(option_value, "bool", this->name);
  }
  void SetDefault() override {
    this->field = this->def;
  }
};

class COptionStringList : public COptionBase {
  std::string * & field; // Reference to the fieldname
  std::string name; // identifier for the option
  unsigned short & size;

public:
  COptionStringList(std::string option_field_name, unsigned short & list_size, std::string * & option_field) : field(option_field), size(list_size) {
    this->name = option_field_name;
  }

  ~COptionStringList() override {};
  std::string SetValue(std::vector<std::string> option_value) override {
    COptionBase::SetValue(option_value);
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if (option_size == 1 && option_value[0].compare("NONE") == 0) {
      this->size = 0;
      return "";
    }
    this->size = option_size;

    // Parse all of the options
    std::string * vals = new std::string[option_size];
    for (unsigned long i  = 0; i < option_size; i++) {
      vals[i].assign(option_value[i]);
    }
    this->field = vals;
    return "";
  }

  void SetDefault() override {
    this->size = 0; // There is no default value for list
  }
};

#endif // OPTION_STRUCTURE_H

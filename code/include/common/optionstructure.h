/*!
 * \file OptionStructure.h
 * \brief Classes for different Options in rtsn
 * \author S. Schotthoefer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#ifndef OPTION_STRUCTURE_H
#define OPTION_STRUCTURE_H

#include <map>
#include <sstream>
#include <string>
#include <vector>

// Base Class for all kinds of options
class OptionBase
{

  private:
    std::vector<std::string> _value;

  public:
    OptionBase() {}

    virtual ~OptionBase() {}

    std::vector<std::string> GetValue();

    virtual std::string SetValue( std::vector<std::string> value );

    virtual void SetDefault() = 0;

    std::string optionCheckMultipleValues( std::vector<std::string>& option_value, std::string type_id, std::string option_name );

    std::string badValue( std::vector<std::string>& option_value, std::string type_id, std::string option_name );
};

class OptionDouble : public OptionBase
{
    double& _field;       // Reference to the fieldname
    double _def;          // Default value
    std::string _name;    // identifier for the option

  public:
    OptionDouble( std::string option_field_name, double& option_field, double default_value );

    ~OptionDouble() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionString : public OptionBase
{
    std::string& _field;    // Reference to the fieldname
    std::string _def;       // Default value
    std::string _name;      // identifier for the option

  public:
    OptionString( std::string option_field_name, std::string& option_field, std::string default_value );

    ~OptionString() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionInt : public OptionBase
{
    int& _field;          // Reference to the feildname
    int _def;             // Default value
    std::string _name;    // identifier for the option

  public:
    OptionInt( std::string option_field_name, int& option_field, int default_value );

    ~OptionInt() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionULong : public OptionBase
{
    unsigned long& _field;    // Reference to the feildname
    unsigned long _def;       // Default value
    std::string _name;        // identifier for the option

  public:
    OptionULong( std::string option_field_name, unsigned long& option_field, unsigned long default_value );

    ~OptionULong() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionUShort : public OptionBase
{
    unsigned short& _field;    // Reference to the feildname
    unsigned short _def;       // Default value
    std::string _name;         // identifier for the option

  public:
    OptionUShort( std::string option_field_name, unsigned short& option_field, unsigned short default_value );

    ~OptionUShort() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionLong : public OptionBase
{
    long& _field;         // Reference to the feildname
    long _def;            // Default value
    std::string _name;    // identifier for the option

  public:
    OptionLong( std::string option_field_name, long& option_field, long default_value );

    ~OptionLong() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionBool : public OptionBase
{
    bool& _field;         // Reference to the feildname
    bool _def;            // Default value
    std::string _name;    // identifier for the option

  public:
    OptionBool( std::string option_field_name, bool& option_field, bool default_value );

    ~OptionBool() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

class OptionStringList : public OptionBase
{
    std::vector<std::string>& _field;    // Reference to the fieldname
    std::string _name;                   // identifier for the option
    unsigned short& _size;

  public:
    OptionStringList( std::string option_field_name, unsigned short& list_size, std::vector<std::string>& option_field );

    ~OptionStringList() override {}

    std::string SetValue( std::vector<std::string> option_value ) override;

    void SetDefault() override;
};

// Template classes (TODO: Find a way to split code to .cpp)

template <class Tenum> class OptionEnum : public OptionBase
{

    std::map<std::string, Tenum> _map;
    Tenum& _field;        // Reference to the fieldname
    Tenum _def;           // Default value
    std::string _name;    // identifier for the option

  public:
    OptionEnum( std::string option_field_name, const std::map<std::string, Tenum> m, Tenum& option_field, Tenum default_value )
        : _field( option_field ) {
        this->_map  = m;
        this->_def  = default_value;
        this->_name = option_field_name;
    }

    ~OptionEnum() override {}

    std::string SetValue( std::vector<std::string> option_value ) override {
        OptionBase::SetValue( option_value );
        // Check if there is more than one string
        std::string out = optionCheckMultipleValues( option_value, "enum", this->_name );
        if( out.compare( "" ) != 0 ) {
            return out;
        }
        // Check to see if the enum value is in the map
        if( this->_map.find( option_value[0] ) == _map.end() ) {
            std::string str;
            str.append( this->_name );
            str.append( ": invalid option value " );
            str.append( option_value[0] );
            str.append( ". Check current RTSN options in config_template.cfg." );
            return str;
        }
        // If it is there, set the option value
        Tenum val    = this->_map[option_value[0]];
        this->_field = val;
        return "";
    }

    void SetDefault() override { this->_field = this->_def; }
};

template <class Tenum> class OptionEnumList : public OptionBase
{

    std::map<std::string, Tenum> _map;
    std::vector<Tenum>& _field;    // Reference to the fieldname
    std::string _name;             // identifier for the option
    unsigned short& _size;

  public:
    OptionEnumList( std::string option_field_name, const std::map<std::string, Tenum> m, std::vector<Tenum>& option_field, unsigned short& list_size )
        : _field( option_field ), _size( list_size ) {
        this->_map  = m;
        this->_name = option_field_name;
    }

    ~OptionEnumList() override {}

    std::string SetValue( std::vector<std::string> option_value ) override {
        OptionBase::SetValue( option_value );
        if( option_value.size() == 1 && option_value[0].compare( "NONE" ) == 0 ) {
            this->_size = 0;
            return "";
        }
        // size is the length of the option list
        this->_size = option_value.size();
        std::vector<Tenum> enums( _size );
        // unsigned short* enums = new unsigned short[_size];
        for( int i = 0; i < this->_size; i++ ) {
            // Check to see if the enum value is in the map
            if( this->_map.find( option_value[i] ) == _map.end() ) {
                std::string str;
                str.append( this->_name );
                str.append( ": invalid option value " );
                str.append( option_value[i] );
                str.append( ". Check current RTSN options in config_template.cfg." );
                return str;
            }
            // If it is there, set the option value
            enums[i] = this->_map[option_value[i]];
        }
        this->_field = enums;
        return "";
    }

    void SetDefault() override {
        // No default to set
        _size = 0;
    }
};

#endif    // OPTION_STRUCTURE_H

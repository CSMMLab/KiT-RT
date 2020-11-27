/*!
 * \file Config.cpp
 * \brief Classes for different Optiontypes in rtsn
 * \author S. Schotthoefer
 *
 * Disclaimer: This class structure was copied and modifed with open source permission from SU2 v7.0.3 https://su2code.github.io/
 */

#include "common/optionstructure.h"

// --- Members of OptionBase ----

OptionBase::~OptionBase() {}

std::vector<std::string> OptionBase::GetValue() { return _value; }

std::string OptionBase::SetValue( std::vector<std::string> value ) {
    this->_value = value;
    return "";
}

std::string OptionBase::OptionCheckMultipleValues( std::vector<std::string>& option_value, std::string type_id, std::string option_name ) {
    if( option_value.size() != 1 ) {
        std::string newString;
        newString.append( option_name );
        newString.append( ": multiple values for type " );
        newString.append( type_id );
        return newString;
    }
    return "";
}

std::string OptionBase::BadValue( std::vector<std::string>& option_value, std::string type_id, std::string option_name ) {
    std::string newString;
    newString.append( option_name );
    newString.append( ": improper option value for type " );
    newString.append( type_id );
    return newString;
}

// ---- Memebers of OptionDouble

OptionDouble::OptionDouble( std::string option_field_name, double& option_field, double default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionDouble::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    // check if there is more than one value
    std::string out = OptionCheckMultipleValues( option_value, "double", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    std::istringstream is( option_value[0] );
    double val;
    if( is >> val ) {
        this->_field = val;
        return "";
    }
    return BadValue( option_value, "double", this->_name );
}

void OptionDouble::SetDefault() { this->_field = this->_def; }

// ---- Members of OptionString

OptionString::OptionString( std::string option_field_name, std::string& option_field, std::string default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionString::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    // check if there is more than one value
    std::string out = OptionCheckMultipleValues( option_value, "double", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    this->_field.assign( option_value[0] );
    return "";
}

void OptionString::SetDefault() { this->_field = this->_def; }

// --- Members of OptionInt

OptionInt::OptionInt( std::string option_field_name, int& option_field, int default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionInt::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    std::string out = OptionCheckMultipleValues( option_value, "int", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    std::istringstream is( option_value[0] );
    int val;
    if( is >> val ) {
        this->_field = val;
        return "";
    }
    return BadValue( option_value, "int", this->_name );
}
void OptionInt::SetDefault() { this->_field = this->_def; }

// ---- Members of OptionULong

OptionULong::OptionULong( std::string option_field_name, unsigned long& option_field, unsigned long default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionULong::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    std::string out = OptionCheckMultipleValues( option_value, "unsigned long", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    std::istringstream is( option_value[0] );
    unsigned long val;
    if( is >> val ) {
        this->_field = val;
        return "";
    }
    return BadValue( option_value, "unsigned long", this->_name );
}

void OptionULong::SetDefault() { this->_field = this->_def; }

// ---- Members of OptionUShort

OptionUShort::OptionUShort( std::string option_field_name, unsigned short& option_field, unsigned short default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionUShort::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    std::string out = OptionCheckMultipleValues( option_value, "unsigned short", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    std::istringstream is( option_value[0] );
    unsigned short val;
    if( is >> val ) {
        this->_field = val;
        return "";
    }
    return BadValue( option_value, "unsigned short", this->_name );
}

void OptionUShort::SetDefault() { this->_field = this->_def; }

// ---- Members of OptionLong

OptionLong::OptionLong( std::string option_field_name, long& option_field, long default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionLong::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    std::string out = OptionCheckMultipleValues( option_value, "long", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    std::istringstream is( option_value[0] );
    long val;
    if( is >> val ) {
        this->_field = val;
        return "";
    }
    return BadValue( option_value, "long", this->_name );
}

void OptionLong::SetDefault() { this->_field = this->_def; }

// ---- Members of OptionBool

OptionBool::OptionBool( std::string option_field_name, bool& option_field, bool default_value ) : _field( option_field ) {
    this->_def  = default_value;
    this->_name = option_field_name;
}

std::string OptionBool::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    // check if there is more than one value
    std::string out = OptionCheckMultipleValues( option_value, "bool", this->_name );
    if( out.compare( "" ) != 0 ) {
        return out;
    }
    if( option_value[0].compare( "YES" ) == 0 ) {
        this->_field = true;
        return "";
    }
    if( option_value[0].compare( "NO" ) == 0 ) {
        this->_field = false;
        return "";
    }
    return BadValue( option_value, "bool", this->_name );
}

void OptionBool::SetDefault() { this->_field = this->_def; }

// --- members of OptionStringList

OptionStringList::OptionStringList( std::string option_field_name, unsigned short& list_size, std::vector<std::string>& option_field )
    : _field( option_field ), _size( list_size ) {
    this->_name = option_field_name;
}

std::string OptionStringList::SetValue( std::vector<std::string> option_value ) {
    OptionBase::SetValue( option_value );
    // The size is the length of option_value
    unsigned short option_size = option_value.size();
    if( option_size == 1 && option_value[0].compare( "NONE" ) == 0 ) {
        this->_size = 0;
        return "";
    }
    this->_size = option_size;

    // Parse all of the options
    this->_field.resize( this->_size );
    for( unsigned long i = 0; i < option_size; i++ ) {
        this->_field.at( i ) = option_value[i];
    }
    return "";
}

void OptionStringList::SetDefault() {
    this->_size = 0;    // There is no default value for list
}

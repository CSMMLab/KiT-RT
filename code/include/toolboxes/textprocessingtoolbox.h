/*!
 * \file textprocessingtoolbox.h
 * \brief File with helper functions for text processing
 */

#ifndef TEXTPROCESSINGTOOLBOX_H
#define TEXTPROCESSINGTOOLBOX_H

#include "common/typedef.h"
#include <iostream>
#include <string>
#include <vector>

namespace TextProcessingToolbox {

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in,out] str - string we want to convert
 */
inline void StringToUpperCase( std::string& str ) { std::transform( str.begin(), str.end(), str.begin(), ::toupper ); }

/*!
 * \brief utility function for splitting strings to at delimiter
 * \param [in] str - string we want to split, delimiter - delimiter character
 *        [out] vector of string tokens that were separated by the delimiter
 */
inline std::vector<std::string> Split( const std::string& s, char delimiter ) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream( s );
    while( std::getline( tokenStream, token, delimiter ) ) {
        tokens.push_back( token );
    }
    return tokens;
}

/*!
 * \brief utility function for getting walltime
 * \param  [in,out] V
 */
// TODO

/*!
 * \brief utility function for printing a VectorVector
 * \param vectorIn VectorVector we want to print
 */
inline void PrintVectorVector( const VectorVector vectorIn ) {
    unsigned dimOuter = vectorIn.size();
    unsigned dimInner = vectorIn[0].size();

    for( unsigned idx_Outer = 0; idx_Outer < dimOuter; idx_Outer++ ) {
        for( unsigned idx_Inner = 0; idx_Inner < dimInner; idx_Inner++ ) {
            if( vectorIn[idx_Outer][idx_Inner] > 0 )
                printf( "|  %.2f ", vectorIn[idx_Outer][idx_Inner] );
            else
                printf( "| %.2f ", vectorIn[idx_Outer][idx_Inner] );
        }
        printf( " |\n" );
    }
}

/*!
 * \brief utility function for returning the last number in a string
 * \param str string to be checked
 */
inline int GetTrailingNumber( std::string const& str ) { return std::stoi( str.substr( str.find_first_of( "0123456789" ), str.length() - 1 ) ); }

/*!
 * \brief utility function for checking if a string has a certain ending
 * \param value string to be checked
 * \param ending string to be checked for
 */
inline bool StringEndsWith( std::string const& value, std::string const& ending ) {
    if( ending.size() > value.size() ) return false;
    return std::equal( ending.rbegin(), ending.rend(), value.rbegin() );
}

}    // namespace TextProcessingToolbox

#endif    // TEXTPROCESSINGTOOLBOX_H

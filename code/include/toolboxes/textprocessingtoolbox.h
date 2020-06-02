/*!
 * \file TextProcessingToolbox.h
 * \brief File with helper functions for text processing
 *
 *
 */

#ifndef TEXTPROCESSINGTOOLBOX_H
#define TEXTPROCESSINGTOOLBOX_H

#include <string>
#include <vector>

#include "typedef.h"

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
 * \param [in] VectorVector we want to print
 */
inline const void PrintVectorVector( const VectorVector vectorIn ) {
    unsigned dimOuter = vectorIn.size();
    unsigned dimInner = vectorIn[0].size();

    for( unsigned idx_Outer = 0; idx_Outer < dimOuter; idx_Outer++ ) {
        for( unsigned idx_Inner = 0; idx_Inner < dimInner; idx_Inner++ ) {
            if( vectorIn[idx_Outer][idx_Inner] > 0 )
                printf( "|  %.2f ", vectorIn[idx_Outer][idx_Inner] );
            else
                printf( "| %.2f ", vectorIn[idx_Outer][idx_Inner] );
        }
        std::cout << " |\n";
    }
}

/*!
 * \brief utility function for printing a Matrix
 * \param [in] matrixIn - Matrix we want to print
 */
// inline const void PrintMatrix( const Matrix matrixIn ) {
//    unsigned dimOuter = vectorIn.size();
//    unsigned dimInner = vectorIn[0].size();
//
//    for( unsigned idx_Outer = 0; idx_Outer < dimOuter; idx_Outer++ ) {
//        for( unsigned idx_Inner = 0; idx_Inner < dimInner; idx_Inner++ ) {
//            if( vectorIn[idx_Outer][idx_Inner] > 0 )
//                printf( "|  %.2f ", vectorIn[idx_Outer][idx_Inner] );
//            else
//                printf( "| %.2f ", vectorIn[idx_Outer][idx_Inner] );
//        }
//        std::cout << " |\n";
//    }
//}
}    // namespace TextProcessingToolbox

#endif    // TEXTPROCESSINGTOOLBOX_H

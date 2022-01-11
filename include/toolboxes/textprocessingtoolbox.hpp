/*!
 * \file textprocessingtoolbox.h
 * \brief File with helper functions for text processing
 */

#ifndef TEXTPROCESSINGTOOLBOX_H
#define TEXTPROCESSINGTOOLBOX_H

#include "common/typedef.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace TextProcessingToolbox {

/*!
 * @brief utility function for converting strings to uppercase
 * @param str - string to be converted
 */
inline void StringToUpperCase( std::string& str ) { std::transform( str.begin(), str.end(), str.begin(), ::toupper ); }

/*!
 * @brief utility function for splitting strings to at delimiter
 * @param s string to be split
 * @param delimiter delimiter character
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
 * @brief utility function for printing a VectorVector
 * @param vectorIn VectorVector we want to print
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

inline void PrintMatrix( const Matrix mat ) { std::cout << mat << std::endl; }

inline void PrintMatrixToFile( const Matrix mat, std::string filename, unsigned matsize ) {
    std::ofstream myfile;
    myfile.open( filename );
    for( unsigned i = 0; i < matsize; i++ ) {
        for( unsigned j = 0; j < matsize - 1; j++ ) {
            myfile << mat( i, j ) << ",";
        }
        myfile << mat( i, matsize - 1 );
        myfile << "\n";
    }
    myfile.close();
}

inline void PrintVectorToFile( const Vector vec, std::string filename, unsigned vecsize ) {
    std::ofstream myfile;
    myfile.open( filename );
    for( unsigned i = 0; i < vecsize; i++ ) {
        myfile << vec[i] << "\n";
    }
    myfile.close();
}

inline void PrintVectorVectorToFile( const VectorVector vecvec, std::string filename, unsigned size_outer, unsigned size_inner ) {
    std::ofstream myfile;
    myfile.open( filename );
    for( unsigned i = 0; i < size_outer; i++ ) {
        for( unsigned j = 0; j < size_inner - 1; j++ ) {
            myfile << vecvec[i][j] << ",";
        }
        myfile << vecvec[i][size_inner - 1];
        myfile << ";\n";
    }
    myfile.close();
}

/*!
 * @brief utility function for returning the last number in a string
 * @param str string to be checked
 */
inline int GetTrailingNumber( std::string const& str ) { return std::stoi( str.substr( str.find_first_of( "0123456789" ), str.length() - 1 ) ); }

/*!
 * @brief utility function for checking if a string has a certain ending
 * @param value string to be checked
 * @param ending string to be checked for
 */
inline bool StringEndsWith( std::string const& value, std::string const& ending ) {
    if( ending.size() > value.size() ) return false;
    return std::equal( ending.rbegin(), ending.rend(), value.rbegin() );
}

}    // namespace TextProcessingToolbox

#endif    // TEXTPROCESSINGTOOLBOX_H

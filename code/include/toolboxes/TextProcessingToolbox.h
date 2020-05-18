/*!
 * \file TextProcessingToolbox.h
 * \brief File with helper functions for text processing
 *
 *
 */
#ifndef TEXTPROCESSINGTOOLBOX_H
#define TEXTPROCESSINGTOOLBOX_H

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in,out] str - string we want to convert
 */
inline void StringToUpperCase( std::string& str ) {
    for( auto& c : str ) c = std::toupper( c );
}

#endif    // TEXTPROCESSINGTOOLBOX_H

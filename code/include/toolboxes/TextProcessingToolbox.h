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
inline void StringToUpperCase(std::string & str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

#endif // TEXTPROCESSINGTOOLBOX_H

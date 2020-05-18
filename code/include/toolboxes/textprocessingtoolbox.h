/*!
 * \file TextProcessingToolbox.h
 * \brief File with helper functions for text processing
 *
 *
 */

#ifndef TEXTPROCESSINGTOOLBOX_H
#define TEXTPROCESSINGTOOLBOX_H

#include <vector>
#include <string>


namespace TextProcessingToolbox {

/*!
 * \brief utility function for converting strings to uppercase
 * \param[in,out] str - string we want to convert
 */
inline void StringToUpperCase(std::string & str) {
  std::transform(str.begin(), str.end(), str.begin(), ::toupper);
}

/*!
 * \brief utility function for splitting strings to at delimiter
 * \param [in,out] str - string we want to split, delimiter - delimiter character
 */
inline std::vector<std::string> Split(const std::string& s, char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter))
  {
    tokens.push_back(token);
  }
  return tokens;
}

/*!
 * \brief utility function for getting walltime
 * \param
 */
//TODO

}

#endif // TEXTPROCESSINGTOOLBOX_H

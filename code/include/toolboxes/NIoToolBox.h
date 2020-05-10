#ifndef NIOTOOLBOX_H
#define NIOTOOLBOX_H

/* Toolbox for IO helperfunctions
 */

#include <vector>
#include <string>

namespace NIoToolBox {

inline std::vector<std::string> split(const std::string& s, char delimiter)
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

}
#endif // NIOTOOLBOX_H

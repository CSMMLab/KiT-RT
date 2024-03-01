#include "common/typedef.hpp"

inline double dot( Vector& vec1, Vector& vec2 ) {
    /* performs dot product of vec1 and vec2. Assumes that vec1 and vec2 are of same size.*/
    double result = 0.;
    for( unsigned i = 0; i < vec1.size(); i++ ) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

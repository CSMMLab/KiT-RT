/* set of linear algebra functions*/
#include <cmath>
#include <vector>

namespace LinAlg {

inline double dot( std::vector<double>& vec1, std::vector<double>& vec2, unsigned len ) {
    double result = 0.0;
    for( unsigned i = 0; i < len; i++ ) {
        result += vec1[i] * vec2[i];
    }
    return result;
}

inline double l2_norm( std::vector<double>& vec1, unsigned len ) {
    double result = 0.0;
    for( unsigned i = 0; i < len; i++ ) {
        result += vec1[i] * vec1[i];
    }
    return sqrt( result );
}

inline double mat_vec( std::vector<double>& mat, std::vector<double>& vec1, std::vector<double>& vec2, unsigned n_row, unsigned n_col ) {
    /* dim(mat)=n_row x n_col
       dim(vec1)=n_row
       dim(vec2)=n_col
    */
    for( unsigned i = 0; i < n_row; i++ ) {
        for( unsigned j = 0; j < n_col; j++ ) {
            vec2[j] += mat[i * n_col + j] * vec1[i];
        }
    }
}

}    // namespace LinAlg
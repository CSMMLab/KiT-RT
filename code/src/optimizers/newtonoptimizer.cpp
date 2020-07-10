
#include "optimizers/newtonoptimizer.h"
#include "settings/config.h"

Vector NewtonOptimizer::Solve( Vector u ) {

    // if we have quadratic entropy, then alpha = u;
    if( _settings->GetEntropyName() == QUADRATIC ) return u;

    return Vector( 1, 0.0 );    // dummy
}

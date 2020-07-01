#include "optimizers/optimizerbase.h"
#include "optimizers/newtonoptimizer.h"

OptimizerBase::OptimizerBase( Config* settings ) {
    _entropy  = EntropyBase::Create( settings );
    _settings = settings;
}

OptimizerBase* OptimizerBase::Create( Config* settings ) {
    switch( settings->GetOptimizerName() ) {
        case NEWTON: return new NewtonOptimizer( settings );
        // extend to other optimizers
        default: return new NewtonOptimizer( settings );
    }
}

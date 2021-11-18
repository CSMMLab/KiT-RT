#include "optimizers/optimizerbase.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/mloptimizer.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "optimizers/regularizednewtonoptimizer.hpp"

OptimizerBase::OptimizerBase( Config* settings ) {
    _entropy  = EntropyBase::Create( settings );
    _settings = settings;
}

OptimizerBase::~OptimizerBase() { delete _entropy; }

OptimizerBase* OptimizerBase::Create( Config* settings ) {
    switch( settings->GetOptimizerName() ) {
        case NEWTON: return new NewtonOptimizer( settings );
        case REGULARIZED_NEWTON: return new RegularizedNewtonOptimizer( settings );
        case ML:
            return new MLOptimizer( settings );

            // extend to other optimizers
        default: return new NewtonOptimizer( settings );
    }
}

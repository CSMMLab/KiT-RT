#include "optimizers/optimizerbase.hpp"
#include "common/config.hpp"
#include "entropies/entropybase.hpp"
#include "optimizers/neuralnetworkoptimizer.hpp"
#include "optimizers/newtonoptimizer.hpp"
#include "optimizers/partregularizednewtonoptimizer.hpp"
#include "optimizers/regularizednewtonoptimizer.hpp"
#include "toolboxes/errormessages.hpp"

OptimizerBase::OptimizerBase( Config* settings ) {
    _entropy  = EntropyBase::Create( settings );
    _settings = settings;
}

OptimizerBase::~OptimizerBase() { delete _entropy; }

OptimizerBase* OptimizerBase::Create( Config* settings ) {
    switch( settings->GetOptimizerName() ) {
        case NEWTON: return new NewtonOptimizer( settings );
        case REGULARIZED_NEWTON: return new RegularizedNewtonOptimizer( settings );
        case PART_REGULARIZED_NEWTON: return new PartRegularizedNewtonOptimizer( settings );
        case ML:
            return new NeuralNetworkOptimizer( settings );
            // extend to other optimizers
        default: ErrorMessages::Error( "Optimizer of choice not implemented.", CURRENT_FUNCTION ); return new NewtonOptimizer( settings );
    }
}

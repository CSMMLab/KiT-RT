#include "optimizers/optimizerbase.h"
#include "optimizers/mloptimizer.h"
#include "common/config.h"
#include "entropies/entropybase.h"
#include "optimizers/newtonoptimizer.h"


OptimizerBase::OptimizerBase( Config* settings ) {
    _entropy  = EntropyBase::Create( settings );
    _settings = settings;
}

OptimizerBase::~OptimizerBase() { delete _entropy; }

OptimizerBase* OptimizerBase::Create( Config* settings ) {
    switch( settings->GetOptimizerName() ) {
        case NEWTON: return new NewtonOptimizer( settings );
        case ML:
            return new MLOptimizer( settings );

            // extend to other optimizers
        default: return new NewtonOptimizer( settings );
    }
}

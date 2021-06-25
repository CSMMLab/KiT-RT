#include "entropies/entropybase.h"
#include "common/config.h"
#include "entropies/maxwellboltzmannentropy.h"
#include "entropies/quadraticentropy.h"

EntropyBase* EntropyBase::Create( Config* settings ) {

    switch( settings->GetEntropyName() ) {
        case QUADRATIC: return new QuadraticEntropy();
        case MAXWELL_BOLTZMANN: return new MaxwellBoltzmannEntropy();
        default: return new QuadraticEntropy();
    }
}

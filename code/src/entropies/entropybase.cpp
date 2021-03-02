#include "common/pch.h"
#include "entropies/maxwellboltzmannentropy.h"    // for MaxwellBoltzmannEntropy
#include "entropies/quadraticentropy.h"           // for QuadraticEntropy

EntropyBase* EntropyBase::Create( Config* settings ) {

    switch( settings->GetEntropyName() ) {
        case QUADRATIC: return new QuadraticEntropy();
        case MAXWELL_BOLTZMANN: return new MaxwellBoltzmannEntropy();
        default: return new QuadraticEntropy();
    }
}

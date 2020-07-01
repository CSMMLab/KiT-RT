#include "entropies/entropybase.h"
#include "entropies/maxwellboltzmannentropy.h"
#include "entropies/quadraticentropy.h"

EntropyBase* EntropyBase::Create( Config* settings ) {

    switch( settings->GetEntropyName() ) {
        case QUADRATIC: return new QuadraticEntropy();
        case MAXWELL_BOLZMANN: return new MaxwellBoltzmannEntropy();
        default: return new QuadraticEntropy();
    }
}

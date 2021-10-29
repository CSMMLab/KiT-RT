#include "entropies/entropybase.hpp"
#include "common/config.hpp"
#include "entropies/maxwellboltzmannentropy.hpp"
#include "entropies/quadraticentropy.hpp"

EntropyBase* EntropyBase::Create( Config* settings ) {

    switch( settings->GetEntropyName() ) {
        case QUADRATIC: return new QuadraticEntropy();
        case MAXWELL_BOLTZMANN: return new MaxwellBoltzmannEntropy();
        default: return new QuadraticEntropy();
    }
}

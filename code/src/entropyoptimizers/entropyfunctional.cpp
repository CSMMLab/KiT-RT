#include "entropyfunctional.h"

EntropyFunctionalBase::CreateEntropyFunctional() {

    return new QuadraticEntropy();
    // switch( name ) TODO: Implement other entropy functionals.
}

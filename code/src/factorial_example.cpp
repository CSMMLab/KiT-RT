#include "factorial_example.h"

unsigned int Factorial( unsigned int number ) { return number > 1 ? Factorial( number - 1 ) * number : 1; }

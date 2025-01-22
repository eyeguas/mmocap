
#include "de_d40_mm.h"

namespace gummocap{
  
  namespace algorithms{
    
    /**
 * For every variable in vars that violates the bounds, this function set the value
 * to the symmetric one with regards to the violated bound, until no violation is encountered.
 */
void makeSymmetric(double* vars, double minbound, double maxbound, int dimension) {

    for (int i = 0; i < dimension; i++) {

        while (vars[i] > maxbound || vars[i] < minbound) { //The symmetric value mignt violate the other bound

            double range = maxbound - minbound;

            if (vars[i] > maxbound) {
                double resto = fmod((vars[i] - maxbound), range);
                vars[i] = maxbound + resto;

                vars[i] = 2 * maxbound - vars[i];
            } else if (vars[i] < minbound) {

                double resto = fmod((minbound - vars[i]), range);
                vars[i] = minbound - resto;

                vars[i] = 2 * minbound - vars[i];
            }
        }
    }
}

/**
 * For every variable in vars that violates the bounds, this function set the value
 * to the symmetric one with regards to the violated bound, until no violation is encountered.
 */
void makeSymmetric(double* vars, double* x, double minbound, double maxbound, int dimension) {

    for (int i = 0; i < dimension; i++) {
	double max=x[i]+maxbound;
	double min=x[i]-minbound;
        while (vars[i] > max || vars[i] < min) { //The symmetric value mignt violate the other bound

            double range = maxbound + minbound;

            if (vars[i] > max) {
                double resto = fmod((vars[i] - max), range);
                vars[i] = max + resto;

                vars[i] = 2 * max - vars[i];
            } else if (vars[i] < min) {

                double resto = fmod((min - vars[i]), range);
                vars[i] = min - resto;

                vars[i] = 2 * min - vars[i];
            }
        }
    }
}

/**
 * This function return value if bounds are not violated. Otherwise, it return the symmetric
 * value with regards to the violated bound.
 */
/*
double makeSymmetric(double value, double minbound, double maxbound) {

    while (value > maxbound || value < minbound) { //The symmetric value mignt violate the other bound

        double range = maxbound - minbound;

        if (value > maxbound) {
            double resto = fmod((value - maxbound), range);
            value = maxbound + resto;

            value = 2 * maxbound - value;
        } else if (value < minbound) {

            double resto = fmod((minbound - value), range); //calcular el resto de esa división
            value = minbound - resto;

            value = 2 * minbound - value;
        }
    }

    return value;
}
*/

    
  }
  
}
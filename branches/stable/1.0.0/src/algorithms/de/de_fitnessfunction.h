
/* 
 * File:   fitnessfunction.h
 * Author: carlos
 *
 * Created on 7 de diciembre de 2009, 12:44
 */

#ifndef _DE_FITNESSFUNCTION_H
#define	_DE_FITNESSFUNCTION_H


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <limits>
#include <values.h>
#include <sys/timeb.h>


namespace gummocap{
  
  namespace algorithms{

    typedef long double tFitness;

    using namespace std;

    
/**
 * This abstract class defines the concept of a fitness function.
 * It incorporates a member to count the number of times it is invoked.
 */
class FitnessFunction {
protected:
    int numFEs;

public:
    FitnessFunction() { numFEs = 0; };
    virtual long double fitness(const double *vars) = 0;
    virtual long double fitness(const double *vars, int dim) = 0;
    virtual int getDim() = 0;
    virtual ostringstream* getName() = 0;
    int getNumFEs(){ return this->numFEs; }

    /**
     * This function defines which is the best fitness value given two as parameters.
     * It is useful to design maximising and minimisin problems.
     */
    virtual long double compare(long double f1, long double f2) = 0;
};


/**
 * This class implements the concept of a Fitness Function that obtains the fitness value
 * of a given solution from a c function given as parameter.
 */
class WrapperFunction : public FitnessFunction{
    static tFitness (*ff)(int, const double*);
    static int dim;
    ostringstream name;

public:

    WrapperFunction(tFitness(*ff)(int, const double*), int dim, char *name=NULL) : FitnessFunction() {
        this->ff = ff;
        this->dim = dim;
        this->name.str("");

        if (name != NULL){
            this->name << name;
            this->name << "_";
            this->name << dim;
        }
    }

    void setFF(tFitness(*ff)(int, const double*), int dim, char *name = NULL) {
        this->ff = ff;
        this->dim = dim;
        this->name.str("");

        if (name != NULL){
            this->name << name;
            this->name << "_";
            this->name << dim;
        }
    }

    static long double newff(const double* vars) {
        return ff(dim, vars);
    }
    
    int getDim(){
        return dim;
    }

    long double fitness(const double *vars, int dimension){
        this->numFEs++;
        return (WrapperFunction::newff((double *)vars));
    }

    long double fitness(const double *vars){
        this->numFEs++;
        return (WrapperFunction::newff((double *)vars));
    }

    ostringstream* getName(){
        return &(this->name);
    }

    /**
     * In this class, the functions are to be minimised.
     */
    long double compare(long double f1, long double f2){

        return (f2 - f1);
    }
};

}
  
}

#endif	/* _FITNESSFUNCTION_H */


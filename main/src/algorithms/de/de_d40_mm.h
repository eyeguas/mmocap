
/* 
 * File:   de_d40_mm.h
 * Author: carlos
 *
 * Created on 7 de julio de 2010, 12:41
 */

#ifndef _DE_D40_MM_H
#define	_DE_D40_MM_H

#include "de_fitnessfunction.h"
#include <string.h>
#include <math.h>
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <algorithm>
#include <string>
#include <values.h>
#include <omp.h>

namespace gummocap{
  
  namespace algorithms{

#define mymin(x,y) (((x)<(y))?(x):(y))
#define mymax(x,y) (((x)>(y))?(x):(y))

using namespace std;

/**
 * For every variable in vars that violates the bounds, this function set the value
 * to the symmetric one with regards to the violated bound, until no violation is encountered.
 */
void makeSymmetric(double* vars, double minbound, double maxbound, int dimension);

/**
 * For every variable in vars that violates the bounds, this function set the value
 * to the symmetric one with regards to the violated bound, until no violation is encountered.
 */
void makeSymmetric(double* vars, double* x, double minbound, double maxbound, int dimension);



/**
 * This class implement a solution to a continuous optimisation problem as a vector of double values.
 * It defines classical constructors, destructor, copy operations, and accessing functions
 */
class Solution_v1 {
    double *vectorVars;
    long double fitnessValue;
    int dimension;

public:

    Solution_v1(int dimension) {
        this->vectorVars = new double[dimension];
        this->dimension = dimension;
    }

    Solution_v1(Solution_v1 &sol) {
        this->vectorVars = new double[sol.dimension];
        this->dimension = sol.dimension;
        memcpy(this->vectorVars, sol.vectorVars, sizeof (double) * sol.dimension);
        this->fitnessValue = sol.fitnessValue;
    }

    Solution_v1 & operator =(Solution_v1 &sol) {
        memcpy(this->vectorVars, sol.vectorVars, sizeof (double) * dimension);
        this->fitnessValue = sol.fitnessValue;
        return (*this);
    }

    virtual Solution_v1* copy(Solution_v1 *sol) {
        memcpy(this->vectorVars, sol->vectorVars, sizeof (double) * dimension);
        this->fitnessValue = sol->fitnessValue;
        return this;
    }

    ~Solution_v1() {
        delete [] vectorVars;
    }

    long double& fitness() {
        return fitnessValue;
    }

    double& operator [](int i) {
        return vectorVars[i];
    }

    double*& vars() {
        return vectorVars;
    }
};


/**
 * This class extend previous one with the parameters to select the leading and correcting vectors
 * according to assortative mating schemes
 *
 * It defines classical constructors, destructor, copy operations and accesing functions.
 */
class Solution_v2 : public Solution_v1 {
    int preferredLeadingNass;
    int preferredCorrectingNass;

public:

    Solution_v2(int dimension) : Solution_v1(dimension) {
        preferredLeadingNass = 0;
        preferredCorrectingNass = 0;
    }

    Solution_v2(Solution_v2& sol) : Solution_v1(sol) {
        preferredLeadingNass = sol.preferredLeadingNass;
        preferredCorrectingNass = sol.preferredCorrectingNass;
    }

    ~Solution_v2() {
    };

    Solution_v2 & operator=(Solution_v2& sol) {
        ((Solution_v1 *)this)->operator =(sol);
        this->preferredLeadingNass = sol.preferredLeadingNass;
        this->preferredCorrectingNass = sol.preferredCorrectingNass;
	return *this;
    }

    Solution_v2* copy(Solution_v2* sol) {
        this->Solution_v1::copy(sol);
        this->preferredLeadingNass = sol->preferredLeadingNass;
        this->preferredCorrectingNass = sol->preferredCorrectingNass;
	return this;
    }

    int& prefLeadingNass() {
        return preferredLeadingNass;
    }

    int & prefCorrectingNass() {
        return preferredCorrectingNass;
    }
};


/**
 * This class defines the algorithm DE-D^{40}+M^m
 */
class DE_D40_MM {
    const gsl_rng_type *T;
    gsl_rng **r;

    vector< Solution_v2* > *pop;
//    int indexBest;
    
    /************************************
     * DOBLE POBLACION PARA PARALELIZAR
     * *******************************/
    vector< Solution_v2* > *newPop;
    Solution_v2 **newSols;

    //int indexBest; //Quitado para parelelizar

    gummocap::algorithms::FitnessFunction *ff;
    vector<FitnessFunction*> ffs;
    double maxBound;
    double minBound;

    int evals;
    int dimension;
    
    /** DE parameters*/
    int popSize;
    double m_CR;
    double m_F;

//    long double bestFitness;


    /**Differential role variables*/
    int numLeadingSols;
    int numReceptorsSols;
    int numPlacingSols;
    int numCorrectingSols;

public:

    /**
     * Constructor
     * @param seed is the initial seed for the gsl random number generator
     */
    DE_D40_MM(gummocap::algorithms::FitnessFunction *ff, int dimension, double minbound, double maxbound, long int seed) {

        //Parameters initialisation
        popSize = 60;
        m_CR = 0.9;	// Crossover constant
        m_F = 0.5;      // Controls the amplification of the differential variation (Mutation)

        numReceptorsSols = popSize;
        numLeadingSols = 40;
        numPlacingSols = 40;
        numCorrectingSols = 40;

        this->pop = new vector< Solution_v2* >();
	this->newPop = new vector< Solution_v2* >();
	
	
		
        gsl_rng_env_setup();
        T = gsl_rng_default;

        r = new gsl_rng*[numReceptorsSols];

        for (int i = 0; i < numReceptorsSols; i++) {
            r[i] = gsl_rng_alloc(T);
        }

        gsl_rng_set(r[0], seed);

        for (int i = 1; i < numReceptorsSols; i++) {
            gsl_rng_set(r[i], gsl_rng_uniform_int(r[i - 1], MAXINT));
        }
        
	this->ff = ff;
        evals = 0;
        this->dimension = dimension;
        this->maxBound = maxbound;
        this->minBound = minbound;
//        this->indexBest = 0;

        //Initialising population
        this->randomFill(popSize, minbound, maxbound);
//        bestFitness = pop->at(indexBest)->fitness();
        sort(pop->begin(), pop->end(), DE_D40_MM::myCmpFunc);
    }
    /**
     * Constructor
     * @param seed is the initial seed for the gsl random number generator
     */
    DE_D40_MM(FitnessFunction *ff, int dimension, double* x, double minbound, double maxbound, long int seed) {

        //Parameters initialisation
        popSize = 60;
        m_CR = 0.9;
        m_F = 0.5;

        numReceptorsSols = popSize;
        numLeadingSols = 40;
        numPlacingSols = 40;
        numCorrectingSols = 40;

        this->pop = new vector< Solution_v2* >();
	this->newPop = new vector< Solution_v2* >();

	gsl_rng_env_setup();
        T = gsl_rng_default;

        r = new gsl_rng*[numReceptorsSols];

        for (int i = 0; i < numReceptorsSols; i++) {
            r[i] = gsl_rng_alloc(T);
        }

        gsl_rng_set(r[0], seed);

        for (int i = 1; i < numReceptorsSols; i++) {
            gsl_rng_set(r[i], gsl_rng_uniform_int(r[i - 1], MAXINT));
        }
	
	this->ff = ff;
        evals = 0;
        this->dimension = dimension;
        this->maxBound = maxbound;
        this->minBound = minbound;
//        this->indexBest = 0;

        //Initialising population
	this->randomFill(popSize, x, minbound, maxbound);
//        bestFitness = pop->at(indexBest)->fitness();
        sort(pop->begin(), pop->end(), DE_D40_MM::myCmpFunc);
    }
    /**
     * Constructor
     * @param seed is the initial seed for the gsl random number generator
     */
    DE_D40_MM(vector<FitnessFunction* > ffs, int dimension, double* x, unsigned int population_size, unsigned int nLPC, double cr, double f, double minbound, double maxbound, long int seed) {

        //Parameters initialisation
        popSize = population_size;
        m_CR = cr;
        m_F = f;

        numReceptorsSols = popSize;
        numLeadingSols = nLPC;
        numPlacingSols = nLPC;
        numCorrectingSols = nLPC;

        this->pop = new vector< Solution_v2* >();
//	this->newPop = new vector< Solution_v2* >();
	
	newSols = new Solution_v2*[this->numReceptorsSols];
	for(int i=0;i<this->numReceptorsSols;i++)
	newSols[i] = new Solution_v2(dimension);

	
        gsl_rng_env_setup();
        T = gsl_rng_default;

        r = new gsl_rng*[numReceptorsSols];

        for (int i = 0; i < numReceptorsSols; i++) {
            r[i] = gsl_rng_alloc(T);
        }

        gsl_rng_set(r[0], seed);

        for (int i = 1; i < numReceptorsSols; i++) {
            gsl_rng_set(r[i], gsl_rng_uniform_int(r[i - 1], MAXINT));
        }
	
	this->ffs = ffs;
	this->ff=ffs[0];
        evals = 0;
        this->dimension = dimension;
        this->maxBound = maxbound;
        this->minBound = minbound;
//        this->indexBest = 0;

        //Initialising population
	this->randomFill(popSize, x, ffs, minbound, maxbound);
//        bestFitness = pop->at(indexBest)->fitness();
        sort(pop->begin(), pop->end(), DE_D40_MM::myCmpFunc);
    }

    /**
     * Destructor
     */
    ~DE_D40_MM() {

        vector< Solution_v2* >::iterator it = pop->begin();
        vector< Solution_v2* >::iterator endIt = pop->end();

        while (it != endIt) {
            delete (*it);
            it++;
        }

        pop->clear();
        delete pop;

        for (int i = 0; i < numReceptorsSols; i++)
            gsl_rng_free(r[i]);

        delete [] r;
    }


    /**
     * This function generate num random solutions and insert them in member pop.
     * Solutions are evaluated by the object ff, which is supposed to implement the fitness function.
     */
    void randomFill(int num, double minBound, double maxBound) {

        double range = maxBound - minBound;

        for (int i = 0; i < num; i++) {
            Solution_v2 *newSol = new Solution_v2(dimension);

            for (int j = 0; j < dimension; j++)
                (*newSol)[j] = gsl_rng_uniform(r[i]) * range + minBound;

            pop->push_back(newSol);
            newSol->fitness() = ff->fitness(newSol->vars());

            //            if (pop->at(indexBest)->fitness() > newSol->fitness()) {  QUITADO PARA PARELELIZAR
            //                indexBest = pop->size() - 1;
            //            }

            evals++;
        }
    }

     /**
     * This function generate num random solutions and insert them in member pop.
     * Solutions are evaluated by the object ff, which is supposed to implement the fitness function.
     */
    void randomFill(int num, double* x, double minBound, double maxBound) {

	// [x-minBound,x+maxBound]
        double range = maxBound + minBound;

        for (int i = 0; i < num; i++) {
            Solution_v2 *newSol = new Solution_v2(dimension);

            for (int j = 0; j < dimension; j++)
                (*newSol)[j] = gsl_rng_uniform(r[i]) * range + (x[j] - minBound);

            pop->push_back(newSol);
            newSol->fitness() = ff->fitness(newSol->vars());

            //if (pop->at(indexBest)->fitness() > newSol->fitness()) {
            //    indexBest = pop->size() - 1;
            //}

            evals++;
        }
      }

/**
     * This function generate num random solutions and insert them in member pop.
     * Solutions are evaluated by the object ff, which is supposed to implement the fitness function.
     */
    void randomFill(int num, double* x, vector<FitnessFunction* > ffs, double minBound, double maxBound) {

	// [x-minBound,x+maxBound]
        double range = maxBound + minBound;

	for (int i = 0; i < num; i++) {
	    Solution_v2 *newSol = new Solution_v2(dimension);

            for (int j = 0; j < dimension; j++)
                (*newSol)[j] = gsl_rng_uniform(r[i]) * range + (x[j] - minBound);

            pop->push_back(newSol);
	    evals++;
        }
        #pragma omp parallel for
        for (int i = 0; i < num; i++) {
	    int threadId=omp_get_thread_num();
            pop->at(i)->fitness() = ffs[threadId]->fitness(pop->at(i)->vars());
        }
        /*
        indexBest=0;
        for (int i = 0; i < num; i++) {
	  if (pop->at(indexBest)->fitness() > pop->at(i)->fitness()) {
                indexBest = i;
            }
	}
	*/
      }

    /**
     * This function creates a permutation of integer numbers (from 0 to dim - 1)
     */
    void createPermutation(int *perm, int dim, int generator) {

        for (int i = 0; i < dim; i++)
            perm[i] = i;

        for (int i = dim - 1; i >= 0; i--) {
            int randI = gsl_rng_uniform_int(r[generator], dim);
            int aux = perm[i];
            perm[i] = perm[randI];
            perm[randI] = aux;
        }
    }

    /**
     * This function implements the linear combination of vectors along with the exponential crossover operation
     *
     * @param s1 is the placing vector
     * @param s2 is the leading vector
     * @param s3 is the correcting vector
     * @param origin is the receiving vector
     * @param offspring is the new trial vector
     */
    void crossover(double *s1, double *s2, double *s3, double *origin, double *offspring, int generator) {

        int *perm = new int[dimension];
        memcpy(offspring, origin, sizeof (double) * dimension);

        int n = gsl_rng_uniform_int(r[generator], dimension);

        for (int i = 0; gsl_rng_uniform(r[generator]) < m_CR && (i < dimension); i++) {
            int k = n;
            offspring[k] = s1[k] + m_F * (s2[k] - s3[k]);
            n = (n + 1) % dimension;
        }

        delete [] perm;
    }

    /**
     * This function compue a distance metric between two double vectors as the sum of the absolute differences
     *
     * @param stopSum: this parameter can be used to stop computing the distance if it is going to be greater than the specified value. That may be useful to reduce computing time.
     */
    double computeDistance(double *s1, double *s2, int dimension, double stopSum = MAXDOUBLE) {

        double sum = 0.;

        for (int i = 0; i < dimension && sum < stopSum; i++) {
            sum += fabs(s1[i] - s2[i]);
        }

        return sum;
    }

    /**
     * This function implements the positive assortative mating selection scheme.
     *
     * @param posibleIndexes is a set of available indexes for the selection of the mating vector.
     * @param nass is the number of random vector considered
     * @param firstParent is the vector with which the distance will be computied. The closest vector is to be selected.
     * @return It returns the index, from a set of available indexes, of the selected vector
     */
    vector<int>::iterator positiveAssMating(vector<int> *posibleIndexes, int nass, double *firstParent, int generator) {

        vector< vector<int>::iterator > *randomIndexes = new vector< vector<int>::iterator>;
        sampleNParents(randomIndexes, posibleIndexes, nass, generator);

        double minDist = computeDistance(pop->at(*(randomIndexes->at(0)))->vars(), firstParent, dimension);
        int indexMinDist = 0;

        for (int i = 1; i < nass; i++) {
            double newDist = computeDistance(pop->at(*(randomIndexes->at(i)))->vars(), firstParent, dimension, minDist);

            if (newDist < minDist) {
                minDist = newDist;
                indexMinDist = i;
            }
        }

        vector<int>::iterator selected = randomIndexes->at(indexMinDist);
        randomIndexes->clear();
        delete randomIndexes;
        return selected;
    }

    /**
     * This function implements the negative assortative mating selection scheme.
     *
     * @param posibleIndexes is a set of available indexes for the selection of the mating vector.
     * @param nass is the number of random vector considered
     * @param firstParent is the vector with which the distance will be computied. The furthest vector is to be selected.
     * @return It returns the index, from a set of available indexes, of the selected vector
     */
    vector<int>::iterator negativeAssMating(vector<int> *posibleIndexes, int nass, double *firstParent, int generator) {

        vector< vector<int>::iterator > *randomIndexes = new vector< vector<int>::iterator>;
        sampleNParents(randomIndexes, posibleIndexes, nass, generator);

        double maxDist = computeDistance(pop->at(*(randomIndexes->at(0)))->vars(), firstParent, dimension);
        int indexMinDist = 0;

        for (int i = 1; i < nass; i++) {
            double newDist = computeDistance(pop->at(*(randomIndexes->at(i)))->vars(), firstParent, dimension);

            if (newDist > maxDist) {
                maxDist = newDist;
                indexMinDist = i;
            }
        }

        vector<int>::iterator selected = randomIndexes->at(indexMinDist);
        randomIndexes->clear();
        delete randomIndexes;
        return selected;
    }

    /**
     * This function selects n vectors by means of positive assortative mating, negative assortative mating, or random mating
     * according to the value of parameter nass.
     * 
     * @param parents selected parents are stored in that vector.
     * @param posibleIndexes is a set of available indexes to be considered in the selection process
     * @param firstParent is a vector with which distance values are measured.
     * @param n is the number of vector to be selected. It comes from previous versions of the code. The actual algorithm always use n=1
     * @param nass is the parameter of the selection operator. If it is negative, then, negative assortative mating is applied with the absolute value for nass. If it is positive, then, positive assortative mating is applied. If it is 0, random mating is applied. Notice that negative or assortative mating with nass=1 is equal to random mating.
     */
    void sampleNParentsAss(vector< double *> *parents, vector< int > *posibleIndexes, double *firstParent, int n, int nass, int generator) {

        if (nass < 0) {

            for (int i = 0; i < n; i++) {
                vector<int>::iterator selected = negativeAssMating(posibleIndexes, -nass, firstParent, generator);
                parents->push_back(pop->at(*selected)->vars());
                posibleIndexes->erase(selected);
            }
        } else if (nass > 0) {

            for (int i = 0; i < n; i++) {
                vector<int>::iterator selected = positiveAssMating(posibleIndexes, nass, firstParent, generator);
                parents->push_back(pop->at(*selected)->vars());
                posibleIndexes->erase(selected);
            }

        } else {

            for (int i = 0; i < n; i++) {
                vector<int>::iterator selected = positiveAssMating(posibleIndexes, 1, firstParent, generator);
                parents->push_back(pop->at(*selected)->vars());
                posibleIndexes->erase(selected);
            }
        }
    }

    /**
     * This function sample n int values (indexes of the population) randomly from a set of posible values
     *
     * @param indexes The function append the selected indexes
     * @param posibleIndexes is the set of available indexes to consider for the selection.
     */
    void sampleNParents(vector< vector<int>::iterator > *indexes, vector<int> *posibleIndexes, int n, int generator) {

        int size = posibleIndexes->size();
        int *perm = new int[size];
        createPermutation(perm, size, generator);

        for (int i = 0; i < n; i++) {

            indexes->push_back(posibleIndexes->begin() + perm[i]);
        }

        delete [] perm;
    }

    /**
     * Dummy function for the STL sort function. It says that solution a should be before solution b if its fitness is smaller.
     */
    static bool myCmpFunc(Solution_v2* a, Solution_v2* b) {

        return (a->fitness() < b->fitness());
    }

    /**
     * This function implements the iteration of the algorithm
     */
    void iterate(long int maxEvals) {

        vector< double *> *parentsLeading = new vector< double *>();
        vector< double *> *parentsCorrecting = new vector< double *>();
        double *offspring = new double[dimension];
        //indexBest = popSize - 1;

        //For each receiving vector in the population
// #pragma omp parallel for	
        for (int i = this->numReceptorsSols - 1; i >= 0 && evals < maxEvals; i--) {

            //Initialize posible indexes for the selection processes
            vector<int> posibleIndexesLeading;
            vector<int> posibleIndexesCorrecting;

            for (int j = 0; j < numLeadingSols; j++)
                posibleIndexesLeading.push_back(j);

            for (int j = popSize - 1, jaux = 0; j >= 0 && jaux < numCorrectingSols; j--, jaux++)
                posibleIndexesCorrecting.push_back(j);


            //Randomly select the placing solution
            int placingSelected = gsl_rng_uniform_int(r[i], numPlacingSols);
            double *current_firstParent = pop->at(placingSelected)->vars();


            //Retrieve the parameters to select the leading and correcting vectors
            int leadingNass = pop->at(placingSelected)->prefLeadingNass();
            int correctingNass = pop->at(placingSelected)->prefCorrectingNass();


            //Randomly perturb previous parameters
            switch (gsl_rng_uniform_int(r[i], 3) - 1) {
                case 1: leadingNass = mymin(4, leadingNass + 1);
                    break;
                case -1: leadingNass = mymax(-4, leadingNass - 1);
                    break;
            }

            switch (gsl_rng_uniform_int(r[i], 3) - 1) {
                case 1: correctingNass = mymin(4, correctingNass + 1);
                    break;
                case -1: correctingNass = mymax(-4, correctingNass - 1);
                    break;
            }


            //Remove the index of the placing vector before selecting the leading vector
            if (placingSelected < numLeadingSols)
                posibleIndexesLeading.erase(posibleIndexesLeading.begin() + placingSelected);

            //Remove the index of the placing vector before selecting the correcting vector.
            if (placingSelected > (popSize - numCorrectingSols))
                posibleIndexesCorrecting.erase(posibleIndexesCorrecting.begin() + (placingSelected - (popSize - numCorrectingSols)));

            //Select the leading vector. Care is paid if the number of available indexes is inferior to the nass parameter value.
            if (leadingNass < 0)
                sampleNParentsAss(parentsLeading, &posibleIndexesLeading, current_firstParent, 1, mymax((double) leadingNass, -1 * ((int) (posibleIndexesLeading.size()) - 1)),i);
            else
                sampleNParentsAss(parentsLeading, &posibleIndexesLeading, current_firstParent, 1, mymin((double) leadingNass, posibleIndexesLeading.size() - 1),i);

            //Select the correcting vector. Care is paid if the number of available indexes is inferior to the nass parameter value.
            if (correctingNass < 0)
                sampleNParentsAss(parentsCorrecting, &posibleIndexesCorrecting, current_firstParent, 1, mymax((double) correctingNass, -1 * ((int) (posibleIndexesCorrecting.size()) - 1)),i);
            else
                sampleNParentsAss(parentsCorrecting, &posibleIndexesCorrecting, current_firstParent, 1, mymin((double) correctingNass, posibleIndexesCorrecting.size() - 1),i);

            //Apply mutation and crossover operation
            crossover(current_firstParent, parentsLeading->at(0), parentsCorrecting->at(0), pop->at(i)->vars(), offspring, i);

            //make the solution to be within the search domain.
            makeSymmetric(offspring, minBound, maxBound, dimension);

            //clear auxiliary estructures
            posibleIndexesLeading.clear();
            posibleIndexesCorrecting.clear();

            //Evaluate the new trial solution.
            long double newFitness = ff->fitness(offspring);
            evals++;
            
            //clear auxiliary estructures
            parentsLeading->clear();
            parentsCorrecting->clear();

            //Apply the selection operation, i.e., replace the target vector if the new trial solution is better.
            if (newFitness < pop->at(i)->fitness()) {
		/*
                pop->at(i)->fitness() = newFitness;
                memcpy(pop->at(i)->vars(), offspring, sizeof (double) * dimension);
                pop->at(i)->prefLeadingNass() = leadingNass;
                pop->at(i)->prefCorrectingNass() = correctingNass;

                if (newFitness < pop->at(indexBest)->fitness()) {
                    indexBest = i;
                }

                if (newFitness < bestFitness) {
                    bestFitness = newFitness;
                }
                */
                Solution_v2 *newSol = new Solution_v2(dimension);
                memcpy(newSol->vars(), offspring, sizeof (double) * dimension);
                newSol->prefCorrectingNass() = correctingNass;
                newSol->prefLeadingNass() = leadingNass;
                newPop->push_back(newSol);
            } 
            else {
		Solution_v2 *newSol = new Solution_v2(*(pop->at(i)));
                newPop->push_back(newSol);
            }
        }

	for (int i = this->numReceptorsSols - 1; i >= 0 && evals < maxEvals; i--) {
            delete pop->at(i);
            pop->at(i) = newPop->at(i);
        }

        newPop->clear();

        //Sort the population
        sort(pop->begin(), pop->begin() + numReceptorsSols, DE_D40_MM::myCmpFunc);

        //Free auxiliary estructures.
        delete parentsLeading;
        delete parentsCorrecting;
        delete [] offspring;
    }

/**
     * This function implements the iteration of the algorithm
     */
/*
    void iterate(double* x, long int maxEvals) {

        vector< double *> *parentsLeading = new vector< double *>();
        vector< double *> *parentsCorrecting = new vector< double *>();
        double *offspring = new double[dimension];
       // indexBest = popSize - 1;

        //For each receiving vector in the population
        for (int i = this->numReceptorsSols - 1; i >= 0 && evals < maxEvals; i--) {

            //Initialize posible indexes for the selection processes
            vector<int> posibleIndexesLeading;
            vector<int> posibleIndexesCorrecting;

            for (int j = 0; j < numLeadingSols; j++)
                posibleIndexesLeading.push_back(j);

            for (int j = popSize - 1, jaux = 0; j >= 0 && jaux < numCorrectingSols; j--, jaux++)
                posibleIndexesCorrecting.push_back(j);


            //Randomly select the placing solution
            int placingSelected = gsl_rng_uniform_int(r[i], numPlacingSols);
            double *current_firstParent = pop->at(placingSelected)->vars();


            //Retrieve the parameters to select the leading and correcting vectors
            int leadingNass = pop->at(placingSelected)->prefLeadingNass();
            int correctingNass = pop->at(placingSelected)->prefCorrectingNass();


            //Randomly perturb previous parameters
            switch (gsl_rng_uniform_int(r[i], 3) - 1) {
                case 1: leadingNass = mymin(4, leadingNass + 1);
                    break;
                case -1: leadingNass = mymax(-4, leadingNass - 1);
                    break;
            }

            switch (gsl_rng_uniform_int(r[i], 3) - 1) {
                case 1: correctingNass = mymin(4, correctingNass + 1);
                    break;
                case -1: correctingNass = mymax(-4, correctingNass - 1);
                    break;
            }


            //Remove the index of the placing vector before selecting the leading vector
            if (placingSelected < numLeadingSols)
                posibleIndexesLeading.erase(posibleIndexesLeading.begin() + placingSelected);

            //Remove the index of the placing vector before selecting the correcting vector.
            if (placingSelected > (popSize - numCorrectingSols))
                posibleIndexesCorrecting.erase(posibleIndexesCorrecting.begin() + (placingSelected - (popSize - numCorrectingSols)));

            //Select the leading vector. Care is paid if the number of available indexes is inferior to the nass parameter value.
            if (leadingNass < 0)
                sampleNParentsAss(parentsLeading, &posibleIndexesLeading, current_firstParent, 1, mymax((double) leadingNass, -1 * ((int) (posibleIndexesLeading.size()) - 1)),i);
            else
                sampleNParentsAss(parentsLeading, &posibleIndexesLeading, current_firstParent, 1, mymin((double) leadingNass, posibleIndexesLeading.size() - 1),i);

            //Select the correcting vector. Care is paid if the number of available indexes is inferior to the nass parameter value.
            if (correctingNass < 0)
                sampleNParentsAss(parentsCorrecting, &posibleIndexesCorrecting, current_firstParent, 1, mymax((double) correctingNass, -1 * ((int) (posibleIndexesCorrecting.size()) - 1)),i);
            else
                sampleNParentsAss(parentsCorrecting, &posibleIndexesCorrecting, current_firstParent, 1, mymin((double) correctingNass, posibleIndexesCorrecting.size() - 1),i);

            //Apply mutation and crossover operation
            crossover(current_firstParent, parentsLeading->at(0), parentsCorrecting->at(0), pop->at(i)->vars(), offspring,i);

            //make the solution to be within the search domain.
            makeSymmetric(offspring, x, minBound, maxBound, dimension);

            //clear auxiliary estructures
            posibleIndexesLeading.clear();
            posibleIndexesCorrecting.clear();

            //Evaluate the new trial solution.
            long double newFitness = ff->fitness(offspring);
	    
            evals++;
            
            //clear auxiliary estructures
            parentsLeading->clear();
            parentsCorrecting->clear();

            //Apply the selection operation, i.e., replace the target vector if the new trial solution is better.
	    if (newFitness < pop->at(i)->fitness()) {
               Solution_v2 *newSol = new Solution_v2(dimension);
                memcpy(newSol->vars(), offspring, sizeof (double) * dimension);
		newSol->fitness() = newFitness;
                newSol->prefCorrectingNass() = correctingNass;
                newSol->prefLeadingNass() = leadingNass;
                newPop->push_back(newSol);
            }
            else
	    {
		Solution_v2 *newSol = new Solution_v2(*(pop->at(i)));
                newPop->push_back(newSol);
		
	    }
        }
        

        for (int i = this->numReceptorsSols - 1; i >= 0 && evals < maxEvals; i--) {
            delete pop->at(i);
            pop->at(i) = newPop->at(i);
        }

        newPop->clear();


        //Sort the population
        sort(pop->begin(), pop->begin() + numReceptorsSols, DE_D40_MM::myCmpFunc);

        //Free auxiliary estructures.
        delete parentsLeading;
        delete parentsCorrecting;
        delete [] offspring;
    }
    */
    /**
     * This function implements the iteration of the algorithm
     */
    void iterate(double* x, long int maxEvals) {

       // indexBest = popSize - 1;

    

	 if(evals<maxEvals){
        //For each receiving vector in the population
                
	    
        #pragma omp parallel for
        //for (int i = this->numReceptorsSols - 1; i >= 0; i--) {
	for (int i = 0; i < this->numReceptorsSols; i++) {

	    vector< double *> *parentsLeading = new vector< double *>();
	    vector< double *> *parentsCorrecting = new vector< double *>();
	    double *offspring = new double[dimension];
	    double *current_firstParent= new double[dimension];

            //Initialize posible indexes for the selection processes
            vector<int> posibleIndexesLeading;
            vector<int> posibleIndexesCorrecting;

            for (int j = 0; j < numLeadingSols; j++)
                posibleIndexesLeading.push_back(j);

            for (int j = popSize - 1, jaux = 0; j >= 0 && jaux < numCorrectingSols; j--, jaux++)
                posibleIndexesCorrecting.push_back(j);


            //Randomly select the placing solution
            int placingSelected = gsl_rng_uniform_int(r[i], numPlacingSols);
            memcpy(current_firstParent, pop->at(placingSelected)->vars(), sizeof (double) * dimension);
		

            //Retrieve the parameters to select the leading and correcting vectors
            int leadingNass = pop->at(placingSelected)->prefLeadingNass();
            int correctingNass = pop->at(placingSelected)->prefCorrectingNass();


            //Randomly perturb previous parameters
            switch (gsl_rng_uniform_int(r[i], 3) - 1) {
                case 1: leadingNass = mymin(4, leadingNass + 1);
                    break;
                case -1: leadingNass = mymax(-4, leadingNass - 1);
                    break;
            }

            switch (gsl_rng_uniform_int(r[i], 3) - 1) {
                case 1: correctingNass = mymin(4, correctingNass + 1);
                    break;
                case -1: correctingNass = mymax(-4, correctingNass - 1);
                    break;
            }

            //Remove the index of the placing vector before selecting the leading vector
            if (placingSelected < numLeadingSols)
                posibleIndexesLeading.erase(posibleIndexesLeading.begin() + placingSelected);

            //Remove the index of the placing vector before selecting the correcting vector.
            if (placingSelected > (popSize - numCorrectingSols))
                posibleIndexesCorrecting.erase(posibleIndexesCorrecting.begin() + (placingSelected - (popSize - numCorrectingSols)));

	    //Select the leading vector. Care is paid if the number of available indexes is inferior to the nass parameter value.
            if (leadingNass < 0)
                sampleNParentsAss(parentsLeading, &posibleIndexesLeading, current_firstParent, 1, mymax((double) leadingNass, -1 * ((int) (posibleIndexesLeading.size()) - 1)),i);
            else
                sampleNParentsAss(parentsLeading, &posibleIndexesLeading, current_firstParent, 1, mymin((double) leadingNass, posibleIndexesLeading.size() - 1),i);

            //Select the correcting vector. Care is paid if the number of available indexes is inferior to the nass parameter value.
            if (correctingNass < 0)
                sampleNParentsAss(parentsCorrecting, &posibleIndexesCorrecting, current_firstParent, 1, mymax((double) correctingNass, -1 * ((int) (posibleIndexesCorrecting.size()) - 1)),i);
            else
                sampleNParentsAss(parentsCorrecting, &posibleIndexesCorrecting, current_firstParent, 1, mymin((double) correctingNass, posibleIndexesCorrecting.size() - 1),i);

            //Apply mutation and crossover operation
            crossover(current_firstParent, parentsLeading->at(0), parentsCorrecting->at(0), pop->at(i)->vars(), offspring,i);

            //make the solution to be within the search domain.
            makeSymmetric(offspring, x, minBound, maxBound, dimension);

            //clear auxiliary estructures
            //posibleIndexesLeading.clear();
            //posibleIndexesCorrecting.clear();

            //Evaluate the new trial solution.
	    int threadId=omp_get_thread_num();
            long double newFitness = ffs[threadId]->fitness(offspring);
	    
            evals++;
            
            //clear auxiliary estructures
//            parentsLeading->clear();
//            parentsCorrecting->clear();

            //Apply the selection operation, i.e., replace the target vector if the new trial solution is better.
	    if (newFitness < pop->at(i)->fitness()) {
                 Solution_v2 *newSol = new Solution_v2(dimension);
                memcpy(newSol->vars(), offspring, sizeof (double) * dimension);
		newSol->fitness() = newFitness;
                newSol->prefCorrectingNass() = correctingNass;
                newSol->prefLeadingNass() = leadingNass;
		*(newSols[i])=*newSol;
            }
            else
	    {
		Solution_v2 *newSol = new Solution_v2(*(pop->at(i)));
		*(newSols[i])=*newSol;
	    }
	    delete parentsLeading;
	    delete parentsCorrecting;
	    
	    delete [] offspring;
	    delete [] current_firstParent;
        }
	}

      if (evals < maxEvals){
	#pragma omp parallel for
        //for (int i = this->numReceptorsSols - 1; i >= 0; i--) {
	  for (int i = 0; i < this->numReceptorsSols; i++) {
            *(pop->at(i)) = *(newSols[i]);
        }
      }


        //Sort the population
        sort(pop->begin(), pop->begin() + numReceptorsSols, DE_D40_MM::myCmpFunc);

        //Free auxiliary estructures.
    }

    /**
     * This function runs the algorithm until a maximum number of evaluations is reached.
     */
    long double run(double* x, double* xout, int maxEvals) {


        while (evals < maxEvals) {
            iterate(x,maxEvals);
        }
        
        double bestFitness = MAXDOUBLE;
	unsigned int indexBest=this->pop->size() - 1;

        for (int i = this->pop->size() - 1; i >= 0; i--) {
            if (bestFitness > pop->at(i)->fitness()){
                bestFitness = pop->at(i)->fitness();
		indexBest=i;
	    }
        }
        
        for(int i=0;i<dimension;i++){
	  xout[i]=pop->at(indexBest)->vars()[i];
	}

        return bestFitness;
    }


    /**
     * This function runs the algorithm until a maximum number of evaluations is reached.
     */
    long double run(int maxEvals) {


        while (evals < maxEvals) {
            iterate(maxEvals);
        }
        
        double bestFitness = MAXDOUBLE;

        for (int i = this->pop->size() - 1; i >= 0; i--) {
            if (bestFitness > pop->at(i)->fitness())
                bestFitness = pop->at(i)->fitness();
        }

        return bestFitness;
    }
};


  }
  
}

#endif	/* _DE_D40_MM_H */


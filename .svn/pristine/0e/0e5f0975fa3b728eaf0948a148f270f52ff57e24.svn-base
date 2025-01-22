/***************************************************************************
 *   Copyright (C) 2004 by salinas                                         *
 *   salinas@decsai.ugr.es                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef gummocap_AnnealedPFParallel_H
#define gummocap_AnnealedPFParallel_H
#include "particle.h"
#include <cstdlib>
#include <omp.h>
#include <gu/guexception.h>
#include <gu/gulinearscale.h>
#include <gu/gumatrix.h> 
namespace gummocap {
namespace algorithms { 
/**
\brief Implmentation of the Annealed Particle filter:"Deutscher J, Reid I. Articulated Body Motion Capture by Stochastic Search. International Journal of Computer Vision. 2005;61(2):185."
\code
\endcode

@author salinas
*/
class AnnealedPFParallel {
public:
    struct Parameters
    {
      Parameters(int NParticles=300,int NLayers=10,double DesiredAlpha=0.5 ){
	nParticles=NParticles;
	nLayers=NLayers;
	desiredAlpha=DesiredAlpha;
       }
      int nParticles;
      int nLayers;
      double desiredAlpha;
    };
    /**Empty constructor
    */
    AnnealedPFParallel();

    /**Destructor
    */
    ~AnnealedPFParallel();
    /**Initializes the object to start the filter
    *@param nparticles number of particles
    *@param noise vector of initial noise employed in the first layer
    *@param initiParticle initial particle
    *@param  M number of layers
    */
    void setParams(Particle initParticle, vector<double> noise,Parameters  P)throw(GUException);


    
    /**Start the iterative process
     * \code
     * for( Alg.start();Alg.end()==0; ) Alg.nextStep();
     * \endcode
     */
     void start()throw(GUException);
     /**Gives a new step of the iterative process, ie., runs the next layer
      */
     void nextStep()throw(GUException);
     /**Indicates if the last layer has been executed.
      */
     bool end()throw(GUException);
     
    /**Calculates the probability of the particle passed
     * @param p particle
     * @param threadId id of the current thread performing the operation.
    */
    virtual double calculateParticleProb(Particle &p,unsigned int threadId)=0;
    /**Sets the seed for the random generator
      */
    void setRandomSeed(unsigned int seed)
    {
        srand(seed);
    }
    /**Debugging info
     */
    void printState()
    {
        for (unsigned int i=0;i<currParticles->size();i++)
            cout <<"i:="<<i<<" "<<*(*currParticles)[i];
    }

  /**Returns the number of parallel threads that will be employed to do calcualtion
     */
    int getNThreads() {
        return omp_get_max_threads();
    }

    /**Calculates an estimation of the state
     * vector using the mean of all the particles using the probability of each one
     */
    Particle getUnimodalEstimate(); 
    /**Sets a mask so that only active elements are employed for the search
    */
    void setMask(vector<bool> activeElements)throw(GUException);
protected:
    /** Coge de forma adecuada los elementos para la siguiente iteraci�
     */
    void samplingAndPropagation();
    /**
     */
    void observe();
 
    /**
     */
    Particle* selectRandomParticle( vector<Particle*>&pv);
    
 
    
    
    void deallocate();
    vector<Particle*> particles,particles2;
    vector<Particle*> *currParticles,*prevParticles,*aux;
    vector<double> _noise,_originalWeights;
    Parameters _params;
    bool _areParms;
 
    //control the iterations
    int _currentLayer;
    float _currentExp;
    
    vector<double> mean,sqmean,var;//employed for calculus of variances
    vector<bool> _activeElements;
    //debug
    float getCurrentSurvivalRate();
    float  getSurvivalRate(float pval);
    double getOptimalExp(double initExp=1);
};

};
};
#endif

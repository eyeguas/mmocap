#include "annealedpfparallel.h"
#include <omp.h>
#include <cmath>
#include <gu/gurandom.h>
namespace gummocap {
namespace algorithms {

#define _ANNPG_RN (double(rand())/double(RAND_MAX))
/**Empty constructor
    */
AnnealedPFParallel::AnnealedPFParallel()
{
    _areParms=false;
}

/**Destructor
*/
AnnealedPFParallel::~AnnealedPFParallel()
{
    deallocate();

}
/**
 */
void AnnealedPFParallel::deallocate()
{
    for (unsigned int i=0;i<particles.size();i++)       delete particles[i];
    for (unsigned int i=0;i<particles2.size();i++)       delete particles2[i];
    particles.clear();
    particles2.clear();

}
/**
*/
void AnnealedPFParallel::setParams(Particle initParticle, vector<double> noise,Parameters  P)throw(gu::Exception)
{
    if (noise.size()!=initParticle.size())
        throw gu::Exception("AnnealedPFParallel::setParams noise has not the same dimensions than initParticle");
    _noise=noise;
    _activeElements.resize(noise.size());
    for(unsigned int i=0;i<_activeElements.size();i++) _activeElements[i]=true;
    _params=P; 
    //delete old partcicles if any
    deallocate();
    currParticles=&particles;
    prevParticles=&particles2;

    //now create the particles an assign, equal prob to all particles
    double initProb=1./double(_params.nParticles );
    _originalWeights.resize(_params.nParticles );
    float accProb=0;
    for (int i=0;i<_params.nParticles;i++) {
        currParticles->push_back(new Particle(initParticle));
        currParticles->back()->prob=initProb;
        accProb+=initProb;
	currParticles->back()->accProb=accProb;	
        prevParticles->push_back(new Particle(*currParticles->back()));	
        _originalWeights[i]=initProb;
    } 
    _areParms=true;
    _currentLayer=-1;
}

 
/**
*
*/
void AnnealedPFParallel::start()   throw(gu::Exception)
{
    if (!_areParms) throw  gu::Exception("AnnealedPFParallel::start no call to setParams()");
    if (_currentLayer!=-1 && _currentLayer<_params.nLayers) throw  gu::Exception("AnnealedPFParallel::start end of layers not reached last time");
    _currentLayer=0;
}

/**
*
*/
void AnnealedPFParallel::nextStep()throw(gu::Exception)
{
    if (!_areParms) throw  gu::Exception("AnnealedPFParallel::nextStep no call to setParams()");
    if (_currentLayer>=_params.nLayers) throw  gu::Exception("AnnealedPFParallel::nextStep  number of iterations exceeded");
    if (_currentLayer==-1) throw  gu::Exception("AnnealedPFParallel::nextStep  no call to start");

    if (_currentLayer==0) {
        _currentExp=1; 
    } 
    samplingAndPropagation();
    observe();
    _currentExp=getOptimalExp(_currentExp);
    _currentLayer++;
}
/**
*
*/
bool AnnealedPFParallel::end()throw(gu::Exception)
{
    if (!_areParms) throw  gu::Exception("AnnealedPFParallel::nextStep no call to setParams()");
    return _currentLayer==_params.nLayers;
}

/** Coge de forma adecuada los elementos para la siguiente iteraci�
 */
void AnnealedPFParallel::samplingAndPropagation() {
    //interchange prev and current particles vectors
    aux=currParticles;
    currParticles=prevParticles;
    prevParticles=aux;

    //resized appropriateadly
    mean.resize((*prevParticles)[0]->size());
    sqmean.resize((*prevParticles)[0]->size());
    var.resize((*prevParticles)[0]->size());
    for (unsigned int i=0;i<mean.size();i++)     mean[i]=sqmean[i]=0;
    
    
    //Start doing propagation based on probability
    Particle *part;
    double sumProb=0;
//     cout<<"NP="<<prevParticles->size()<<endl;
    for (unsigned int p=0;p<prevParticles->size();p++)
    {
        //selecciona una particula
        part=selectRandomParticle(*prevParticles);
        //copy the particle
        *(*currParticles)[p]=*part;
        //update mean and sqmean
        for (unsigned int j=0;j<mean.size();j++) {
            mean[j]+=(*part)[j] * part->prob;
            sqmean[j]+=(*part)[j]*(*part)[j] * part->prob;
	    sumProb+= part->prob;
        }
    }


    //now, calculate variances
//     double invN=1./float(mean.size());
    double invN=1./float(sumProb);
    for (unsigned int i=0;i<mean.size();i++) {
        mean[i]*=invN;
        sqmean[i]*=invN;
        var[i]=  10* sqmean[i] - (mean[i]*mean[i]);	//variance
        if (var[i]<0) {
            var[i]=0;    //it might happen due to precision errors
        }
         else if (var[i]>_noise[i]) var[i]=_noise[i];
    }
    
    //finally, lets propagate particles according to this value
    if (_currentLayer==0)//first time, user defined noise
    {
        for (unsigned int i=0;i<currParticles->size();i++)
        {
            Particle *p=(*currParticles)[i];
            for (unsigned int j=0;j<p->size();j++){
	      if (_activeElements[j])//variate only if allowed by the mask
		  (*p)[j]+=gu::Random::nextGaussian() * _noise[j];
	    }
        }
    }
    else {//next, use variance information 
        
/* for (unsigned int i=0;i<(*currParticles)[0] ->size();i++)
   cout<< var[i]<<",";
 cout<<endl;*/
        for (unsigned int i=0;i<currParticles->size();i++)
        {
            Particle *p=(*currParticles)[i];
            for (unsigned int j=0;j<p->size();j++){ 
	      if (_activeElements[j])//variate only if allowed by the mask
                (*p)[j]+=gu::Random::nextGaussian() * var[j];
	    }
        }
    }
}
//////////////////////////////////
//
//
//////////////////////////////////
void AnnealedPFParallel::observe()
{
    //Calculate accumulated probability and do observation
    double probSum=0;
#pragma omp parallel for reduction(+:probSum)
    for (unsigned int i=0;i<currParticles->size();i++)
    {
        _originalWeights[i]=calculateParticleProb(*(*currParticles)[i],omp_get_thread_num());
        //applis the expent to sharpen of smooth the functoin accoring to survival rate
        (*currParticles)[i]->prob=std::pow(_originalWeights[i],_currentExp);
        probSum+=(*currParticles)[i]->prob;
//         cout<<*(*currParticles)[i]<<"->"<<_originalWeights[i]<<",probSum="<<probSum<<", exp="<<_currentExp<<endl;

    }
    //ahora normalizamos las probabilidades y calculamos la acumulada
    double accProb=0;
    probSum=1./probSum;
    for (unsigned int i=0;i<currParticles->size();i++) {
        (*currParticles)[i]->prob*=probSum;
        accProb+=(*currParticles)[i]->prob;
        (*currParticles)[i]->accProb=accProb;
    }

}

/**Calculates an estimation of the state
 * vector using the mean of all the particles using the probability of each one
 */
Particle AnnealedPFParallel::getUnimodalEstimate()
{
    Particle res(*currParticles->back());
    //se inicia a cero caad elemento del vector de la particula
    for (unsigned int i=0;i<res.size();i++) res[i]=0;
    //se comienza calacula la media ponderada por la probabilidad de cada particula
    for (unsigned int i=0;i<currParticles->size();i++) {
        for (unsigned int j=0;j<res.size();j++) {
            res[j]+=(*(*currParticles)[i])[j]*(*currParticles)[i]->prob;
        }
    }
    return res;
}

///////////////////////////////////////////
//
//
///////////////////////////////////////////

Particle* AnnealedPFParallel::selectRandomParticle( vector<Particle*>&pv)
{
    //genera un numero aleatorio
    double random=_ANNPG_RN;
    int init=0,end=pv.size()-1,middle;
    if (random<pv.front()->accProb) {
        return  pv.front();
    }
    //ahora encontramos por busqueda binaria la partiucla con dicha
    //probabilidad acumulada
    while (true) {
        if (init==end || init+1==end) return pv[init];
        middle=(end+init)/2;
        if (pv[middle]->accProb<=random && random<pv[middle+1]->accProb) {
            return pv[middle+1];
        }
        else if (random > pv[middle]->accProb)
            init=middle;
        else  end=middle;
    }
}

float AnnealedPFParallel::getCurrentSurvivalRate()
{
    double D=0;
    for (unsigned int i=0;i<currParticles->size();i++)
    {

        D+=(*currParticles)[i]->prob *(*currParticles)[i]->prob;
    }
    double alpha=1./(D *float(currParticles->size()));
    cout<<"alpha="<<alpha<<" OF "<<1./float(currParticles->size())<<endl;
    return alpha;

}




double AnnealedPFParallel::getOptimalExp(double initExp)
{
    double bestExp=1,currExp=initExp,increment=10;
    double maxErrorAllowed=0.1,minError=9999;
    double incSig=0;
    int nIter=0;
    //We seach for the exp that makes
    //_alphaDesired= B(exp)/D

    //decide first the initial direction

    double aplhaposError= fabs( getSurvivalRate( currExp+increment ) -_params.desiredAlpha);
    double aplhaFirstError= fabs( getSurvivalRate( currExp ) - _params.desiredAlpha);

    if (aplhaFirstError<aplhaposError) return currExp;
    else {
        incSig=1;
        currExp+=increment;
        minError=aplhaposError;
    }
    bestExp=currExp;

//where D is the survival Diagnosis
    //now, do a binary search to find the best exp
//     cout<<"seaching optimal exp..."<<currExp<<" "<<minError<<endl;

    bool  someImprovement=true;
    while (minError>maxErrorAllowed && nIter< 50  && someImprovement)
    {
        currExp+=incSig*increment;
        float currSR=getSurvivalRate( currExp );
        double currError=fabs( currSR-_params.desiredAlpha);
        if (currError>minError) {
            incSig*=-1;
            increment*=0.5;

        }
        else if (currError<minError) {
            minError=currError;
            bestExp=currExp;
        }
//       else someImprovement=false;
        nIter++;
//       cout<<"seaching optimal exp..."<<currExp<<" "<<minError<<" iter="<<nIter<<" SR="<<currSR<<" CE="<<currError<<endl;
    }
    return bestExp;
}

float AnnealedPFParallel::getSurvivalRate(float pval)
{
    double  probSum=0;
    vector<double> powWeights(_originalWeights.size());
    for (unsigned int i=0;i<_originalWeights.size();i++)
        probSum+=powWeights[i]=pow(_originalWeights[i],pval);
    //normalizamos
    float invSum=1./probSum;
    double D=0;
    for (unsigned int i=0;i<_originalWeights.size();i++) {

        powWeights[i]*=invSum;
        D+=powWeights[i]*powWeights[i];
    }
    //now, calculate

    return 1./(D *float(powWeights.size()));
}

/////////////////////////////////
//
//
/////////////////////////////////
void AnnealedPFParallel::setMask(vector<bool> activeElements)throw(gu::Exception)
{
if (  activeElements.size()!=_activeElements.size()) throw gu::Exception("AnnealedPFParallel::setMask the mask must have the same size than variable in the search space");
_activeElements=activeElements;
}

}
}

#include "pso.h"

namespace gummocap
{
namespace algorithms
{
  
  
 /**
 *
 *
 */
  PSO::PSO(){
      _areParams=false;
      _xx=NULL;
      _vx=NULL;
      _tx=NULL;
      _pbestx=NULL;
      _pbest=NULL;
      _pf=NULL;
      _besteverx=NULL;
  }
 /**
 *
 *
 */ 
  PSO::~PSO(){
    if(_areParams){
      for(unsigned int i=0;i<_params.number_of_agents;i++){
	delete[] _xx[i];
	delete[] _vx[i];
	delete[] _tx[i];
	delete[] _pbestx[i];
      }
      delete[] _xx;
      delete[] _vx;
      delete[] _tx;
      delete[] _pbestx;
      delete[] _pbest;
      delete[] _pf;
      delete[] _besteverx;
      
    }
  }
 /**
 *
 *
 */
  void PSO::setParams(Params  p){
    
    _params=p; 
    _xx=new double*[_params.number_of_agents];
    _vx=new double*[_params.number_of_agents];
    _tx=new double*[_params.number_of_agents];
    _pbestx=new double*[_params.number_of_agents];
    _pbest=new double[_params.number_of_agents];
    _pf=new double[_params.number_of_agents];
    _besteverx=new double[_params.dimension];
    
    for(unsigned int i=0;i<_params.dimension;i++){
      _besteverx[i]=0.0;
    }
    
    for(unsigned int i=0;i<_params.number_of_agents;i++){
      _xx[i]=new double[_params.dimension];
      _vx[i]=new double[_params.dimension];
      _tx[i]=new double[_params.dimension];
      _pbestx[i]=new double[_params.dimension];
    }
    _besteverf=INFINITE;
    _areParams=true;
  }
  /**
 *
 *
 */
  void PSO::initialize(){
    	
	for (unsigned int a=0;a<_params.number_of_agents;a++)
		{
			for (unsigned int b=0;b<_params.dimension;b++)
			{
				_xx[a][b] = (double) GURandom::nextFloat(_params.irang_l,_params.irang_r);
				//_xx[a][b] = (double) ((_params.irang_r - _params.irang_l)*(GURandom::nextGaussian()) + _params.irang_l);
				_pbestx[a][b]=_xx[a][b];
				_vx[a][b] = _params.maxv*GURandom::nextFloat(1.0);

				if ((GURandom::nextFloat(1.0)) > 0.5) _vx[a][b]=-_vx[a][b];
				_tx[a][b]=_xx[a][b]+_vx[a][b];
			}
			_pbest[a]=INFINITE;
			_pf[a]=INFINITE;
		}
	_weight_up=_params.initial_weight;
	_C=_params.stop_maxeval/_params.number_of_agents;
	_incr_c=log(10.0*_params.initial_weight)/_C;
	_c=0;
	
	_current_it=1;
	_current_nevals=0;
	_gbest=0;
    
  }
  /**
 *
 *
 */
  void PSO::initialize(double* x){
    	
	for (unsigned int a=0;a<_params.number_of_agents;a++)
		{
			for (unsigned int b=0;b<_params.dimension;b++)
			{
				//_xx[a][b] = (double) GURandom::nextFloat(x[b]-_params.irang_l,x[b]+_params.irang_r);
				_xx[a][b] = (double) ((_params.irang_l)*(GURandom::nextGaussian()) + x[b]);
				_pbestx[a][b]=_xx[a][b];
				_vx[a][b] = _params.maxv*GURandom::nextFloat(1.0);

				if ((GURandom::nextFloat(1.0)) > 0.5) _vx[a][b]=-_vx[a][b];
				_tx[a][b]=_xx[a][b]+_vx[a][b];
			}
			_pbest[a]=INFINITE;
			_pf[a]=INFINITE;
		}
	_weight_up=_params.initial_weight;
	_C=_params.stop_maxeval/_params.number_of_agents;
	_incr_c=log(10.0*_params.initial_weight)/_C;
	_c=0;
	
	_current_it=1;
	_current_nevals=0;
	_gbest=0;
    
  }
  /**
 *
 *
 */
  void PSO::updateInertiaWeight(){
		
		//_weight_up = (_params.initial_weight-0.4) * (_params.stop_maxiter - _current_it) /_params.stop_maxiter +0.4;    //time variant weight, linear from weight to 0.4
		//_weight_up=_params.initial_weight;		//constant inertia weight
		
		// Inertia weight rule (Markerless human articulated tracking using hierarchical particle swarm optimization)
		_weight_up = _params.initial_weight/exp(_c);
		_c+=_incr_c;
		
  }
  /**
 *
 *
 */ 
  void PSO::updateBest(){
    
    for (unsigned int a=0;a<_params.number_of_agents;a++)
    {
	_current_nevals++;
	if (_pf[a] < _pbest[a])
	{
	  _pbest[a]=_pf[a];
	  for (unsigned int b=0;b<_params.dimension;b++) 
	    _pbestx[a][b]=_xx[a][b];
	  if (_pbest[a] < _pbest[_gbest])
	  {
	    _gbest=a;
	  }
	}

    }
			
  }
   /**
 *
 *
 */ 
  void PSO::updateVelocityAndPosition(){
    
    for (unsigned int a=0;a<_params.number_of_agents;a++)
    {
	/* Asynchronous version */
	for (unsigned int b=0;b<_params.dimension;b++)
	{
		_vx[a][b] = _weight_up*_vx[a][b] + _params.c1*(GURandom::nextFloat(1.0))*(_pbestx[a][b]-_xx[a][b]) +
					_params.c2*(GURandom::nextFloat(1.0))*(_pbestx[_gbest][b]-_xx[a][b]);
		if (_vx[a][b]>_params.maxv)
			_vx[a][b]=_params.maxv;
		else if (_vx[a][b]<-_params.maxv)
			_vx[a][b]=-_params.maxv;
		_tx[a][b]=_xx[a][b]+_vx[a][b];
		// Update positions. Define new coordinates
		_xx[a][b] =_tx[a][b];
	}
    }
  }
   /**
 *
 *
 */ 
  bool PSO::checkTerminationCriteria(){
    
    if ((_pbest[_gbest] <= _params.stop_cutoff) || (_current_it >= _params.stop_maxiter) || (_current_nevals >= _params.stop_maxeval))
    {
      return true;
    }
    else
      return false;
  }
    /**
 *
 *
 */ 
  void PSO::nextIteration(){
    
    updateInertiaWeight();
    evaluate();
    updateBest();
    updateVelocityAndPosition();
    _current_it++;
    
  }
   /**
 *
 *
 */ 
   void PSO::resetBestEver(){
     _besteverf=INFINITE;
   }
   /**
 *
 *
 */ 
  void PSO::run(double* x){
    initialize(x);
    while(!checkTerminationCriteria())
      nextIteration();
    if(_pbest[_gbest]<_besteverf){
      _besteverf=_pbest[_gbest];
      for(unsigned int i=0;i<_params.dimension;i++)
	_besteverx[i]=_pbestx[_gbest][i];
    }
  }
   /**
 *
 *
 */ 
  void PSO::run(){
    initialize();
    while(!checkTerminationCriteria())
      nextIteration();
    if(_pbest[_gbest]<_besteverf){
      _besteverf=_pbest[_gbest];
      for(unsigned int i=0;i<_params.dimension;i++)
	_besteverx[i]=_pbestx[_gbest][i];
    }
  }
   /**
 *
 *
 */ 
  void PSO::execute(){
    for(unsigned int i=0;i<_params.no_runs;i++){
      run();
    }
  }
  /**
 *
 *
 */ 
  void PSO::execute(double* x){
    for(unsigned int i=0;i<_params.no_runs;i++){
      run(x);
    }
  }
    /**
 *
 *
 */ 
  double* PSO::getBest(){
    return _besteverx;
  }
    /**
 *
 *
 */ 
  double PSO::getBestFitness(){
    return _besteverf;
  }
}

}
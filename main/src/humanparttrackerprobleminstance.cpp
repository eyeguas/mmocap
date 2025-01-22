#include "humanparttrackerprobleminstance.h"

namespace gummocap {
  
  namespace humantracker{

   HumanPartTrackerProblemInstance::HumanPartTrackerProblemInstance(Random *random, unsigned int dim) {
    m_random = random;
    m_ndim = dim;
    m_init = false;
  }
  /*
  *
  *
  *
  */
  HumanPartTrackerProblemInstance::~HumanPartTrackerProblemInstance(void) {
  }
  /*
  *
  *
  *
  */
  void HumanPartTrackerProblemInstance::init(tChromosomeReal* xglobal, unsigned int inipart, unsigned int endpart, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func) {
    _bodyModel=sk;
    _TheViewSet=VSet;
    _fitness_function=fitness_func;
    _xglobal=xglobal;
    _inipart=inipart;
    _endpart=endpart;
    m_init = true;
  }
  /*
  *
  *
  *
  */
  ProblemPtr HumanPartTrackerProblemInstance::get(unsigned int maxeval){
    //string name;
    //double min=-0.1, max=0.1;

    ProblemPtr prob (new HumanPartTrackerProblem(_xglobal,_inipart,_endpart,_bodyModel,_TheViewSet,_fitness_function));
    prob->setDimension(m_ndim);
    /*
    for (unsigned i = 0; i < m_ndim; ++i) {
	prob->setDomainValues(i, min, max, false);
    }
    */
    // Define the optimum criteria
    prob->setOptimize(0., 1e-8);
    //prob->setMaxEval(10000*m_ndim);
    prob->setMaxEval(maxeval);
    prob->setMinimize();

    return prob;
  }
    
  }
  
  
}

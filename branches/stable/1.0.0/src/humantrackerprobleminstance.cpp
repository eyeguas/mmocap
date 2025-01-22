#include "humantrackerprobleminstance.h"

namespace gummocap {
  
  namespace humantracker{

   HumanTrackerProblemInstance::HumanTrackerProblemInstance(Random *random, unsigned int dim) {
    m_random = random;
    m_ndim = dim;
    m_init = false;
  }
  /*
  *
  *
  *
  */
  HumanTrackerProblemInstance::~HumanTrackerProblemInstance(void) {
  }
  /*
  *
  *
  *
  */
  void HumanTrackerProblemInstance::init(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func) {
    _bodyModel=sk;
    _TheViewSet=VSet;
    _fitness_function=fitness_func;
    m_init = true;
  }
  /*
  *
  *
  *
  */
  void HumanTrackerProblemInstance::setViewSet(guscene::view::ViewSet *VSet){
     _TheViewSet=VSet;
     
   }
  /*
  *
  *
  *
  */
  ProblemPtr HumanTrackerProblemInstance::get(unsigned int maxeval){
    //string name;
    //double min=-0.1, max=0.1;

    
    ProblemPtr pproblem(new HumanTrackerProblem(_bodyModel,_TheViewSet,_fitness_function));
    pproblem->setDimension(m_ndim);
    /*
    for (unsigned i = 0; i < m_ndim; ++i) {
	prob->setDomainValues(i, min, max, false);
    }
    */
    // Define the optimum criteria
    pproblem->setOptimize(0., 1e-8);
    //prob->setMaxEval(10000*m_ndim);
    pproblem->setMaxEval(maxeval);
    pproblem->setMinimize();
    
    return pproblem;
  }
    
  }
  
  
}
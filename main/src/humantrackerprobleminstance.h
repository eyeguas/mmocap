#ifndef __HUMANTRACKER_PROBLEMINSTANCE__
#define __HUMANTRACKER_PROBLEMINSTANCE__

#include <gurealea/realea/problemfactory.h>
#include "humantrackerproblem.h"
#include <gurealea/realea/random.h>
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include "poseevaluators.h"

#include <math.h>

using namespace std;
using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace gummocap::poseevaluators;

namespace gummocap {
  
  namespace humantracker{

class HumanTrackerProblemInstance : public ProblemFactory {
   public:
   
   HumanTrackerProblemInstance(Random *random, unsigned int dim);
  /*
  *
  *
  *
  */
  ~HumanTrackerProblemInstance(void);
  /*
  *
  *
  *
  */
  void init(BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func);
  /*
  *
  *
  *
  */
  void setViewSet(guscene::view::ViewSet *VSet);
  /*
  *
  *
  *
  */ 
  ProblemPtr get(unsigned int maxeval);
    private:
        Random *m_random;
	unsigned int m_ndim;
	bool m_init;
	PoseEvaluator* _fitness_function;
	BodyModel<HumanSkeleton> _bodyModel;
	guscene::view::ViewSet *_TheViewSet;
	
};
  
  }
}

#endif

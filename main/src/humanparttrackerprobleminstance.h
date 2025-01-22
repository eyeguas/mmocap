#ifndef __HUMANPARTTRACKER_PROBLEMINSTANCE__
#define __HUMANPARTTRACKER_PROBLEMINSTANCE__

#include <gurealea/realea/problemfactory.h>
#include "humanparttrackerproblem.h"
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

class HumanPartTrackerProblemInstance : public ProblemFactory {
   public:
   
   HumanPartTrackerProblemInstance(Random *random, unsigned int dim);
  /*
  *
  *
  *
  */
  ~HumanPartTrackerProblemInstance(void);
  /*
  *
  *
  *
  */
  void init(tChromosomeReal* xglobal, unsigned int inipart, unsigned int endpart, BodyModel<HumanSkeleton> &sk,guscene::view::ViewSet *VSet, PoseEvaluator* fitness_func);
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
	tChromosomeReal* _xglobal;
	unsigned int _inipart;
	unsigned int _endpart;
};
  
  }
}

#endif

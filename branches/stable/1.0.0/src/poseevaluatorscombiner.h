#ifndef _gummocap_PoseEvaluatorCombiner_
#define _gummocap_PoseEvaluatorCombiner_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>
#include "poseevaluators.h"
namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace guscene::view;
/**
 * \brief Combines multiples PoseEvaluators providing a single value
 */
class PoseEvaluatorCombiner: public PoseEvaluator
{
public:
      /**Modes of combining the info from multiple evaluators
     */
    enum CombMode  {ADD,AVRG,MULT};
    /**
     */
    PoseEvaluatorCombiner();
    /**
     */
    PoseEvaluatorCombiner( const PoseEvaluatorCombiner &X);
    /**
     */
    
    ~PoseEvaluatorCombiner();

    /**Sets the parameter
     */
    
    void setParams(CombMode mode=AVRG);
    /**Adds a new pose evaluator
     */
    void add(PoseEvaluator *pe);
    
    /**
     * Evaluates the solution given
     */
     double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;

     /**Makes a copy of this
     */
    PoseEvaluator * makeCopy();
    /**Returns a string identifying the evaluator
     */
    string getId();

    /**Call when a new frame is obtained and before a call to evaluate. Some evaluators
     * might require this in order to intialize data
     * @param vs next frame to be analyzed
     * @param frameId identifies the frame (for instance the frame number)
     */
    void startFrame(ViewSet *vs,int frameId);
    /**Call when a frame has being processed passing the best estimation obtained
     */
    void frameFinished(BodyModel<HumanSkeleton> *body );
    /**Number of evaluations
     */
    unsigned int numberOfEvaluations(){return _nEvals;}

private:
  vector<PoseEvaluator*> _eval;
  bool _areParams,_destroyMem;
  CombMode _mode;
  int _nEvals;
 
  
};

}
}

#endif






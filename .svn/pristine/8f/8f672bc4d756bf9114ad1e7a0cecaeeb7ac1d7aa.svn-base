#ifndef _gummocap_POSEEVALUATORS_
#define _gummocap_POSEEVALUATORS_
#include <gumocap/bodymodel.h>
#include <gumocap/humanskeleton.h>
#include <guscene/viewset.h>

namespace gummocap
{

namespace poseevaluators
{

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;
using namespace guscene::view;

class PoseEvaluator
{
public:
    /**Destructor
     */
    virtual ~PoseEvaluator() {}
    /**
     * Evaluates the solution given
     */
    virtual double evaluate ( BodyModel<HumanSkeleton> *body,ViewSet *vs ) =0;
    /**
     * Evaluates the solution given
     */
    virtual double evaluate ( BodyModel<HumanSkeleton> *body,ViewSet *vs,__attribute__ ( ( unused ) ) unsigned int n_inipart,__attribute__ ( ( unused ) ) unsigned int n_endpart ) {return evaluate ( body,vs );}
    /**Makes a copy of this
     */
    virtual PoseEvaluator * makeCopy() =0;
    /**Returns a string identifying the evaluator
     */
    virtual string getId() =0;
    /**Call when a new frame is obtained and before a call to evaluate. Some evaluators
     * might require this in order to intialize data
     * @param vs next frame to be analyzed
     * @param frameId identifies the frame (for instance the frame number)
     */
    virtual void startFrame ( ViewSet *vs,int frameId ) {}
    /**Call when a frame has being processed passing the best estimation obtained
     */
    virtual void frameFinished ( BodyModel<HumanSkeleton> *body ) {}
    /**Number of evaluations
     */
    virtual unsigned int numberOfEvaluations() =0;
    /**Returns the individual fitness value of in evaluated camera of the last evaluation
     */
    virtual vector<float> getViewFitness()throw(GUException){throw GUException(" getViewEvaluations not implemented");}
    /**Set a mask to select these views to be employed
     */
    virtual void setValidViews(vector<bool> vviews)throw(GUException){throw GUException(" PoseEvaluator::getViewEvaluations not implemented");}
};

}
}

#endif






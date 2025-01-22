#ifndef _gummocap_FogVH_POSEEVALUATORS_
#define _gummocap_FogVH_POSEEVALUATORS_

#include "poseevaluators.h"
#include <omp.h>

namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class FogVH: public PoseEvaluator
            /**\brief This class performs visual hull inersection of forground points analyzed
            * For evaluating of a pose estimate, the shape of human subjects is evaluated by silhouette comparison of the model projection masks and
            * binary foreground masks extracted using a background substraction technique. We use the normalized sum of overlapped pixels in the foreground
            * by the model as error measure. The objective will be to minimize this error.
            */
{
public:
    /**
     */
    FogVH();
    /**
     */
     FogVH(const  FogVH &CF);
    /**
     */
    ~FogVH();
    /**
     */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs )throw(gu::Exception);
    /**
     */
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs);
    /**
     */
    PoseEvaluator*makeCopy();
    /**
     */
    string getId(){return "fogvh";}
    /**
     */
    unsigned int numberOfEvaluations(){return _nEvals;}
private: 
    unsigned int _nEvals;
    svector<int> _invalid_v;
    svector<float> _projections; 
    guscene::view::ViewSet *_TheViewSet;
    bool _areParams;
    svector<float> _intersections;//vector with the intersection results
};


}
}

#endif

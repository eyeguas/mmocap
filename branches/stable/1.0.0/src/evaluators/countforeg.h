#ifndef _gummocap_COUNTFOREG_POSEEVALUATORS_
#define _gummocap_COUNTFOREG_POSEEVALUATORS_

#include "poseevaluators.h"
#include <omp.h>

namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class CountForeground: public PoseEvaluator
            /**\brief This class represents the Xor fitness function.
            * For evaluating of a pose estimate, the shape of human subjects is evaluated by silhouette comparison of the model projection masks and
            * binary foreground masks extracted using a background substraction technique. We use the normalized sum of overlapped pixels in the foreground
            * by the model as error measure. The objective will be to minimize this error.
            */
{
public:
    /**
     */
    CountForeground();
    /**
     */
     CountForeground(const  CountForeground &CF);
    /**
     */
    ~CountForeground();
    /**
     */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs,bool parallel=true);
    /**
     */
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs);
    /**
     */
    PoseEvaluator*makeCopy();
    /**
     */
    string getId(){return "fog";}
    /**
     */
    unsigned int numberOfEvaluations(){return _nEvals;}
private:
    unsigned int _nEvals;
    bool _parallel;  
    vector<svector<int> >_invalid_v;
    vector<svector<float> >_projections;
    vector<CameraParameters *>_camPrms;
    guscene::view::ViewSet *_TheViewSet;
    bool _areParams;
};


}
}

#endif

#ifndef _gummocap_STEV_POSEEVALUATORS_
#define _gummocap_STEV_POSEEVALUATORS_

#include "poseevaluators.h"
#include <omp.h>
#include <gu/gutimestatistics.h>

namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class StEv: public PoseEvaluator
            /**\brief This class represents a fitness function based exclusively on stereo evaluation
            */
{
public:
    /**
     */
    StEv();
    /**
     */
    StEv(const  StEv&X);

    /**
     */
    ~StEv();
    /**Set the required params.
     * @params body to be evaluate
     * @param vs view set from which evaluation is done
     * @param parallel whether multiples threads must be employed
    */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs )throw(GUException);
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;

    PoseEvaluator * makeCopy();
        string getId(){return "stxor";}

    /**
     */
    unsigned int numberOfEvaluations(){return _nEvals;}
private: 
    unsigned int _nEvals;
    guscene::view::ViewSet *_TheViewSet;
    BodyModel<HumanSkeleton>* _bodyModel;
    bool _areParams;

    svector<int>  _invalid_v;
    svector<float> _projections;
    svector<float>  _zvertices;
    ///does xor, counts the non zero pixels, and set to zero the first image
    int xor_special(IplImage *im1,IplImage *im2);
    void getZverticesRelativeToCamera(gumocap::core::svector<Point3D<float> > &vertices,svector<float> &vz,CameraParameters *cp);
    //debug
    GUTimeStatistics Timer;
//     void writeVRML(string path ,int view);
};


}
}


#endif

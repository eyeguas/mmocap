#ifndef _gummocap_STXOR_POSEEVALUATORS_
#define _gummocap_STXOR_POSEEVALUATORS_

#include "poseevaluators.h"
#include <omp.h>
#include <gu/gutimestatistics.h>

namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class StXor: public PoseEvaluator
            /**\brief This class represents the Xor fitness function combined with stereo evaluation
            * For evaluating of a pose estimate, the shape of human subjects is evaluated by silhouette comparison of the model projection masks and
            * binary foreground masks extracted using a background substraction technique. We use the sum of pixelwise logical symmetric difference
            * as error measure. The objective will be to minimize this error.
            */
{
public:
    /**
     */
    StXor();
    /**
     */
    StXor(const  StXor&X);

    /**
     */
    ~StXor();
    /**Set the required params.
     * @params body to be evaluate
     * @param vs view set from which evaluation is done
     * @param parallel whether multiples threads must be employed
    */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs )throw(gu::Exception);
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;

    PoseEvaluator * makeCopy();
        string getId(){return "stxor";}

    /**
     */
    unsigned int numberOfEvaluations(){return _nEvals;}
private: 
    unsigned int _nEvals;
    IplImage *_proj,*_xorimg;
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
    gu::TimeStatistics Timer;
//     void writeVRML(string path ,int view);
};


}
}


#endif

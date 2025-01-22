#ifndef _gummocap_XORGrad_POSEEVALUATORS_
#define _gummocap_XORGrad_POSEEVALUATORS_

#include "poseevaluators.h"
#include <omp.h>
namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class XorGrad: public PoseEvaluator
            /**\brief This class represents the Xor  fitness function with gradient information in the borders
            * For evaluating of a pose estimate, the shape of human subjects is evaluated by silhouette comparison of the model projection masks and
            * binary foreground masks extracted using a background substraction technique. We use the sum of pixelwise logical symmetric difference
            * as error measure. The objective will be to minimize this error.
            */
{
  public:
    /**
     */
    XorGrad();
    /**
     */
    XorGrad(const  XorGrad&X);

    /**
     */
    ~XorGrad();
    /**Set the required params.
     * @params body to be evaluate
     * @param vs view set from which evaluation is done
     * @param dotPThres threshold employed to consider a point as part of the border (it has to be a negative value).
     * The smaller it is, the more points will be included
    */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs,float  dotPThres=0.15 )throw(GUException);
    /**
    */
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;
    /**
    */
    PoseEvaluator * makeCopy();
    string getId() {
        return "xorgrad";
    }
    /**
     */
    unsigned int numberOfEvaluations() {
        return _nEvals;
    }
private:
    unsigned int _nEvals;   
    IplImage * _proj;
    guscene::view::ViewSet *_TheViewSet;
    BodyModel<HumanSkeleton>* _bodyModel;

    bool _areParams;

    svector<int>_invalid_v;
    svector<float> _projections;
    svector<float> _dotproducts;
    //threshold employed to consider a point as part of the border
    float _dotPThres;
    ///does xor, counts the non zero pixels, and set to zero the first image
    int xor_special(IplImage *im1,IplImage *im2);


};


}
}

#endif

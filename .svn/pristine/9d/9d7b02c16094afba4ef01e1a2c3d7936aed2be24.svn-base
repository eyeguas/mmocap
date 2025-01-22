#ifndef _gummocap_XORColor_POSEEVALUATORS_
#define _gummocap_XORColor_POSEEVALUATORS_

#include "poseevaluators.h"
#include <omp.h>
#include "rgcolorhistogram.h"
namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class XorColor: public PoseEvaluator
            /**\brief This class represents the XorColor fitness function with color information
            * For evaluating of a pose estimate, the shape of human subjects is evaluated by silhouette comparison of the model projection masks and
            * binary foreground masks extracted using a background substraction technique. We use the sum of pixelwise logical symmetric difference
            * as error measure. The objective will be to minimize this error.
            */
{
public:
    /**
     */
    XorColor();
    /**
     */
    XorColor(const  XorColor&X);

    /**
     */
    ~XorColor();
    /**Set the required params.
     * @params body to be evaluate
     * @param vs view set from which evaluation is done
     * @param parallel whether multiples threads must be employed
    */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs )throw(GUException);
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;

    PoseEvaluator * makeCopy();
       string getId(){return "xorcolor";}
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
    svector<RGColorHistogram> _ColorModel; //histogram of each body part
    svector<RGColorHistogram> _auxColorModel; //histogram of each body part
    svector<RGColorHistogram> *_currColorModel;//pointer to the color model where it is going to be written


    ///does xor, counts the non zero pixels, and set to zero the first image
    int xor_special(IplImage *im1,IplImage *im2);

    ///sort body parts according their distance to the camera
    void sortBodyPart(vector<pair<unsigned int,float> > &indexDist,CameraParameters *cp);

};


}
}

#endif

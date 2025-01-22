#ifndef _gummocap_XOR_POSEEVALUATORS_
#define _gummocap_XOR_POSEEVALUATORS_

#include "poseevaluators.h"

#include <gu/gutimestatistics.h>
#include <omp.h>
#include <guimage/pointprojector.h>
namespace gummocap
{

namespace poseevaluators {

using namespace gumocap;
using namespace gumocap::models;
using namespace gumocap::skeleton;

class Xor: public PoseEvaluator
            /**\brief This class represents the Xor fitness function.
            * For evaluating of a pose estimate, the shape of human subjects is evaluated by silhouette comparison of the model projection masks and
            * binary foreground masks extracted using a background substraction technique. We use the sum of pixelwise logical symmetric difference
            * as error measure. The objective will be to minimize this error.
            */
{
public:
    /**
     */
    Xor();
    /**
     */
    Xor(const  Xor&X);

    /**
     */
    ~Xor();
    /**Set the required params.
     * @params body to be evaluate
     * @param vs view set from which evaluation is done
     * @param parallel whether multiples threads must be employed
    */
    void setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs);
    double evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) ;


    PoseEvaluator * makeCopy();
    string getId() {
        return "xor";
    }
    /**
     */
    unsigned int numberOfEvaluations() {
        return _nEvals;
    }
    
    vector<float> getViewFitness()throw(gu::Exception){
      return _viewEval;
    }
     void setValidViews(vector<bool> vviews)throw(gu::Exception){_validViews=vviews; }
private:
    unsigned int _nEvals;
    vector<float> _viewEval; //fitness in the last evaluation of each view
    vector<bool> _validViews;
    IplImage * _proj ;
    guscene::view::ViewSet *_TheViewSet;
    BodyModel<HumanSkeleton>* _bodyModel;
    bool _areParams;
 
    vector<float> _meshAvgSize;
    guimage::PointProjector ThePProj;
    ///does xor, counts the non zero pixels, and set to zero the first image
    int xor_special(IplImage *im1,IplImage *im2);
    int xor_special(IplImage *im1,IplImage *im2,int minX,int minY,int maxX,int maxY);
    gu::TimeStatistics Timer,Timer2;
     gu::TimeStatistics STTimer;
};


}
}

#endif


#ifndef _HumanTracker_H_
#define _HumanTracker_H_

#include <guscene/viewset.h>
#include <gumocap/skeleton.h>

namespace gummocap
{
namespace algorithms
{

/**
 * \brief Base class for all human trackers
 */
class HumanTracker
{
public:

    /**
    */
    ~HumanTracker() {}
    /**
    */
    virtual void  step ( guscene::view::ViewSet *VSet ) throw ( GUException ) =0;
    /**
     */
    virtual gumocap::skeleton::SkeletonPose getCurrentEstimation() =0;
    /**
     */
    virtual void enableDebug ( bool enable ) =0;
    /**Set a mask to select these views to be employed
      */
    virtual void setValidViews ( vector<bool> vviews ) throw ( GUException ) {throw GUException ( " HumanTracker::setValidViews not implemented for current algoritm" );}

};
};
};

#endif


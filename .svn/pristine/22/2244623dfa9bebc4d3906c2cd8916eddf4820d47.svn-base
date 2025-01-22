#ifndef _gummocap_OptFlowPointMatcher_H_
#define _gummocap_OptFlowPointMatcher_H_
#include <guscene/viewset.h>

namespace gummocap {
namespace poseevaluators {
  using namespace guscene::view;
/**Performs optical flow point matching of good tracking corners
 */
class  OptFlowPm
{
public:
    /**
     */
    OptFlowPm();
    /**
     */
    ~OptFlowPm();
    /**Sets the initial params
     */
    void setParams(ViewSet *vs,bool useSubPixRefinement=false);
    /**Find new features to track
    */
    void findFeatures(ViewSet *vs,vector< vector<CvPoint2D32f> > &features)throw(gu::Exception);
    /**find matches of previous image in the passed
     * @param newLocations new locations of the points and (-1,-1) if it was not located
    */
    void findMatches(ViewSet *vs,vector< vector<CvPoint2D32f> >&newLocations)throw(gu::Exception);

private:
    bool _findFeaturesCall;
    vector<CvPoint2D32f*> _points1,_points2;
    vector<CvPoint2D32f*>*_curPoints,*_prevPoints;
    vector<int> _nPoints;
    vector<IplImage *> _grey1,_grey2;
    vector<IplImage *> *_currGrey,*_prevGrey;
    vector<IplImage *> _eig,_tmp;
    vector<IplImage *> _pyramid1,_pyramid2;
    vector<IplImage *> *_curPyramid,*_prevPyramid;
    vector<char *>_status;
    int _flags ;
    bool _useSubPix;
    
};
};
};
#endif

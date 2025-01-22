#include "metaevaluator.h"
namespace gummocap
{
namespace poseevaluators
{
/**
 *
 *
 */
MetaEvaluator::MetaEvaluator()
{

}
/**
*
*
*/
vector<MetaEvaluator::CameraConfidence> MetaEvaluator::evaluate ( PoseEvaluator *PE,BodyModel<HumanSkeleton> *body,ViewSet *vs )
{
    positionConfidence.resize ( vs->size() );
    //calculate the confidence based on the position
    vector<float> camConf ( vs->size() );
    //calcualte the cos of the angles between the cameras
    //Person center
    Point3Df pPos=body->getSkptr()->get ( HumanSkeleton::Root ).xyz;

    for ( unsigned int i=0;i<vs->size();i++ )
    {
        camConf[i]=0;
        Point3Df ciPos,cjPos;
        vs->get ( i )->getCameraParameters()->getCameraLocation ( ciPos[0],ciPos[1],ciPos[2] );
        for ( unsigned int j=0;j<vs->size();j++ )
        {
            vs->get ( j )->getCameraParameters()->getCameraLocation ( cjPos[0],cjPos[1],cjPos[2] );
            camConf[i]+=getSin ( ciPos,cjPos,pPos );

        }
    }
    //now, determine the normalized values
    double confSum=0;
    for ( unsigned int i=0;i<vs->size();i++ ) confSum+=camConf[i];
    for ( unsigned int i=0;i<vs->size();i++ ) camConf[i]/=confSum;

    vector<MetaEvaluator::CameraConfidence> retCC ( vs->size() );
    for ( unsigned int i=0;i<vs->size();i++ )
    {
        retCC[i].camIndex=i;
        retCC[i].conf=camConf[i];
    }
    //now evaluate and adapt the value accordingly
    PE->evaluate ( body,vs );
    vector<float> vfitness=PE->getViewFitness();

    for ( unsigned int i=0;i<vs->size();i++ )
    {
        retCC[i].camIndex=i;
        retCC[i].conf*=1-vfitness[i];
        cout<<"cam :"<<i<<" Fitness="<<1-vfitness[i]<<endl;
    }
    sort ( retCC.begin(),retCC.end(),MetaEvaluator::CameraConfidence::Greater() );

    return retCC;
}


/** @brief  Calculate the cosine of the angle between the two cameras and the point
  *@param  xp, yp, zp : Position of the point in the world coordinates
 */
float MetaEvaluator::getSin ( Point3Df cam1, Point3Df cam2, Point3Df point )
{
    float x1p = cam1[0] - point[0];
    float y1p = cam1[1] - point[1];
    float z1p = cam1[2] - point[2];
    float x2p = cam2[0] - point[0];
    float y2p = cam2[1] - point[1];
    float z2p = cam2[2] - point[2];
    float cosAngle = ( x1p*x2p + y1p*y2p + z1p*z2p ) / ( sqrt ( x1p*x1p + y1p*y1p + z1p*z1p ) *
                     sqrt ( x2p*x2p + y2p*y2p + z2p*z2p ) );
    return 1-fabs ( cosAngle );
}
/**
 */
vector<bool> MetaEvaluator::getViewMaskConf ( vector<MetaEvaluator::CameraConfidence> &conf,float confLevel )
{
    assert(confLevel>=0 &&confLevel<=1);
    double sum=0;
    vector<bool> validViews ( conf.size() );
    for ( unsigned int i=0;i<conf.size();i++ )
    {
        if ( sum<confLevel ) validViews[conf[i].camIndex]=true;
        else validViews[conf[i].camIndex]=false;
        sum+=conf[i].conf;
        cout<<" cam:"<<conf[i].camIndex<<" "<<conf[i].conf<<" selected="<<validViews[conf[i].camIndex]<<" sum="<<sum<<endl;
    }
    return validViews;
}

/**
 */
vector<bool> MetaEvaluator::getViewMasknCams ( vector<MetaEvaluator::CameraConfidence> &conf,int n )
{
    vector<bool> validViews ( conf.size() );
    for ( unsigned int i=0;i<validViews.size();i++ )  validViews[i]=false;

    if (n>=0) {//use tne n best cameras
        if (n>conf.size()) n=conf.size();
        for ( unsigned int i=0;i<n;i++ )  validViews[conf[i].camIndex]=true;
    }
    else {
        n=-n;
        if (n>conf.size()) n=conf.size();
        for ( int i=conf.size()-n;i<conf.size(); i++ ) {cout<<"i"<<i<<endl; validViews[conf[i].camIndex]=true;}
    }


//     double sum=0;
//     cout<<"Using cams"<<endl;
//     for ( unsigned int i=0;i<conf.size();i++ )
//     {
//         sum+=conf[i].conf;
//         cout<<" cam:"<<conf[i].camIndex<<" "<<conf[i].conf<<" selected="<<validViews[conf[i].camIndex]<<" sum="<<sum<<endl;
//     }


    return validViews;
}
}
};

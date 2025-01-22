#include "countforeg.h"

namespace gummocap
{

namespace poseevaluators {

/*
*
*
*/
CountForeground::CountForeground() {
    _areParams=false;
    _parallel=true;
    
}
/*
*
*
*/
CountForeground::~CountForeground() {
  
 }
/*
*
*
*/
CountForeground::CountForeground(const  CountForeground &CF)
{
    _invalid_v=CF._invalid_v;
    _projections=CF._projections;
    _camPrms=CF._camPrms;
    _parallel=CF._parallel;
    _areParams=CF._areParams;
    _nEvals=CF._nEvals;
}

/*
*
*
*/
PoseEvaluator*CountForeground::makeCopy()
{
       return  new CountForeground(*this);
}
/*
*
*
*/
void CountForeground::setParams(BodyModel<HumanSkeleton> *bodyModel,ViewSet *TheViewSet,bool parallel) {
    _parallel=parallel;
    _areParams=true;
    _nEvals=0;
    //allocate auxiliar vectors with the projection of the points
    if ( _projections.size()==0) {
        //Determine the maximum number of points in a bodypart
        unsigned int maxP=0;
        for (unsigned int i=0;i<bodyModel->size();i++) {
            if (maxP<bodyModel->getMesh(i).getVertices().size())
                maxP=bodyModel->getMesh(i).getVertices().size();
        }
        _invalid_v.resize(TheViewSet->size());
        _projections.resize(TheViewSet->size());

        for (unsigned int i=0;i<TheViewSet->size();i++) {
            _invalid_v[i].resize(maxP*2);
            _projections[i].resize(maxP*2);
        }
    }
    _camPrms.resize(TheViewSet->size());
    for (unsigned int v=0;v<TheViewSet->size();v++) {
        _camPrms[v]=TheViewSet->get(v)->getCameraParameters  ();
    }
}
/*
*
*
*/
double CountForeground::evaluate(BodyModel<HumanSkeleton> *bodyModel,ViewSet *TheViewSet) {

  _nEvals++;
    double result=1.; // 1. --> HIGH VALUE FOR FITNESS FUNCTION
    double sumAvrg=0;

vector<IplImage*> masks(TheViewSet->size());

    for (unsigned int v=0;v<TheViewSet->size();v++) {
        masks[v]=TheViewSet->get(v)->getMaskImage();
    }

//Here the main loop, the point analysis
    int totalPointsAnalysed=0;
    int totalCount=0;

    if (_parallel==false)
        omp_set_num_threads(1);

     #pragma omp parallel for   reduction(+:totalPointsAnalysed,totalCount)
    for (unsigned int v=0;v<TheViewSet->size();v++) {
        for (unsigned part=0;part<bodyModel->size();part++) {
            gumocap::core::svector<Point3D<float> > &vertices=bodyModel->getMesh(part).getVertices();
            //project the body part
            _camPrms[v]->project3DPoints_4bytes_nodist((float*) vertices._data ,_projections[v]._data,vertices._size,_invalid_v[v]._data );
            //ok, now for each point calcuale the value
            int index=0;
            for (unsigned point=0;point<vertices.size() ;point++,index+=2) {
                if (!_invalid_v[v][index]) { //the points has projected
                    //note : _projections[v][index] is the x projection and _projections[v][index+1] the y projection in the image
                    unsigned char *mask_pix=(unsigned char *)(masks[v]->imageData + masks[v]->widthStep*int(_projections[v][index+1])+ int(_projections[v][index]));
                    if (*mask_pix==0) totalCount++;
                    totalPointsAnalysed++;
                }
            }
        }
    }


    if (_parallel==false)
        omp_set_num_threads(omp_get_max_threads() );
    
   return float(totalCount)/float(totalPointsAnalysed);
}
/*
*
*
*/

// double CountForeground::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs, unsigned int n_inipart, unsigned int n_endpart) {
//     //GUTimeMark Timer;
//     //Timer.init();
//     double result=1.; // 1. --> HIGH VALUE FOR FITNESS FUNCTION
// 
//     TheViewSet=vs;
//     bodyModel=body;
// 
//     //if(checkFeasibility(_auxBodyModel,bodyModel)){
//     //int t1=Timer.end();
//     //cout<<"t1="<<t1<<endl;
//     double sumAvrg=0;
// #pragma omp parallel for reduction(+:sumAvrg)
//     for (unsigned int v=0;v<TheViewSet->size();v++) {
//         CameraParameters *cp=TheViewSet->get(v)->getCameraParameters  ();
//         IplImage *mask=TheViewSet->get(v)->getMaskImage();
//         for (unsigned part=0;part<bodyModel->size();part++) {
//             int count=0;
//             gumocap::core::svector<Point3D<float> > &vertices=bodyModel->getMesh(part).getVertices();
//             for (unsigned point=0;point<vertices.size() ;point++) {
//                 Point3Df Tpoint=vertices[point];
//                 int xo,yo;
//                 if ( cp->project3DPoint  (Tpoint.xyzh[0],Tpoint.xyzh[1], Tpoint.xyzh[2], xo,yo)) {
//                     unsigned char *pix=(unsigned char *)  mask->imageData+ mask->widthStep * yo;
//                     if (pix[xo]>125) count++;
//                 }
//             }
//             if (vertices.size()>0)
//                 sumAvrg+= float(count) / float(vertices.size());
//         }
//     }
// 
//     //int t2=Timer.end();
//     //cout<<"time="<<t1<<" "<<t2<<" "<<sumAvrg<<endl;
//     if (sumAvrg>0)
//         result=1./sumAvrg;
//     //}
// 
//     return result;
// 
// }


}

}

#include "fogvh.h"

namespace gummocap
{

namespace poseevaluators {

/*
*
*
*/
FogVH::FogVH() {
    _areParams=false;

}
/*
*
*
*/
FogVH::~FogVH() {

}
/*
*
*
*/
FogVH::FogVH(const  FogVH &CF)
{
    _invalid_v=CF._invalid_v;
    _projections=CF._projections; 
    _areParams=CF._areParams;
    _nEvals=CF._nEvals;
}

/*
*
*
*/
PoseEvaluator*FogVH::makeCopy()
{
    return  new FogVH(*this);
}
/*
*
*
*/
void FogVH::setParams(BodyModel<HumanSkeleton> *bodyModel,ViewSet *TheViewSet ) throw(GUException)
{
    _nEvals=0;
    _areParams=true;
    int totalNumberOfPoints=0;
    //check that image intergrals are valid
    for(unsigned int i=0;i<TheViewSet->size();i++)
      if(TheViewSet->get(i)->getMaskIntegral()==NULL) throw GUException("FogVH::setParams integral masks must be calculated");
    //allocate auxiliar vectors with the projection of the points
        //Determine the maximum number of points in a bodypart
        unsigned int maxP=0;
        for (unsigned int i=0;i<bodyModel->size();i++) {
	  totalNumberOfPoints+=bodyModel->getMesh(i).getVertices().size();
            if (maxP<bodyModel->getMesh(i).getVertices().size())
                maxP=bodyModel->getMesh(i).getVertices().size();
	}
        for (unsigned int i=0;i<TheViewSet->size();i++) {
            _invalid_v.resize(maxP*2);
            _projections.resize(maxP*2);
        }

 
    _intersections.resize(totalNumberOfPoints );
}
/*
*
*
*/
double FogVH::evaluate(BodyModel<HumanSkeleton> *bodyModel,ViewSet *TheViewSet) {
  _nEvals++;
    //set the result to 1(all occupied)
    _intersections=1;
    
    for (unsigned int v=0;v<TheViewSet->size();v++) {
	    IplImage *maskImg=TheViewSet->get(v)->getMaskIntegral();
	    int wsize=2;
	    float invFac=1./(255.*(2*wsize+1)*(2*wsize+1));
	    CameraParameters *cp=TheViewSet->get(v)->getCameraParameters();
	    int intsPoint=0;
        for (unsigned part=0;part<bodyModel->size();part++) {
            gumocap::core::svector<Point3D<float> > &vertices=bodyModel->getMesh(part).getVertices();
            //project the body part
            cp->project3DPoints_4bytes_nodist((float*) vertices._data ,_projections._data,vertices._size,_invalid_v._data );
            //ok, now for each point calcuale the value
            int index=0;
            for (unsigned point=0;point<vertices.size() ;point++,index+=2,intsPoint++) {
                if (!_invalid_v[index] /*&& _intersections[intsPoint] */) { //the points has projected
		    //get the information form the integral image
		    int xo=_projections[index];
		    int yo=_projections[index+1];		    
		    int *y1=(int *)(maskImg->imageData + maskImg->widthStep*(yo-wsize));
		    int *y2=(int *)(maskImg->imageData + maskImg->widthStep*(yo+wsize));		    		    		    
		    float val= invFac*float(y2[xo+wsize] - y2[xo-wsize] - y1[xo+wsize] +y1[xo-wsize]);
// 		    if (val>1) cout<<"ERR="<<val<<" "<<float(y2[xo+wsize] - y2[xo-wsize] - y1[xo+wsize] +y1[xo-wsize])<<endl;
                    _intersections[intsPoint]= min( _intersections[intsPoint], val );
                }
            }
        }
    }

    //now, count the non-empty points
    float sumEmpty=0;
    for(unsigned int i=0;i<_intersections.size();i++)
	sumEmpty+=_intersections[i];
    float success=float(sumEmpty)/float(_intersections.size());
    
    //the error
	float error=1- success;
    return error;
}
/*
*
*
*/

// double FogVH::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs, unsigned int n_inipart, unsigned int n_endpart) {
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

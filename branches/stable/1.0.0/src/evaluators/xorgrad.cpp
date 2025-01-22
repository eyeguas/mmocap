#include "xorgrad.h"
#include <gu/gutimestatistics.h>
#include <opencv/cv.h>
#include <opencv/highgui.h>
#include <pmmintrin.h>
#include <mmintrin.h>
#include <emmintrin.h>
#include <gu/gualloc.h>
#include <gu/gurandom.h>

// namespace gummocap
// {
// 
// namespace poseevaluators {
// 
// /*
// *
// *
// */
// XorGrad::XorGrad() {
//     _proj=NULL;
//     _areParams=false;
// }
// 
// /*
// *
// *
// */
// XorGrad::XorGrad(const  XorGrad&X)
// {
//     _areParams=X._areParams;
//     _TheViewSet=X._TheViewSet;
//     _bodyModel=X._bodyModel;
//     _invalid_v=X._invalid_v;
//     _projections=X._projections;
//     _dotproducts=X._dotproducts;
//     _dotPThres=X._dotPThres;
//     _nEvals=X._nEvals;
//     if (X._proj!=NULL)
//         _proj=cvCloneImage(X._proj );
//     else
//         _proj=NULL;
// 
// }
// /*
// *
// *
// */
// 
// PoseEvaluator * XorGrad::makeCopy()
// {
//     return new XorGrad(*this);
// }
// /*
// *
// *
// */
// XorGrad::~XorGrad() {
//     if (_proj!=NULL)
//         cvReleaseImage(&(_proj));
// 
// }
// /*
// *
// *
// */
// void XorGrad::setParams(BodyModel<HumanSkeleton> *body,ViewSet *vs,float  dotPThres ) throw(GUException)
// {
//     _TheViewSet=vs;
//     _bodyModel=body;
//     _dotPThres=dotPThres; 
//     _nEvals=0;
//     //check that the view set has the extension of grandients
//     for (unsigned int i=0;i<vs->size();i++) {
// 	gummocap::view::ViewExtension  *VE= dynamic_cast<gummocap::view::ViewExtension *>( vs->get(i));
//         if (VE==NULL) throw GUException("XorGrad::setParams input view set is not of ViewExtension type. So gradient information is not available");
//     }
// 
//     //Determine the maximum number of points in a bodypart and free previous possible memory
//     unsigned int maxP=0;
//     for (unsigned int i=0;i<body->size();i++) {
//         if (maxP<body->getMesh(i).getVertices().size())
//             maxP=body->getMesh(i).getVertices().size();
//     }
//     //allocation of memory
//     //check all images are of same sime and allocate for only one
//     IplImage *mask0=_TheViewSet->get(0)->getMaskImage();
//     for (unsigned int v=1;v<_TheViewSet->size();v++) {
//         IplImage *mask=_TheViewSet->get(v)->getMaskImage();
//         if (mask->width!=mask0->width || mask->height!=mask0->height)
//             throw GUException("XorGrad::setParams images of different sizes not supported yet");
//     }
//     //allocate only for one image
//     _proj=cvCreateImage(cvSize(mask0->width,mask0->height),mask0->depth,mask0->nChannels);
// 
//     //allocate auxiliar vectors with the projection of the points
//     _invalid_v.resize(maxP*2);
//     _projections.resize(maxP*2);
//     _dotproducts.resize(maxP);
//     cvSetZero(_proj);
//     _areParams=true;
// 
// }
// 
// 
// 
// /*
// *
// *
// */
// double XorGrad::evaluate(BodyModel<HumanSkeleton> *body,ViewSet *vs) {
// 
//     assert(_areParams);
//     _nEvals++;
//     //    double result=1.; // 1. --> HIGH VALUE FOR FITNESS FUNCTION
//     double sumAvrgXor=0,sumAvrgGrad=0;
//     _TheViewSet=vs;
//     _bodyModel=body;
//     Point3Df camLoc;
//     int width=_TheViewSet->get(0)->getMaskImage()->width;
//     int height=_TheViewSet->get(0)->getMaskImage()->height;
//     gummocap::view::ViewExtension *viewE;
//     //Here the main loop, the point analysis
//     
//     for (unsigned int v=0;v<_TheViewSet->size();v++) {
//         CameraParameters *cp=_TheViewSet->get(v)->getCameraParameters  ();
//         viewE=(gummocap::view::ViewExtension *)_TheViewSet->get(v);
// 	IplImage *gradImg=viewE->getGradientImage();
// 	int sumGrad=0,nPointsGrad=0;
//         for (unsigned  part=0;part<_bodyModel->size();part++) {
// 	  
//             gumocap::core::svector<Point3D<float> > &vertices=_bodyModel->getMesh(part).getVertices();
//             gumocap::core::svector<Point3D<float> > &normals=_bodyModel->getMesh(part).getNormals();
//             //project the body part
//             cp->project3DPoints_4bytes_nodist((float*) vertices._data ,_projections._data,vertices._size,_invalid_v._data );
//             //calculate the points of the border
//             cp->getCameraLocation( camLoc.xyzh[0],camLoc.xyzh[1],camLoc.xyzh[2]);
// 
//             for (unsigned int n=0;n<normals.size();n++)
//                 _dotproducts[n]=normals[n].dot( (vertices[n]-camLoc).normalize()  );
// 
//             //ok, now for each point calcuale the value
//             int index=0;
//             for (unsigned point=0;point<vertices.size() ;point++,index+=2) {
//                 if (_invalid_v[index]==0   ) { //the point is  projected and visible
//                     //note : _projections[v][index] is the x projection and _projections[v][index+1] the y projection in the image
//                     int yo=_projections[index+1];
//                     int xo=_projections[index];
//                     //do only if not in a border
//                     unsigned char *pix=(unsigned char *)  _proj ->imageData+ _proj->widthStep * yo;
//                     if (xo>=1 && xo+1<width && yo>=1 && yo+1<height) {
//                         pix[xo]=125;
//                         pix[xo-1]=125;
//                         pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * (yo-1);
//                         pix[xo]=125;
//                         pix=(unsigned char *)  _proj->imageData+ _proj->widthStep * (yo+1);
//                         pix[xo]=125;
//                     }
//                     //now, gradient information around borders 
//                     if (  -_dotPThres<_dotproducts[point] &&  _dotproducts[point]<_dotPThres) {
//                         pix=(unsigned char *)  gradImg ->imageData+ gradImg->widthStep * yo;
// 			sumGrad+=pix[xo];
// 			nPointsGrad++;
//                     }
// 
//                 }
//             }
//         }
//  
// 
// 
//         int count=xor_special(_proj,_TheViewSet->get(v)->getMaskImage());
//         sumAvrgXor+=float(count)/float( width* height);
// 	sumAvrgGrad+=1.-(float(sumGrad)/float(255.*nPointsGrad));
// // 	cout<<"grad="<<float(sumGrad)/float(255.*nPointsGrad)<<endl;
// // 	       cvNamedWindow("proj");
// //         cvNamedWindow("grad");
// //         cvShowImage("proj", _proj);
// //         cvShowImage("grad", viewE->getGradientImage());
// //         cvWaitKey(0);
//     }
// 
// 
//     return (sumAvrgXor+sumAvrgGrad)/float(2*_TheViewSet->size())   ;
// 
// }
// 
// 
// /**
// */
// int XorGrad::xor_special(IplImage *im1,IplImage *im2)
// {
// //    if( !gu_isaligned( im1) ||!gu_isaligned( im2) )
// //    {
// //     cerr<<" XorGrad::xor_special not aligned images"<<__FILE__<<":"<<__LINE__<<endl;
// //     exit(0);
// //    }
// //    //yet, only for some sizes
// //    if( im1->width%16!=0)
// //    {
// //     cerr<<" XorGrad::xor_special not appropriated size for images"<<__FILE__<<":"<<__LINE__<<endl;
// //     exit(0);
// //    }
// 
// 
//     __m128i *end=(__m128i *)(im1->imageData+im1->widthStep*im1->height);
//     __m128i *ptr1=(__m128i *)im1->imageData;
//     __m128i *ptr2=(__m128i *)im2->imageData;
//     __m128i zeros=_mm_set_epi8( 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
//     __m128i ones=_mm_set_epi8( 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);
// 
//     unsigned char _aux[16] __attribute__((aligned(16)));
//     short int *_siaux=(short int *)&_aux;
//     __m128i *iaux=(__m128i *)&_aux;
//     int sum=0;
// 
//     for ( ; ptr1<end;ptr1++,ptr2++) {
//         //first load
// 
//         __m128i r1= _mm_load_si128(ptr1);
//         __m128i r2= _mm_load_si128(ptr2);
//         __m128i raux;
//         //does xor
//         raux=_mm_xor_si128(r1,r2);
//         //afte the two following instructions, the result is raux[i]==0? 0:1
//         raux=_mm_cmpeq_epi8(zeros,raux); // raux[i]==0? 0xfff:0x000
// 
//         raux=_mm_andnot_si128(raux,ones); //
//         //use this instruction to add the value of the 8 first and the second 8 uchar and convert them into 16 bits
//         raux= _mm_sad_epu8(raux,zeros);
//         _mm_store_si128(iaux,raux);//store and sum
//         sum=sum+_siaux[0]+_siaux[4];
//         _mm_store_si128(ptr1,zeros);//Set zero in the im1
//     }
//     return sum;
// 
// }
// }
// 
// }

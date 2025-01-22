#ifndef _GUPVM_RGColorHistogram_H
#define _GUPVM_RGColorHistogram_H

#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>
#include <gu/guexception.h>
#include <fstream>

using namespace std;
namespace gummocap
{
namespace poseevaluators
{
  #define _gummocap_RGCH_NDIV 6
#define _gummocap_RGCH_Size 36
#define _def_RG_ADD_RGB(this,R,G,B,val)   if (int(R+G+B)!=0){\
		    float nr=float(R)/255.;\
		    float ng=float(G)/255.;\
		    float nb=float(B)/255.;\
		    float Y=0.299*nr+0.587*ng+0.114*nb;\
		    float Cr=(nr-Y)+0.5;\
		    float Cb=(nb-Y)+0.5;\
		    int indexR=(  Cr*float(_gummocap_RGCH_NDIV-1)  ) +0.5;\
		    int indexB=(  Cb*float(_gummocap_RGCH_NDIV-1)  )+0.5;\
		    int indexF=indexR*_gummocap_RGCH_NDIV+indexB;\
		    this._data[indexF]+=val;\
		    }\
		    else this._data[0]+=val;

/** \brief
*This class represents a   color histogram in the YCrCb  normalized color space
*/
class RGColorHistogram
{
public:

 

    /**Types of comparison that can be made between histograms. Intersection, Bhattacharyya (needs the histograms to be normalized)
    */
    enum ComparisonType {Intersection,Bhattacharyya};
    /**Empty constructor
    */
    RGColorHistogram();


    /**Creates the histogram indicating the number of divisions
    */
    RGColorHistogram ( int nDivision );
    /**
    */
    ~RGColorHistogram();
    /**Copy constructor
    */
    RGColorHistogram ( const RGColorHistogram &CH );

    /**
      */
    RGColorHistogram & operator= ( const   RGColorHistogram &CH );



    /**Adds a color to the histogram.
    *@param r,g,b
    *@param val value to be added in the histogram
     */
    void addRGB ( int r,int g,int b ,float val=1 );
    /**Adds a color line to the histogram.
    *@param rgb pointer to scan line
    *@param ncolors number of pixel to add
    *@param val value to be added in the histogram
     */
    void addRGB ( unsigned char *rgb,int ncolors,float val=1 );

    /**Resets the information in the object
    */
    void reset();

    /**Number of subdivisions of the histogram in each axis of the color space
     */
    int getNDivisions() {
        return _gummocap_RGCH_NDIV;
    }


    /**Compares this histogram with the one passed using the method indicated.
    * Intersection: sum I  min(H1(I),H2(I)) => 1 equal, 0 differnt
    * Bhattacharyya: sqrt(1-sum I (sqrt(H1(I)*H2(I))))=>0 equal, 1 different
    */
    double compare ( RGColorHistogram &Model,ComparisonType type=Intersection );

    ///DEbug
    void print()
    {
        for ( int i=0;i<size();i++ ) cout<< _data [i]<<" ";
        cout<<endl;
    }

    /**Normalizes the histogram
    */
    void normalize();
    /**Updates this color histogram using the one passed as this=(1-innovation)*thia+innovation*CH
     *@param CH the new color histogram
     *@param innovation
     */
    void update ( RGColorHistogram &CH,double alpha );

    /**Adds the infromation of the  histogram passed to this
    */
    void add(const RGColorHistogram &hist);

    /** return a pointer to the data
     */
    float *getData() {
        return _data;
    }


    int size()const  {
        return _gummocap_RGCH_Size;
    }
// 		private:
    float _data[_gummocap_RGCH_Size];
};
};
}
#endif

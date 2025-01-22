/***************************************************************************
 *   Copyright (C) 2004 by salinas                                         *
 *   salinas@decsai.ugr.es                                                 *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef gummocap_PARTICLE_H
#define gummocap_PARTICLE_H
#include <vector>
#include <iostream>
using namespace std;
/**\brief This namespace contains the tools to implement condensation algorithm
*/
namespace gummocap{
namespace algorithms {

/**
\brief Represents a particle of the particle filter.
@author salinas
*/

class Particle:public vector<double>{
public:
   /** Constructor
    */
    Particle():vector<double>()
    {}
     /** Constructor
    */
    Particle(int n):vector<double>(n)
    {}
    

   /** Destructor
    */
    ~Particle()
    {
    }
    /**Copy constructor
    */
    Particle(const Particle&P):vector<double>(P)
    {
    prob=P.prob;
    accProb=P.accProb;
    }
    
    /**Assign operator
    */
    Particle &operator=(const Particle&P)
    {
    resize(P.size());
    copy(P.begin(),P.end(),begin());
    prob=P.prob;
    accProb=P.accProb;
    return *this;
    }


    /**Ostream Operator
    */
    friend ostream & operator<<(ostream &str,const Particle&P )
    {
     cout <<"|";
    	for(unsigned int i=0;i<P.size();i++)
    		str<<P[i]<<" ";
	str<<" prob="<<P.prob<<" acc="<<P.accProb<<endl;
	return str;
    }
    //vector<float> state;//state vector
    double prob; //probability of the particle
    double accProb;//accumulated probility
};

};
};
#endif

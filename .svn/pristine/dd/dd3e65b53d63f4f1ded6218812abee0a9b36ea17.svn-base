#include "humanparttracker.h"

namespace gummocap
{
namespace algorithms
{
  
/**
*/
HumanPartTracker::HumanPartTracker ()
{
    _enabledParams=false;
    
}
/**
*/
HumanPartTracker::~HumanPartTracker ()
{
  if(_enabledParams){
  }
}
/**Sets the required params and initializes
*/
void HumanPartTracker::setParams(unsigned int nBodyParts, unsigned int nbodyjoints, unsigned int nlarmjoints, unsigned int nrarmjoints, unsigned int nllegjoints, unsigned int nrlegjoints)
{
    _initrans=0;
    _endtrans=2; //Translation
    _inirot=3;
    _endrot=5;   //Rotation
    _inibody=6;
    _endbody=_inibody+(3*nbodyjoints)-1;
    _inilarm=_endbody+1;
    _endlarm=_inilarm+(3*nlarmjoints)-1;
    _inirarm=_endlarm+1;
    _endrarm=_inirarm+(3*nrarmjoints)-1;
    _inilleg=_endrarm+1;
    _endlleg=_inilleg+(3*nllegjoints)-1;
    _inirleg=_endlleg+1;
    _endrleg=_inirleg+(3*nrlegjoints)-1;

    _BODY_PARTS=nBodyParts;
    if(nBodyParts==6)
      setIniEndBodyParts6();
    else
      setIniEndBodyParts12();
    
    _enabledParams=true;
}
/**
  *
  *
  */
void HumanPartTracker::setIniEndBodyParts6(){

    _ini_end_body_parts.push_back(_initrans);
    _ini_end_body_parts.push_back(_endrot);
    _ini_end_body_parts.push_back(_inibody);
    _ini_end_body_parts.push_back(_endbody);
    _ini_end_body_parts.push_back(_inilarm);
    _ini_end_body_parts.push_back(_endlarm);
    _ini_end_body_parts.push_back(_inirarm);
    _ini_end_body_parts.push_back(_endrarm);
    _ini_end_body_parts.push_back(_inilleg);
    _ini_end_body_parts.push_back(_endlleg);
    _ini_end_body_parts.push_back(_inirleg);
    _ini_end_body_parts.push_back(_endrleg);


}
/**
  *
  *
  */
void HumanPartTracker::setIniEndBodyParts12(){
  
  // BODY_P
     _ini_end_body_parts.push_back(_initrans);
     _ini_end_body_parts.push_back(_endtrans);
     // BODY_O
     _ini_end_body_parts.push_back(_inirot);
     _ini_end_body_parts.push_back(_endrot);
     // TORSO
     // * 2 BODY JOINTS: ONLY SPINE_BOTTOM
     // * 3 BODY JOINTS: SPINE_BOTTOM + SPINE_UP
     // * Not NECK
     _ini_end_body_parts.push_back(_inibody);
     _ini_end_body_parts.push_back(_endbody-3);
     // LU_ARM (LEFT UPPER ARM)
     // * LEFT SHOULDER
     _ini_end_body_parts.push_back(_inilarm);
     _ini_end_body_parts.push_back(_inilarm+2);	
     // LL_ARM (LEFT LOWER ARM)
     // * 2 LEFT-ARM JOINTS: ONLY LEFT ELBOW
     // * 3 LEFT-ARM JOINTS: LEFT ELBOW + LEFT WRIST
     _ini_end_body_parts.push_back(_inilarm+3);
     _ini_end_body_parts.push_back(_endlarm);	
     // RU_ARM (RIGHT UPPER ARM)
     // * RIGHT SHOULDER
     _ini_end_body_parts.push_back(_inirarm);
     _ini_end_body_parts.push_back(_inirarm+2);	
     // RL_ARM (RIGHT LOWER ARM)
     // * 2 RIGHT-ARM JOINTS: ONLY RIGHT ELBOW
     // * 3 RIGHT-ARM JOINTS: RIGHT ELBOW + RIGHT WRIST
     _ini_end_body_parts.push_back(_inirarm+3);
     _ini_end_body_parts.push_back(_endrarm);	
     // HEAD
     // * NECK JOINT
     _ini_end_body_parts.push_back(_endbody-2);
     _ini_end_body_parts.push_back(_endbody);
     // LU_LEG (LEFT UPPER LEG)
     // * LEFT HIP
     _ini_end_body_parts.push_back(_inilleg);
     _ini_end_body_parts.push_back(_inilleg+2);	
     // LL_LEG (LEFT LOWER LEG)
     // * 2 LEFT-LEG JOINTS: ONLY LEFT KNEE
     // * 3 LEFT-LEG JOINTS: LEFT KNEE + LEFT ANKLE
     _ini_end_body_parts.push_back(_inilleg+3);
     _ini_end_body_parts.push_back(_endlleg);	
     // RU_LEG (RIGHT UPPER LEG)
     // * RIGHT HIP
     _ini_end_body_parts.push_back(_inirleg);
     _ini_end_body_parts.push_back(_inirleg+2);	
     // RL_LEG (RIGHT LOWER LEG)
     // * 2 RIGHT-LEG JOINTS: ONLY RIGHT KNEE
     // * 3 RIGHT-LEG JOINTS: RIGHT KNEE + RIGHT ANKLE
     _ini_end_body_parts.push_back(_inirleg+3);
     _ini_end_body_parts.push_back(_endrleg);	
     
}
 
 
 
}

}





//////////////////////////////////////////
// Read in the list of 34 ternary forms //
//////////////////////////////////////////

void Read34Ternaries(long QQ[34][6]);


////////////////////////////////////////////////////////////////////////
// Makes any of the 34 ternary theta functions with precision <= 10^7 //
////////////////////////////////////////////////////////////////////////

void MakeOneThetaFrom34Ternaries();


// Look for the exceptions < 10^7 for the 34 ternaries
void FindTernaryExceptions_by_theta();


// Look for the exceptions < 10^7 for the 34 ternaries
void FindTernaryExceptions_New(); 


//////////////////////////////////////////////////////////////////
// Identifies the upper-left ternary in each of our 6560 forms, //
// and describes those whose ternary is not regular.            //
//////////////////////////////////////////////////////////////////

void Find_4var_with_Irreg_Ternaries(); 



// ---------------------------------------------------

#if !defined(SCRAPS3_H)
#define SCRAPS3_H

#include "scraps3.cc"

#endif

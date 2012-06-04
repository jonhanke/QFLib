

#include <gmp.h>
#include <gmpxx.h>

#include <valarray>

//#include "vec_mat.h"
//#include "vectors.h"


// Notes:
//    - Sets will mean vectors assumed to be in increasing order, with 
//        no duplicates (i.e. strictly increasing).
//    - Routines ending in the word "Ordered" assume that the input is an 
//        ordered (i.e. strictly increasing) vector, also referred to as a set.  
//        The entries should also be positive integers, but this may not be 
//        enforced (I need to look).
//
//    - Types of returns:
//      -----------------
//        Tuples:
//            VectorAppend - append as tuples 
//            MakeVector - makes a tuple
//            VectorRemove - makes a tuple
//        Sets:
//            OrderVector - makes a set from a tuple
//            VectorUnion - makes a set
//            VectorIntersection - makes a set
//            VectorComplement - makes a set
//        Boolean (long):
//            IsDisjointOrdered - 0 or 1
//            CheckVectorIntersection - 0 or 1
//            CheckVectorOrdered - 0 or 1




//////////////////////////////////////////////////////////////////////
// Checks if two ordered (strictly increasing) vectors are disjoint //
//////////////////////////////////////////////////////////////////////

bool IsDisjointOrdered(const valarray<size_t> & v, const valarray<size_t> & w);


/////////////////////////////////////////////////
// Gives a vector by an arithmetic progression //
/////////////////////////////////////////////////

valarray<size_t> MakeVector(size_t length, size_t start, size_t increment);


//////////////////////////////////////////////////////////////////////////////////
// Reorders the vector v in strictly increasing order (so no repetitions occur) //
//////////////////////////////////////////////////////////////////////////////////

valarray<size_t> OrderVector(const valarray<size_t> & v);


//////////////////////////////////////////////
// Creates the union (v + w) of two vectors //
//////////////////////////////////////////////

valarray<size_t> VectorAppend(const valarray<size_t> & v, const valarray<size_t> & w);


//////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) union of two vectors //
//////////////////////////////////////////////////////////////////
valarray<size_t> VectorUnion(const valarray<size_t> & v, const valarray<size_t> & w);


/////////////////////////////////////////////////////////////////////////////////////
// Finds their ordered union, assuming they are both ordered (strictly increasing) //
/////////////////////////////////////////////////////////////////////////////////////
valarray<size_t> VectorUnionOrdered(const valarray<size_t> & v, const valarray<size_t> & w);


/////////////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) intersection of two vectors //
/////////////////////////////////////////////////////////////////////////
valarray<size_t> VectorIntersection(const valarray<size_t> & v, const valarray<size_t> & w);


/////////////////////////////////////////////////////////////////////////////////
// Finds the ordered (strictly increasing) intersection of two ordered vectors //
/////////////////////////////////////////////////////////////////////////////////
valarray<size_t> VectorIntersectionOrdered(const valarray<size_t> & v, const valarray<size_t> & w);


////////////////////////////////////
// Remove the entries in w from v //
////////////////////////////////////
valarray<size_t> VectorRemove(const valarray<size_t> & v, const valarray<size_t> & w);


//////////////////////////////////////////////////////////////
// Returns an ordered vector with the entries of v not in w //
//////////////////////////////////////////////////////////////
valarray<size_t> VectorComplement(const valarray<size_t> & v, const valarray<size_t> & w);


///////////////////////////////////////////////////////////////
// Returns an ordered vector with the entries of v not in w, //
// assuming both v and w are ordered (strictly increasing)   //
///////////////////////////////////////////////////////////////
valarray<size_t> VectorComplementOrdered(const valarray<size_t> & v, const valarray<size_t> & w);


/////////////////////////////////////////
// Checks if the two vectors intersect //
/////////////////////////////////////////
bool CheckVectorIntersection(const valarray<size_t> & v, const valarray<size_t> & w);


///////////////////////////////////////////////////
// Checks if the two (ordered) vectors intersect //
///////////////////////////////////////////////////
bool CheckVectorIntersectionOrdered(const valarray<size_t> & v, const valarray<size_t> & w);


/////////////////////////////////////////////////
// Checks if v is a strictly increasing vector //
/////////////////////////////////////////////////
bool CheckVectorOrdered(const valarray<size_t> & v);



//////////////////////////////////////////
// Finds the GCD of a vector of numbers //
//////////////////////////////////////////
mpz_class GCD(const valarray<mpz_class> & v);




// ---------------------------------------------

#if !defined(VECTORS_H)
#define VECTORS_H


#include "valarray_extras.h"  


#include "vectors.cc"


#endif




// inversewishart.cpp: implementation of the InverseWishart class.
//
//////////////////////////////////////////////////////////////////////
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines

#include "newran.h"
#include "tnt.h"
using namespace TNT;
#include "stdafx.h"

#include "inversewishart.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

InverseWishart::InverseWishart()
{
	mt0 = 0;		//t0
	ms = 0;			//number of factors
	mpS0 = NULL;	//a pointer to the initial s0 matrix
	mpThetaSquare = NULL;	//a pointer to the theta * transpose(theta) matrix
}

//clearn up function
InverseWishart::~InverseWishart()
{
	if (mpS0 != NULL) {
		delete mpS0; mpS0 = NULL;
	}
	if (mpThetaSquare != NULL) {
		delete mpThetaSquare; mpThetaSquare = NULL;
	}
}

//initialize the InverseWishart sampler
//t0: the t0 parameter in W(t0,S0) settings
//s: number of distintictive theta's
//k: the number of factors
//S0: initial S0 matrix, S0 = b/a * t0 * I_k
void InverseWishart::Init(int t0, int s, int k, SymmetricMatrix & S0, SymmetricMatrix & ThetaSquare) {
	mt0 = t0;
	ms = s;
	mk = k;
	//calculate mB such that mB * mBt = inverse((s0 + theta * transpose(theta)))
	mpS0 = new SymmetricMatrix(1); *mpS0 = S0;		
	mpThetaSquare = new SymmetricMatrix(1); *mpThetaSquare = ThetaSquare; 
	SymmetricMatrix BSSum = *mpS0 + *mpThetaSquare;	
	SymmetricMatrix BSSumInv = BSSum.i();
	mB = Cholesky(BSSumInv);
	mBt = mB.t();
}

//sample an IW matix from W(t0,S0)
void InverseWishart::Next(Array2D<double> & RT) {
	int i,j;

	LowerTriangularMatrix U(mk);
	for (i = 0; i < mk; i++) {
		ChiSq cs(mt0 + ms - i);
		U[i][i] = sqrt(cs.Next());
		for ( j = 0; j < i; j++) {
			U[i][j] = mNormal.Next(); 
		}
	}
	Matrix K = mB * U * U.t() * mBt;			
	Matrix T = K.i();
	//RT is returned
	for (i = 0; i < mk; i++) {
		for (j = 0; j< mk; j++) {
			RT[i][j] = T[i][j];
		}
	}
	//std::cout << "det = " << T.Determinant() << endl;
}


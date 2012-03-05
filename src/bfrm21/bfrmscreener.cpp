// bfrmscreener.cpp: implementation of the BfrmScreener class.
//
//////////////////////////////////////////////////////////////////////

#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines

#include "newran.h"
#include "tnt.h"
using namespace TNT;

#include <math.h>
#include "stdafx.h"
#include "mdp.h"
#include "Model.h"
#include "bfrmrawdata.h"
#include "bfrmprior.h"

#include "bfrmresult.h"
#include "bfrm.h"
#include "util.h"


#include "bfrmscreener.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BfrmScreener::BfrmScreener()
{

}

BfrmScreener::~BfrmScreener()
{

}

void BfrmScreener::LoadData(Model& model) {
	std::cout << "BfrmScreener::LoadData: This function should not be called" << endl;
	exit (1);
}

void BfrmScreener::Run(Model& model) {
	//std::cout << "Entering variable selection mode" << endl;
	int NonGaussian = model.mnNonGaussianFactors;
	int nExtra = NonGaussian?50:0;
	model.mnNonGaussianFactors = 0; //force to be gaussian first
	for (int it= 1; it <= model.mnMCSelect + model.mnBurninSelect + nExtra; it++) {
		if (NonGaussian && it > nExtra ) {
			model.mnNonGaussianFactors = NonGaussian;
		}
		Iterate(model.mnNLatentFactors, model.mnNDesignControls, model.mnPriorPiStandard, 
			1,model.mdDPAlpha, model.mnNonGaussianFactors, model.mnInverseWishart, it, model.mnPriorG);
		if(it% model.mnPrintIteration == 0 || model.mnPrintIteration < 2) { 
			std::cout << "MCMC iteration " << it << endl;
			std::cout << "Sparsity " << mdSparsity << endl;
		}
		if (it>model.mnBurninSelect + nExtra) {
			Accumulate(model.mnNLatentFactors, model.mnNDesignControls, model.mnPriorPiStandard);
		}
	}
	Summarise(model.mnNLatentFactors, model.mnNDesignControls ,model.mnMCSelect, model.mnPriorPiStandard,1 );
}


int BfrmScreener::SetData(Array2D<double>& AllData, Array2D<double>& AllMask, Array1D<int>& indicator,
			 Array2D<double>& FH, Array1D<double>& Weight, int yfactors, 
			 Array1D<int>& VariablesIn, int nVarIn, bool bHasXMask) {
	int i,j;

	maOriginalIndex.clear();
	//X first
	mnVariables = indicator.dim() + yfactors;
	mbHasXMask = bHasXMask;
	mnSampleSize = FH.dim2();
	mX = Array2D<double>(mnVariables, mnSampleSize);
	if (mbHasXMask) {
		mXMask = Array2D<double>(mnVariables - yfactors, mnSampleSize);
	}
	for (j = 0; j < mnSampleSize; j++) {
		for (i = 0; i < yfactors; i++) {
			mX[i][j] = AllData[i][j];
		}
		for (i = 0; i < nVarIn; i++) {
			mX[i+yfactors][j] = AllData[yfactors + VariablesIn[i]][j];
		}
		if (mbHasXMask) {
			for (i = 0; i < nVarIn; i++) {
				mXMask[i][j] = AllMask[VariablesIn[i]][j];
			}
		}
	}
	for (i = 0; i < nVarIn; i++) {
		maOriginalIndex.push_back(VariablesIn[i]);
	}
	
	int count = nVarIn + yfactors;
	for (i = 0; i < indicator.dim(); i++) {
		if (indicator[i]) {
			for (j = 0; j < mnSampleSize; j++) {
				mX[count][j] = AllData[yfactors + i][j];
			}
			if (mbHasXMask) {
				for (j = 0; j < mnSampleSize; j++) {
					mXMask[count-yfactors][j] = AllMask[i][j];
				}
			}
			count++;
			maOriginalIndex.push_back(i);
		}
	}

	//mH
	mH = Array2D<double>(FH.dim2(), FH.dim1());
	for (i = 0; i < FH.dim1(); i++) {
		for (j = 0; j < FH.dim2(); j++) {
			mH[j][i] = FH[i][j];
		}
	}
	//mWeight
	mnLeaveout = 0;
	mWeight = Array1D<double>(Weight.dim());
	for (i = 0; i < mWeight.dim(); i++) {
		mWeight[i] = Weight[i];
		if (mWeight[i] == 0) { mnLeaveout++; }
	}
	return 1;
}

void BfrmScreener::Initialise(Model& model, Bfrm &fitted)
{
	int i, j;
	int design = model.mnNDesignControls;
	int nTotalFactor = design + model.mnNLatentFactors + mnYFactors;
	mPrior.Initialise(model, mnYFactors, mnVariables);
	mPrior.InitialisePsi(mIter.mPsi, mnVariables, mnYFactors - mnVarSurvival - mnVarContinuous, mnYFactors - mnVarContinuous );
	mDP.init(model.mnNObservations,model.mdDPAlpha,model.mdDPAlphaX, model.mdaPriorAlpha[0], model.mdaPriorAlpha[1]);
	//carry over Psi's
	for (i = 0; i < fitted.mMeanResult.mPsi.dim(); i++) {
		mIter.mPsi[i] = fitted.mMeanResult.mPsi[i];
	}
	mPrior.InitialiseTau(mIter.maTau, nTotalFactor);
	mPrior.InitialiseRho(mIter.maRho, nTotalFactor);
	
	if (model.mnPriorPiStandard) {
		mPrior.InitialisePibIJ(mIter.maPi, nTotalFactor,mnVariables);
		//carry over pib2's
		for (i = 0; i < fitted.mMeanResult.maPi.dim1(); i++) {
			for (j = 0; j < fitted.mMeanResult.maPi.dim2(); j++) {
				mIter.maPi[i][j] = fitted.mMeanResult.maPi[i][j];
			}
		}
	}

	mIter.mPostPib = Array2D<double>(mnVariables, nTotalFactor); mIter.mPostPib = 0;
	if (model.mnNDesignControls) {			//only do this if there is an intercept factor
		for (i = 0; i < mnVariables; i++) {
			mIter.mPostPib[i][0] = 1.0;			
		}
	}
	//carry over mPostPib's
	for (i = 0; i < fitted.mMeanResult.mPostPib.dim1(); i++) {
		for (j = 0; j < fitted.mMeanResult.mPostPib.dim2(); j++) {
			mIter.mPostPib[i][j] = fitted.mMeanResult.mPostPib[i][j];
		}
	}

	if (model.mnInverseWishart) {
		mPrior.InitInverseWishart(nTotalFactor - model.mnNDesignControls, model.mdInverseWishartT0);
		mPrior.InitInverseWishart(nTotalFactor - model.mnNDesignControls, mIter.mT);
	}

	SetBMasks(model.mnShapeOfB, 0);
	
	InitiliseX();		
	//we also borrow the imputed values from the fitted model
	if (mbHasXMask) {			//fill the missings with fitted
		for ( i = 0; i < fitted.mMeanResult.mX.dim1(); i++) { 
			for (j = 0; j < fitted.mMeanResult.mX.dim2(); j++) {
				if (mXMask[i][j] > 0) {
					mX[i + mnYFactors][j] = fitted.mMeanResult.mX[i][j];
				}
			}
		}
	}

	Array2D<double> Mean = RowMeanXY();
	//initilize B matrix
	int ntotalzeros = 0;
	mIter.mB = Array2D<double>(mnVariables, nTotalFactor); mIter.mB = 0;
	int nVar = fitted.mMeanResult.mB.dim1();
	for (i = 0; i < mnVariables; i++) {
		for (j = 0; j < nTotalFactor; j++) {
			if (i < nVar) {
				mIter.mB[i][j] = fitted.mMeanResult.mB[i][j];
			} else {
				if (design && j == 0) {
					mIter.mB[i][0] = Mean[i][0];
				} else {
					int nMask = GetBMask(i,j,model.mnNDesignControls, 0);
					if (nMask == 0) {
						mIter.mB[i][j] = 0;
						ntotalzeros++;
					} else if (nMask == 1) {
						//July 29 2005, use carlos's initialisation instead 
						mIter.mB[i][j] = 1.0; //mPrior.SampleRestrictedNormal(0,1,1,0);
					} else {
						//July 29 2005, use carlos's initialisation instead
						mIter.mB[i][j] = 0.0;
					}
				}
			}
		}
	}
	
	//initilise  F matrix
	mIter.mF = Array2D<double>(nTotalFactor, mnSampleSize);
	//design part
	for (i = 0; i < nTotalFactor; i++) {		//not design here
		for (j = 0; j < mnSampleSize; j++) {
			mIter.mF[i][j] = mWeight[j] ? mH[j][i] : 0;			
		}
	}
	mIter.mAlpha = Array1D<double>(1);
	//initilise residules, part 1: X 
	mResX = mX - matmult(mIter.mB,mIter.mF);
	if (mnLeaveout != 0) {			//if there are leave-out samples
		for (j = 0; j < mnSampleSize; j++) {
			if (mWeight[j] == 0) {
				for (i = 0; i < mnVariables; i++) {
					mResX[i][j] = 0;
				}
			}
		}
	}
	
	//estimate tau
	mIter.maTau = BfrmData::Var(mIter.mF);
	for (i = 0; i < nTotalFactor; i++) {
		if (!mIter.maTau[i]) { mIter.maTau[i] = 1e-9; }
	}

	if (model.mnInverseWishart) {
		for (i = design; i < nTotalFactor; i++) {
			mIter.maTau[i] = 1.0;
		}
	}
	//addjust F and B so F ~ N(0,1)
	Array2D<double> T(nTotalFactor,nTotalFactor); T = 0; T[0][0] = 1.0;
	for (i = 1; i < nTotalFactor; i++) {
		T[i][i] = sqrt(mIter.maTau[i]);
	}
	mIter.mB = matmult(mIter.mB, T);
	
	for (i = 1; i < nTotalFactor; i++) {
		T[i][i] = 1.0 / sqrt(mIter.maTau[i]);
	}
	mIter.mF = matmult(T, mIter.mF);
	
	//now restore mpib & mtau
	for (i = 0; i < nTotalFactor; i++) {
		mIter.maTau[i] = fitted.mMeanResult.maTau[i];
	}
	if (model.mnInverseWishart) {
		for (i = design; i < nTotalFactor; i++) {
			mIter.maTau[i] = 1.0;
		}
	}

	ManageHelperVariables(true,model.mnNLatentFactors);

	mResult.Init(mnVariables, nTotalFactor, mnSampleSize, mnYFactors, model.mnPriorPiStandard,mbHasXMask, model.mnNDesignControls);
	mMeanResult.Init(mnVariables, nTotalFactor, mnSampleSize, mnYFactors, model.mnPriorPiStandard,mbHasXMask, model.mnNDesignControls);
	if (mbHasXMask) {
		mIter.mX = Array2D<double>(mnVariables - mnYFactors, mnSampleSize); mIter.mX = 0;
	}
}


bool BfrmScreener::InitResponse(int nbinary, int nCategorical, int nSurvival, int nContinuous, int nExtra) {
	mnExtraVariables = nExtra;
	mnVarBinary = nbinary;
	mnVarCategorical = nCategorical;
	mnVarSurvival = nSurvival;
	mnVarContinuous = nContinuous;
	mnYFactors = mnVarBinary + mnVarCategorical + mnExtraVariables + mnVarSurvival + mnVarContinuous;
	return true;
}

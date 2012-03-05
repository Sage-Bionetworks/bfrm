// bfrmscreener.h: interface for the BfrmScreener class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BFRMSCREENER_H__41008F2C_2B05_4E9A_BA9F_0DCFC9F7AB05__INCLUDED_)
#define AFX_BFRMSCREENER_H__41008F2C_2B05_4E9A_BA9F_0DCFC9F7AB05__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "bfrm.h"

class BfrmScreener : public Bfrm  
{
public:
	BfrmScreener();
	virtual ~BfrmScreener();
	virtual void LoadData(Model& model);
	virtual void Initialise(Model& model, Bfrm &fitted);
	int SetData(Array2D<double>& AllData, Array2D<double>& AllMask, Array1D<int>& indicator,
		Array2D<double>& FH, Array1D<double>& Weight, int yfactors, Array1D<int>& VariablesIn, 
		int nVarIn, bool bHasXMask);
	void Run(Model& model);
	bool InitResponse(int nbinary, int nCategorical, int nSurvival, int nContinuous, int nExtra);	
};

#endif // !defined(AFX_BFRMSCREENER_H__41008F2C_2B05_4E9A_BA9F_0DCFC9F7AB05__INCLUDED_)

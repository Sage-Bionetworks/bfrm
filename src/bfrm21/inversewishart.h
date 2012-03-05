// inversewishart.h: interface for the InverseWishart class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_INVERSEWISHART_H__980F9A29_6C48_4332_84EA_CBAC788C8B65__INCLUDED_)
#define AFX_INVERSEWISHART_H__980F9A29_6C48_4332_84EA_CBAC788C8B65__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class InverseWishart  
{
public:
	InverseWishart();
	virtual ~InverseWishart();
	int mt0;
	int ms;
	int mk;
	SymmetricMatrix* mpS0;
	SymmetricMatrix* mpThetaSquare;
	LowerTriangularMatrix mB;
	UpperTriangularMatrix mBt;
	Normal mNormal;
	void Init(int t0, int s, int k, SymmetricMatrix & S0, SymmetricMatrix & ThetaSquare);
	void Next(Array2D<double> & RT);
};

#endif // !defined(AFX_INVERSEWISHART_H__980F9A29_6C48_4332_84EA_CBAC788C8B65__INCLUDED_)

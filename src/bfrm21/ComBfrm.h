// ComBfrm.h : Declaration of the CComBfrm

#ifndef __COMBFRM_H_
#define __COMBFRM_H_
#include "newran.h"
#include "tnt.h"
#include "jama_svd.h"
#include "jama_cholesky.h"
using namespace TNT;
using namespace JAMA;
#include "Model.h"
#include "bfrm.h"	// Added by ClassView
#include "resource.h"       // main symbols
#include "bfrmcp.h"

/////////////////////////////////////////////////////////////////////////////
// CComBfrm
class Bfrm;
class ATL_NO_VTABLE CComBfrm : 
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CComBfrm, &CLSID_ComBfrm>,
	public IDispatchImpl<IComBfrm, &IID_IComBfrm, &LIBID_Bfrm>,
	public CProxy_IBfrmConnectorEvents< CComBfrm >,
	public IConnectionPointContainerImpl<CComBfrm>
{
public:
	CComBfrm()
	{
		mnDoneInit = 0;
		mbCheckPoint = false;
		mbStop = false;
//		mbDoneSummary = false;
		mbUseWeight = false;
		mbSilent = false;
	}

DECLARE_REGISTRY_RESOURCEID(IDR_COMBFRM)

DECLARE_PROTECT_FINAL_CONSTRUCT()

BEGIN_COM_MAP(CComBfrm)
	COM_INTERFACE_ENTRY(IComBfrm)
	COM_INTERFACE_ENTRY(IDispatch)
	COM_INTERFACE_ENTRY_IMPL(IConnectionPointContainer)
END_COM_MAP()

// IComBfrm
public:
	STDMETHOD(get_WeightFile)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_WeightFile)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_DesignFile)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_DesignFile)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_DataFile)(/*[out, retval]*/ BSTR *pVal);
	STDMETHOD(put_DataFile)(/*[in]*/ BSTR newVal);
	STDMETHOD(get_PriorTaub3)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorTaub3)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorTaua3)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorTaua3)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorTaub2)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorTaub2)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorTaua2)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorTaua2)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorTaub1)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorTaub1)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorTaua1)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorTaua1)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorPiMultiplier)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorPiMultiplier)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorPiProportion)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorPiProportion)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorPsib)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorPsib)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorPsia)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorPsia)(/*[in]*/ double newVal);
	STDMETHOD(get_PriorMeanValue)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_PriorMeanValue)(/*[in]*/ double newVal);
	STDMETHOD(get_CheckPoint)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_CheckPoint)(/*[in]*/ long newVal);
	STDMETHOD(get_Designs)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Designs)(/*[in]*/ long newVal);
	STDMETHOD(SetWeight)(VARIANT data);
	STDMETHOD(SetDesign)(VARIANT data);
	STDMETHOD(LoadParameters)(BSTR FileName);
	STDMETHOD(get_F)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_Bz)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_B)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(LoadData)();
	STDMETHOD(Run)();
	STDMETHOD(Summarise)();
	STDMETHOD(Accumulate)();
	STDMETHOD(Iterate)();
	STDMETHOD(Initialise)();
	STDMETHOD(SetX)(VARIANT data);
	STDMETHOD(AddHint)(BSTR key, BSTR value);
	STDMETHOD(get_Burnin)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Burnin)(/*[in]*/ long newVal);
	STDMETHOD(get_Factors)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Factors)(/*[in]*/ long newVal);
	STDMETHOD(get_Samples)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Samples)(/*[in]*/ long newVal);
	STDMETHOD(get_Variables)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Variables)(/*[in]*/ long newVal);
	STDMETHOD(get_Iterations)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Iterations)(/*[in]*/ long newVal);
protected:
	//int mnmc;
	//int mnit;
	int mnDoneInit;
	int mbCheckPoint;
	int mbStop;
//	int mbDoneSummary;
	int mbUseWeight;
	int mbSilent;
private:
	Bfrm T;
	Model model;
	BSTR toBSTR(const char* msg);
	
public :
	STDMETHOD(get_Mask)(/*[out, retval]*/ BOOL *pVal);
	STDMETHOD(put_Mask)(/*[in]*/ BOOL newVal);
	STDMETHOD(get_BayesFactorSelect)(double threshold, /*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_Sparsity)(/*[out, retval]*/ double *pVal);
	STDMETHOD(put_Sparsity)(/*[in]*/ double newVal);
	STDMETHOD(get_Silent)(/*[out, retval]*/ BOOL *pVal);
	STDMETHOD(put_Silent)(/*[in]*/ BOOL newVal);
	STDMETHOD(Resume)();
	STDMETHOD(get_Print)(/*[out, retval]*/ long *pVal);
	STDMETHOD(put_Print)(/*[in]*/ long newVal);
	STDMETHOD(get_StopAt)(/*[out, retval]*/ BOOL *pVal);
	STDMETHOD(put_StopAt)(/*[in]*/ BOOL newVal);
	STDMETHOD(get_UseCheckPoint)(/*[out, retval]*/ BOOL *pVal);
	STDMETHOD(put_UseCheckPoint)(/*[in]*/ BOOL newVal);
	STDMETHOD(get_PostPib)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_Pib)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_Tau)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(get_Psi)(/*[out, retval]*/ VARIANT *pVal);
	STDMETHOD(SaveParameters)(BSTR FileName);

BEGIN_CONNECTION_POINT_MAP(CComBfrm)
	CONNECTION_POINT_ENTRY(DIID__IBfrmConnectorEvents)
END_CONNECTION_POINT_MAP()

};

#endif //__COMBFRM_H_

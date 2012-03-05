// ComBfrm.cpp : Implementation of CComBfrm
#include "stdafx.h"
#include "bfrm_i.h"
#include "Bfrmc.h"
#include "ComBfrm.h"

/////////////////////////////////////////////////////////////////////////////
// CComBfrm

STDMETHODIMP CComBfrm::get_Iterations(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.nmc;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_Iterations(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	model.nmc = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Variables(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.NVARIABLES ;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_Variables(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	model.NVARIABLES = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Samples(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.NOBSERVATIONS;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_Samples(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	if (model.NOBSERVATIONS != newVal) {
		model.NOBSERVATIONS = newVal;
		T.InitialiseWeight("", model.NOBSERVATIONS);
		if (mbUseWeight) {
			AfxMessageBox("Sample size has changed and its weight indicator needs to be reset!");
		}
	}
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Factors(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.NFACTORS;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_Factors(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	model.NFACTORS = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Burnin(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.burnin;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_Burnin(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	model.burnin = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::AddHint(BSTR key, BSTR value)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	// TODO: Add your implementation code here
	return S_OK;
}

BSTR CComBfrm::toBSTR(const char* msg) {
	return CString(msg).AllocSysString();
}
STDMETHODIMP CComBfrm::LoadParameters(BSTR FileName)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	CString s =FileName;
	string ss = s;
	if (ss == "") {
		ss = "parameters.txt";
		Fire_Progress(toBSTR("default parameters file parameters.txt is used\n"));
	}
	if (!model.Load(ss)) {
		Fire_Progress(toBSTR(model.GetErrorMessage().c_str()));
		if (!mbSilent) {
			AfxMessageBox(model.GetErrorMessage().c_str());
		}
		return S_FALSE;
	}
	T.InitialiseH("", 1, 0);   //dummy way to initalise H
	T.InitialiseWeight("",model.NOBSERVATIONS);
	return S_OK;
}

STDMETHODIMP CComBfrm::SetX(VARIANT data)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	SAFEARRAY *psa = data.parray;
    if (psa->cDims != 2) {
		return S_OK;
	}
	int nSize1 = psa->rgsabound[1].cElements;
	int nSize2 = psa->rgsabound[0].cElements;
	long index[2];
	double *d = new double[nSize1 * nSize2];
	for (long j =0; j<nSize2; j++)
	{
		index[1] = j;
		for (long i = 0;  i< nSize1; i++) {
			index[0] = i;
			SafeArrayGetElement(psa, index, &d[j * nSize1 + i]);
		}
	}
	
	if (!T.LoadData(d,model.NVARIABLES,model.NOBSERVATIONS)) {
		Fire_Progress(toBSTR("Loading X matrix failed!\n"));
		if (!mbSilent) {
			AfxMessageBox("Loading X matrix failed!");
		}
	}
	delete [] d;
	return S_OK;
}

STDMETHODIMP CComBfrm::SetDesign(VARIANT data)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	SAFEARRAY *psa = data.parray;
    if (psa->cDims != 2) {
		return S_OK;
	}
	int nSize1 = psa->rgsabound[1].cElements;
	int nSize2 = psa->rgsabound[0].cElements;
	long index[2];
	double *d = new double[nSize1 * nSize2];
	for (long j =0; j<nSize2; j++)
	{
		index[1] = j;
		for (long i = 0;  i< nSize1; i++) {
			index[0] = i;
			SafeArrayGetElement(psa, index, &d[j * nSize1 + i]);
		}
	}
	if (!T.InitialiseH(d, model.NOBSERVATIONS, model.NDESIGNS)) {
		Fire_Progress(toBSTR("Loading design matrix failed!\n"));
		if (!mbSilent) {
			AfxMessageBox("Loading design matrix failed!");
		}
	}
	delete [] d;
	return S_OK;
}

STDMETHODIMP CComBfrm::SetWeight(VARIANT data)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	SAFEARRAY *psa = data.parray;
	if (psa->cDims != 2) {
		return S_OK;
	}
	int nSize1 = psa->rgsabound[1].cElements;
	int nSize2 = psa->rgsabound[0].cElements;
	long index[2];
	double *d = new double[nSize1 * nSize2];
	for (long j =0; j<nSize2; j++)
	{
		index[1] = j;
		for (long i = 0;  i< nSize1; i++) {
			index[0] = i;
			SafeArrayGetElement(psa, index, &d[j * nSize1 + i]);
		}
	}

	if (!T.InitialiseWeight(d,model.NOBSERVATIONS)) {
		Fire_Progress(toBSTR("Loading weight matrix failed!\n"));
		if (!mbSilent) {
			AfxMessageBox("Loading weight matrix failed!");
		}
	}
	mbUseWeight = true;
	return S_OK;
}

STDMETHODIMP CComBfrm::LoadData()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	
	//read data
	if (!T.InitialiseH(model.HFILE,model.NOBSERVATIONS, model.NDESIGNS)) {
		Fire_Progress(toBSTR("Loading design matrix failed!\n"));
		if (!mbSilent) {
			AfxMessageBox("Loading design matrix failed!");
		}
		return S_FALSE;
	}
	
	if (!T.InitialiseWeight(model.WEIGHTFILE,model.NOBSERVATIONS)) {
		Fire_Progress(toBSTR("Loading weight matrix failed!\n"));
		if (!mbSilent) {
			AfxMessageBox("Loading weight matrix failed!");
		}
		return S_FALSE;
	}
	
	if (!T.LoadData(model.DATAFILE,model.NVARIABLES,model.NOBSERVATIONS)) {
		Fire_Progress(toBSTR("Loading X matrix failed!\n"));
		if (!mbSilent) {
			AfxMessageBox("Loading X matrix failed!");
		}
		return S_FALSE;
	}
	return S_OK;
}

STDMETHODIMP CComBfrm::Initialise()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	if (model.NDESIGNS == 0) {
		T.SubtractSampleMeans();
		Fire_Progress(toBSTR("Standardized to zero mean!\n"));
	}
	T.InitialisePrior(model);
	T.Initialise(model.NFACTORS, model.NDESIGNS);
	T.InitialiseIterResult(model.NFACTORS, model.NDESIGNS);
	Fire_Progress(toBSTR("Done Init\n"));
	mnDoneInit = 1;
	return S_OK;
}

STDMETHODIMP CComBfrm::Iterate()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	T.SampleFactorModel(model.NFACTORS, model.NDESIGNS, model.MeanValue);
	return S_OK;
}

STDMETHODIMP CComBfrm::Accumulate()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	T.SaveSamples(model.NFACTORS, model.NDESIGNS);
	return S_OK;
}

STDMETHODIMP CComBfrm::Summarise()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	if (model.CheckPoint > 0) {
		T.CheckPoint(model.nmc + model.burnin, "checkpoint_last.chp",true);
		Fire_Progress(toBSTR("The last iteration was saved in checkpoint_last.chp\n"));
	}
	T.Summarise(model.NFACTORS, model.NDESIGNS ,model.nmc);
	return S_OK;
}

STDMETHODIMP CComBfrm::Run()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	mbStop = 0;
	Initialise();
	int nDoneIter = 0;
	if (mbCheckPoint) {
		nDoneIter = T.CheckPoint(nDoneIter,"checkpoint.chp",false);
		if (!nDoneIter) {
			Fire_Progress(toBSTR("Reading checkpoint.chp failed!\n"));
			if (!mbSilent) {
				AfxMessageBox("Reading checkpoint.chp failed!");
			}
			return S_FALSE;
		}
	}

	Stopwatch sw;
	sw.start();
	for (int it=1 + nDoneIter; it <= model.nmc + model.burnin; it++) {
		if(it% model.noisy == 0 || model.noisy < 2) { 
			Fire_DoneIterate(it);
		}
		Iterate();
		if (it>model.burnin) {
			Accumulate();
		}
		if ((model.CheckPoint > 0) && ((model.CheckPoint == 1) || (it% model.CheckPoint == 0))) {
			T.CheckPoint(it,"checkpoint.chp",true);
			Fire_Progress(toBSTR("Saveing Check Point done\n"));
		}
	
		if (mbStop == -1 || mbStop == it) {
			mbStop = it;		 
			Fire_Progress(toBSTR("Iteration stoped by request\n"));
			if (mbStop - model.burnin > 0) {
				T.Summarise(model.NFACTORS, model.NDESIGNS ,mbStop - model.burnin);
			} else {
				Fire_Progress(toBSTR("Burnin stage is not done yet. The result is not valid!\n"));
			}
			return S_OK;
		}
	}
	mbStop = model.nmc + model.burnin;
	Summarise();
	//cout << sw.read() << endl;
	sw.stop();
	return S_OK;
}


STDMETHODIMP CComBfrm::get_B(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = model.NVARIABLES;
	int nSize2 = model.NFACTORS + model.NDESIGNS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			//double dv = (double)(i * nSize2 + j);
			double dv = mnDoneInit?T.GetResult(i,j,0):0;
			SafeArrayPutElement(psa, index, &dv);
		}
		//Copy the SAFEARRAY data into the returning VARIANT variable
		V_VT(pVal) = VT_ARRAY | VT_R8;
		V_ARRAY(pVal) = psa;
	}
	
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Bz(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = model.NVARIABLES;
	int nSize2 = model.NFACTORS + model.NDESIGNS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			//double dv = (double)(i * nSize2 + j);
			double dv = mnDoneInit?T.GetResult(i,j,1):0;
			SafeArrayPutElement(psa, index, &dv);
		}
	}
	
	//Copy the SAFEARRAY data into the returning VARIANT variable
	V_VT(pVal) = VT_ARRAY | VT_R8;
	V_ARRAY(pVal) = psa;


	return S_OK;
}

STDMETHODIMP CComBfrm::get_F(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = model.NFACTORS + model.NDESIGNS;
	int nSize2 = model.NOBSERVATIONS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			//double dv = (double)(i * nSize2 + j);
			double dv = mnDoneInit?T.GetResult(i,j,2):0;
			SafeArrayPutElement(psa, index, &dv);
		}
	}
	
	//Copy the SAFEARRAY data into the returning VARIANT variable
	V_VT(pVal) = VT_ARRAY | VT_R8;
	V_ARRAY(pVal) = psa;

	return S_OK;
}



STDMETHODIMP CComBfrm::get_Designs(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.NDESIGNS;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_Designs(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.NDESIGNS = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_CheckPoint(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.CheckPoint;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_CheckPoint(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.CheckPoint = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorMeanValue(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.MeanValue; 
	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorMeanValue(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	model.MeanValue = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorPsia(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorPsia;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorPsia(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorPsia = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorPsib(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.PriorPsib;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorPsib(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	model.PriorPsib = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorPiProportion(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorPiProportion;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorPiProportion(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorPiProportion = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorPiMultiplier(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	*pVal = model.PriorPiMultiplier;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorPiMultiplier(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorPiMultiplier = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorTaua1(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorTaua1;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorTaua1(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorTaua1 = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorTaub1(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorTaub1;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorTaub1(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorTaub1 = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorTaua2(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorTaua2;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorTaua2(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorTaua2 = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorTaub2(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorTaub2;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorTaub2(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorTaub2 = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorTaua3(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorTaua3;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorTaua3(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorTaua3 = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PriorTaub3(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.PriorTaub3;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_PriorTaub3(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.PriorTaub3 = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_DataFile(BSTR *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	CString s = model.DATAFILE.c_str();
	*pVal=s.AllocSysString();

	return S_OK;
}

STDMETHODIMP CComBfrm::put_DataFile(BSTR newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	CString s =newVal;
	string ss = s;
	model.DATAFILE = ss;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_DesignFile(BSTR *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	CString s = model.HFILE.c_str();
	*pVal=s.AllocSysString();

	return S_OK;
}

STDMETHODIMP CComBfrm::put_DesignFile(BSTR newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	CString s =newVal;
	string ss = s;
	model.HFILE = ss;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_WeightFile(BSTR *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	CString s = model.WEIGHTFILE.c_str();
	*pVal=s.AllocSysString();
	return S_OK;
}

STDMETHODIMP CComBfrm::put_WeightFile(BSTR newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	CString s =newVal;
	string ss = s;
	model.WEIGHTFILE = ss;
	return S_OK;
}

STDMETHODIMP CComBfrm::SaveParameters(BSTR FileName)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	CString s =FileName;
	string ss = s;
	if (ss == "") {
		ss = "parameters.txt";
	}
	if (!model.Save(ss)) {
		return S_FALSE;
	}
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Psi(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = model.NVARIABLES;
	int nSize2 = 1;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			double dv = mnDoneInit?T.GetResult(i,j,3):0;
			SafeArrayPutElement(psa, index, &dv);
		}
		//Copy the SAFEARRAY data into the returning VARIANT variable
		V_VT(pVal) = VT_ARRAY | VT_R8;
		V_ARRAY(pVal) = psa;
	}

	return S_OK;
}

STDMETHODIMP CComBfrm::get_Tau(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = 1;
	int nSize2 = model.NFACTORS + model.NDESIGNS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			double dv = mnDoneInit?T.GetResult(j,i,4):0;
			SafeArrayPutElement(psa, index, &dv);
		}
		//Copy the SAFEARRAY data into the returning VARIANT variable
		V_VT(pVal) = VT_ARRAY | VT_R8;
		V_ARRAY(pVal) = psa;
	}

	return S_OK;
}

STDMETHODIMP CComBfrm::get_Pib(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = 1;
	int nSize2 = model.NFACTORS + model.NDESIGNS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			double dv = mnDoneInit?T.GetResult(j,i,5):0;
			SafeArrayPutElement(psa, index, &dv);
		}
		//Copy the SAFEARRAY data into the returning VARIANT variable
		V_VT(pVal) = VT_ARRAY | VT_R8;
		V_ARRAY(pVal) = psa;
	}

	return S_OK;
}

STDMETHODIMP CComBfrm::get_PostPib(VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	//Create a SAFEARRAY to hold the objects list
	int nSize1 = model.NVARIABLES;
	int nSize2 = model.NFACTORS + model.NDESIGNS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			double dv = mnDoneInit?T.GetResult(i,j,6):0;
			SafeArrayPutElement(psa, index, &dv);
		}
		//Copy the SAFEARRAY data into the returning VARIANT variable
		V_VT(pVal) = VT_ARRAY | VT_R8;
		V_ARRAY(pVal) = psa;
	}

	return S_OK;
}

STDMETHODIMP CComBfrm::get_UseCheckPoint(BOOL *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = mbCheckPoint;
	return S_OK;
}

STDMETHODIMP CComBfrm::put_UseCheckPoint(BOOL newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	mbCheckPoint = newVal;
	return S_OK;
}

STDMETHODIMP CComBfrm::get_StopAt(BOOL *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = mbStop;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_StopAt(BOOL newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	mbStop = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::get_Print(long *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.noisy;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_Print(long newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.noisy = newVal;

	return S_OK;
}

STDMETHODIMP CComBfrm::Resume()
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	//if (mbStop > 0 && !mbDoneSummary) {
	if (mbStop > 0) {
		try {	
			Stopwatch sw;
			sw.start();
			int itStart = mbStop + 1; 
			//Fire_Progress(toBSTR("Iteration resumed by request\n"));
			for (int it=itStart; it <= model.nmc + model.burnin; it++) {
				if(model.noisy < 2 || it% model.noisy == 0) { 
					Fire_DoneIterate(it);
				}
				Iterate();
				if (it>model.burnin) {
					Accumulate();
				}
				if ((model.CheckPoint > 0) && ((model.CheckPoint == 1) || (it% model.CheckPoint == 0))) {
					T.CheckPoint(it,"checkpoint.chp",true);
				}
				if (mbStop == -1 || mbStop == it) {
					mbStop = it; 
					Fire_Progress(toBSTR("Iteration stoped by request\n"));
					if (mbStop - model.burnin > 0) {
						T.Summarise(model.NFACTORS, model.NDESIGNS ,mbStop - model.burnin);
					} else {
						Fire_Progress(toBSTR("Burnin stage is not done yet. The result is not valid!\n"));
					}
					return S_OK;
				}
			}
			mbStop = model.nmc + model.burnin;
			Summarise();
		} catch (...) {
			return S_FALSE;
		}
	} else {
		AfxMessageBox("Can not resume from here!");
	}
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Silent(BOOL *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = mbSilent;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_Silent(BOOL newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	mbSilent = newVal;
	
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Sparsity(double *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = T.mSparsity;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_Sparsity(double newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	// do nothing, only to get matlab inspector showing this property

	return S_OK;
}

STDMETHODIMP CComBfrm::get_BayesFactorSelect(double threshold, VARIANT *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())
	//Create a SAFEARRAY to hold the objects list
	int nSize1 = model.NVARIABLES;
	int nSize2 = model.NFACTORS + model.NDESIGNS;

	SAFEARRAY FAR* psa;
	SAFEARRAYBOUND rgsabound[2];
    rgsabound[0].lLbound = 0;
	rgsabound[0].cElements = nSize1;
	rgsabound[1].lLbound = 0;
	rgsabound[1].cElements = nSize2;
	psa = SafeArrayCreate(VT_R8, 2, rgsabound);

	long index[2];
	
	//Copy the objects list into the SAFEARRAY
	for (long i=0; i<nSize1; i++)
	{
		for (long j = 0; j < nSize2; j++) {
			index[0] = i; index[1] = j;
			double pib = T.GetResult(j,0,5);
			double postpib = T.GetResult(i,j,6);
			double dv;
			if (pib == 0) {
				dv = 0;
			} else if (postpib == 1) {
				dv = 1;
			} else {
				double ratio = postpib * (1-pib) / (pib * (1-postpib));
				if (ratio > threshold) {
					dv = 1;
				} else {
					dv = 0;
				}
			}
			SafeArrayPutElement(psa, index, &dv);
		}
		//Copy the SAFEARRAY data into the returning VARIANT variable
		V_VT(pVal) = VT_ARRAY | VT_R8;
		V_ARRAY(pVal) = psa;
	}
	return S_OK;
}

STDMETHODIMP CComBfrm::get_Mask(BOOL *pVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	*pVal = model.Mask;

	return S_OK;
}

STDMETHODIMP CComBfrm::put_Mask(BOOL newVal)
{
	AFX_MANAGE_STATE(AfxGetStaticModuleState())

	model.Mask = newVal;

	return S_OK;
}

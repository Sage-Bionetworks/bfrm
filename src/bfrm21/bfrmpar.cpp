/*
AUTHOR: QUANLI WANG (QUANLI@STAT.DUKE.EDU)
*/
#include <stdlib.h>
#include "/usr/include/mpi.h"
#include "newran.h"
#include "tnt.h"
#include "jama_svd.h"
#include "jama_cholesky.h"
using namespace TNT;
using namespace JAMA;
#include "stdafx.h"
#include "bfrm.h"
#include "Model.h"
#include "BfrmMPI.h"

Bfrm T;
Model model;

void SendMessageFindPredictors(int* work, int jobsRunning) {
	MPI_Send(work,
		1,
		MPI_INT,
		jobsRunning,
		FIND_PREDICTORS,
		MPI_COMM_WORLD);
}

void master()
{
	/*
	int jobsRunning;
	int ntasks;

	MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
*/
	for (int it=1; it <= model.nmc + model.nit; it++) {
		if(it% model.noisy == 0) { 
			std::cout << "MCMC iteration " << it << endl;
			std::cout << "Sparsity " << T.mSparsity << endl;
		}
		T.SampleFactorModel(model.NFACTORS);

		if (it>model.nit) {
			//update mF mPsi mD mpib mBz mB using B,Bz,F,D,Psi,pib
			T.SaveSamples(model.NFACTORS);
		}
	}
	T.Summarise(model.NFACTORS,model.nmc);
/*
	jobsRunning = 1;
	for(gene=0;gene<Data.NumberOfGenes;gene++) {
		if (!tracer.IsDone(gene)) {
			work = gene;
			printf("Working on gene :: %d\n",gene+1);
			if(jobsRunning<ntasks) {
				SendMessageFindPredictors(&work, jobsRunning); //send a new work request
				jobsRunning++;
			} else {
				RecvMessage(workresults, status);
				AppendSamplingResult(buffer, workresults); //save the models in a file
				tracer.WriteTrace(workresults);
				SendMessageFindPredictors(&work, status.MPI_SOURCE);
			}
		}
	}
	for(gene=1;gene<jobsRunning;gene++)	{
		RecvMessage(workresults, status);
		AppendSamplingResult(buffer, workresults); //save the models in a file
		tracer.WriteTrace(workresults);
	}
	//Tell all the slaves to exit
	SendMessageKillSlave(ntasks);
*/	
}

void slave(int myrank)
{
	int i,j;
	int npred;
	int work;
	int workresults[Data.NumberOfGenes];
	MPI_Status status;

	//read the dataset
    Data.ReadData();
	Data.Allocate();
	GS.Allocate(Data);
	//seed the random number generator
	//the seed is different for each slave
	GS.unifOpen.seed((unsigned int)(randomseed+100*myrank));

	int notDone = 1;
	set<int> varFreq;
	while(notDone) {
        int lenwork = 1; 
		SlaveRecvMessage(&work, lenwork, status);
		int targetGene = work;
		switch(status.MPI_TAG) {
			case FIND_PREDICTORS:
				GS.gibbsSamplingReg(targetGene,varFreq, Data);
				npred = varFreq.size();
				workresults[0] = targetGene;
				workresults[1] = npred;
				j = 0;
				for (set<int>::iterator it = varFreq.begin(); it!= varFreq.end(); it++) {
					workresults[2+j] = *it;
					j++;
				}
				SlaveSendMesage(workresults, npred+2);
				break;
       
			case DIETAG:
				notDone = 0;
				break;
		}
	}
	//clean memory
	Data.Cleanup();
	GS.Cleanup();
}
/*
int main(int argc,char* argv[])
{
	

	//////////////////////////////////
	//the process finds its own rank//
	//////////////////////////////////  
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank); 
	if(0==myrank) {
		master();
	}  else {
		slave(myrank);
	}
	//close the MPI session
	MPI_Finalize();
	return(0);
}
*/
int main()
{
	int myrank;
	MPI_Init(&argc, &argv);
	if (argc != 2) {
	  cout << "usage: " <<  argv[0] << " " << "parameters_file.txt" << endl;
	  MPI_Finalize();
	  exit(1);
	}
	
	if (!model.Load(argv[1])) {
	  cout << model.GetErrorMessage() << endl;
	  MPI_Finalize();
	  exit(1);
	}

	//read data
	T.LoadData(model.DATAFILE,model.NVARIABLES,model.NOBSERVATIONS);
	Stopwatch sw;
	//Initialise
	T.SubtractSampleMeans();
	T.Initialise(model.NFACTORS, model.INITMETHOD, 
		model.BFILE,model.FFILE, model.DFILE);
	T.InitialisePrior(model.NFACTORS, model.BetaA);
	T.InitialiseIterResult(model.NFACTORS);
	//std::cout<<"Done init" << endl;
	sw.start();
	
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank); 
	if(0==myrank) {
		master();
	}  else {
		slave(myrank);
	}
	cout << sw.read() << endl;
	sw.stop();
	MPI_Finalize();
	return 0;
}


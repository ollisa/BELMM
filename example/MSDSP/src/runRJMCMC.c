/*
DATE: 2023-08-13 
NAME: runRJMCMC.c 
DESC: Bash Generated Script to Run Parallel RJMCMC Algorithm

    Run make sure to properly 
    Number of process is NMODES*CHAINS+ + CHAINS+ 1
 
    mpicc -o <name> <file>
    mpirun -np <nProc>
    Make sure to 
    export TMPDIR=/tmp

*/

/* Required Headers */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*MPI Header */
#include <mpi.h>
#include <string.h>

/*GNU Scientif Library Packages */
//Factorial 
//#include <gsl/gsl_sf_fact.h>
//Rando Number
#include <gsl/gsl_rng.h>





/* Constants */
//Simulation Info
#define CHAINS 1
#define CHAINLEN 4000 
#define NMODEL 3
#define COLLEN NMODEL*CHAINS
#define MAXPARAM 205

#define STRINGLEN (35)*sizeof(int)

//Differnt Keneral Density Estiamte (kde) implemnted
#define UNIFORM 0 
#define GAUSSIAN 1


#define PI 3.14159265358979323846

/* ------------------------------------------ */
/*FUNCTION: proposedModel
INPUT: current model
OUTPUT: Calculates proposed move
DESC: Get propsoed move.  */
int proposedModel(int k){
  double pU; 
  int pMove; 
  
  pU = ((double)rand()/ (double)RAND_MAX);

  if (pU < .333333){
    pMove = 0;
  }

  else if (pU < .666666){
    pMove = 1;
  }
  else{
    pMove = 2;
  }
  return pMove;
}




/* ------------------------------------------ */
/*FUNCTION: initialModel
INPUT: 
OUTPUT: Initial Model
DESC: intialize model  */
int initModel(){
  double uInitModel; 
  int initModel;  

  uInitModel = ((double)rand()/ (double)RAND_MAX);
  

  if (uInitModel < .333333){
    initModel = 0;
  }

  else if (uInitModel < .666666){
    initModel = 1;
  }
  else{
    initModel = 2;
  }
  return initModel;
}


/*
AUTHOR: John Chavis
NAME: bodypRJMCMC.c 
DESCRIPTION: Body of the Phase Two of CU-MSDSp


MIT License

Copyright (c) 2020 John T. Chavis III

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

/* ------------------------------------------ */
/*FUNCTION: normpdf
INPUT:  shape, rate
OUTPUT: norm pdf
DESC:   */
double normPdf(float x, float mean, float sd) {
  double val; 

  val = 1/sqrt(2*PI*pow(sd,2)) * (exp(-0.5*(pow(x-mean/sd,2)))); 

  return val; 
}

/* ------------------------------------------ */
/*FUNCTION: sort
INPUT:  data, TOTDATA
OUTPUT:
DESC:   */
//Sort function for finding interquartile range 
//Return order of index in sorted list 
//Bubble sort algorithm  

double * sort(double *x, int runTime){
    int i;
    int j;
    double temp;
    double * sortData; 

    sortData = (double *)malloc(runTime*sizeof(double));

    //Copy data so I dont ruin original list 
    for (i =0; i <runTime; i ++){
      *(sortData  +i) = *(x + i); 
    }

    for (i = 0; i < (runTime - 1); ++i){
        for (j = 0; j < runTime - 1 - i; ++j ){
            if (*(sortData + j) > *(sortData + j+1)){
                temp = *(sortData + j+1);
                *(sortData + j+1) = *(sortData +j);
                *(sortData + j) = temp;
            }
        }
    }

    return sortData; 
}


/* ------------------------------------------ */
/*FUNCTION: IQR
INPUT:  data, TOTDATA
OUTPUT:
DESC:   */
//Index of IQR
//From https://www.hackerrank.com/challenges/s10-interquartile-range/problem
//https://www.geeksforgeeks.org/interquartile-range-iqr/
double IQR(double *x, int totX) { 
  double q1; 
  double q3; 
  double iqr; 
  //Index bounds
  int q1IndxLower; 
  int q1IndxUpper; 
  int q3IndxLower; 
  int q3IndxUpper;


  int medianIndx; 
  //Find median index
  //If even then split list evenly
  if (totX %2 == 0){ 
    medianIndx = totX/2;
  }
  //The total data is an odd number shoudl also be easier
  else{
    medianIndx = (totX -1)/2; 

  }

  //See if the median index inteself splits the data evenly
  if (medianIndx %2 == 0){
    //Have to find average 
    //q1 index 
    q1IndxUpper = medianIndx/2; 
    q1IndxLower =  medianIndx/2 -1; 

    //Calcualte average 
    q1 = (*(x + q1IndxUpper) + *(x+ q1IndxLower))/2.0;


    //Q3 index 
    q3IndxUpper = medianIndx + medianIndx/2; 
    q3IndxLower =  medianIndx + medianIndx/2 + 1; 

    q3 = (*(x + q3IndxUpper) + *(x + q3IndxLower))/2.0;
  }
  //Otehrwise we haev an unevent number, should be easier- we should never get here
  else{
    q1 = *(x + (medianIndx-1)/2); 
    q3 = *(x + medianIndx + (medianIndx-1)/2);
  }

  //printf("Iqr %f %f\n", q1, q3);
  iqr = q3-q1; 
    return iqr; 
} 


//KDe esitmates
/* ------------------------------------------ */
/*FUNCTION: sampleStd
INPUT:  data, TOTDATA
OUTPUT: sampleStd
DESC:   */
double sampleStd(double *x, int totX) {
  double val; 
  double mean; 
  int i; 

  double sd; 

  //Find eman of data 
  val = 0; 

  for (i=0; i< totX; i++){
    val = val + *(x + i); 
  }


  mean = val/totX; 

  val = 0; 

  for (i =0; i < totX; i++){
    val = val + pow((*(x + i) - mean),2); 
  }

  sd = sqrt(val/(totX-1)); 

  return sd; 
}


/* ------------------------------------------ */
/*FUNCTION: kde
INPUT:  data,, type
OUTPUT: kde estimates compared for each data point 
DESC:   */
double * kde(double *subData, int type, int * nParamVec){
  int i; 
  int j; 

  double *kdeEstimate;

  double h; 
  double sd; 

  //Interquartile 
  double iqr; 

  double a;  
  double tmp; 
  double tmpMult; 

  double val; 

  //Model 
  int k;
  int nParam; 
  int param; 

  
  //Create kde esitmate
  kdeEstimate = (double *)malloc(CHAINLEN*NMODEL*sizeof(double));

  if (type == GAUSSIAN){
    //USe Silverman method h d is dimension 
    //n = is total data
    //Sigma is standard deviatoin 
    //h = (4/(d+2))^{1/(d+4)} * n ^{-1/(d+4)} \sigma 
    //I can also ue more robust interquartile method 
    //Either h = 1.06\sd N^{−1/5} 
    // Or h = 0.9AN^{-1/5} where A = min (sd, IQR/1.3) IQR - interquartile range
    //The latter is the more robust measrure better at handing multimodel densitiyies
    //LEt's try sinterquartle method 
    for (k= 0; k <NMODEL; k++){
      nParam = *(nParamVec + k); 

      //Find H vector 
      double *hVec = (double *)malloc(nParam*sizeof(double));

      //Collect the data and h vector 
      for (param = 0; param < nParam; param++){
        //Build vector to collect h?...Same with data from subset 
        double *tmpData = (double *)malloc(CHAINLEN*sizeof(double));

        for (i =0; i < CHAINLEN; i ++){
          *((tmpData + i)) = *((subData + k*CHAINLEN*MAXPARAM + i*MAXPARAM) + param); 

          //printf("sub = %lf\n", *((tmpData + i)));
        }

        //Sort the data 

        double * sortData = (double *)malloc(CHAINLEN*sizeof(double));
        sortData = sort(tmpData, CHAINLEN); 

        //Get interquartilerange 
        iqr = IQR(sortData, CHAINLEN); 

        //Find sample standard deviation 
        sd = sampleStd(sortData, CHAINLEN);

         //Great! Now we have two options to find "optiamal h"
  
        if (sd < iqr/1.34){
          a = sd; 
        }
        else{
          a = iqr/1.34; 
        }

        h = 0.9*a*pow(CHAINLEN, -1/5); 
          
        *((hVec + param)) = h; 


        //printf("H = %lf\n", h); 
        free(tmpData); 
        free(sortData); 
      }
      
      //Great now caluclate Kernel desntiy 
      //p(x) = 1/n \sum_{i=1}^{n} (1/h*...h_{d}} \Pi_{d=1}^{D}K(x-xd/h))
      for (i= 0; i < CHAINLEN; i++){
        tmp = 0; 
        tmpMult = 1; 

        for (j = 0; j < CHAINLEN; j ++){
          tmpMult = 1; 
          for (param = 0; param < nParam; param++){
            h =  *((hVec + param)); 

            //Kde part (x- x_i)/h
            val = normPdf((*((subData + k*CHAINLEN*MAXPARAM + i*MAXPARAM) + param) -  *((subData + k*CHAINLEN*MAXPARAM + j*MAXPARAM) + param))/h,0,1);
            //printf("What about here\n");
            //printf("Val = %lf\n", *((subData + k*CHAINLEN*MAXPARAM + j*MAXPARAM) + param)); 
           // printf("val = %lf\n", val/h); 
            tmpMult = tmpMult*val/h; 
          }
          
          tmp = tmp + tmpMult; 
        }

        //printf("%d %d %d %d %d \n", i*NMODEL + k, runTime*NMODEL, i, runTime, k); 
        *((kdeEstimate + i*NMODEL) + k) = log(tmp/CHAINLEN); 
        //printf("Tmp = %lf\n", log(tmp/CHAINLEN)); 
        //printf("Made it out\n"); 
      }
      free(hVec); 
    }
  }
  
  //printf("Please make it\n:");
  return kdeEstimate; 
}




/* ------------------------------------------ */
/*FUNCTION: acceptCralwer(.)
INPUT: data, kC, kP, cParam, pParam, cTime,  pTimme, kdeEstimates
OUTPUT: acceptance probability
DESC: Calcualte acceptance probability given priors adn proposla move */
double acceptCrawler(int kC, int kP, double cLikelihood, double pLikelihood, double * kdeEstimates, int cTime, int pTime){
  //Calculates the Likelihood*prior*proposal for current model 
  double tmpVal; 

  double modelAcceptProb; 

  double cPosterior; 
  double pPosterior; 




  //Calulate the posteiror estimates for each model 
  //cPrposeed corresponds to whaever distribution you pull kC 
  cPosterior = *((kdeEstimates + cTime *NMODEL) + kC); 
  pPosterior = *((kdeEstimates + pTime *NMODEL) + kP); 

  tmpVal = (pLikelihood+ cPosterior) - (cLikelihood + pPosterior); 

  //printf("kC %d cLike %lf cPost %lf kP %d pLike %lf pPost %lf\n", kC, cLikelihood, cPosterior, kP, pLikelihood, pPosterior); 
  //Get minimum z value for accepance//If acceptrance is infinit return accept as 0 i.e reject
  if (isinf(tmpVal) || isnan(tmpVal)){
    modelAcceptProb = 0.1; 
  }
  else if (tmpVal >0){
    //This is the exp()
    modelAcceptProb = 1;
  }
  else{
    modelAcceptProb = exp(tmpVal); 
  }

  return modelAcceptProb;
} 




/* ------------------------------------------------------------------ */
/*                    RJMCMC                                            */
/* ------------------------------------------------------------------ */
//Name: rjmcmc Crawler
//INPUT: 
//Output: Model samples
//Description: peforms the fixed mcmc 
double * rjmcmcCrawler(int chain, double *A, int * nParamVec, int myId){
  int kC; 
  int kP; 
  double acceptCount; 

  double acceptProb;  

  //Random number generation 
  double u;  

  int link;

  int chainPos; 
  int indx;
  int chainPosP; 
  int indxP; 



  //Create  paramter vector for this 
  double cLikelihood;
  double pLikelihood; 

  
  //Index for pulling a radnom sample 
  int curTime; 
  int propTime; 

  int randIndx;
  //Create gold stand chain
  //CHAIN LEN x  |modelIndex, Acceptance| 
  double * crawlRes = (double *)malloc((CHAINLEN*2)*sizeof(double));


  //Kde estimate [Poisson, NB]
  //Try to generalize later.....
  double * kdeEstimates = (double *)malloc((CHAINLEN*NMODEL)*sizeof(double));
  double *subData = (double *)malloc((CHAINLEN*MAXPARAM*NMODEL)*sizeof(double)); 

  //Model 
  int k;
  int nParam; 
  int param; 
  int i; 
  
  //For Random Number
  //For Setting the seeds
  time_t tme = time(NULL);
  const unsigned m = 1664525u;
  const unsigned c = 1013904223u;
  unsigned t_hashed = (unsigned) tme;
  t_hashed = m*t_hashed + c;

  gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);
  gsl_rng_set(r, t_hashed ^ myId);

  //Copy Sub data over (should I do this instead of sending A all the way over....memory issue)
  for (k = 0; k < NMODEL; k++){
    nParam = *(nParamVec + k); 

    //Calculate chain pois 
    chainPos = k*CHAINS + chain; 

    for (i=0; i < CHAINLEN; i++){
      indx = chainPos*CHAINLEN*(MAXPARAM+1) + i*(MAXPARAM+1); 

      for (param = 1; param <= nParam; param ++){
        *((subData + k*CHAINLEN*MAXPARAM + i*MAXPARAM) + (param-1) ) = *(A + indx + param);
        //printf("val %lf\n", *(A + indx + param));
      }
    }
  }

  //Calculate Kernel density estimate 
  //This will be tricky we need to figure out 
  kdeEstimates = kde(subData, GAUSSIAN, nParamVec);

  //Initialize 
  acceptCount = 0; 

  kC  = initModel();

  *((crawlRes + 0*2) + 0) = kC; 
  *((crawlRes + 0*2) + 1) = acceptCount; 



  //Simulation 
  for (link = 1; link < CHAINLEN; link++){

    //Get current model 
    kC = *((crawlRes + (link-1)*2) + 0);

    //Find chain position (column)
    chainPos = kC*CHAINS + chain; 

    //Upper intt -lower int +1 + lower 
    //curTime = (rand() % CHAINLEN ) -1;
    //curTime =  (gsl_rng_get(r) % CHAINLEN) -1;
    curTime = gsl_rng_uniform_int(r, CHAINLEN); 
    indx = chainPos*CHAINLEN*(MAXPARAM+1) + (curTime)*(MAXPARAM+1); 

    //Grab likelihood from our dataStructure 
    cLikelihood = *(A + indx + 0); 



    //Propose a new model 
    kP = proposedModel(kC); 

    //Find chain position (column)
    chainPosP = kP*CHAINS + chain; 

    //Upper intt -lower int +1 + lower 
    //propTime = (rand() % CHAINLEN ) -1;
    //propTime = (gsl_rng_get(r) % CHAINLEN) -1; 
    propTime = gsl_rng_uniform_int(r, CHAINLEN); 
    
    indxP = chainPosP*CHAINLEN*(MAXPARAM+1) + (propTime)*(MAXPARAM+1); 

    pLikelihood = *(A + indxP + 0); 

    //printf("cLike %lf pLike %lf \n", cLikelihood, pLikelihood);   
    //Acceptane Prob 
    //u = ((double)rand()/ (double)RAND_MAX);
    u = gsl_rng_uniform(r); 
    acceptProb = acceptCrawler(kC, kP, cLikelihood, pLikelihood, kdeEstimates, curTime, propTime); 

    if (u < acceptProb){ 
      *((crawlRes + (link)*2) + 0) = kP;
      acceptCount = acceptCount + 1; 
      curTime = propTime;
    }
    else{
      *((crawlRes + (link)*2) + 0) = kC;
    }
  
    *((crawlRes + (link)*2) + 1) = acceptCount/link; 


  }

  free(subData); 
  free(kdeEstimates);
  //printf("Issue here"); 
  return crawlRes;
}





/* ============================== MAIN FUNCTION  ==================================== */
//
//                  MAIN FUNCTION 
////
/* ============================== MAIN FUNCTION  ==================================== */
int main(int argc, char **argv){
  //For timing 
  clock_t begin;
  clock_t end; 
  double timeSpent; 


  //Create result matrix. This the 3D array of pointers
  //Initialize variables ............
  // Read data file
  int i; 
  int j; 
  int h; 

  //For opening up data 
  
  double *data;
  int totData;

  //For keep track of chains, params time
  int t; 
  int param; 
  int chain; 
  // Model indicator 
  int k; 
  int numParam; 

  int col; 
  int colCount; 

  double acceptProb; 

  //For saving gold standard chains 
  char resName[STRINGLEN]; 
  const char *resSumName; 




  int totEleCount; 
  double *A; 
  double * crawlerRes; 



 
  totEleCount = CHAINLEN*COLLEN*(MAXPARAM+1); 
  A = (double*) malloc(totEleCount*sizeof(double));

  //For collcting paramter s
  int *nParamVec; 

  nParamVec = (int*) malloc(NMODEL*sizeof(int)); 

  //For determinging it to run rjmcmc crawler
  //To keep track of convergence 
  int convTime; 
  int runTime; 

  /*---- MPI INITIZIALIZATION ----- */
  //Determine number of processor to use
  int numProcs;

  //Number of workers 
  int numWorkers; 

  //id of each process 
  int myId;

  //For sending messages
  int source;
  int dest;

  //KEep track of any offesning and which model weer workign with 
  int offset; 
  int cols; 


  //Initialize MPI and create group of works, Make one worker to be main processor
  //Will create random numbers and send them 
  MPI_Init(&argc, &argv);

  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myId);

  //Last processor is server for random numbers
  numWorkers = numProcs -1; 

  FILE *f; 

  int row; 
  double x; 

  //char *fileName;   
  char fileName[200];
  char buffer[50]; 
  
  //Creat the A matrix 

  //Start clock
  //begin = clock();
  /*----------------------------------- MAIN SERVER-----------------------*/
  if (myId == 0){
    
   
   

    //Send data, lengths and a column of 
    //cols = COLLEN/numWorkers;
    //offset = 0; 

     /*------------  READ DATA ------------- */ 
    //Reading the paramters and logliklihoods into A data structure
    //Read for each model and chain 
    printf("Reading in Gold Standard Chains from CmdSTAN....\n");
    for (k=0; k < NMODEL; k++){

      for (chain=0; chain < CHAINS; chain++){
        
        //printf("Model %d Chain %d \n", k, chain);

        //snprintf(fileName, STRINGLEN, "goldStandardChains/model%d/gsc%d.csv",k,chain+1);
        //snprintf(fileName, "goldStandardChains/model%d/gsc%d.csv",k,chain+1);
        strcpy(fileName, "goldStandardChains/model");
        sprintf(buffer, "%d", k); 
        strcat(fileName, buffer); 
        strcat(fileName, "/gsc"); 
        
        sprintf(buffer, "%d", chain+1);
        strcat(fileName, buffer); 
        strcat(fileName,".csv"); 

        //Hard coded for right now
        
        
        //fileName = "goldStandardChains/model0/gsc1.csv";       
        //printf("Strin names %s\n", fileName);
        //printf("Strin names %s\n", string);
        f  = fopen(fileName, "r");

        //printf("We made it here\n");
        //Read row and colum nsize 
        fscanf(f, "%d", &row); 
        fscanf(f, "%d", &col); 

        //printf("HOw about here\n"); 
        //printf("%d %d\n", row, col);
        
        //Collect paramters for this model 
        *(nParamVec + k) = col-1; 

        for (i =0; i < row; i++){
          for (j =0; j <col; j++){
            fscanf(f, "%lf", &x);
            //printf("%lf\n",x);
            *(A + i*(MAXPARAM+1) +(k*CHAINS+chain)*(MAXPARAM+1)*CHAINLEN + j) = x; 
          }
        }
      }
    }

  
    printf("CmdSTAN CHAINS Read...\n");
    

    //A should be correct.....
    //SEND chains to crawler 
    //Runtime is 1/2 off convtime for chain sampling purpose
    crawlerRes = (double*)malloc((CHAINLEN)*CHAINS *sizeof(double));
    chain = 0; 
    /*-------------- RJMCMC CRawler -------- */ 
    for (dest = 1; dest <= CHAINS; dest++){

      MPI_Send(&chain, 1, MPI_INT, dest, 1, MPI_COMM_WORLD); 
      MPI_Send(&totEleCount, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
      MPI_Send(A, totEleCount, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD); 
      MPI_Send(nParamVec, NMODEL, MPI_INT, dest, 1, MPI_COMM_WORLD); 

      printf("RJMCMC Chain %d/%d Sent by Master\n", chain+1, CHAINS);
         //printf("%zd \n", dataSize);
  
      chain = chain + 1; 
      
      //printf("PLEase WORK\n");
    }


    //Collect Data 

    for (source = 1; source <= CHAINS; source++){
      //GEt request 
     // MPI_Recv(&request, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
      double *tmpRes = (double *)malloc(((CHAINLEN)*2)*sizeof(double));
      MPI_Recv(&chain, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status); 
      MPI_Recv(tmpRes, CHAINLEN*2, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &status); 
  


      //Model Chain 
      snprintf(resName, STRINGLEN, "modelSelection/modelIndexChain%d.csv",chain+1);
      f = fopen(resName,"w");
      for (t=0; t < CHAINLEN; t++){
        k = *((tmpRes + t*2) + 0);
        *(crawlerRes + t*CHAINS + chain)  = k;
        fprintf(f, "%d", k);
        fprintf(f,"\n"); 
          
      }
      fclose(f); 
 
     
      //Acceptance 
      snprintf(resName, STRINGLEN, "modelSelection/acceptChain%d.csv",chain+1);

      f = fopen(resName,"w");

      for (t=0; t < CHAINLEN; t++){
        //Get model index
        acceptProb = *((tmpRes + t*2) + 1); 
        fprintf(f, "%lf", acceptProb);
        fprintf(f,"\n"); 
          
      }
      fclose(f); 


      free(tmpRes); 
      printf("RJMCMC Chain %d/%d Received by Master\n", chain+1, CHAINS);
    }

   
    free(A); 
  }


  //-------WORKERS-------------
  else if (myId <= CHAINS){
    /*------------- Crawler ----------- */
    //REceive work for the cralwer 
    //MPI_Recv(&request, 1, MPI_INT, source,1, MPI_COMM_WORLD, &status); 
    //Set seed 

    source = 0; 

    MPI_Recv(&chain, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status); 
    MPI_Recv(&totEleCount, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);

    A = (double*) malloc(totEleCount*sizeof(double));
    MPI_Recv(A, totEleCount, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &status); 

    nParamVec = (int*) malloc(NMODEL*sizeof(int)); 
    MPI_Recv(nParamVec, NMODEL, MPI_INT, source, 1, MPI_COMM_WORLD, &status);


     //Pefform the crawler 
    double *tmpRes = (double *)malloc(((CHAINLEN*2))*sizeof(double));
    tmpRes = rjmcmcCrawler(chain,A, nParamVec, myId);
  

    MPI_Send(&chain, 1, MPI_INT, source, 2, MPI_COMM_WORLD);
    MPI_Send(tmpRes, CHAINLEN*2, MPI_DOUBLE, source, 2, MPI_COMM_WORLD);
    //printf("We make it here\n");
    free(tmpRes);
    printf("RJMCMC Chain %d/%d Completed by Worker %d\n", chain+1, CHAINS, myId); 
  }

  MPI_Finalize(); 
  return 0;
}


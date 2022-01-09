#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include<iostream>
#include <string.h>
#include <time.h>
using namespace std;

extern void InitSort( void );

class TSort 
{
 public:
  TSort();
  ~TSort();
  void Index( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd );
  void Index( int* Arg, int numOfArg, int* indexOrderd, int numOfOrd );
  void Index_B( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd );
  void Index_B( int* Arg, int numOfArg, int* indexOrderd, int numOfOrd );
  void Sort( int* Arg, int numOfArg );
};

extern TSort* tSort;
class TIndi 
{
public:
  TIndi();
  ~TIndi();
  void Define( int N );
  TIndi& operator = ( const TIndi& src );   /* Copy */
  bool operator == (  const TIndi& indi2 ); /* Return true if two tours are the same, false otherwise */

  int fN;                 /* Number of cities */
  int** fLink;            /* fLink[i][]: two vertices adjacent to i */
  int fEvaluationValue;   /* Tour length of */
};

class TEvaluator {
 public:
  TEvaluator();
  ~TEvaluator();
  void SetInstance( char filename[] );       /* Set the instance */
  void DoIt( TIndi& indi );                  /* Set the value of indi.fEvaluationValue */
  void WriteTo( FILE* fp, TIndi& indi );     /* Write an tour */
  bool ReadFrom( FILE* fp, TIndi& indi );    /* Read an tour */
  bool CheckValid( int* array, int value ); /* Check an tour */ 

  int fNearNumMax;  /* Maximum number of k (see below) */
  int **fNearCity;  /* NearCity[i][k]: k-th nearest city from */
  int  **fEdgeDis;  /* EdgeDis[i][j]: distance between i and j */
  int Ncity;        /* Number of cities */
  double *x;        /* x[i]: x-coordinate of */
  double *y;        /* y[i]: x-coordinate of */
};


class TCross {
public:
  TCross( int N );
  ~TCross();
  void SetParents( const TIndi& tPa1, const TIndi& tPa2,     /* Set information of the parent tours */
		   int flagC[ 10 ], int numOfKids ); 
  void DoIt( TIndi& tKid, TIndi& tPa2, int numOfKids,        /* Main procedure of EAX */
	     int flagP, int flagC[ 10 ], int** fEdgeFreq ); 
  void SetABcycle( const TIndi& parent1, const TIndi& parent2, /* Step 2 of EAX */
		   int flagC[ 10 ], int numOfKids );
  void FormABcycle();                                   /* Store an AB-cycle found */
  void Swap(int &a,int &b);                             /* Swap */ 
  void ChangeSol( TIndi& tKid, int ABnum, int type );   /* Apply an AB-cycle to an intermediate solution */
  void MakeCompleteSol( TIndi& tKid );                  /* Step 5 of EAX */
  void MakeUnit();                                      /* Step 5-1 of EAX */ 
  void BackToPa1( TIndi& tKid );                        /* Undo the parent p_A */
  void GoToBest( TIndi& tKid );                         /* Modify tKid to the best offspring solution */
  void IncrementEdgeFreq( int **fEdgeFreq );            /* Increment fEdgeFreq[][] */
  int Cal_ADP_Loss( int **fEdgeFreq );                  /* Compute the difference in the averate distance */
  double Cal_ENT_Loss( int **fEdgeFreq );               /* Compute the difference in the edge entropy */

  void SetWeight( const TIndi& parent1, const TIndi& parent2 ); /* Block2 */
  int Cal_C_Naive();                                            /* Block2 */
  void Search_Eset( int num );                                  /* Block2 */
  void Add_AB( int AB_num );                                    /* Block2 */
  void Delete_AB( int AB_num );                                 /* Block2 */

  void CheckValid( TIndi& indi );                               /* For debug */


  int fNumOfGeneratedCh;
  TEvaluator* eval;			 
  int fNumOfPop;
  
private:
  int fFlagImp;         
  int fN;
  TIndi tBestTmp;
  int r,exam;
  int exam_flag;
  int **near_data;
  int *koritsu, *bunki, *kori_inv, *bun_inv;
  int koritsu_many,bunki_many;
  int st,ci,pr,stock,st_appear;
  int *check_koritsu;
  int *fRoute;
  int flag_st,flag_circle,pr_type;
  int ch_dis;
  int **fABcycle;
  int *fPermu;
  int fEvalType;
  int fEsetType;
  int fNumOfABcycleInESet;
  int fNumOfABcycle;
  int fPosiCurr;
  int fMaxNumOfABcycle;

  int *fC;
  int *fJun; 
  int *fOrd1, *fOrd2; 

  // Speed Up Start
  int *fOrder;    
  int *fInv;      
  int **fSegment; 
  int *fSegUnit;  
  		       
  int fNumOfUnit; 
  int fNumOfSeg;  
  int *fSegPosiList;
  int fNumOfSPL;    
  int *LinkAPosi;   
  int **LinkBPosi;  
  int *fPosiSeg;    
  int *fNumOfElementInUnit; 
  int *fCenterUnit;         
  int fNumOfElementInCU;    
  int *fListOfCenterUnit;   
  int fNumOfSegForCenter;   
  int *fSegForCenter;       

  int *fGainAB;             
  int fGainModi;            
  int fNumOfModiEdge;				 
  int fNumOfBestModiEdge;				 
  int **fModiEdge;				 
  int **fBestModiEdge;			
  int fNumOfAppliedCycle;
  int fNumOfBestAppliedCycle;
  int *fAppliedCylce;
  int *fBestAppliedCylce;
  // Speed Up End

  // Block2
  int *fNumOfElementINAB;    
  int **fInEffectNode; 
  int **fWeight_RR; 
  int *fWeight_SR;  
  int *fWeight_C;   
  int *fUsedAB, fNumOfUsedAB;
  int fNum_C, fNum_E;
  int fTmax, fMaxStag;
  int *fMoved_AB;
  int fNumOfABcycleInEset;
  int *fABcycleInEset;
  int fDis_AB;     
  int fBest_Num_C, fBest_Num_E;
};

class TKopt {
public:
  TKopt( int N );
  ~TKopt();
  void SetInvNearList();
  void TransIndiToTree( TIndi& indi );
  void TransTreeToIndi( TIndi& indi );

  void DoIt( TIndi& tIndi );             /* Apply a local search with the 2-opt neighborhood */
  void Sub();
  int GetNext( int t );
  int GetPrev( int t );
  void IncrementImp( int flagRev );
  void CombineSeg( int segL, int segS );

  void CheckDetail();
  void CheckValid();

  void Swap(int &a,int &b);
  int Turn( int &orient );

  void MakeRandSol( TIndi& indi );      /* Set a random tour */


  TEvaluator* eval;

private:
  int fN;

  int fFixNumOfSeg;
  int fNumOfSeg;   
  int **fLink;     
  int *fSegCity;   
  int *fOrdCity;   

  int *fOrdSeg;    
  int *fOrient;    
  int **fLinkSeg;  
  int *fSizeSeg;   
  int **fCitySeg;  

  int *fT;
  int fFlagRev;  
  int fTourLength;

  int *fActiveV;
  int **fInvNearList; 
  int *fNumOfINL;     
  
  int *fArray;
  int *fCheckN;
  int *fGene;
  int *fB;
};

 
class TEnvironment {
public:
  TEnvironment(); 
  ~TEnvironment();
  void Define();                         /* Define the variables */

  void DoIt();                           /* Main procedure of the GA */
  void Init();                           /* Initialization of the GA */
  bool TerminationCondition();           /* Decide whether to proceed to next stage 
					    (or treminate the GA) */
  void SetAverageBest();                 /* Compute average and best tour lengths of the population */
  void InitPop();                        /* Create an initial population */
  void SelectForMating();                /* Determine a set of pairs of parents at each generation */
  void SelectForSurvival( int s );       /* Not do anything */
  void GenerateKids( int s );            /* Generate offspring solutions from a selected pair of  
                                            parents. Selection for survival is also performed here. */
  void GetEdgeFreq();                    /* Compute the frequency of the edges of the population */

  void PrintOn( int n, char* dstFile );  /* Display and write summary of results */
  void WriteBest( char* dstFile );       /* Write the best tour */
  void WritePop( int n, char* dstFile ); /* Write the population */
  void ReadPop( char* fileName );        /* Read an initial population */


  TEvaluator* fEvaluator;                /* Distance of the edges */
  TCross* tCross;                        /* Eede assembly crossover */
  TKopt* tKopt;                          /* Local search with the 2-opt neighborhood */
  char *fFileNameTSP;                    /* File name of an TSP instance */
  char *fFileNameInitPop;                /* File name of an initial population */

  int fNumOfPop;                         /* Number of population members (N_pop in the paper) */
  int fNumOfKids;                        /* Number of offspring solutions (N_ch in the paper) */
  TIndi* tCurPop;                        /* Current population members */
  TIndi tBest;                           /* Best solution in the current population */
  int fCurNumOfGen;                      /* The current number of generations */
  long int fAccumurateNumCh;             /* The accumulated number of offspring solutions */
  int fBestNumOfGen;                     /* The number of generations at which the current best 
                                            solution was found */ 
  long int fBestAccumeratedNumCh;        /* The accumulated number of offspring solutions at which 
                                            the current best solution was found */
  int **fEdgeFreq;                       /* The frequency of the edges of the population */
  double fAverageValue;                  /* The average tour lengths of the population */
  int fBestValue;                        /* The tour lenght of the best tour in the population */
  int fBestIndex;                        /* Index of the best population member */
  int* fIndexForMating;                  /* Mating list (r[] in the paper) */

  int fStagBest;                         /* The number of generations during which no improvement  
                                            is found in the best tour */
  int fFlagC[ 10 ];                      /* Specify configurations of EAX and selection strategy */
  int fStage;                            /* Current stage */
  int fMaxStagBest;                      /* If fStagBest = fMaxStagBest, proceed to the next stage */
  int fCurNumOfGen1;                     /* Number of generations at which Stage I is terminated */

  clock_t fTimeStart, fTimeInit, fTimeEnd;  /* Use them to measure the execution time */
};




TSort* tSort = NULL;

extern void InitURandom( int );    
extern void InitURandom( void );   

class TRandom {
 public:
  TRandom();
  ~TRandom();
  int Integer( int minNumber, int maxNumber ); 
  double Double( double minNumber, double maxNumber );
  void Permutation( int *array, int numOfelement, int numOfSample );
  double NormalDistribution( double mu, double sigma );
  void Shuffle( int *array, int numOfElement );
};

extern TRandom* tRand;





void InitSort( void )
{
  tSort = new TSort();
}

TSort::TSort()
{
}

TSort::~TSort()
{
}


void TSort::Index( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest = 0;
  double valueBest;
  int *checked;
  checked = new int [ numOfArg ];

  assert( Arg[0] < 99999999999.9 );

  for( int i = 0 ; i < numOfArg ; ++i ) 
    checked[ i ] = 0;
  
  for( int i = 0; i < numOfOrd; ++i )
  {
    valueBest = 99999999999.9;
    for( int j = 0; j < numOfArg; ++j )
    {
      if( ( Arg[j] < valueBest ) && checked[j]==0){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    indexOrderd[ i ]=indexBest;
    checked[ indexBest ]=1;
  }

  delete [] checked;
}


void TSort::Index_B( double* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest = 0; 
  double valueBest;
  int *checked;
  checked = new int [ numOfArg ];

  assert( Arg[0] > -99999999999.9 );

  for( int i = 0 ; i < numOfArg ; ++i ) 
    checked[ i ] = 0;
  
  for( int i = 0; i < numOfOrd; ++i )
  {
    valueBest = -99999999999.9;
    for( int j = 0; j < numOfArg; ++j )
    {
      if( ( Arg[j] > valueBest ) && checked[j]==0){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    indexOrderd[ i ]=indexBest;
    checked[ indexBest ]=1;
  }

  delete [] checked;
}


void TSort::Index( int* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest = 0;
  int valueBest;
  int *checked;
  checked = new int [ numOfArg ];

  assert( Arg[0] < 99999999 );

  for( int i = 0 ; i < numOfArg ; ++i ) 
    checked[ i ] = 0;
  
  for( int i = 0; i < numOfOrd; ++i )
  {
    valueBest = 99999999;
    for( int j = 0; j < numOfArg; ++j )
    {
      if( ( Arg[j] < valueBest ) && checked[j]==0){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    indexOrderd[ i ]=indexBest;
    checked[ indexBest ]=1;
  }

  delete [] checked;
}


void TSort::Index_B( int* Arg, int numOfArg, int* indexOrderd, int numOfOrd )
{
  int indexBest = 0;
  int valueBest;
  int *checked;
  checked = new int [ numOfArg ];

  assert( Arg[0] > -999999999 );

  for( int i = 0 ; i < numOfArg ; ++i ) 
    checked[ i ] = 0;
  
  for( int i = 0; i < numOfOrd; ++i )
  {
    valueBest = -999999999;
    for( int j = 0; j < numOfArg; ++j )
    {
      if( ( Arg[j] > valueBest ) && checked[j]==0){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    indexOrderd[ i ]=indexBest;
    checked[ indexBest ]=1;
  }

  delete [] checked;
}


void TSort::Sort( int* Arg, int numOfArg )
{
  int indexBest = 0;
  int valueBest;
  int stock;

  assert( Arg[0] < 99999999 );

  for( int i = 0; i < numOfArg; ++i )
  {
    valueBest = 99999999;
    for( int j = i; j < numOfArg; ++j )
    {
      if( ( Arg[j] < valueBest ) ){
	valueBest = Arg[j];
	indexBest = j;
      }
    }
    stock = Arg[ i ];
    Arg[ i ] = Arg[ indexBest ];
    Arg[ indexBest ] = stock;
  }
}

TKopt::TKopt( int N )
{
  fN = N;

  fLink = new int* [ fN ];
  for( int i = 0; i < fN; ++i ) 
    fLink[ i ] = new int [ 2 ];
  fOrdCity = new int [ fN ];
  fOrdSeg = new int [ fN ];
  fSegCity = new int [ fN ];
  fOrient = new int [ fN ];
  fLinkSeg = new int* [ fN ];
  for( int i = 0; i < fN; ++i ) 
    fLinkSeg[ i ] = new int [ 2 ];
  fSizeSeg = new int [ fN ];
  fCitySeg = new int* [ fN ];
  for( int i = 0; i < fN; ++i ) 
    fCitySeg[ i ] = new int [ 2 ];

  fT = new int [ 5 ];

  fActiveV = new int [ fN ];
  fInvNearList = new int* [ fN ];
  for( int i = 0; i < fN; ++i ) 
    fInvNearList[ i ] = new int [ 500 ];
  fNumOfINL = new int [ fN ];

  fArray = new int [ fN+2 ]; 
  fCheckN = new int [ fN ];
  fB = new int [ fN ];
  fGene = new int [ fN ];
}

TKopt::~TKopt()
{
  for( int i = 0; i < fN; ++i ) 
    delete [] fLink[ i ];
  delete [] fLink;

  for( int i = 0; i < fN; ++i ) 
    delete [] fLinkSeg[ i ];
  delete [] fLinkSeg;

  for( int i = 0; i < fN; ++i ) 
    delete [] fCitySeg[ i ];
  delete [] fCitySeg;

  for( int i = 0; i < fN; ++i ) 
    delete [] fInvNearList[ i ];
  delete [] fInvNearList;

  delete [] fOrdCity;
  delete [] fOrdSeg;
  delete [] fSegCity;
  delete [] fOrient;
  delete [] fSizeSeg;
  delete [] fT;
  delete [] fActiveV;
  delete [] fNumOfINL;
  delete [] fArray;
  delete [] fCheckN;
  delete [] fB;
  delete [] fGene;
}


void TKopt::SetInvNearList()
{
  assert( eval->fNearNumMax == 50 );

  int c;
  for( int i = 0; i < fN; ++i ) 
    fNumOfINL[ i ] = 0;

  for( int i = 0; i < fN; ++i ){ 
    for( int k = 0; k < 50; ++k ){ 
      c = eval->fNearCity[i][k];
      if( fNumOfINL[c] < 500 )
	fInvNearList[ c ][ fNumOfINL[c]++ ] = i;  
      else{
	printf( "Check fNumOfINL[c] < 500 ) in kopt.cpp \n" );
	fflush( stdout );
      }
    }
  } 
}


void TKopt::TransIndiToTree( TIndi& indi )
{
  int num;
  int size;
  int orient;

  fArray[1] = 0; 
  for( int i = 2; i <= fN; ++i )
    fArray[i] = indi.fLink[ fArray[i-1] ][ 1 ]; 
  fArray[0] = fArray[fN]; 
  fArray[fN+1] = fArray[1]; 

  num = 1;
  fNumOfSeg = 0;

  while(1){
    orient = 1;
    size = 0;

    fOrient[ fNumOfSeg ] = orient;
    fOrdSeg[ fNumOfSeg ] = fNumOfSeg;

    fLink[ fArray[ num ] ][ 0 ] = -1;
    fLink[ fArray[ num ] ][ 1 ] = fArray[ num+1 ];
    fOrdCity[ fArray[ num ] ] = size;
    fSegCity[ fArray[ num ] ] = fNumOfSeg;
    fCitySeg[ fNumOfSeg ][ this->Turn(orient) ] = fArray[ num ];
    ++num;
    ++size;

    for( int i = 0; i < (int)sqrt( fN )-1; ++i )
    {
      if( num == fN )
	break;
      fLink[ fArray[ num ] ][ 0 ] = fArray[ num-1 ];
      fLink[ fArray[ num ] ][ 1 ] = fArray[ num+1 ];
      fOrdCity[ fArray[ num ] ] = size;
      fSegCity[ fArray[ num ] ] = fNumOfSeg;
      ++num;
      ++size;
    }

    if( num == fN-1 ){
      fLink[ fArray[ num ] ][ 0 ] = fArray[ num-1 ];
      fLink[ fArray[ num ] ][ 1 ] = fArray[ num+1 ];
      fOrdCity[ fArray[ num ] ] = size;
      fSegCity[ fArray[ num ] ] = fNumOfSeg;
      ++num;
      ++size;
    }
    
    fLink[ fArray[ num ] ][ 0 ] = fArray[ num-1 ];
    fLink[ fArray[ num ] ][ 1 ] = -1;
    fOrdCity[ fArray[ num ] ] = size;
    fSegCity[ fArray[ num ] ] = fNumOfSeg;
    fCitySeg[ fNumOfSeg ][ orient ] = fArray[ num ];
    ++num;
    ++size;

    fSizeSeg[ fNumOfSeg ] = size;
    ++fNumOfSeg; 

    if( num == fN+1 )
      break;
  }


  for( int s = 1; s < fNumOfSeg-1; ++s ){
    fLinkSeg[ s ][ 0 ] = s-1;
    fLinkSeg[ s ][ 1 ] = s+1;
  }
  fLinkSeg[ 0 ][ 0 ] = fNumOfSeg-1;
  fLinkSeg[ 0 ][ 1 ] = 1;
  fLinkSeg[ fNumOfSeg-1 ][ 0 ] = fNumOfSeg-2;
  fLinkSeg[ fNumOfSeg-1 ][ 1 ] = 0;

  fTourLength = indi.fEvaluationValue;
  fFixNumOfSeg = fNumOfSeg;
}


void TKopt::TransTreeToIndi( TIndi& indi )
{
  int t_p, t_n;
  for( int t = 0; t < fN; ++t )
  {
    t_p = this->GetPrev( t );
    t_n = this->GetNext( t );
    
    indi.fLink[ t ][ 0 ] = t_p;
    indi.fLink[ t ][ 1 ] = t_n;
  }
  eval->DoIt( indi );
}


void TKopt::DoIt( TIndi& tIndi )
{ 
  this->TransIndiToTree( tIndi );
  //  this->CheckDetail();           // Check
  //  this->CheckValid();            // Check
  this->Sub();
  this->TransTreeToIndi( tIndi );
}


void TKopt::Sub()
{
  int t1_st; 
  int t_next;
  int dis1, dis2;
  
  for( int t = 0; t < fN; ++t ) 
    fActiveV[ t ] = 1;
  
LLL1: t1_st = rand()%fN;
  fT[1] = t1_st;

  while(1)   // t1's loop
  {
    fT[1] = this->GetNext( fT[1] );
    if( fActiveV[ fT[1] ] == 0 )
      goto EEE;
    
    ////
    fFlagRev = 0;
    fT[2] = this->GetPrev( fT[1] );
    for( int num1 = 1; num1 < 50; ++num1 )
    {
      fT[4] = eval->fNearCity[ fT[1] ][ num1 ]; 
      fT[3] = this->GetPrev( fT[4] );
      dis1 = eval->fEdgeDis[fT[1]][fT[2]] - eval->fEdgeDis[fT[1]][fT[4]];

      if( dis1 > 0 ){
	dis2 = dis1 + 
               eval->fEdgeDis[fT[3]][fT[4]] - eval->fEdgeDis[fT[3]][fT[2]];
 
	if( dis2 > 0 ){
	  this->IncrementImp( fFlagRev );

	  for( int a = 1; a <= 4; ++a )
	    for( int k = 0; k < fNumOfINL[fT[a]]; ++k )
	      fActiveV[ this->fInvNearList[fT[a]][k] ] = 1;
	  
	  goto LLL1;
	}
      }
      else break;
    }

    ////
    fFlagRev = 1;
    fT[2] = this->GetNext( fT[1] );
    for( int num1 = 1; num1 < 50; ++num1 )
    {
      fT[4] = eval->fNearCity[ fT[1] ][ num1 ]; 
      fT[3] = this->GetNext( fT[4] );
      dis1 = eval->fEdgeDis[fT[1]][fT[2]] - eval->fEdgeDis[fT[1]][fT[4]];

      if( dis1 > 0 ){
	dis2 = dis1 + 
               eval->fEdgeDis[fT[3]][fT[4]] - eval->fEdgeDis[fT[3]][fT[2]];
 
	if( dis2 > 0 ){
	  this->IncrementImp( fFlagRev );
	  
	  for( int a = 1; a <= 4; ++a )
	    for( int k = 0; k < fNumOfINL[fT[a]]; ++k )
	      fActiveV[ this->fInvNearList[fT[a]][k] ] = 1;

	  goto LLL1;
	}
      }
      else break;
    }

    fActiveV[ fT[1] ] = 0;
  EEE:;
    if( fT[1] == t1_st ) 
      break;
  }
}


int TKopt::GetNext( int t )
{
  int t_n, seg, orient;

  seg = fSegCity[ t ];
  orient = fOrient[ seg ];

  t_n = fLink[ t ][ orient ];
  if( t_n == -1 ){
    seg = fLinkSeg[ seg ][ orient ];
    orient = Turn( fOrient[ seg ] );
    t_n = fCitySeg[ seg ][ orient ];
  }
  return t_n;
}


int TKopt::GetPrev( int t )
{
  int t_p, seg, orient;

  seg = fSegCity[ t ];
  orient = fOrient[ seg ];

  t_p = fLink[ t ][ this->Turn( orient ) ];
  if( t_p == -1 ){
    seg = fLinkSeg[ seg ][ Turn(orient) ];
    orient = fOrient[ seg ];
    t_p = fCitySeg[ seg ][ orient ];
  }
  return t_p;
}

void TKopt::Swap(int &a,int &b)
{
  int s;
  s=a;
  a=b;
  b=s;
}

int TKopt::Turn( int &orient )
{
  assert( orient == 0 || orient == 1 );
  if( orient == 0 )
    return 1;
  else if( orient == 1 )
    return 0;
  else 
    assert( 1 == 2 );
}

void TKopt::IncrementImp( int flagRev )
{
  int t1_s, t1_e, t2_s, t2_e;
  int seg_t1_s, seg_t1_e, seg_t2_s, seg_t2_e;
  int ordSeg_t1_s, ordSeg_t1_e, ordSeg_t2_s, ordSeg_t2_e;
  int orient_t1_s, orient_t1_e, orient_t2_s, orient_t2_e;
  int numOfSeg1, numOfSeg2;
  int curr;
  int ord;

  int flag_t2e_t1s;    
  int flag_t2s_t1e;    
  int length_t1s_seg;  
  int length_t1e_seg;
  int seg;

  // Seg1: b->d path
  // Seg2: c->a path

  if( fFlagRev == 0 ){
    t1_s = fT[1];
    t1_e = fT[3];
    t2_s = fT[4];
    t2_e = fT[2];
  }
  else if( fFlagRev == 1 ){
    t1_s = fT[2];
    t1_e = fT[4];
    t2_s = fT[3];
    t2_e = fT[1];
  }
  
  seg_t1_s = fSegCity[ t1_s ];
  ordSeg_t1_s = fOrdSeg[ seg_t1_s ];
  orient_t1_s = fOrient[ seg_t1_s ];
  seg_t1_e = fSegCity[ t1_e ];
  ordSeg_t1_e = fOrdSeg[ seg_t1_e ];
  orient_t1_e = fOrient[ seg_t1_e ];
  seg_t2_s = fSegCity[ t2_s ];
  ordSeg_t2_s = fOrdSeg[ seg_t2_s ];
  orient_t2_s = fOrient[ seg_t2_s ];
  seg_t2_e = fSegCity[ t2_e ];
  ordSeg_t2_e = fOrdSeg[ seg_t2_e ];
  orient_t2_e = fOrient[ seg_t2_e ];
  
  //////////////////// Type1 ////////////////////////
  if( ( seg_t1_s == seg_t1_e ) && ( seg_t1_s == seg_t2_s ) && ( seg_t1_s == seg_t2_e ) ){ 

    if( (fOrient[seg_t1_s] == 1 && (fOrdCity[ t1_s ] > fOrdCity[ t1_e ])) || 
        (fOrient[seg_t1_s] == 0 && (fOrdCity[ t1_s ] < fOrdCity[ t1_e ]))){
      this->Swap( t1_s, t2_s );
      this->Swap( t1_e, t2_e );
      this->Swap( seg_t1_s, seg_t2_s );
      this->Swap( seg_t1_e, seg_t2_e );
      this->Swap( ordSeg_t1_s, ordSeg_t2_s );
      this->Swap( ordSeg_t1_e, ordSeg_t2_e );
      this->Swap( orient_t1_s, orient_t2_s );
      this->Swap( orient_t1_e, orient_t2_e );
    }

    curr = t1_s;
    ord = fOrdCity[ t1_e ];

    while(1)
    {
      this->Swap( fLink[curr][0], fLink[curr][1] );
      fOrdCity[ curr ] = ord;
      if( curr == t1_e )
	break;
      curr = fLink[curr][Turn(orient_t1_s)];
      if( orient_t1_s == 0 )
	++ord;
      else 
	--ord;
    }

    fLink[t2_e][orient_t1_s] = t1_e;
    fLink[t2_s][Turn(orient_t1_s)] = t1_s;
    fLink[t1_s][orient_t1_s] = t2_s;
    fLink[t1_e][Turn(orient_t1_s)] = t2_e;

    //    this->CheckDetail();              // Check
    //    this->CheckValid();               // Check
    return;
  }
  //////////////////// Type1 ///////////////////////


  if( ordSeg_t1_e >= ordSeg_t1_s )
    numOfSeg1 = ordSeg_t1_e - ordSeg_t1_s + 1;
  else
    numOfSeg1 = ordSeg_t1_e - ordSeg_t1_s + 1 + fNumOfSeg;
  if( ordSeg_t2_e >= ordSeg_t2_s )
    numOfSeg2 = ordSeg_t2_e - ordSeg_t2_s + 1;
  else 
    numOfSeg2 = ordSeg_t2_e - ordSeg_t2_s + 1 + fNumOfSeg;

  if( numOfSeg1 > numOfSeg2 ){
    this->Swap( numOfSeg1, numOfSeg2 );
    this->Swap( t1_s, t2_s );
    this->Swap( t1_e, t2_e );
    this->Swap( seg_t1_s, seg_t2_s );
    this->Swap( seg_t1_e, seg_t2_e );
    this->Swap( ordSeg_t1_s, ordSeg_t2_s );
    this->Swap( ordSeg_t1_e, ordSeg_t2_e );
    this->Swap( orient_t1_s, orient_t2_s );
    this->Swap( orient_t1_e, orient_t2_e );
  }

  if( fLink[ t2_e ][ orient_t2_e ] == -1 ) flag_t2e_t1s = 1;
  else flag_t2e_t1s = 0;
  if( fLink[ t2_s ][ this->Turn(orient_t2_s) ] == -1 ) flag_t2s_t1e = 1;
  else flag_t2s_t1e = 0;

  length_t1s_seg = abs( fOrdCity[ t2_e ] 
                      - fOrdCity[ fCitySeg[ seg_t2_e ][ orient_t2_e ] ] );
  length_t1e_seg = abs( fOrdCity[ t2_s ] 
          - fOrdCity[ fCitySeg[ seg_t2_s ][ this->Turn(orient_t2_s) ] ] );
  
  ///////////////////// Type2 /////////////////
  if( seg_t1_s == seg_t1_e )
  {
    if( flag_t2e_t1s == 1 && flag_t2s_t1e == 1 )
    {
      orient_t1_s = Turn( fOrient[ seg_t1_s ] ); 
      fOrient[ seg_t1_s ] = orient_t1_s;
      fCitySeg[ seg_t1_s ][ orient_t1_s ] = t1_s;
      fCitySeg[ seg_t1_s ][ Turn(orient_t1_s) ] = t1_e;
      fLinkSeg[ seg_t1_s ][ orient_t1_s ] = seg_t2_s;
      fLinkSeg[ seg_t1_s ][ Turn(orient_t1_s) ] = seg_t2_e;

      //      this->CheckDetail();              // Check
      //      this->CheckValid();               // Check
      return;
    }
    
    if( flag_t2e_t1s == 0 && flag_t2s_t1e == 1 )
    {
      curr = t1_e;
      ord = fOrdCity[ t1_s ];
      while(1)
      {
	this->Swap( fLink[curr][0], fLink[curr][1] );
	fOrdCity[ curr ] = ord;
	if( curr == t1_s )
	  break;
	curr = fLink[curr][orient_t2_e];
	if( orient_t2_e == 0 )
	  --ord;
	else
	  ++ord;
      }

      fLink[t2_e][orient_t2_e] = t1_e;
      fLink[t1_s][orient_t2_e] = -1;
      fLink[t1_e][Turn(orient_t2_e)] = t2_e;
      fCitySeg[seg_t2_e][orient_t2_e] = t1_s;

      //      this->CheckDetail();              // Check
      //      this->CheckValid();               // Check

      return;
    }

    if( flag_t2e_t1s == 1 && flag_t2s_t1e == 0 )
    {
      curr = t1_s;
      ord = fOrdCity[ t1_e ];
      while(1)
      {
	this->Swap( fLink[curr][0], fLink[curr][1] );
	fOrdCity[ curr ] = ord;
	if( curr == t1_e )  
	  break;
	curr = fLink[curr][Turn(orient_t2_s)];
	if( orient_t2_s == 0 )
	  ++ord;
	else
	  --ord;
      }

      fLink[t2_s][Turn(orient_t2_s)] = t1_s;
      fLink[t1_e][Turn(orient_t2_s)] = -1;
      fLink[t1_s][orient_t2_s] = t2_s;
      fCitySeg[seg_t2_s][Turn(orient_t2_s)] = t1_e;

      //      this->CheckDetail();              // Check
      //      this->CheckValid();               // Check

      return;
    }
  }


  ///////////////////// Type3 /////////////////

  if( flag_t2e_t1s == 1 ){
    fLinkSeg[seg_t1_s][Turn(orient_t1_s)] = seg_t2_s;
  }
  else
  {
    seg_t1_s = fNumOfSeg++;
    orient_t1_s = orient_t2_e;
    fLink[ t1_s ][Turn(orient_t1_s)] = -1;
    fLink[ fCitySeg[seg_t2_e][orient_t2_e]][orient_t1_s] = -1;
    fOrient[seg_t1_s] = orient_t1_s;
    fSizeSeg[seg_t1_s] = length_t1s_seg;
    fCitySeg[seg_t1_s][Turn(orient_t1_s)] = t1_s;   
    fCitySeg[seg_t1_s][orient_t1_s] = fCitySeg[seg_t2_e][orient_t2_e];
    fLinkSeg[seg_t1_s][Turn(orient_t1_s)] = seg_t2_s;
    fLinkSeg[seg_t1_s][orient_t1_s] = fLinkSeg[seg_t2_e][orient_t2_e];
    seg = fLinkSeg[seg_t2_e][orient_t2_e];
    fLinkSeg[seg][Turn(fOrient[seg])] = seg_t1_s;
  }

  if( flag_t2s_t1e == 1 ){
    fLinkSeg[seg_t1_e][orient_t1_e] = seg_t2_e;
  }
  else
  {
    seg_t1_e = fNumOfSeg++;
    orient_t1_e = orient_t2_s;
    fLink[ t1_e ][orient_t1_e] = -1;
    fLink[ fCitySeg[seg_t2_s][Turn(orient_t2_s)] ][Turn(orient_t1_e)] = -1;
    fOrient[seg_t1_e] = orient_t1_e;
    fSizeSeg[seg_t1_e] = length_t1e_seg;
    fCitySeg[seg_t1_e][orient_t1_e] = t1_e;   
    fCitySeg[seg_t1_e][Turn(orient_t1_e)] = fCitySeg[seg_t2_s][Turn(orient_t2_s)];
    fLinkSeg[seg_t1_e][orient_t1_e] = seg_t2_e;
    fLinkSeg[seg_t1_e][Turn(orient_t1_e)] = fLinkSeg[seg_t2_s][Turn(orient_t2_s)];
    seg = fLinkSeg[seg_t2_s][Turn(orient_t2_s)];
    fLinkSeg[seg][fOrient[seg]] = seg_t1_e;
  }

  fLink[t2_e][orient_t2_e] = -1;
  fSizeSeg[seg_t2_e] -= length_t1s_seg;
  fCitySeg[seg_t2_e][orient_t2_e] = t2_e;
  fLinkSeg[seg_t2_e][orient_t2_e] = seg_t1_e;
  fLink[t2_s][Turn(orient_t2_s)] = -1;
  fSizeSeg[seg_t2_s] -= length_t1e_seg;
  fCitySeg[seg_t2_s][Turn(orient_t2_s)] = t2_s;
  fLinkSeg[seg_t2_s][Turn(orient_t2_s)] = seg_t1_s;

  seg = seg_t1_e;
  while(1)
  {
    fOrient[seg] = Turn(fOrient[seg]); 
    if( seg == seg_t1_s )
      break;
    seg = fLinkSeg[seg][fOrient[seg]];
  }
  

  if( fSizeSeg[seg_t2_e] < length_t1s_seg )
  {  
    seg = fLinkSeg[seg_t2_e][Turn(fOrient[seg_t2_e])];
    fLinkSeg[seg][fOrient[seg]] = seg_t1_s;
    seg = fLinkSeg[seg_t2_e][fOrient[seg_t2_e]];
    fLinkSeg[seg][Turn(fOrient[seg])] = seg_t1_s;
    seg = fLinkSeg[seg_t1_s][Turn(fOrient[seg_t1_s])];
    fLinkSeg[seg][fOrient[seg]] = seg_t2_e;
    seg = fLinkSeg[seg_t1_s][fOrient[seg_t1_s]];
    fLinkSeg[seg][Turn(fOrient[seg])] = seg_t2_e;
    
    this->Swap( fOrient[seg_t2_e], fOrient[seg_t1_s] );
    this->Swap( fSizeSeg[seg_t2_e], fSizeSeg[seg_t1_s] );
    this->Swap( fCitySeg[seg_t2_e][0], fCitySeg[seg_t1_s][0] );
    this->Swap( fCitySeg[seg_t2_e][1], fCitySeg[seg_t1_s][1] );
    this->Swap( fLinkSeg[seg_t2_e][0], fLinkSeg[seg_t1_s][0] );
    this->Swap( fLinkSeg[seg_t2_e][1], fLinkSeg[seg_t1_s][1] );
    this->Swap( seg_t2_e, seg_t1_s );
  }

  if( fSizeSeg[seg_t2_s] < length_t1e_seg )
  {  
    seg = fLinkSeg[seg_t2_s][Turn(fOrient[seg_t2_s])];
    fLinkSeg[seg][fOrient[seg]] = seg_t1_e;
    seg = fLinkSeg[seg_t2_s][fOrient[seg_t2_s]];
    fLinkSeg[seg][Turn(fOrient[seg])] = seg_t1_e;
    seg = fLinkSeg[seg_t1_e][Turn(fOrient[seg_t1_e])];
    fLinkSeg[seg][fOrient[seg]] = seg_t2_s;
    seg = fLinkSeg[seg_t1_e][fOrient[seg_t1_e]];
    fLinkSeg[seg][Turn(fOrient[seg])] = seg_t2_s;

    this->Swap( fOrient[seg_t2_s], fOrient[seg_t1_e] );
    this->Swap( fSizeSeg[seg_t2_s], fSizeSeg[seg_t1_e] );
    this->Swap( fCitySeg[seg_t2_s][0], fCitySeg[seg_t1_e][0] );
    this->Swap( fCitySeg[seg_t2_s][1], fCitySeg[seg_t1_e][1] );
    this->Swap( fLinkSeg[seg_t2_s][0], fLinkSeg[seg_t1_e][0] );
    this->Swap( fLinkSeg[seg_t2_s][1], fLinkSeg[seg_t1_e][1] );
    this->Swap( seg_t2_s, seg_t1_e );
  }


  while( fNumOfSeg > fFixNumOfSeg )
  {
    if( fSizeSeg[ fLinkSeg[fNumOfSeg-1][0] ] < 
        fSizeSeg[ fLinkSeg[fNumOfSeg-1][1] ] )  
      this->CombineSeg( fLinkSeg[fNumOfSeg-1][0], fNumOfSeg-1 );
    else 
      this->CombineSeg( fLinkSeg[fNumOfSeg-1][1], fNumOfSeg-1 );
  }

  int ordSeg = 0;
  seg = 0;

  while(1)
  {
    fOrdSeg[seg] = ordSeg;
    ++ordSeg;

    seg = fLinkSeg[seg][ fOrient[seg] ];
    if( seg == 0 )
      break;
  }

  // this->CheckDetail();              // Check
  // this->CheckValid();               // Check
  return;
}


void TKopt::CombineSeg( int segL, int segS )
{
  int seg;
  int t_s, t_e, direction; t_s = 0; t_e = 0; direction = 0;
  int ord; ord = 0;
  int increment; increment = 0;
  int curr, next;

  if( fLinkSeg[segL][fOrient[segL]] == segS ){ 
    fLink[fCitySeg[segL][fOrient[segL]]][fOrient[segL]] = 
    fCitySeg[segS][Turn(fOrient[segS])]; 

    fLink[fCitySeg[segS][Turn(fOrient[segS])]][Turn(fOrient[segS])] = fCitySeg[segL][fOrient[segL]];
    ord = fOrdCity[fCitySeg[segL][fOrient[segL]]]; 

    fCitySeg[segL][fOrient[segL]] = fCitySeg[segS][fOrient[segS]]; 
    fLinkSeg[segL][fOrient[segL]] = fLinkSeg[segS][fOrient[segS]]; 
    seg = fLinkSeg[segS][fOrient[segS]]; 
    fLinkSeg[seg][Turn(fOrient[seg])] = segL;

    t_s = fCitySeg[segS][Turn(fOrient[segS])];
    t_e = fCitySeg[segS][fOrient[segS]];
    direction = fOrient[segS];


    if( fOrient[segL] == 1 )
      increment = 1;
    else 
      increment = -1;
  }
  else if( fLinkSeg[segL][Turn(fOrient[segL])] == segS ){ 
    fLink[fCitySeg[segL][Turn(fOrient[segL])]][Turn(fOrient[segL])] = 
    fCitySeg[segS][fOrient[segS]]; 

    fLink[fCitySeg[segS][fOrient[segS]]][fOrient[segS]] = fCitySeg[segL][Turn(fOrient[segL])];
    ord = fOrdCity[fCitySeg[segL][Turn(fOrient[segL])]]; 

    fCitySeg[segL][Turn(fOrient[segL])] = fCitySeg[segS][Turn(fOrient[segS])]; 
    fLinkSeg[segL][Turn(fOrient[segL])] = fLinkSeg[segS][Turn(fOrient[segS])]; 
    seg = fLinkSeg[segS][Turn(fOrient[segS])]; 
    fLinkSeg[seg][fOrient[seg]] = segL;

    t_s = fCitySeg[segS][fOrient[segS]];
    t_e = fCitySeg[segS][Turn(fOrient[segS])];
    direction = Turn(fOrient[segS]);


    if( fOrient[segL] == 1 )
      increment = -1;
    else 
      increment = 1;
  }
  curr = t_s;
  ord = ord + increment;

  while(1)
  {
    fSegCity[curr] = segL;
    fOrdCity[curr] = ord;

    next = fLink[curr][direction];
    if( fOrient[segL] != fOrient[segS] )
      this->Swap( fLink[curr][0], fLink[curr][1] ); 

    if( curr == t_e )
      break;
    curr = next;
    ord += increment;
  }
  fSizeSeg[segL] += fSizeSeg[segS];
  --fNumOfSeg;
}



void TKopt::CheckDetail()
{
  int seg, seg_p, seg_n;
  int ord, ord_p, ord_n;
  int orient;
  int curr;

  seg = 0;

  for( int s = 0; s < fNumOfSeg; ++s ){

    seg = s;
    orient = fOrient[ seg ];
    seg_p = fLinkSeg[ seg ][ this->Turn(orient) ];
    seg_n = fLinkSeg[ seg ][ orient ];

    ord = fOrdSeg[ seg ];
    ord_p = ord - 1 ;
    if( ord_p < 0 ) 
      ord_p = fNumOfSeg - 1;
    ord_n = ord + 1;
    if( ord_n >= fNumOfSeg ) 
      ord_n = 0; 

    assert( ord_p == fOrdSeg[seg_p] ); 
    assert( ord_n == fOrdSeg[seg_n] );

    curr = fCitySeg[ s ][ 0 ];
    int count = 0;

    while(1)
    {
      ++count; 
      if( curr == fCitySeg[ s ][1] )
	break;
      curr = fLink[curr][1];
      assert( curr != -1 );
    }    
    assert( count == fSizeSeg[s] ); 
  }

  
  int t;
  int t_n, t_p, t_s, t_e;

  for( t = 0; t < fN; ++t )
  {
    seg = fSegCity[ t ];
    orient = fOrient[ seg ];
    t_s = fCitySeg[ seg ][ 0 ]; 
    t_e = fCitySeg[ seg ][ 1 ]; 

    t_p = fLink[ t ][ 0 ];      
    t_n = fLink[ t ][ 1 ];

    if( t == t_s ){
      assert( t_p == -1 ); 
    }
    else{
      assert( t_p != -1 );       
      assert( t == fLink[ t_p ][ 1 ] );      
      assert( seg == fSegCity[ t_p ] );      
      assert( fOrdCity[ t ] == fOrdCity[ t_p ] + 1 );
    } 

    if( t == t_e ){
      assert( t_n == -1 ); 
    }
    else{
      assert( t_n != -1 );       
      assert( t == fLink[ t_n ][ 0 ] );      
      assert( seg == fSegCity[ t_n ] );      
      assert( fOrdCity[ t ] == fOrdCity[ t_n ] - 1 );
    } 
  }
}


void TKopt::CheckValid()
{
  int t_st, t_c, t_n, st;
  int count;
  int seg, orient;
  int Invalid = 0;

  for( int i = 0; i < fN; ++i )  
    fCheckN[ i ] = 0;

  t_st = rand() % fN;
  t_n = t_st;
  
  count = 0;
  while(1)
  {
    t_c = t_n;
    fCheckN[ t_c ] = 1;
    ++count;

    seg = fSegCity[ t_c ];
    orient = fOrient[ seg ];

    t_n = this->GetNext( t_c );

    if( t_n == t_st )
      break;

    if( count == fN+1 ){
      Invalid = 1;
      break;
    }
  }

  for( int i = 0; i < fN; ++i )  
    if( fCheckN[ i ] == 0 )
      Invalid = 1;


  if( Invalid == 1 ){
    printf( "Invalid \n" ); fflush( stdout );
    assert( 1 == 2 );
  }
  
}


void TKopt::MakeRandSol( TIndi& indi )
{
  int k, r;

  for( int j = 0; j < fN; ++j ) 
    fB[j] = j;

  for( int i = 0; i < fN; ++i )
  {  
    r = rand() % (fN-i);
    fGene[i] = fB[r];
    fB[r] = fB[fN-i-1];
  }
   
  for( int j2 = 1 ; j2 < fN-1; ++j2 )
  {
    indi.fLink[fGene[j2]][0] = fGene[j2-1];
    indi.fLink[fGene[j2]][1] = fGene[j2+1];
  }
  indi.fLink[fGene[0]][0] = fGene[fN-1];
  indi.fLink[fGene[0]][1] = fGene[1];  
  indi.fLink[fGene[fN-1]][0] = fGene[fN-2];
  indi.fLink[fGene[fN-1]][1] = fGene[0]; 

  eval->DoIt( indi );
}

TRandom* tRand = NULL;

void InitURandom()
{
  int seed;
  unsigned short seed16v[3];

  seed = 1111;     

  seed16v[0] = 100;
  seed16v[1] = 200;
  seed16v[2] = seed;

  seed48( seed16v );
  tRand = new TRandom();
  srand( seed );

  // srand((unsigned int)time(NULL));
}


void InitURandom( int dd )
{
  int seed;
  unsigned short seed16v[3];

  seed = dd;       

  seed16v[0] = 100;
  seed16v[1] = 200;
  seed16v[2] = seed;

  seed48( seed16v );
  tRand = new TRandom();
  srand( seed );
}


TRandom::TRandom()
{
}


TRandom::~TRandom()
{
}


int TRandom::Integer( int minNumber, int maxNumber )
{
  return minNumber + (int)(drand48() * (double)(maxNumber - minNumber + 1));
}


double TRandom::Double( double minNumber, double maxNumber )
{
  return minNumber + drand48() * (maxNumber - minNumber);
}


void TRandom::Permutation( int *array, int numOfElement, int numOfSample )
{
  int i,j,k,r;
  int *b;

  if( numOfElement <= 0 )
    return;

  b= new int[numOfElement];

  for(j=0;j<numOfElement;j++) b[j]=0;
  for(i=0;i<numOfSample;i++)
  {  
    r=rand()%(numOfElement-i);
    k=0;
    for(j=0;j<=r;j++)
    {
      while(b[k]==1)
      {
	k++;
      }
      k++;
    }
    array[i]=k-1;
    b[k-1]=1;
  }
  delete [] b;
}


double TRandom::NormalDistribution( double mu, double sigma )
{
  double U1,U2,X;
  double PI = 3.1415926;
  
  while( 1 ){
    U1 = this->Double( 0.0, 1.0 );
    if( U1 != 0.0 ) break;
  }
  U2 = this->Double( 0.0, 1.0 );

  X = sqrt(-2.0*log(U1)) * cos(2*PI*U2);
  return( mu + sigma*X );
}


void TRandom::Shuffle( int *array, int numOfElement )
{
  int *a;
  int *b;

  a = new int[numOfElement];
  b = new int[numOfElement];

  this->Permutation( b, numOfElement, numOfElement ); 

  for( int i = 0; i < numOfElement; ++i )
    a[ i ] = array[ i ]; 
  for( int i = 0; i < numOfElement; ++i )
    array[ i ] = a[ b[ i ] ]; 

  delete [] a;
  delete [] b;
}



TEvaluator::TEvaluator()
{
  fEdgeDis = NULL;
  fNearCity = NULL;
  Ncity = 0;
  fNearNumMax = 50;  
}

TEvaluator::~TEvaluator()
{
  for ( int i = 0; i < Ncity; ++i ) 
    delete[] fEdgeDis[ i ];
  delete[] fEdgeDis;
  for ( int i = 0; i < Ncity; ++i ) 
    delete[] fNearCity[ i ];
  delete[] fNearCity;
  
  delete [] x;
  delete [] y;
}

void TEvaluator::SetInstance( char filename[] )
{
  FILE* fp;
  int n;
  int *check;
  char word[ 80 ], type[ 80 ];

  fp = fopen( filename, "r" );


  ////// read instance //////
  while( 1 ){
    if( fscanf( fp, "%s", word ) == EOF )
      break;

    if( strcmp( word, "DIMENSION" ) == 0 ){
      fscanf( fp, "%s", word ); 
      assert( strcmp( word, ":" ) == 0 );
      fscanf( fp, "%d", &Ncity ); 
    } 

    if( strcmp( word, "EDGE_WEIGHT_TYPE" ) == 0 ){
      fscanf( fp, "%s", word ); 
      assert( strcmp( word, ":" ) == 0 );
      fscanf( fp, "%s", type ); 
    } 

    if( strcmp( word, "NODE_COORD_SECTION" ) == 0 ) 
      break;
      

  }
  if( strcmp( word, "NODE_COORD_SECTION" ) != 0 ){
    printf( "Error in reading the instance\n" );
    exit(0);
  }

  x = new double [ Ncity ]; 
  y = new double [ Ncity ]; 
  int checkedN[ Ncity ];

  int xi, yi; 
  for( int i = 0; i < Ncity; ++i ) 
  {
    fscanf( fp, "%d", &n );
    assert( i+1 == n ); 
    fscanf( fp, "%s", word ); 
    x[ i ] = atof( word );
    fscanf( fp, "%s", word ); 
    y[ i ] = atof( word );
  }

  fclose(fp);
  //////////////////////////

  fEdgeDis = new int* [ Ncity ];
  for( int i = 0; i < Ncity; ++i ) fEdgeDis[ i ] = new int [ Ncity ];
  fNearCity = new int* [ Ncity ];
  for( int i = 0; i < Ncity; ++i ) fNearCity[ i ] = new int [ fNearNumMax+1 ];

  /////////////// Remove duplicated nodes /////////////
  ////////////////////////////////////////////////////


  if( strcmp( type, "EUC_2D" ) == 0  ) {
    for( int i = 0; i < Ncity ; ++i )
    {
      for( int j = 0; j < Ncity ; ++j )
      {
	fEdgeDis[ i ][ j ]=(int)(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))+0.5);
      }
    }
  }
  else if( strcmp( type, "ATT" ) == 0  ) { 
    for( int i = 0; i < Ncity; ++i )
    {
      for( int j = 0; j < Ncity; ++j ) 
      {
	double r = (sqrt(((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]))/10.0));
	int t = (int)r;
	if( (double)t < r ) fEdgeDis[ i ][ j ] = t+1;
	else fEdgeDis[ i ][ j ] = t; 
      }
    }
  }
  else if( strcmp( type, "CEIL_2D" ) == 0  ) {  
    for( int i = 0; i < Ncity ; ++i )
    {
      for( int j = 0; j < Ncity ; ++j )
      {
	fEdgeDis[ i ][ j ]=(int)ceil(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])));
      }
    }
  }
  else{
    printf( "EDGE_WEIGHT_TYPE is not supported\n" );
    exit( 1 );
  }

  int ci;
  int j1 ,j2 ,j3;
  int city_num = 0;
  int min_dis;

  for( ci = 0; ci < Ncity; ++ci )
  {
    for( j3 = 0; j3 < Ncity; ++j3 ) checkedN[ j3 ] = 0;
    checkedN[ ci ] = 1;
    fNearCity[ ci ][ 0 ] = ci;
    for( j1 = 1; j1 <= fNearNumMax; ++j1 ) 
    {
      min_dis = 100000000;
      for( j2 = 0; j2 < Ncity; ++j2 )
      { 
	if( fEdgeDis[ ci ][ j2 ] <= min_dis && checkedN[ j2 ] == 0 )
	{
	  city_num = j2;
          min_dis = fEdgeDis[ ci ][ j2 ];
	}
      }
      fNearCity[ ci ][ j1 ] = city_num;
      checkedN[ city_num ] = 1;
    }
  }
}


void TEvaluator::DoIt( TIndi& indi )
{
  int d;
  d = 0;  
  for( int i = 0; i < Ncity; ++i )
  {  
    d = d + fEdgeDis[ i ][ indi.fLink[i][0] ];
    d = d + fEdgeDis[ i ][ indi.fLink[i][1] ];
  }
  indi.fEvaluationValue = d/2;
}


void TEvaluator::WriteTo( FILE* fp, TIndi& indi ) 
{
  assert( Ncity == indi.fN );
  int Array[ Ncity ];
  int curr, next, pre, st, count;

  count = 0;
  pre = -1;
  curr = 0;
  st = 0;
  while( 1 )
  {
    Array[ count++ ] = curr + 1;

    if( count > Ncity ){
      printf( "Invalid\n" );
      return;
    } 
 
    if( indi.fLink[ curr ][ 0 ] == pre )
      next = indi.fLink[ curr ][ 1 ];
    else 
      next = indi.fLink[ curr ][ 0 ];

    pre = curr;
    curr = next;
    if( curr == st )
      break;
  }

  if( this->CheckValid( Array, indi.fEvaluationValue ) == false ){
    printf( "Individual is invalid \n" );
  }

  fprintf( fp, "%d %d\n", indi.fN, indi.fEvaluationValue );
  for( int i = 0; i < indi.fN; ++i )
    fprintf( fp, "%d ", Array[ i ] );
  fprintf( fp, "\n" ); 
}


bool TEvaluator::ReadFrom( FILE* fp, TIndi& indi )
{
  assert( Ncity == indi.fN );
  int Array[ Ncity ];
  int curr, next, pre, st, count;
  int N, value;

  if( fscanf( fp, "%d %d", &N, &value ) == EOF ) 
    return false;
  assert( N == Ncity );
  indi.fN = N;
  indi.fEvaluationValue = value;

  for( int i = 0; i < Ncity; ++i ){ 
    if( fscanf( fp, "%d", &Array[ i ] ) == EOF )
      return false;
  }

  if( this->CheckValid( Array, indi.fEvaluationValue ) == false ){
    printf( "Individual is invalid \n" );
    return false;
  }

  for( int i = 0; i < Ncity; ++i ){ 
    Array[ i ] -= 1; 
  }

  for( int i = 1; i < Ncity-1; ++i ){ 
    indi.fLink[ Array[ i ] ][ 0 ] = Array[ i-1 ]; 
    indi.fLink[ Array[ i ] ][ 1 ] = Array[ i+1 ]; 
  }
  indi.fLink[ Array[ 0 ] ][ 0 ] = Array[ Ncity-1 ]; 
  indi.fLink[ Array[ 0 ] ][ 1 ] = Array[ 1 ]; 
  indi.fLink[ Array[ Ncity-1 ] ][ 0 ] = Array[ Ncity-2 ]; 
  indi.fLink[ Array[ Ncity-1 ] ][ 1 ] = Array[ 0 ]; 

  return true;
}    


bool TEvaluator::CheckValid( int* array, int value )
{
  int check[ Ncity ];
  int distance;

  for( int i = 0; i < Ncity; ++i ){
    check[ i ] = 0;
  }

  for( int i = 0; i < Ncity; ++i )
    ++check[ array[ i ]-1 ];

  for( int i = 0; i < Ncity; ++i ){
    if( check[ i ] != 1 ){
      return false;
    }
  }
    
  distance = 0;  
  for( int i = 0; i < Ncity-1; ++i ){
    distance += fEdgeDis[ array[ i ]-1 ][ array[ i+1 ]-1 ];
  }
  distance += fEdgeDis[ array[ Ncity-1 ]-1 ][ array[ 0 ]-1 ];
  if( distance != value ){
    return false;
  }
  return true;
}



TIndi::TIndi()
{                
  fN = 0;
  fLink = NULL;
  fEvaluationValue = 0;
}
 

TIndi::~TIndi()
{
  for ( int i = 0; i < fN; ++i ) 
    delete[] fLink[ i ];
  delete[] fLink;
}


void TIndi::Define( int N )
{
  fN = N;
  
  fLink = new int* [ fN ];
  for( int i = 0; i < fN; ++i ) 
    fLink[ i ] = new int [ 2 ];
} 


TIndi& TIndi::operator = ( const TIndi& src )
{
  fN = src.fN;

  for ( int i = 0; i < fN; ++i ) 
    for ( int j = 0; j < 2; ++j ) 
      fLink[i][j] = src.fLink[i][j];
  fEvaluationValue = src.fEvaluationValue;

  return *this;
}


bool TIndi::operator == ( const TIndi& src )
{
  int curr,next,pre;
  int flag_identify;

  if( fN != src.fN )  
    return false;
  if( fEvaluationValue != src.fEvaluationValue )  
    return false;
  
  curr = 0;
  pre = -1;
  for( int i = 0; i < fN; ++i )
  {
    if( fLink[curr][0] == pre ) 
      next = fLink[curr][1];
    else 
      next = fLink[curr][0];
	
    if( src.fLink[curr][0] != next && src.fLink[curr][1] != next ) 
    {
      return false;
    }

    pre = curr;    
    curr = next; 
  }

  return true;
}

 


TCross::TCross( int N )
{
  fMaxNumOfABcycle = 2000; /* Set an appropriate value (2000 is usually enough) */

  fN = N;
  tBestTmp.Define( fN );

  near_data = new int* [ fN ];
  for ( int j = 0; j < fN; ++j ) 
    near_data[j] = new int [ 5 ];

  fABcycle = new int* [ fMaxNumOfABcycle ];
  for ( int j = 0; j < fMaxNumOfABcycle; ++j ) 
    fABcycle[j] = new int [ 2*fN + 4 ];

  koritsu = new int [ fN ];
  bunki = new int [ fN ];
  kori_inv = new int [ fN ];
  bun_inv = new int [ fN ];
  check_koritsu = new int [ fN ];
  fRoute = new int [ 2*fN + 1 ];
  fPermu = new int [ fMaxNumOfABcycle ];

  fC = new int [ 2*fN+4 ];
  fJun = new int[ fN+ 1 ];
  fOrd1 = new int [ fN ];
  fOrd2 = new int [ fN ];

  // Speed Up Start
  fOrder = new int [ fN ];
  fInv = new int [ fN ];
  fSegment = new int* [ fN ];
  for ( int j = 0; j < fN; ++j ) 
    fSegment[ j ] = new int [ 2 ];
  fSegUnit = new int [ fN ]; 
  fSegPosiList = new int[ fN ];
  LinkAPosi = new int [ fN ];
  LinkBPosi = new int* [ fN ];
  for ( int j = 0; j < fN; ++j ) 
    LinkBPosi[ j ] = new int [ 2 ];
  fPosiSeg = new int [ fN ];
  fNumOfElementInUnit = new int [ fN ]; 
  fCenterUnit = new int [ fN ]; 
  for ( int j = 0; j < fN; ++j ) 
    fCenterUnit[ j ] = 0;
  fListOfCenterUnit = new int [ fN+2 ]; 
  fSegForCenter = new int [ fN ]; 
  fGainAB = new int [ fN ]; 
  fModiEdge = new int* [ fN ]; 				 
  for ( int j = 0; j < fN; ++j ) 
    fModiEdge[ j ] = new int [ 4 ]; 				 
  fBestModiEdge = new int* [ fN ]; 				 
  for ( int j = 0; j < fN; ++j ) 
    fBestModiEdge[ j ] = new int [ 4 ]; 				 
  fAppliedCylce = new int [ fN ];
  fBestAppliedCylce = new int [ fN ];
  // Speed Up End

  // Block2
  fNumOfElementINAB = new int [ fMaxNumOfABcycle ];
  fInEffectNode = new int* [ fN ];
  for( int i = 0; i < fN; ++i )
    fInEffectNode[ i ] = new int [ 2 ];
  fWeight_RR = new int* [ fMaxNumOfABcycle ];
  for( int i = 0; i < fMaxNumOfABcycle; ++i )
    fWeight_RR[ i ] = new int [ fMaxNumOfABcycle ];
  fWeight_SR = new int [ fMaxNumOfABcycle ];
  fWeight_C = new int [ fMaxNumOfABcycle ];
  fUsedAB = new int [ fN ];
  fMoved_AB = new int [ fN ];
  fABcycleInEset = new int [ fMaxNumOfABcycle ];
}

TCross::~TCross()
{
  delete [] koritsu;
  delete [] bunki;
  delete [] kori_inv;
  delete [] bun_inv;
  delete [] check_koritsu;
  delete [] fRoute;
  delete [] fPermu;

  for ( int j = 0; j < fN; ++j ) 
    delete[] near_data[ j ];
  delete[] near_data;

  for ( int j = 0; j < fMaxNumOfABcycle; ++j ) 
    delete[] fABcycle[ j ];
  delete[] fABcycle;

  delete [] fC;
  delete [] fJun; 
  delete [] fOrd1; 
  delete [] fOrd2; 


  // Speed Up Start
  delete [] fOrder;
  delete [] fInv;

  for ( int j = 0; j < fN; ++j ) 
    delete[] fSegment[ j ];
  delete[] fSegment;
  delete[] fSegUnit;
  delete [] fSegPosiList;
  delete [] LinkAPosi;
  for ( int j = 0; j < fN; ++j ) 
    delete[] LinkBPosi[ j ];
  delete [] LinkBPosi;
  delete [] fPosiSeg;
  delete [] fNumOfElementInUnit; 
  delete [] fCenterUnit;
  delete [] fListOfCenterUnit;
  delete [] fSegForCenter;
  delete [] fGainAB;

  for ( int j = 0; j < fN; ++j ) 
    delete[] fModiEdge[ j ];
  delete [] fModiEdge;
  for ( int j = 0; j < fN; ++j ) 
    delete[] fBestModiEdge[ j ];
  delete [] fBestModiEdge;
  delete [] fAppliedCylce;
  delete [] fBestAppliedCylce;
  // Speed Up End
  

  // Block2
  delete [] fNumOfElementINAB;
  for ( int j = 0; j < fN; ++j ) 
    delete [] fInEffectNode[ j ];
  delete [] fInEffectNode;
  for( int i = 0; i < fMaxNumOfABcycle; ++i )
    delete [] fWeight_RR[ i ];
  delete [] fWeight_SR;
  delete [] fWeight_C;
  delete [] fUsedAB;
  delete [] fMoved_AB;
  delete [] fABcycleInEset;
}

void TCross::SetParents( const TIndi& tPa1, const TIndi& tPa2, int flagC[ 10 ], int numOfKids )
{
  this->SetABcycle( tPa1, tPa2, flagC, numOfKids ); 

  fDis_AB = 0;   

  int curr, next, st, pre;
  st = 0;
  curr=-1;
  next = st;
  for( int i = 0; i < fN; ++i )
  {
    pre=curr;
    curr=next;
    if( tPa1.fLink[curr][0] != pre ) 
      next = tPa1.fLink[ curr ][ 0 ];
    else 
      next=tPa1.fLink[curr][1];
    
    if( tPa2.fLink[ curr ][ 0 ] != next && tPa2.fLink[ curr ][ 1 ] != next )  
      ++fDis_AB; 

    fOrder[ i ] = curr;
    fInv[ curr ] = i;
  }

  assert( next == st );

  if( flagC[ 1 ] == 2 ){           /* Block2 */
    fTmax = 10;                    /* Block2 */
    fMaxStag = 20;                 /* Block2 (1:Greedy LS, 20:Tabu Search) */
    this->SetWeight( tPa1, tPa2 ); /* Block2 */
  }
}


void TCross::DoIt( TIndi& tKid, TIndi& tPa2, int numOfKids, int flagP, int flagC[ 10 ], int **fEdgeFreq )
{
  int Num;     
  int jnum, centerAB; 
  int gain;
  int BestGain;  
  double pointMax, point;
  double DLoss;

  fEvalType = flagC[ 0 ];              /* 1:Greedy, 2:---, 3:Distance, 4:Entropy */
  fEsetType = flagC[ 1 ];              /* 1:Single-AB, 2:Block2 */

  assert( fEvalType == 1 || fEvalType == 3 || fEvalType == 4 );
  assert( fEsetType == 1 || fEsetType == 2 );

  if ( numOfKids <= fNumOfABcycle ) 
    Num = numOfKids;
  else 
    Num = fNumOfABcycle;

  if( fEsetType == 1 ){         /* Single-AB */
    tRand->Permutation( fPermu, fNumOfABcycle, fNumOfABcycle ); 
  }
  else if( fEsetType == 2 ){    /* Block2 */
    for( int k =0; k< fNumOfABcycle; ++k )
      fNumOfElementINAB[ k ] = fABcycle[ k ][ 0 ];
    tSort->Index_B( fNumOfElementINAB, fNumOfABcycle, fPermu, fNumOfABcycle );
  }


  fNumOfGeneratedCh = 0;
  pointMax = 0.0;
  BestGain = 0;       
  fFlagImp = 0;  

  for( int j =0; j < Num; ++j )
  { 
    fNumOfABcycleInEset = 0;
    if( fEsetType == 1 ){         /* Single-AB */
      jnum = fPermu[ j ];
      fABcycleInEset[ fNumOfABcycleInEset++ ] = jnum; 
    }
    else if( fEsetType == 2 ){    /* Block2 */
      jnum = fPermu[ j ];
      centerAB = jnum;
      for( int s = 0; s < fNumOfABcycle; ++s ){ 
	if( s == centerAB )
	  fABcycleInEset[ fNumOfABcycleInEset++ ] = s; 
	else{
	  if( fWeight_RR[ centerAB ][ s ] > 0 && fABcycle[ s ][ 0 ] < fABcycle[ centerAB ][ 0 ] ){
	    if( rand() %2 == 0 )
	      fABcycleInEset[ fNumOfABcycleInEset++ ] = s; 
	  }
	}
      }
      this->Search_Eset( centerAB );    
    }


    fNumOfSPL = 0;          
    gain = 0;               
    fNumOfAppliedCycle = 0; 
    fNumOfModiEdge = 0;	    

    fNumOfAppliedCycle = fNumOfABcycleInEset; 
    for( int k = 0; k < fNumOfAppliedCycle; ++k ){   
      fAppliedCylce[ k ] = fABcycleInEset[ k ]; 
      jnum = fAppliedCylce[ k ];
      this->ChangeSol( tKid, jnum, flagP );
      gain += fGainAB[ jnum ];                       
    }

    this->MakeUnit();                                   
    this->MakeCompleteSol( tKid );                      
    gain += fGainModi;                                  
    
    ++fNumOfGeneratedCh;

    if( fEvalType == 1 )       /* Greedy */
      DLoss = 1.0;
    else if( fEvalType == 2 ){  
      assert( 1 == 2 );
    }
    else if( fEvalType == 3 )  /* Distance preservation */
      DLoss = this->Cal_ADP_Loss( fEdgeFreq );     
    else if( fEvalType == 4 )  /* Entropy preservation */
      DLoss = this->Cal_ENT_Loss( fEdgeFreq );             

    if( DLoss <= 0.0 ) DLoss = 0.00000001; 

    point = (double)gain / DLoss;                         
    tKid.fEvaluationValue = tKid.fEvaluationValue - gain; 
    
    // if( pointMax < point ){
    if( pointMax < point && (2 * fBest_Num_E < fDis_AB || 
			     tKid.fEvaluationValue != tPa2.fEvaluationValue ) ){   
      pointMax = point;
      BestGain = gain;        
      fFlagImp = 1;  

      fNumOfBestAppliedCycle = fNumOfAppliedCycle;
      for( int s = 0; s < fNumOfBestAppliedCycle; ++s )
	fBestAppliedCylce[ s ] = fAppliedCylce[ s ];
      
      fNumOfBestModiEdge = fNumOfModiEdge;	  	
      for( int s = 0; s < fNumOfBestModiEdge; ++s ){
	fBestModiEdge[ s ][ 0 ] = fModiEdge[ s ][ 0 ];
	fBestModiEdge[ s ][ 1 ] = fModiEdge[ s ][ 1 ];
	fBestModiEdge[ s ][ 2 ] = fModiEdge[ s ][ 2 ];
	fBestModiEdge[ s ][ 3 ] = fModiEdge[ s ][ 3 ];
      }	

    }

    this->BackToPa1( tKid ); 
    tKid.fEvaluationValue = tKid.fEvaluationValue + gain;
  }

  if( fFlagImp == 1 ){           
    this->GoToBest( tKid ); 
    tKid.fEvaluationValue = tKid.fEvaluationValue - BestGain;
    this->IncrementEdgeFreq( fEdgeFreq );
  }
}


void TCross::SetABcycle( const TIndi& tPa1, const TIndi& tPa2, int flagC[ 10 ], int numOfKids )
{
  bunki_many=0; koritsu_many=0;
  for( int j = 0; j < fN ; ++j )
  {
    near_data[j][1]=tPa1.fLink[j][0];
    near_data[j][3]=tPa1.fLink[j][1];

    near_data[j][0] = 2;
    
    koritsu[koritsu_many]=j;
    koritsu_many++;

    near_data[j][2]=tPa2.fLink[j][0];
    near_data[j][4]=tPa2.fLink[j][1];
  }
  for(int j = 0; j < fN; ++j ) 
  {
    check_koritsu[j]=-1;
    kori_inv[koritsu[j]]=j;
  }

  /**************************************************/

  fNumOfABcycle=0; 
  flag_st=1;                   
  while(koritsu_many!=0)
  {                                                               
    if(flag_st==1)          
    {
      fPosiCurr=0;
      r=rand()%koritsu_many;
      st=koritsu[r];    
      check_koritsu[st]=fPosiCurr;
      fRoute[fPosiCurr]=st;
      ci=st;
      pr_type=2;
    }
    else if(flag_st==0)    
    {
      ci=fRoute[fPosiCurr];   
    }
                        
    flag_circle=0;
    while(flag_circle==0)
    {
      fPosiCurr++;
      pr=ci;
      
      switch(pr_type)
      {
      case 1:                 
	ci=near_data[pr][fPosiCurr%2+1];
	break;
      case 2:   
	r=rand()%2;
	ci=near_data[pr][fPosiCurr%2+1+2*r];
	if(r==0) this->Swap(near_data[pr][fPosiCurr%2+1],near_data[pr][fPosiCurr%2+3]);
	break;
      case 3:   
	ci=near_data[pr][fPosiCurr%2+3];
      }

      fRoute[fPosiCurr]=ci;
      
      if(near_data[ci][0]==2) 
      {   
	if(ci==st)            
	{        
	  if(check_koritsu[st]==0) 
	  {        
	    if((fPosiCurr-check_koritsu[st])%2==0)  
	    {                  
	      if(near_data[st][fPosiCurr%2+1]==pr)
	      {
		this->Swap(near_data[ci][fPosiCurr%2+1],near_data[ci][fPosiCurr%2+3]); 
	      }
	      st_appear = 1;
	      this->FormABcycle();
	      if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto LLL;
	      if( fNumOfABcycle == fMaxNumOfABcycle ) goto LLL;

	      flag_st=0;
	      flag_circle=1;
	      pr_type=1; 
	    }
	    else
	    {
	      this->Swap(near_data[ci][fPosiCurr%2+1],near_data[ci][fPosiCurr%2+3]); 
	      pr_type=2;
	    }
	    check_koritsu[st]=fPosiCurr;
	  } 
	  else                     
	  {         
	    st_appear = 2;
	    this->FormABcycle();
	    if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto LLL;
	    if( fNumOfABcycle == fMaxNumOfABcycle ) goto LLL;

	    flag_st=1;
	    flag_circle=1;
	  }
	}
	else if(check_koritsu[ci]==-1) 
	{
	  check_koritsu[ci]=fPosiCurr;
	  if(near_data[ci][fPosiCurr%2+1]==pr)
	  {
	    this->Swap(near_data[ci][fPosiCurr%2+1],near_data[ci][fPosiCurr%2+3]); 
	  }
	  pr_type=2;
	}
	else if(check_koritsu[ci]>0)   
	{
	  this->Swap(near_data[ci][fPosiCurr%2+1],near_data[ci][fPosiCurr%2+3]); 
	  if((fPosiCurr-check_koritsu[ci])%2==0)  
	  {
	    st_appear = 1;
	    this->FormABcycle();
	    if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto LLL;
	      if( fNumOfABcycle == fMaxNumOfABcycle ) goto LLL;
	      
	    flag_st=0;
	    flag_circle=1;
	    pr_type=1;
	  }
	  else
	  {
	    this->Swap(near_data[ci][(fPosiCurr+1)%2+1],near_data[ci][(fPosiCurr+1)%2+3]); 
	    pr_type=3;
	  }  
	}
      }
      else if(near_data[ci][0]==1)    
      {
	if(ci==st)                    
        {
	  st_appear = 1;
	  this->FormABcycle();
	  if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto LLL;
	  if( fNumOfABcycle == fMaxNumOfABcycle ) goto LLL;

	  flag_st=1;
	  flag_circle=1;
	}
	else pr_type=1;
      }
    }
  }
                                       
  while(bunki_many!=0)
  {            
    fPosiCurr=0;   
    r=rand()%bunki_many;
    st=bunki[r];
    fRoute[fPosiCurr]=st;
    ci=st;
    
    flag_circle=0;
    while(flag_circle==0)
    { 
      pr=ci; 
      fPosiCurr++;
      ci=near_data[pr][fPosiCurr%2+1]; 
      fRoute[fPosiCurr]=ci;
      if(ci==st)                       
      {
	st_appear = 1;
	this->FormABcycle();
	if( flagC[ 1 ] == 1 && fNumOfABcycle == numOfKids ) goto LLL;
	if( fNumOfABcycle == fMaxNumOfABcycle ) goto LLL;
	
	flag_circle=1;
      }
    }
  }

LLL: ;

  if( fNumOfABcycle == fMaxNumOfABcycle ){
    printf( "fMaxNumOfABcycle(%d) must be increased\n", fMaxNumOfABcycle );
    exit( 1 );
  }
}


void TCross::FormABcycle()
{
  int j;
  int st_count;
  int edge_type;
  int st,ci, stock;
  int cem;                   
  int diff;
 
  if(fPosiCurr%2==0) edge_type=1; 
  else edge_type=2;               
  st=fRoute[fPosiCurr];
  cem=0;
  fC[cem]=st;    

  st_count=0;
  while(1)
  {
    cem++;
    fPosiCurr--;
    ci=fRoute[fPosiCurr];
    if(near_data[ci][0]==2)
    {
      koritsu[kori_inv[ci]]=koritsu[koritsu_many-1];
      kori_inv[koritsu[koritsu_many-1]]=kori_inv[ci];
      koritsu_many--;
      bunki[bunki_many]=ci;
      bun_inv[ci]=bunki_many;
      bunki_many++;
    }
    else if(near_data[ci][0]==1)
    {
      bunki[bun_inv[ci]]=bunki[bunki_many-1];
      bun_inv[bunki[bunki_many-1]]=bun_inv[ci];
      bunki_many--;
    }
             
    near_data[ci][0]--;
    if(ci==st) st_count++;
    if(st_count==st_appear) break;
    fC[cem]=ci;  
  }

  if(cem==2)
    return;

  fABcycle[fNumOfABcycle][0]=cem;    

  if(edge_type==2)
  {
    stock=fC[0];
    for( int j=0;j<cem-1;j++) fC[j]=fC[j+1];
    fC[cem-1]=stock;
  }
  
  for( int j=0;j<cem;j++) 
    fABcycle[fNumOfABcycle][j+2]=fC[j];
  fABcycle[fNumOfABcycle][1]=fC[cem-1];
  fABcycle[fNumOfABcycle][cem+2]=fC[0];
  fABcycle[fNumOfABcycle][cem+3]=fC[1];

  fC[ cem ] = fC[ 0 ]; 
  fC[ cem+1 ] = fC[ 1 ]; 
  diff = 0;
  for( j = 0; j < cem/2; ++j ) 
  {
    diff = diff + eval->fEdgeDis[fC[2*j]][fC[1+2*j]]
                - eval->fEdgeDis[fC[1+2*j]][fC[2+2*j]];
  }
  fGainAB[fNumOfABcycle] = diff;
  ++fNumOfABcycle;
}


void TCross::Swap(int &a,int &b)
{
  int s;
  s=a;
  a=b;
  b=s;
}


void TCross::ChangeSol( TIndi& tKid, int ABnum, int type )
{
  int j;
  int cem,r1,r2,b1,b2;
  int po_r1, po_r2, po_b1, po_b2; 

  cem=fABcycle[ABnum][0];  
  fC[0]=fABcycle[ABnum][0];

  if(type==2)   
  {
    for(j=0;j<cem+3;j++) fC[cem+3-j]=fABcycle[ABnum][j+1];
  }
  else for(j=1;j<=cem+3;j++) fC[j]=fABcycle[ABnum][j];

  for(j=0;j<cem/2;j++)
  {                           
    r1=fC[2+2*j];r2=fC[3+2*j];
    b1=fC[1+2*j];b2=fC[4+2*j];

    if(tKid.fLink[r1][0]==r2)
      tKid.fLink[r1][0]=b1;
    else 
      tKid.fLink[r1][1]=b1;
    if(tKid.fLink[r2][0]==r1) 
      tKid.fLink[r2][0]=b2;
    else
      tKid.fLink[r2][1]=b2;   

    po_r1 = fInv[ r1 ]; 
    po_r2 = fInv[ r2 ]; 
    po_b1 = fInv[ b1 ]; 
    po_b2 = fInv[ b2 ]; 
    
    if( po_r1 == 0 && po_r2 == fN-1 )
      fSegPosiList[ fNumOfSPL++ ] = po_r1;
    else if( po_r1 == fN-1 && po_r2 == 0 )
      fSegPosiList[ fNumOfSPL++ ] = po_r2;
    else if( po_r1 < po_r2 )
      fSegPosiList[ fNumOfSPL++ ] = po_r2;
    else if( po_r2 < po_r1 )
      fSegPosiList[ fNumOfSPL++ ] = po_r1;
    else
      assert( 1 == 2 );
    
    LinkBPosi[ po_r1 ][ 1 ] = LinkBPosi[ po_r1 ][ 0 ];
    LinkBPosi[ po_r2 ][ 1 ] = LinkBPosi[ po_r2 ][ 0 ];
    LinkBPosi[ po_r1 ][ 0 ] = po_b1; 
    LinkBPosi[ po_r2 ][ 0 ] = po_b2; 
  }
}


void TCross::MakeCompleteSol( TIndi& tKid )
{
  int j,j1,j2,j3;
  int st,ci,pre,curr,next,a,b,c,d,aa,bb,a1,b1;
  int city_many;
  int remain_unit_many;
  int ucm;
  int unit_num;
  int min_unit_city; 
  int near_num;
  int unit_many;               
  int center_un;               
  int select_un;               
  int diff,max_diff;
  int count;      
  int nearMax;

  fGainModi = 0;         

  while( fNumOfUnit != 1 )
  {    
    min_unit_city = fN + 12345;
    for( int u = 0; u < fNumOfUnit; ++u ) 
    {
      if( fNumOfElementInUnit[ u ] < min_unit_city )
      {
	center_un = u;
        min_unit_city = fNumOfElementInUnit[ u ];
      }
    }  

    st = -1;
    fNumOfSegForCenter = 0;   
    for( int s = 0; s < fNumOfSeg; ++s ){
      if( fSegUnit[ s ] == center_un ){
	int posi = fSegment[ s ][ 0 ];
	st = fOrder[ posi ];    
	fSegForCenter[  fNumOfSegForCenter++ ] = s; 
      }
    } 
    assert( st != -1 );

    curr = -1;
    next = st;
    fNumOfElementInCU = 0;
    while(1){ 
      pre = curr;
      curr = next;
      fCenterUnit[ curr ] = 1;     
      fListOfCenterUnit[ fNumOfElementInCU ] = curr;
      ++fNumOfElementInCU;

      if( tKid.fLink[ curr ][ 0 ] != pre )
	next = tKid.fLink[ curr ][ 0 ];
      else 
	next = tKid.fLink[ curr ][ 1 ]; 

      if( next == st ) break;
    }       
    fListOfCenterUnit[ fNumOfElementInCU ] = fListOfCenterUnit[ 0 ];
    fListOfCenterUnit[ fNumOfElementInCU+1 ] = fListOfCenterUnit[ 1 ];

    assert( fNumOfElementInCU == fNumOfElementInUnit[ center_un ] );

    max_diff = -999999999;
    a1 = -1; b1 = -1;
    nearMax = 10; /* N_near (see Step 5.3 in Section 2.2 of the Online Supplement) */
    /* nearMax must be smaller than or equal to eva->fNearNumMax (kopt.cpp ) */

  RESTART:;
    for( int s = 1; s <= fNumOfElementInCU; ++s )  
    { 
      a = fListOfCenterUnit[ s ];

      for( near_num = 1; near_num <= nearMax; ++near_num )   
      {
	c = eval->fNearCity[ a ][ near_num ];
	if( fCenterUnit[ c ] == 0 )   
	{
	  for( j1 = 0; j1 < 2; ++j1 )
	  {
	    b = fListOfCenterUnit[ s-1+2*j1 ];
            for( j2 = 0; j2 < 2; ++j2 )
	    {
	      d = tKid.fLink[ c ][ j2 ];
	      diff = eval->fEdgeDis[a][b] + eval->fEdgeDis[c][d] -
                     eval->fEdgeDis[a][c] - eval->fEdgeDis[b][d];
	      if( diff > max_diff ) 
	      { 
	        aa = a; bb = b; a1 = c; b1 = d;
	        max_diff = diff;
	      }
	      diff = eval->fEdgeDis[a][b] + eval->fEdgeDis[d][c] - 
		     eval->fEdgeDis[a][d] - eval->fEdgeDis[b][c];
	      if( diff > max_diff ) 
	      {
	        aa = a; bb = b; a1 = d; b1 = c;
	        max_diff = diff;
	      } 
	    }
	  }
	}
      }
    }

    if( a1 == -1 && nearMax == 10 ){  /* This value must also be changed if nearMax is chenged above */
      nearMax = 50;
      goto RESTART;
    }    
    else if( a1 == -1 && nearMax == 50  )
    {       
      int r = rand() % ( fNumOfElementInCU - 1 );
      a = fListOfCenterUnit[ r ];
      b = fListOfCenterUnit[ r+1 ];
      for( j = 0; j < fN; ++j )
      {
	if( fCenterUnit[ j ] == 0 )
        {
	  aa = a; bb = b;
	  a1 = j;
	  b1 = tKid.fLink[ j ][ 0 ];
	  break;
	}
      }
      max_diff = eval->fEdgeDis[aa][bb] + eval->fEdgeDis[a1][b1] -
         	 eval->fEdgeDis[a][a1] - eval->fEdgeDis[b][b1];
    }  

    if( tKid.fLink[aa][0] == bb ) tKid.fLink[aa][0]=a1;
    else tKid.fLink[aa][1] = a1;
    if( tKid.fLink[bb][0] == aa ) tKid.fLink[bb][0] = b1;
    else tKid.fLink[bb][1] = b1;   
    if( tKid.fLink[a1][0] == b1 ) tKid.fLink[a1][0] = aa;
    else tKid.fLink[a1][1] = aa;
    if( tKid.fLink[b1][0] == a1 ) tKid.fLink[b1][0] = bb;
    else tKid.fLink[b1][1] = bb; 


    fModiEdge[ fNumOfModiEdge ][ 0 ] = aa;
    fModiEdge[ fNumOfModiEdge ][ 1 ] = bb;
    fModiEdge[ fNumOfModiEdge ][ 2 ] = a1;
    fModiEdge[ fNumOfModiEdge ][ 3 ] = b1;
    ++fNumOfModiEdge;


    fGainModi += max_diff;
    

    int posi_a1 = fInv[ a1 ];  
    select_un = -1;
    for( int s = 0; s < fNumOfSeg; ++s ){
      if( fSegment[ s ][ 0 ] <= posi_a1 && posi_a1 <=  fSegment[ s ][ 1 ] ){
	select_un = fSegUnit[ s ];       
	break;
      }
    } 
    assert( select_un != -1 );

    for( int s = 0; s < fNumOfSeg; ++s ){
      if( fSegUnit[ s ] == select_un )
	fSegUnit[ s ] = center_un;
    }
    fNumOfElementInUnit[ center_un ] += fNumOfElementInUnit[ select_un ];
    
    for( int s = 0; s < fNumOfSeg; ++s ){
      if( fSegUnit[ s ] == fNumOfUnit - 1 )
	fSegUnit[ s ] = select_un;
    }
    fNumOfElementInUnit[ select_un ] = fNumOfElementInUnit[ fNumOfUnit - 1 ];
    --fNumOfUnit;

    for( int s = 0; s < fNumOfElementInCU; ++s ){
      c = fListOfCenterUnit[ s ];
      fCenterUnit[ c ] = 0;
    }
  }
}  


void TCross::MakeUnit()                    
{
  int flag = 1; 
  for( int s = 0; s < fNumOfSPL; ++s ){
    if( fSegPosiList[ s ] == 0 ){
      flag = 0;
      break;
    }
  }
  if( flag == 1 ) 
  {
    fSegPosiList[ fNumOfSPL++ ] = 0;

    LinkBPosi[ fN-1 ][ 1 ]  = LinkBPosi[ fN-1 ][ 0 ];
    LinkBPosi[ 0 ][ 1 ] = LinkBPosi[ 0 ][ 0 ];
    LinkBPosi[ fN-1 ][ 0 ] = 0; 
    LinkBPosi[ 0 ][ 0 ] = fN-1;

  }

  tSort->Sort( fSegPosiList, fNumOfSPL );     


  fNumOfSeg = fNumOfSPL;
  for( int s = 0; s < fNumOfSeg-1; ++s ){
    fSegment[ s ][ 0 ] = fSegPosiList[ s ];
    fSegment[ s ][ 1 ] = fSegPosiList[ s+1 ]-1;
  }

  fSegment[ fNumOfSeg-1 ][ 0 ] = fSegPosiList[ fNumOfSeg-1 ];
  fSegment[ fNumOfSeg-1 ][ 1 ] = fN - 1;


  for( int s = 0; s < fNumOfSeg; ++s ){
    LinkAPosi[ fSegment[ s ][ 0 ] ] = fSegment[ s ][ 1 ];
    LinkAPosi[ fSegment[ s ][ 1 ] ] = fSegment[ s ][ 0 ];
    fPosiSeg[ fSegment[ s ][ 0 ] ] = s;
    fPosiSeg[ fSegment[ s ][ 1 ] ] = s;
  }

  for( int s = 0; s < fNumOfSeg; ++s )
    fSegUnit[ s ] = -1;
  fNumOfUnit = 0; 

  int p_st, p1, p2, p_next, p_pre; 
  int segNum;

  while(1)
  {
    flag = 0;
    for( int s = 0; s < fNumOfSeg; ++s ){
      if( fSegUnit[ s ] == -1 ){
	p_st = fSegment[ s ][ 0 ]; 
	p_pre = -1;
	p1 = p_st;
	flag = 1;
	break;
      }
    }
    if( flag == 0 )
      break;
    
    while(1)
    {
      segNum = fPosiSeg[ p1 ];
      fSegUnit[ segNum ] = fNumOfUnit;

      p2 = LinkAPosi[ p1 ];
      p_next = LinkBPosi[ p2 ][ 0 ];
      if( p1 == p2 ){
	if( p_next == p_pre )
	  p_next = LinkBPosi[ p2 ][ 1 ];
      } 
      
      if( p_next == p_st ){
	++fNumOfUnit;
	break;
      }

      p_pre = p2;
      p1 = p_next;
    }
  }
  

  for( int s = 0; s < fNumOfUnit; ++s )
    fNumOfElementInUnit[ s ] = 0; 
  
  int unitNum = -1;
  int tmpNumOfSeg = -1;
  for( int s = 0; s < fNumOfSeg; ++s ){
    if( fSegUnit[ s ] != unitNum ){
      ++tmpNumOfSeg;
      fSegment[ tmpNumOfSeg ][ 0 ] = fSegment[ s ][ 0 ];
      fSegment[ tmpNumOfSeg ][ 1 ] = fSegment[ s ][ 1 ];
      unitNum = fSegUnit[ s ];
      fSegUnit[ tmpNumOfSeg ] = unitNum;
      fNumOfElementInUnit[ unitNum ] += 
	fSegment[ s ][ 1 ] - fSegment[ s ][ 0 ] + 1;
    }
    else
    {
      fSegment[ tmpNumOfSeg ][ 1 ] = fSegment[ s ][ 1 ];
      fNumOfElementInUnit[ unitNum ] += 
	fSegment[ s ][ 1 ] - fSegment[ s ][ 0 ] + 1;
    }
  }
  fNumOfSeg = tmpNumOfSeg + 1;  
}


void TCross::BackToPa1( TIndi& tKid )
{
  int aa, bb, a1, b1; 
  int jnum;

  for( int s = fNumOfModiEdge -1; s >= 0; --s ){ 
    aa = fModiEdge[ s ][ 0 ];
    a1 = fModiEdge[ s ][ 1 ];   // ここを変更に注意 
    bb = fModiEdge[ s ][ 2 ];   // ここを変更に注意 
    b1 = fModiEdge[ s ][ 3 ];

    if( tKid.fLink[aa][0] == bb ) tKid.fLink[aa][0] = a1;
    else tKid.fLink[aa][1] = a1;
    if( tKid.fLink[b1][0] == a1 ) tKid.fLink[b1][0] = bb;
    else tKid.fLink[b1][1] = bb; 
    if( tKid.fLink[bb][0] == aa ) tKid.fLink[bb][0] = b1;
    else tKid.fLink[bb][1] = b1;   
    if( tKid.fLink[a1][0] == b1 ) tKid.fLink[a1][0] = aa;
    else tKid.fLink[a1][1] = aa;
  }
  
  for( int s = 0; s < fNumOfAppliedCycle; ++s ){
    jnum = fAppliedCylce[ s ];
    this->ChangeSol( tKid, jnum, 2 );
  }
}

void TCross::GoToBest( TIndi& tKid )
{
  int aa, bb, a1, b1; 
  int jnum;

  for( int s = 0; s < fNumOfBestAppliedCycle; ++s ){
    jnum = fBestAppliedCylce[ s ];
    this->ChangeSol( tKid, jnum, 1 );
  }

  for( int s = 0; s < fNumOfBestModiEdge; ++s )
  { 
    aa = fBestModiEdge[ s ][ 0 ];
    bb = fBestModiEdge[ s ][ 1 ];   
    a1 = fBestModiEdge[ s ][ 2 ];   
    b1 = fBestModiEdge[ s ][ 3 ];

    if( tKid.fLink[aa][0] == bb ) tKid.fLink[aa][0]=a1;
    else tKid.fLink[aa][1] = a1;
    if( tKid.fLink[bb][0] == aa ) tKid.fLink[bb][0] = b1;
    else tKid.fLink[bb][1] = b1;   
    if( tKid.fLink[a1][0] == b1 ) tKid.fLink[a1][0] = aa;
    else tKid.fLink[a1][1] = aa;
    if( tKid.fLink[b1][0] == a1 ) tKid.fLink[b1][0] = bb;
    else tKid.fLink[b1][1] = bb; 
  }
}


void TCross::IncrementEdgeFreq( int **fEdgeFreq )
{
  int j, jnum, cem;
  int r1, r2, b1, b2;
  int aa, bb, a1;
  
  for( int s = 0; s < fNumOfBestAppliedCycle; ++s ){
    jnum = fBestAppliedCylce[ s ];
    
    cem = fABcycle[ jnum ][ 0 ];  
    fC[ 0 ] = fABcycle[ jnum ][ 0 ];

    for( j = 1; j <= cem+3; ++j ) 
      fC[ j ] = fABcycle[ jnum ][ j ];

    for( j = 0; j <cem/2; ++j )
    {                           
      r1 = fC[2+2*j]; r2 = fC[3+2*j]; 
      b1 = fC[1+2*j]; b2 = fC[4+2*j]; 

      // r1 - b1 add    
      // r1 - r2 remove
      // r2 - r1 remove
      // r2 - b2 add

      ++fEdgeFreq[ r1 ][ b1 ];
      --fEdgeFreq[ r1 ][ r2 ];
      --fEdgeFreq[ r2 ][ r1 ];
      ++fEdgeFreq[ r2 ][ b2 ];

    }
  }

  for( int s = 0; s < fNumOfBestModiEdge; ++s )
  { 
    aa = fBestModiEdge[ s ][ 0 ];
    bb = fBestModiEdge[ s ][ 1 ];   
    a1 = fBestModiEdge[ s ][ 2 ];   
    b1 = fBestModiEdge[ s ][ 3 ];

    --fEdgeFreq[ aa ][ bb ];
    --fEdgeFreq[ a1 ][ b1 ];
    ++fEdgeFreq[ aa ][ a1 ];
    ++fEdgeFreq[ bb ][ b1 ];
    --fEdgeFreq[ bb ][ aa ];
    --fEdgeFreq[ b1 ][ a1 ];
    ++fEdgeFreq[ a1 ][ aa ];
    ++fEdgeFreq[ b1 ][ bb ];
  }
}


int TCross::Cal_ADP_Loss( int **fEdgeFreq )
{
  int j, jnum, cem;
  int r1, r2, b1, b2;
  int aa, bb, a1;
  double DLoss; 
  double h1, h2;

  
  DLoss = 0;
  for( int s = 0; s < fNumOfAppliedCycle; ++s ){
    jnum = fAppliedCylce[ s ];
    
    cem = fABcycle[ jnum ][ 0 ];  
    fC[ 0 ] = fABcycle[ jnum ][ 0 ];

    for( j = 1; j <= cem+3; ++j ) 
      fC[ j ] = fABcycle[ jnum ][ j ];

    for( j = 0; j <cem/2; ++j )
    {                           
      r1 = fC[2+2*j]; r2 = fC[3+2*j]; 
      b1 = fC[1+2*j]; b2 = fC[4+2*j]; 

      // r1 - b1 add 
      // r1 - r2 remove
      // r2 - r1 remove
      // r2 - b2 add

      DLoss -= (fEdgeFreq[ r1 ][ r2 ]-1);
      DLoss -= (fEdgeFreq[ r2 ][ r1 ]-1);
      DLoss += fEdgeFreq[ r2 ][ b2 ];
      DLoss += fEdgeFreq[ b2 ][ r2 ];

      // Remove
      --fEdgeFreq[ r1 ][ r2 ]; 
      --fEdgeFreq[ r2 ][ r1 ]; 

      // Add
      ++fEdgeFreq[ r2 ][ b2 ]; 
      ++fEdgeFreq[ b2 ][ r2 ]; 
    }
  }


  for( int s = 0; s < fNumOfModiEdge; ++s )
  { 
    aa = fModiEdge[ s ][ 0 ];
    bb = fModiEdge[ s ][ 1 ];   
    a1 = fModiEdge[ s ][ 2 ];   
    b1 = fModiEdge[ s ][ 3 ];

    DLoss -= (fEdgeFreq[ aa ][ bb ]-1);
    DLoss -= (fEdgeFreq[ bb ][ aa ]-1);
    DLoss -= (fEdgeFreq[ a1 ][ b1 ]-1);
    DLoss -= (fEdgeFreq[ b1 ][ a1 ]-1);

    DLoss += fEdgeFreq[ aa ][ a1 ];
    DLoss += fEdgeFreq[ a1 ][ aa ];
    DLoss += fEdgeFreq[ bb ][ b1 ];
    DLoss += fEdgeFreq[ b1 ][ bb ];

    // Remove
    --fEdgeFreq[ aa ][ bb ];
    --fEdgeFreq[ bb ][ aa ];
    --fEdgeFreq[ a1 ][ b1 ];
    --fEdgeFreq[ b1 ][ a1 ];

    // Add
    ++fEdgeFreq[ aa ][ a1 ];
    ++fEdgeFreq[ a1 ][ aa ];
    ++fEdgeFreq[ bb ][ b1 ];
    ++fEdgeFreq[ b1 ][ bb ];
  }

  
  for( int s = 0; s < fNumOfAppliedCycle; ++s ){
    jnum = fAppliedCylce[ s ];
    
    cem = fABcycle[ jnum ][ 0 ];  
    fC[ 0 ] = fABcycle[ jnum ][ 0 ];

    for( j = 1; j <= cem+3; ++j ) 
      fC[ j ] = fABcycle[ jnum ][ j ];

    for( j = 0; j <cem/2; ++j )
    {                           
      r1 = fC[2+2*j]; r2 = fC[3+2*j]; 
      b1 = fC[1+2*j]; b2 = fC[4+2*j]; 

      ++fEdgeFreq[ r1 ][ r2 ]; 
      ++fEdgeFreq[ r2 ][ r1 ]; 
      --fEdgeFreq[ r2 ][ b2 ]; 
      --fEdgeFreq[ b2 ][ r2 ]; 
    }
  }

  // Modification
  for( int s = 0; s < fNumOfModiEdge; ++s )
  { 
    aa = fModiEdge[ s ][ 0 ];
    bb = fModiEdge[ s ][ 1 ];   
    a1 = fModiEdge[ s ][ 2 ];   
    b1 = fModiEdge[ s ][ 3 ];

    // Remove
    ++fEdgeFreq[ aa ][ bb ];
    ++fEdgeFreq[ bb ][ aa ];

    ++fEdgeFreq[ a1 ][ b1 ];
    ++fEdgeFreq[ b1 ][ a1 ];

    --fEdgeFreq[ aa ][ a1 ];
    --fEdgeFreq[ a1 ][ aa ];

    --fEdgeFreq[ bb ][ b1 ];
    --fEdgeFreq[ b1 ][ bb ];

  }

  return int(DLoss / 2);
}


double TCross::Cal_ENT_Loss( int **fEdgeFreq )
{
  int j, jnum, cem;
  int r1, r2, b1, b2;
  int aa, bb, a1;
  double DLoss; 
  double h1, h2;

  
  DLoss = 0;
// AB-cycle
  for( int s = 0; s < fNumOfAppliedCycle; ++s ){
    jnum = fAppliedCylce[ s ];
    
    cem = fABcycle[ jnum ][ 0 ];  
    fC[ 0 ] = fABcycle[ jnum ][ 0 ];

    for( j = 1; j <= cem+3; ++j ) 
      fC[ j ] = fABcycle[ jnum ][ j ];

    for( j = 0; j <cem/2; ++j )
    {                           
      r1 = fC[2+2*j]; r2 = fC[3+2*j]; 
      b1 = fC[1+2*j]; b2 = fC[4+2*j]; 

      // r1 - b1 add    
      // r1 - r2 remove
      // r2 - r1 remove
      // r2 - b2 add

      // Remove
      h1 = (double)( fEdgeFreq[ r1 ][ r2 ] - 1 )/(double)fNumOfPop;
      h2 = (double)( fEdgeFreq[ r1 ][ r2 ] )/(double)fNumOfPop;
      if( fEdgeFreq[ r1 ][ r2 ] - 1 != 0 )
	DLoss -= h1 * log( h1 );
      DLoss += h2 * log( h2 );
      --fEdgeFreq[ r1 ][ r2 ]; 
      --fEdgeFreq[ r2 ][ r1 ]; 

      // Add
      h1 = (double)( fEdgeFreq[ r2 ][ b2 ] + 1 )/(double)fNumOfPop;
      h2 = (double)( fEdgeFreq[ r2 ][ b2 ])/(double)fNumOfPop;
      DLoss -= h1 * log( h1 );
      if( fEdgeFreq[ r2 ][ b2 ] != 0 )
	DLoss += h2 * log( h2 );
      ++fEdgeFreq[ r2 ][ b2 ]; 
      ++fEdgeFreq[ b2 ][ r2 ]; 
    }
  }

  // Modification
  for( int s = 0; s < fNumOfModiEdge; ++s )
  { 
    aa = fModiEdge[ s ][ 0 ];
    bb = fModiEdge[ s ][ 1 ];   
    a1 = fModiEdge[ s ][ 2 ];   
    b1 = fModiEdge[ s ][ 3 ];

    // Remove
    h1 = (double)( fEdgeFreq[ aa ][ bb ] - 1 )/(double)fNumOfPop;
    h2 = (double)( fEdgeFreq[ aa ][ bb ] )/(double)fNumOfPop;
    if( fEdgeFreq[ aa ][ bb ] - 1 != 0 )
      DLoss -= h1 * log( h1 );
    DLoss += h2 * log( h2 );
    --fEdgeFreq[ aa ][ bb ];
    --fEdgeFreq[ bb ][ aa ];

    h1 = (double)( fEdgeFreq[ a1 ][ b1 ] - 1 )/(double)fNumOfPop;
    h2 = (double)( fEdgeFreq[ a1 ][ b1 ] )/(double)fNumOfPop;
    if( fEdgeFreq[ a1 ][ b1 ] - 1 != 0 )
      DLoss -= h1 * log( h1 );
    DLoss += h2 * log( h2 );
    --fEdgeFreq[ a1 ][ b1 ];
    --fEdgeFreq[ b1 ][ a1 ];

    // Add
    h1 = (double)( fEdgeFreq[ aa ][ a1 ] + 1 )/(double)fNumOfPop;
    h2 = (double)( fEdgeFreq[ aa ][ a1 ])/(double)fNumOfPop;
    DLoss -= h1 * log( h1 );
    if( fEdgeFreq[ aa ][ a1 ] != 0 )
      DLoss += h2 * log( h2 );
    ++fEdgeFreq[ aa ][ a1 ];
    ++fEdgeFreq[ a1 ][ aa ];

    h1 = (double)( fEdgeFreq[ bb ][ b1 ] + 1 )/(double)fNumOfPop;
    h2 = (double)( fEdgeFreq[ bb ][ b1 ])/(double)fNumOfPop;
    DLoss -= h1 * log( h1 );
    if( fEdgeFreq[ bb ][ b1 ] != 0 )
      DLoss += h2 * log( h2 );
    ++fEdgeFreq[ bb ][ b1 ];
    ++fEdgeFreq[ b1 ][ bb ];
  }
  DLoss = -DLoss;  

  // restore EdgeFreq
  for( int s = 0; s < fNumOfAppliedCycle; ++s ){
    jnum = fAppliedCylce[ s ];
    
    cem = fABcycle[ jnum ][ 0 ];  
    fC[ 0 ] = fABcycle[ jnum ][ 0 ];

    for( j = 1; j <= cem+3; ++j ) 
      fC[ j ] = fABcycle[ jnum ][ j ];

    for( j = 0; j <cem/2; ++j )
    {                           
      r1 = fC[2+2*j]; r2 = fC[3+2*j]; 
      b1 = fC[1+2*j]; b2 = fC[4+2*j]; 

      ++fEdgeFreq[ r1 ][ r2 ]; 
      ++fEdgeFreq[ r2 ][ r1 ]; 
      --fEdgeFreq[ r2 ][ b2 ]; 
      --fEdgeFreq[ b2 ][ r2 ]; 
    }
  }

  for( int s = 0; s < fNumOfModiEdge; ++s )
  { 
    aa = fModiEdge[ s ][ 0 ];
    bb = fModiEdge[ s ][ 1 ];   
    a1 = fModiEdge[ s ][ 2 ];   
    b1 = fModiEdge[ s ][ 3 ];

    ++fEdgeFreq[ aa ][ bb ];
    ++fEdgeFreq[ bb ][ aa ];

    ++fEdgeFreq[ a1 ][ b1 ];
    ++fEdgeFreq[ b1 ][ a1 ];

    --fEdgeFreq[ aa ][ a1 ];
    --fEdgeFreq[ a1 ][ aa ];

    --fEdgeFreq[ bb ][ b1 ];
    --fEdgeFreq[ b1 ][ bb ];

  }

  return DLoss;
}


void TCross::SetWeight( const TIndi& tPa1, const TIndi& tPa2 ) 
{
  int cem;
  int r1, r2, v1, v2, v_p;
  int AB_num;

  for( int i = 0; i < fN; ++i ){
    fInEffectNode[ i ][ 0 ] = -1;  
    fInEffectNode[ i ][ 1 ] = -1;
  }

  // Step 1:
  for( int s = 0; s < fNumOfABcycle; ++s ){
    cem = fABcycle[ s ][ 0 ];  
    for( int j = 0; j < cem/2; ++j ){
      r1 = fABcycle[ s ][ 2*j+2 ];  // red edge
      r2 = fABcycle[ s ][ 2*j+3 ]; 

      if( fInEffectNode[ r1 ][ 0 ] == -1 ) fInEffectNode[ r1 ][ 0 ] = s;
      else if ( fInEffectNode[ r1 ][ 1 ] == -1 ) fInEffectNode[ r1 ][ 1 ] = s;
      else assert( 1 == 2 );

      if( fInEffectNode[ r2 ][ 0 ] == -1 ) fInEffectNode[ r2 ][ 0 ] = s;
      else if ( fInEffectNode[ r2 ][ 1 ] == -1 ) fInEffectNode[ r2 ][ 1 ] = s;
      else assert( 1 == 2 );
    }
  }
  
  // Step 2:
  for( int i = 0; i < fN; ++i ){
    if( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] == -1 ){ 
      AB_num = fInEffectNode[ i ][ 0 ];
      v1 = i;

      if( tPa1.fLink[ v1 ][ 0 ] != tPa2.fLink[ v1 ][ 0 ] && tPa1.fLink[ v1 ][ 0 ] != tPa2.fLink[ v1 ][ 1 ] )
	v_p = tPa1.fLink[ v1 ][ 0 ];
      else if( tPa1.fLink[ v1 ][ 1 ] != tPa2.fLink[ v1 ][ 0 ] && tPa1.fLink[ v1 ][ 1 ] != tPa2.fLink[ v1 ][ 1 ] )
	v_p = tPa1.fLink[ v1 ][ 1 ];
      else
	assert( 1 == 2 );

      while( 1 ){
	assert( fInEffectNode[ v1 ][ 0 ] != -1 );
	assert( fInEffectNode[ v1 ][ 1 ] == -1 );
	fInEffectNode[ v1 ][ 1 ] = AB_num;

	if( tPa1.fLink[ v1 ][ 0 ] != v_p )
	  v2 = tPa1.fLink[ v1 ][ 0 ];
	else if( tPa1.fLink[ v1 ][ 1 ] != v_p )
	  v2 = tPa1.fLink[ v1 ][ 1 ];
	else 
	  assert( 1 == 2 );

	if( fInEffectNode[ v2 ][ 0 ] == -1 )
	  fInEffectNode[ v2 ][ 0 ] = AB_num;
	else if( fInEffectNode[ v2 ][ 1 ] == -1 )
	  fInEffectNode[ v2 ][ 1 ] = AB_num;
	else 
	  assert( 1 == 2 );
	
	if( fInEffectNode[ v2 ][ 1 ] != -1 )
	  break;

	v_p = v1;
	v1 = v2;
      }
    }
  }

  // Step 3:
  assert( fNumOfABcycle < fMaxNumOfABcycle );
  for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
    fWeight_C[ s1 ] = 0;
    for( int s2 = 0; s2 < fNumOfABcycle; ++s2 ){
      fWeight_RR[ s1 ][ s2 ] = 0;
    }
  }
  
  for( int i = 0; i < fN; ++i ){
    assert( (fInEffectNode[ i ][ 0 ] == -1 && fInEffectNode[ i ][ 1 ] == -1) ||
	    (fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] != -1) );

    if( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] != -1 ){
      ++fWeight_RR[ fInEffectNode[ i ][ 0 ] ][ fInEffectNode[ i ][ 1 ] ];
      ++fWeight_RR[ fInEffectNode[ i ][ 1 ] ][ fInEffectNode[ i ][ 0 ] ];
    }
    if( fInEffectNode[ i ][ 0 ] != fInEffectNode[ i ][ 1 ] ){
      ++fWeight_C[ fInEffectNode[ i ][ 0 ] ];
      ++fWeight_C[ fInEffectNode[ i ][ 1 ] ];
    }
    
  }
  for( int s1 = 0; s1 < fNumOfABcycle; ++s1 )
    fWeight_RR[ s1 ][ s1 ] = 0;
  

  for( int i = 0; i < fN; ++i ){
    assert( ( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] != -1 ) ||
	    ( fInEffectNode[ i ][ 0 ] == -1 && fInEffectNode[ i ][ 1 ] == -1 ) );
  }

}


int TCross::Cal_C_Naive() 
{
  int count_C;
  int tt;

  count_C = 0;

  for( int i = 0; i < fN; ++i ){
    if( fInEffectNode[ i ][ 0 ] != -1 && fInEffectNode[ i ][ 1 ] != -1 ){
      tt = 0;
      if( fUsedAB[ fInEffectNode[ i ][ 0 ] ] == 1 )
	++tt;
      if( fUsedAB[ fInEffectNode[ i ][ 1 ] ] == 1 )
	++tt;
      if( tt == 1 )
	++count_C;
    }
  }
  return count_C;
}

void TCross::Search_Eset( int centerAB ) 
{
  int nIter, stagImp;
  int delta_weight, min_delta_weight_nt;
  int flag_AddDelete, flag_AddDelete_nt;
  int selected_AB, selected_AB_nt;
  int t_max;
  int jnum;

  fNum_C = 0;  // Number of C nodes in E-set
  fNum_E = 0;  // Number of Edges in E-set 

  fNumOfUsedAB = 0;
  for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
    fUsedAB[ s1 ] = 0;
    fWeight_SR[ s1 ] = 0;
    fMoved_AB[ s1 ] = 0;
  }

  for( int s = 0; s < fNumOfABcycleInEset; ++s )   
  {
    jnum = fABcycleInEset[ s ];
    this->Add_AB( jnum );
  }
  fBest_Num_C = fNum_C;
  fBest_Num_E = fNum_E;
  
  stagImp = 0;
  nIter = 0;
  while( 1 )
  { 
    ++nIter;

    min_delta_weight_nt = 99999999;  
    flag_AddDelete = 0;
    flag_AddDelete_nt = 0;
    for( int s1 = 0; s1 < fNumOfABcycle; ++s1 )
    {
      if( fUsedAB[ s1 ] == 0 && fWeight_SR[ s1 ] > 0 )
      {
	delta_weight = fWeight_C[ s1 ] - 2 * fWeight_SR[ s1 ];   
	if( fNum_C + delta_weight < fBest_Num_C ){
	  selected_AB = s1;
	  flag_AddDelete = 1;
	  fBest_Num_C = fNum_C + delta_weight;
	}
	if( delta_weight < min_delta_weight_nt && nIter > fMoved_AB[ s1 ] ){
	  selected_AB_nt = s1;
	  flag_AddDelete_nt = 1;
	  min_delta_weight_nt = delta_weight;
	}
      }
      else if( fUsedAB[ s1 ] == 1 && s1 != centerAB )
      {
	delta_weight = - fWeight_C[ s1 ] + 2 * fWeight_SR[ s1 ];   
	if( fNum_C + delta_weight < fBest_Num_C ){
	  selected_AB = s1;
	  flag_AddDelete = -1;
	  fBest_Num_C = fNum_C + delta_weight;
	}
	if( delta_weight < min_delta_weight_nt && nIter > fMoved_AB[ s1 ] ){
	  selected_AB_nt = s1;
	  flag_AddDelete_nt = -1;
	  min_delta_weight_nt = delta_weight;
	}
      }
    }
      
    if( flag_AddDelete != 0 ){
      if( flag_AddDelete == 1 ){
	this->Add_AB( selected_AB );
      }
      else if( flag_AddDelete == -1 )
	this->Delete_AB( selected_AB );
      
      fMoved_AB[ selected_AB ] = nIter + tRand->Integer( 1, fTmax ); 
      assert( fBest_Num_C == fNum_C );
      fBest_Num_E = fNum_E;

      fNumOfABcycleInEset = 0;
      for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
	if( fUsedAB[ s1 ] == 1 )
	  fABcycleInEset[ fNumOfABcycleInEset++ ] = s1;
      }
      assert( fNumOfABcycleInEset == fNumOfUsedAB );      
      stagImp = 0;
    }
    else if( flag_AddDelete_nt != 0 ) {
      if( flag_AddDelete_nt == 1 ){
	this->Add_AB( selected_AB_nt );
      }
      else if( flag_AddDelete_nt == -1 )
	this->Delete_AB( selected_AB_nt );

      fMoved_AB[ selected_AB_nt ] = nIter + tRand->Integer( 1, fTmax ); 
    } 

    if( flag_AddDelete == 0 )
      ++stagImp;
    if( stagImp == fMaxStag )
      break;
  }
}


void TCross::Add_AB( int AB_num )  
{
  fNum_C += fWeight_C[ AB_num ] - 2 * fWeight_SR[ AB_num ];   
  fNum_E += fABcycle[ AB_num ][ 0 ] / 2;  

  assert( fUsedAB[ AB_num ] == 0 );
  fUsedAB[ AB_num ] = 1;
  ++fNumOfUsedAB;

  for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
    fWeight_SR[ s1 ] += fWeight_RR[ s1 ][ AB_num ];
  }
}


void TCross::Delete_AB( int AB_num )  
{
  fNum_C -= fWeight_C[ AB_num ] - 2 * fWeight_SR[ AB_num ];   
  fNum_E -= fABcycle[ AB_num ][ 0 ] / 2;  

  assert( fUsedAB[ AB_num ] == 1 );
  fUsedAB[ AB_num ] = 0;
  --fNumOfUsedAB;

  for( int s1 = 0; s1 < fNumOfABcycle; ++s1 ){
    fWeight_SR[ s1 ] -= fWeight_RR[ s1 ][ AB_num ];
  }
}


void TCross::CheckValid( TIndi& indi )
{
  int curr, pre, next, st;
  int count;

  st = 0;
  curr = -1;
  next = st;

  count = 0;
  while(1){ 
    pre = curr;
    curr = next;
    ++count;
    if( indi.fLink[ curr ][ 0 ] != pre )
      next = indi.fLink[ curr ][ 0 ];
    else 
      next = indi.fLink[ curr ][ 1 ]; 
    
    if( next == st ) break;

    if( count > fN ){
      printf( "Invalid = %d\n", count );
      break;
    }
  }       
  if( count != fN )
      printf( "Invalid = %d\n", count );
}



     
void MakeRandSol( TEvaluator* eval , TIndi& indi );
void Make2optSol( TEvaluator* eval , TIndi& indi );

TEnvironment::TEnvironment()
{
  fEvaluator = new TEvaluator();
}


TEnvironment::~TEnvironment()
{
  delete [] fIndexForMating;
  delete [] tCurPop;
  delete fEvaluator;
  delete tCross;

  int N = fEvaluator->Ncity;
  for( int i = 0; i < N; ++i ) 
    delete [] fEdgeFreq[ i ];
  delete [] fEdgeFreq;
}


void TEnvironment::Define()
{
  fEvaluator->SetInstance( fFileNameTSP );
  int N = fEvaluator->Ncity;

  fIndexForMating = new int [ fNumOfPop + 1 ];  

  tCurPop = new TIndi [ fNumOfPop ];
  for ( int i = 0; i < fNumOfPop; ++i )
    tCurPop[i].Define( N );

  tBest.Define( N );

  tCross = new TCross( N );
  tCross->eval = fEvaluator;                 
  tCross->fNumOfPop = fNumOfPop;             

  tKopt = new TKopt( N );
  tKopt->eval = fEvaluator;
  tKopt->SetInvNearList();

  fEdgeFreq = new int* [ N ]; 
  for( int i = 0; i < N; ++i ) 
    fEdgeFreq[ i ] = new int [ N ]; 
}


void TEnvironment::DoIt()
{
  this->fTimeStart = clock();   

  if( fFileNameInitPop == NULL )
    this->InitPop();                       
  else
    this->ReadPop( fFileNameInitPop );     

  this->fTimeInit = clock();    

  this->Init();
  this->GetEdgeFreq();

  while( 1 )
  {
    this->SetAverageBest();
    printf( "%d: %d %lf\n", fCurNumOfGen, fBestValue, fAverageValue );

    if( this->TerminationCondition() ) break;

    this->SelectForMating();

    for( int s =0; s < fNumOfPop; ++s )
    {
      this->GenerateKids( s );     
      this->SelectForSurvival( s ); 
    }
    ++fCurNumOfGen;
  }

  this->fTimeEnd = clock();   
}
 

void TEnvironment::Init()
{
  fAccumurateNumCh = 0;
  fCurNumOfGen = 0;
  fStagBest = 0;
  fMaxStagBest = 0;
  fStage = 1;          /* Stage I */
  fFlagC[ 0 ] = 4;     /* Diversity preservation: 1:Greedy, 2:--- , 3:Distance, 4:Entropy (see Section 4) */
  fFlagC[ 1 ] = 1;     /* Eset Type: 1:Single-AB, 2:Block2 (see Section 3) */ 
} 


bool TEnvironment::TerminationCondition()
{
  if ( fAverageValue - fBestValue < 0.001 )  
    return true;

  if( fStage == 1 ) /* Stage I */      
  {
    if( fStagBest == int(1500/fNumOfKids) && fMaxStagBest == 0 ){ /* 1500/N_ch (See Section 2.2) */
      fMaxStagBest =int( fCurNumOfGen / 10 );                 /* fMaxStagBest = G/10 (See Section 2.2) */   
    } 
    else if( fMaxStagBest != 0 && fMaxStagBest <= fStagBest ){ /* Terminate Stage I (proceed to Stage II) */
      fStagBest = 0;
      fMaxStagBest = 0;
      fCurNumOfGen1 = fCurNumOfGen;
      fFlagC[ 1 ] = 2; 
      fStage = 2;      
    }
    return false;
  }

  if( fStage == 2 ){ /* Stage II */
    if( fStagBest == int(1500/fNumOfKids) && fMaxStagBest == 0 ){ /* 1500/N_ch */
      fMaxStagBest = int( (fCurNumOfGen - fCurNumOfGen1) / 10 ); /* fMaxStagBest = G/10 (See Section 2.2) */
    } 
    else if( fMaxStagBest != 0 && fMaxStagBest <= fStagBest ){ /* Terminate Stage II and GA */
      return true;
    }

    return false;
  }
  return false;
}


void TEnvironment::SetAverageBest() 
{
  int stockBest = tBest.fEvaluationValue;
  
  fAverageValue = 0.0;
  fBestIndex = 0;
  fBestValue = tCurPop[0].fEvaluationValue;
  
  for(int i = 0; i < fNumOfPop; ++i ){
    fAverageValue += tCurPop[i].fEvaluationValue;
    if( tCurPop[i].fEvaluationValue < fBestValue ){
      fBestIndex = i;
      fBestValue = tCurPop[i].fEvaluationValue;
    }
  }
  
  tBest = tCurPop[ fBestIndex ];
  fAverageValue /= (double)fNumOfPop;

  if( tBest.fEvaluationValue < stockBest ){
    fStagBest = 0;
    fBestNumOfGen = fCurNumOfGen;
    fBestAccumeratedNumCh = fAccumurateNumCh;
  }
  else ++fStagBest;
}


void TEnvironment::InitPop()
{
  for ( int i = 0; i < fNumOfPop; ++i ){ 
    tKopt->MakeRandSol( tCurPop[ i ] );    /* Make a random tour */
    tKopt->DoIt( tCurPop[ i ] );           /* Apply the local search with the 2-opt neighborhood */ 
  }
}


void TEnvironment::SelectForMating()
{
  /* fIndexForMating[] <-- a random permutation of 0, ..., fNumOfPop-1 */
  tRand->Permutation( fIndexForMating, fNumOfPop, fNumOfPop ); 
  fIndexForMating[ fNumOfPop ] = fIndexForMating[ 0 ];
}

void TEnvironment::SelectForSurvival( int s )
{
}


void TEnvironment::GenerateKids( int s )
{
  tCross->SetParents( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], fFlagC, fNumOfKids );  

  tCross->DoIt( tCurPop[fIndexForMating[s]], tCurPop[fIndexForMating[s+1]], fNumOfKids, 1, fFlagC, fEdgeFreq );

  fAccumurateNumCh += tCross->fNumOfGeneratedCh;
}


void TEnvironment::GetEdgeFreq()
{
  int N = fEvaluator->Ncity;
  int k0, k1;
  
  for( int j1 = 0; j1 < N; ++j1 )
    for( int j2 = 0; j2 < N; ++j2 ) 
      fEdgeFreq[ j1 ][ j2 ] = 0;

  
  for( int i = 0; i < fNumOfPop; ++i )
  {
    for(int j = 0; j < N; ++j )
    {
      k0 = tCurPop[ i ].fLink[ j ][ 0 ];
      k1 = tCurPop[ i ].fLink[ j ][ 1 ];
      ++fEdgeFreq[ j ][ k0 ];
      ++fEdgeFreq[ j ][ k1 ];
    }
  }
}


void TEnvironment::PrintOn( int n, char* dstFile ) 
{
  printf( "n = %d val = %d Gen = %d \n" , 
	  n, 
	  tBest.fEvaluationValue, 
	  fCurNumOfGen );
	  
  fflush(stdout);

  FILE *fp;
  char filename[ 80 ];
  sprintf( filename, "%s_Result", dstFile );
  fp = fopen( filename, "a");
  
  fprintf( fp, "%d %d %d %d %d\n" , 
	   n, 
	   tBest.fEvaluationValue, 
	   fCurNumOfGen, 
	   (int)((double)(this->fTimeInit - this->fTimeStart)/(double)CLOCKS_PER_SEC), 
	   (int)((double)(this->fTimeEnd - this->fTimeStart)/(double)CLOCKS_PER_SEC) );
  
  fclose( fp );
}


void TEnvironment::WriteBest( char* dstFile ) 
{
  FILE *fp;
  char filename[ 80 ];
  sprintf( filename, "%s_BestSol", dstFile );
  fp = fopen( filename, "a");
  
  fEvaluator->WriteTo( fp, tBest );

  fclose( fp );
}


void TEnvironment::WritePop( int n, char* dstFile ) 
{
  FILE *fp;
  char filename[ 80 ];
  sprintf( filename, "%s_POP_%d", dstFile, n );
  fp = fopen( filename, "w");

  for( int s = 0; s < fNumOfPop; ++s )
    fEvaluator->WriteTo( fp, tCurPop[ s ] );

  fclose( fp );
}


void TEnvironment::ReadPop( char* fileName )
{
  FILE* fp;

  if( ( fp = fopen( fileName, "r" ) ) == NULL ){
    printf( "Read Error1\n"); 
    fflush( stdout );
    exit( 1 );
  }

  for ( int i = 0; i < fNumOfPop; ++i ){ 
    if( fEvaluator->ReadFrom( fp, tCurPop[ i ] ) == false ){
      printf( "Read Error2\n"); 
      fflush( stdout );
      exit( 1 );
    }
  }
  fclose( fp );
}





int main( int argc, char* argv[] )
{
  int maxNumOfTrial;
  
  sscanf( argv[1], "%d", &maxNumOfTrial );
  char* dstFile = argv[2];
  
  TEnvironment* gEnv = NULL;
  gEnv = new TEnvironment();
  InitURandom();  

  int d;
  sscanf( argv[3], "%d", &d );
  gEnv->fNumOfPop = d;
  sscanf( argv[4], "%d", &d );
  gEnv->fNumOfKids = d;
  gEnv->fFileNameTSP = argv[5];
  gEnv->fFileNameInitPop = NULL;
  if( argc == 7 )
    gEnv->fFileNameInitPop = argv[6];

  gEnv->Define();
  float startTicks=(float)clock();//to only account for time spent in solving the problem without storing input values
  for( int n = 0; n < maxNumOfTrial; ++n )
  { 
    gEnv->DoIt();

    gEnv->PrintOn( n, dstFile );       
    gEnv->WriteBest( dstFile );  
    // gEnv->WritePop( n, dstFile );
  }
  cout << "time taken : " << ((float)clock()-startTicks) / CLOCKS_PER_SEC << " secs " << endl;
  return 0;
}
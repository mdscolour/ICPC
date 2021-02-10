#ifndef _MyLib_
#define _MyLib_
/**************************************************************************
// Define "RNG_NAME" for random number, in linux using drand48(), in windows using rand(),
// Define SeedByTime() respectively
**************************************************************************/

#define M_PI       3.14159265358979323846
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<iostream>
#include<vector>
#include<fstream>
#include<string>
#include<iterator>
#include<sstream>
#include<time.h>
#include<algorithm>
//#include"../../tnt/tnt.h"
//using namespace TNT;
using namespace std;

// Define "RNG_NAME" for random number, in linux using drand48(), in windows using rand(),
// the initialization of the seed is in the constructor of class Walk
#ifdef __APPLE__
#include <sys/time.h>
#define RNG_NAME rand
static inline void SeedByTime()
{
	struct timeval tpstart;
	gettimeofday(&tpstart, NULL);
	srand(tpstart.tv_usec);
	//srand48(floor(mytime()));
}
#endif

#ifdef linux
#include <sys/time.h>
#define RNG_NAME drand48
static inline void SeedByTime()
{
	struct timeval tpstart;
	gettimeofday(&tpstart, NULL);
	srand48(tpstart.tv_usec);
	//srand48(floor(mytime()));
}
#endif

#ifdef _WIN32
#include <time.h>
static inline double GetRand()
{
	//srand(time(NULL));
	return (rand()/(double)(RAND_MAX + 1));
}
#define RNG_NAME GetRand
static inline void SeedByTime()
{
	srand(unsigned int(time(NULL)));
}
#endif
////// above is basic definition

//////////////////////////////////////////////////////////////////////////
// Parameters of the SAW starts from here
//////////////////////////////////////////////////////////////////////////
#define PRECISION 1e-6

static inline int Random_integer_uniform(int a,int b){return a+int((b-a)*RNG_NAME());}
static inline int Random_integer_log(int a,int b)//not efficiency
{
	if(b==a)return 0;
	int t,i;
	for (i=0;i<1000;i++)
	{
		t = a+int((b-a)*RNG_NAME());
		if(RNG_NAME()<=log(1+1.0/(t-a+1))/log(b-a+1.0)) break;
	}
	if(i==999)printf("error(Random_integer_log): 1000 trial not success but still in domain.\n");
	return t;
}

static inline int twopowof(int power){int temp = 1; for(int i=0;i<power;i++) temp*=2; return temp;}
static inline int factorial(int n){if(n==1|| n==0)return 1;else return n*factorial(n-1);}

static inline double GaussianPDF(double x,double miu,double sigma){
    return 1.0/(sigma*sqrt(2*M_PI))*exp(-0.5*pow((x-miu)/sigma,2));}

static inline void VectorMultipliedByScalar(vector<double> &v, double k)
{int numV=int(v.size());for(int i=0;i<numV;i++)v[i]*=k;}

static inline void VectorDividedByScalar(vector<double> &v, double k)
{int numV=int(v.size());for(int i=0;i<numV;i++)v[i]/=k;}

static inline void VectorAddedByVector(vector<double> &v, vector<double> v2)
{int numV=int(v.size());for(int i=0;i<numV;i++)v[i]+=v2[i];}

static inline void printfVector(vector<double> &v)
{
    std::cout<<'[';
    for(int i=0;i<v.size();i++) {std::cout << v[i] << " ";}
    std::cout<<']'<<endl;
}

static inline void write2DVector(const char* name,vector<double> vecX,vector<double> vecY)
{
    FILE *fptr;fptr = fopen(name, "w");
    int numV=int(vecX.size());
    for(int i=0;i<numV;i++)fprintf(fptr, "%.6lf %.6lf\n",vecX[i],vecY[i]);
    fclose(fptr);
}

//======================================================================
// Returns interpolated value at x from parallel arrays ( xData, yData )
//   Assumes that xData has at least two elements, is sorted and is strictly monotonic increasing
//   boolean argument extrapolate determines behaviour beyond ends of array (if needed)
static inline double interpolate( vector<double> &xData, vector<double> &yData, double x, bool extrapolate )

{
   
    int size = xData.size();

   
    int i = 0;                                                                  
    // find left end of interval for interpolation
   
    if ( x >= xData[size - 2] )                                                 
    // special case: beyond right end
   
    {
      
        i = size - 2;
   
        
    }
   
    else
   
    {
      
        while ( x > xData[i+1] ) i++;
   
        
    }
   
    double xL = xData[i], yL = yData[i], xR = xData[i+1], yR = yData[i+1];      
    // points on either side (unless beyond ends)
  
    if ( !extrapolate )                                                         
        // if beyond ends of array and not extrapolating
   
    {
      
        if ( x < xL ) yR = yL;
      
        if ( x > xR ) yL = yR;
   
    }

   
    double dydx = ( yR - yL ) / ( xR - xL );                                    
    // gradient

   
    return yL + dydx * ( x - xL );                                              
    // linear interpolation
}

// Create a vector of evenly spaced numbers.
static inline vector<double> linspace(double min, double max, int N)
{
    
    vector<double> range;
    
    double delta = (max-min)/double(N-1);
    
    for(int i=0; i<N; i++) {
        range.push_back(min + i*delta);
    }
    
    return range;
}
static inline vector<double> arange(double start, double stop, double step = 1) 
{
    
    vector<double> values;
    
    for(double value = start; value < stop; value += step)
        
        values.push_back(value);
    
    return values;
}
static inline double fPeriodRestrict(double x,double begin,double end)
{
    double siz = end-begin;
    if(x>end)x-=siz;
    if(x<begin)x+=siz;
    return x;
}
//Converted from MATLAB script at http://billauer.co.il/peakdet.html
//numbers must in the scale under 1e99
//return arguments for the peaks
static inline vector<int> argpeakdet(vector<double> v,double delta)
{
    vector<int> maxtab;vector<int> mintab;
    int numV = int(v.size());
    
    vector<int> x;    
    for(int value = 0; value < numV; value += 1)x.push_back(value);
    
    if(delta<=0)throw range_error("delta should be positive");
    int mnpos;int mxpos;
    double mn=1e99;double mx=-1e99;
    bool lookformax=true;
    for(int i=0;i<numV;i++)
    {
        if(v[i]>mx){mx=v[i];mxpos=x[i];}
        if(v[i]<mn){mn=v[i];mnpos=x[i];}
        if(lookformax)
        {
            if(v[i]<mx-delta)
            {
                maxtab.push_back(mxpos);
                mn=v[i];
                mnpos=x[i];
                lookformax=false;
            }
        }
        else
        {
            if(v[i]>mn+delta)
            {
                mintab.push_back(mnpos);
                mx=v[i];
                mnpos=x[i];
                lookformax=true;
            }
        }
    }
    return maxtab;
}
//double potentialDoubleRMC(double r);
////// above is basic function
#endif 

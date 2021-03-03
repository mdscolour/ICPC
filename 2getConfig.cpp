//#include "/home/li/bin/coreMolecularMC.h"
//#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
#include "coreMolecularMC.h"//only one instance


double LJPotential(double r,double epsilon=1,double sigma=1) {
    return 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));}
double QuarticBond(double r,double k0,double k1,double R){
    return 0.5*k0*pow((r-R),2)+0.25*k1*pow(r-R,4);}
double FeneBond(double r,double k,double drm,double r0){
    return -0.5*k*drm*drm*log(1-pow((r-r0)/drm,2));}
//double potential(double r){return LJPotential(r,1,1);}

//global variables
double global_LJpara;
double global_scaling;
double global_boxl;
double global_boxr;
int global_npart;
double global_minvol;
double global_temperature;
int global_num_configs;
int high_global_num_configs;
int global_num_steps;
int high_global_num_steps;
int global_init_step;
double global_timediff;
double high_global_timediff;

double potential(double r){
    //return FeneBond(r,1,30,3)+QuarticBond(r,1,1,10)+FeneBond(r,-1,30,8);
    return LJPotential(r,3,global_LJpara);
}

void findStepSize()
{
    double scaling,successrate;
    int each_step,successcount;
    ///// find step size
    scaling=0.00;
    //if(argc>=2) scaling=atof(argv[1]);
    while(true)
    {
        scaling+=5;
        // function of potential, temperature, step size
        coreMolecularMC MMC=coreMolecularMC(&potential,global_temperature,scaling,global_minvol);
        MMC.setBox(global_boxl,global_boxr);    
        MMC.addAllPartEvenly(global_npart);
        MMC.readyEverything();
        
        each_step = 1000;
        successcount = MMC.run(each_step);
        //cout<<count;
        
        successrate = float(successcount)/float(each_step*global_npart);
        printf("%d trials, scaling: %lf, successrate : %lf \n",each_step,scaling,successrate);
        if(successrate<0.5)break;
    }    
    printf("final result: scaling: %lf, successrate : %lf \n",scaling,successrate);
    global_scaling = scaling;
}

void computeMCError(double crit1,double crit2)
{
    const char* saveGrName="standardLJTest.gr";

    //vector<int> list1 = {1,10,100,200,500,1000,5000};
    //vector<int> list2 = {100,500,800};
    vector<int> list1 = {10};
    vector<int> list2 = {500,1000,10000};
    vector<vector<double>> recAll={};
    int repeat=1;
    
    int lgr = int(global_LJpara*8);
    double dx=0.01;
    int ngr = int(double(lgr)/dx);
 
    clock_t begin_time = clock();
    double timediff = float(clock() - begin_time)/CLOCKS_PER_SEC;
    vector<double> recTimeDiff ={};
    //// save all the gr
    for(int num_configs: list2)
    for(int num_steps: list1){
    begin_time = clock();
    for(int nrep=0;nrep<repeat;nrep++)
    {
        coreMolecularMC MMC=coreMolecularMC(&potential,global_temperature,
                                            global_scaling,global_minvol);
        MMC.setBox(global_boxl,global_boxr);    
        MMC.addAllPartEvenly(global_npart);
        MMC.readyEverything();
        
        //int init_step = 5000;
        MMC.run(global_init_step);
       
        vector<double> gr(ngr,0.0);//=MMC.getGr(ngr,dx);
    
        for(int i=0;i<num_configs;i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr,MMC.getGr(ngr,dx));
        }
        VectorDividedByScalar(gr,num_configs);
        recAll.push_back(gr);
        printf("%d repeat: num_configs: %d, num_steps : %d \n",nrep,num_configs,num_steps);
    }
    timediff = float(clock() - begin_time)/CLOCKS_PER_SEC;
    recTimeDiff.push_back(timediff/repeat);
    }

    //// compute and save standardLJTest.gr
    vector<double> stdgr = {};
    int numgr = int(recAll.size());
    for(int j=0;j<ngr;j++)
    {
    double tmp = 0;
    for(int i=0;i<numgr;i++){ tmp += recAll[i][j];}
    stdgr.push_back(tmp/=numgr);
    }

    write2DVector(saveGrName,arange(0,lgr,dx),stdgr);
    printf("standardLJTest.gr generated.\n");

    //// compute std
    vector<double> Error = {};
    int iumgr = 0;
    while(iumgr<numgr)
    {
        double tmp = 0;
        for(int i=iumgr;i<iumgr+repeat;i++)
            for(int j=0;j<ngr;j++)
                tmp += pow(recAll[i][j]-stdgr[j],2);
        Error.push_back(sqrt(tmp/(repeat*ngr)));
        iumgr += repeat;
    }

    //// print and write to global
    int tmpcount;
    bool tmpFlag;
    double errtmp;
    //printfVector(Error);
    tmpcount = 0;
    for(int num_configs: list2)
    for(int num_steps: list1)
    {
        errtmp = Error[tmpcount];
        printf("%d x %d, std=%lf, time=%lf \n",num_configs,num_steps,errtmp,recTimeDiff[tmpcount]);
        tmpcount++;
    }

    tmpFlag=true;
    tmpcount = 0;
    for(int num_configs: list2)
    for(int num_steps: list1)
    {
        if(tmpFlag==false)break;
        errtmp = Error[tmpcount];
        if(errtmp<crit1)
        {
            global_num_configs = num_configs;
            global_num_steps = num_steps;
            global_timediff = recTimeDiff[tmpcount];
            tmpFlag=false;///can not go to outside because of double loop
            break;
        }
        tmpcount++;//cout<<tmpcount<<" "<<int(list2.size()*list1.size())<<endl;
        //if(tmpcount==int(list2.size()*list1.size())) throw domain_error("no satisfying step size (fast).");
        if(tmpcount==int(list2.size()*list1.size())){global_num_configs=500;global_num_steps=10;}
    }
    tmpFlag=true;
    tmpcount = 0;
    for(int num_configs: list2)
    for(int num_steps: list1)
    {
        if(tmpFlag==false)break;
        errtmp = Error[tmpcount];
        if(errtmp<crit2)
        {
            high_global_num_configs = num_configs;
            high_global_num_steps = num_steps;
            high_global_timediff = recTimeDiff[tmpcount];
            tmpFlag=false;
            break;
        }      
        tmpcount++;
        //if(tmpcount==int(list2.size()*list1.size())) throw domain_error("no satisfying step size (slow).");
        if(tmpcount==int(list2.size()*list1.size())){high_global_num_configs=10000;high_global_num_steps=10;}
    }
    //printf("final result: num_configs: %d, num_steps : %d \n",global_num_configs,global_num_steps);
    //printf("final result: num_configs: %d, num_steps : %d \n",high_global_num_configs,high_global_num_steps);
    printf("computation of error done!\n");
}


void computeMCSteps()
{
    coreDirectRun::loadTarget("standardGr.gr");
    vector<int> list1 = {1,10,100,200,500,1000,5000};
    vector<int> list2 = {100,500,800};
    //vector<int> list1 = {1,10,20,50,100};
    //vector<int> list2 = {1,2,3};
    int nlist1 = int(list1.size());
    vector<int>  arrNewpOrd = {};
    
    for(int num_configs: list2)
    for(int num_steps: list1)
    {
        coreMolecularMC MMC=coreMolecularMC(&potential,global_temperature,
                                            global_scaling,global_minvol);
        MMC.setBox(global_boxl,global_boxr);    
        MMC.addAllPartEvenly(global_npart);
        MMC.readyEverything();
        
        //int num_configs=500;
        int init_step = 5000;
        MMC.run(init_step);
        
        int ngr = int(coreDirectRun::targetX.size());
        double dx = coreDirectRun::dx;
        int lgr = int(double(ngr)*dx);
        vector<double> gr(ngr,0.0);
        
        for(int i=0;i<num_configs;i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr,MMC.getGr(ngr,dx));
        }
        VectorDividedByScalar(gr,num_configs);
        double newp = coreDirectRun::getDiff(gr);
        arrNewpOrd.push_back(floor(log10(newp)));
        printf("%d x %d, MSD=%lf \n",num_configs,num_steps,newp);
    }
    /////find the result
    //for(int n : arrNewpOrd) {std::cout << n << " ";}std::cout<<endl;
    int minNewpOrd = 99;int argNewp;
    for(int i=0;i<arrNewpOrd.size();i++)
        if(arrNewpOrd[i]<minNewpOrd)
        {
            argNewp=i;//first minimum order
            minNewpOrd=arrNewpOrd[i];
        }
    //printf("%d \n",argNewp);
    global_num_configs = list2[int(argNewp/nlist1)];
    global_num_steps = list1[int(argNewp%nlist1)];
    printf("final result: num_configs: %d, num_steps : %d \n",global_num_configs,global_num_steps);
}
void computeInitSteps()
{
    int num_configs=50;
    int num_steps=1;

    coreDirectRun::loadTarget("standardLJTest.gr");
    vector<int> list1 = {1,5,10,20,50,100,200,500};
    vector<int>  arrNewpOrd = {};
    for(int init_step: list1)
    {
        coreMolecularMC MMC=coreMolecularMC(&potential,global_temperature,
                                            global_scaling,global_minvol);
        MMC.setBox(global_boxl,global_boxr);    
        MMC.addAllPartEvenly(global_npart);
        MMC.readyEverything();
        
        //int num_configs=500;
        //int init_step = 5000;
        MMC.run(init_step);
        
        int ngr = int(coreDirectRun::targetX.size());
        double dx = coreDirectRun::dx;
        int lgr = int(double(ngr)*dx);
        vector<double> gr(ngr,0.0);
        
        for(int i=0;i<num_configs;i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr,MMC.getGr(ngr,dx));
        }
        VectorDividedByScalar(gr,num_configs);
        double newp = coreDirectRun::getDiff(gr);
        printf("%d, MSD=%lf \n",init_step,newp);
        arrNewpOrd.push_back(floor(log10(newp)));
    }
    int minNewpOrd = 99;int argNewp;
        for(int i=0;i<arrNewpOrd.size();i++)
            if(arrNewpOrd[i]<minNewpOrd)
            {
                argNewp=i;//first minimum order
                minNewpOrd=arrNewpOrd[i];
            }
    printf("%d\n",argNewp);
    global_init_step = list1[argNewp];
    if(global_init_step<500)global_init_step=500;
    printf("final result: init_step: %d \n",global_init_step);
}
double runWithGlobal()
{
    coreDirectRun::loadTarget("standardGr.gr");
    coreMolecularMC MMC=coreMolecularMC(&potential,global_temperature,
                                        global_scaling,global_minvol);
    MMC.setBox(global_boxl,global_boxr);    
    MMC.addAllPartEvenly(global_npart);
    MMC.readyEverything();
    
    int num_configs = global_num_configs;
    int num_steps = global_num_steps;
    int init_step = global_init_step;
    MMC.run(init_step);
    
    int ngr = int(coreDirectRun::targetX.size());
    double dx = coreDirectRun::dx;
    int lgr = int(double(ngr)*dx);
    vector<double> gr(ngr,0.0);
    
    for(int i=0;i<num_configs;i++)
    {
        MMC.run(num_steps);
        VectorAddedByVector(gr,MMC.getGr(ngr,dx));
    }
    VectorDividedByScalar(gr,num_configs);
    double newp = coreDirectRun::getDiff(gr);
    return newp;
}
void generateStandardGr()
{
    //cout<<scaling<<endl;
    const char* saveGrName="standardGr.gr";
    coreMolecularMC MMC=coreMolecularMC(&potential,global_temperature,
                                        global_scaling,global_minvol);
    MMC.setBox(global_boxl,global_boxr);    
    MMC.addAllPartEvenly(global_npart);
    MMC.readyEverything();
    
    int num_configs=5000;
    int num_steps=100;
    int init_step = 5000;
    MMC.run(init_step);
    
    int npart = MMC.getNpart();
    int nsteps = 0;
    
    //FILE *fptrinit;
    //fptrinit = fopen("inittest.dat", "w");
    
    int lgr = global_LJpara*2.5;
    double dx=0.01;
    int ngr = int(double(lgr)/dx);
    vector<double> gr(ngr,0.0);//=MMC.getGr(ngr,dx);
    //printfVector(gr);
    
    for(int i=0;i<num_configs;i++)
    {
        MMC.run(num_steps);
        nsteps+=num_steps;
        VectorAddedByVector(gr,MMC.getGr(ngr,dx));
        //MMC.run(num_steps);
    }
    VectorDividedByScalar(gr,num_configs);
    write2DVector(saveGrName,arange(0,lgr,dx),gr);
    printf("standardGr generated.\n");
}

int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

char prename[50];
char lowname[50];char midname[50];char highname[50];
FILE *fptr;
char* itarea;
if(argc>=2) itarea=(argv[1]);

sprintf(prename, "chr2-%s.preconfig", itarea);
//sprintf(prename, "chr2-0.preconfig");
sprintf(lowname, "chr2-%s.lowconfig", itarea);
//sprintf(lowname, "%s.lowconfig", itarea);
sprintf(midname, "chr2-%s.midconfig", itarea);
//sprintf(lowname, "%s.lowconfig", itarea);
sprintf(highname, "chr2-%s.highconfig", itarea);
//sprintf(highname, "%s.highconfig", itarea);

fptr = fopen(prename, "r");
fscanf(fptr, "global_boxl:%lf\n", &global_boxl);
fscanf(fptr, "global_boxr:%lf\n", &global_boxr);
fscanf(fptr, "global_npart:%d\n", &global_npart);
fscanf(fptr, "global_minvol:%lf\n", &global_minvol);
fclose(fptr);

//global_boxl=0;
//global_boxr=15000;
//global_npart = 74;
//global_minvol=100;

global_temperature=1;
global_LJpara = 0.8*(global_boxr-global_boxl)/float(global_npart);
global_LJpara = round(global_LJpara*1)/1.0;
printf("LJpara is %lf\n",global_LJpara);

findStepSize();
//global_scaling=162;
//generateStandardGr();
//computeMCSteps();
//global_init_step=500;
//computeMCError(1.0,0.2);
//computeInitSteps();

//global_boxl=0;
//global_boxr=15000;
//global_npart=74;
//global_minvol=100.000000;
//global_temperature=1.000000;
//global_scaling=40.000000;
//global_num_configs=500;
//global_num_steps=10;
//global_init_step=500;

//const clock_t begin_time = clock();
//runWithGlobal();
//double timediff = float(clock() - begin_time)/CLOCKS_PER_SEC;

//printf("Recommended (fast): scaling: %lf, num_configs: %d, num_steps: %d, init_step: %d , loop time: %lf \n",global_scaling,global_num_configs,global_num_steps,global_init_step,global_timediff);
//printf("Recommended (slow): scaling: %lf, num_configs: %d, num_steps: %d, init_step: %d , loop time: %lf \n",global_scaling,high_global_num_configs,high_global_num_steps,global_init_step,high_global_timediff);

//FILE *fptr;
fptr = fopen(lowname, "w");
fprintf(fptr, "global_boxl:%lf\n", global_boxl);
fprintf(fptr, "global_boxr:%lf\n", global_boxr);
fprintf(fptr, "global_npart:%d\n", int(global_npart));
fprintf(fptr, "global_minvol:%lf\n", global_minvol);
fprintf(fptr, "global_temperature:%lf\n", global_temperature);
fprintf(fptr, "global_scaling:%lf\n", global_scaling);
fprintf(fptr, "global_num_configs:%d\n", 500);
fprintf(fptr, "global_num_steps:%d\n", 10);
fprintf(fptr, "global_init_step:%d\n", 500);
fprintf(fptr, "run time(s) for one loop:%lf\n",-1);
fprintf(fptr, "global_para:%lf\n",global_LJpara);
fclose(fptr);

fptr = fopen(midname, "w");
fprintf(fptr, "global_boxl:%lf\n", global_boxl);
fprintf(fptr, "global_boxr:%lf\n", global_boxr);
fprintf(fptr, "global_npart:%d\n", global_npart);
fprintf(fptr, "global_minvol:%lf\n", global_minvol);
fprintf(fptr, "global_temperature:%lf\n", global_temperature);
fprintf(fptr, "global_scaling:%lf\n", global_scaling);
fprintf(fptr, "global_num_configs:%d\n", 5000);
fprintf(fptr, "global_num_steps:%d\n", 10);
fprintf(fptr, "global_init_step:%d\n", 500);
fprintf(fptr, "run time(s) for one loop:%lf\n",-1);
fprintf(fptr, "global_para:%lf\n",global_LJpara);
fclose(fptr);

fptr = fopen(highname, "w");
fprintf(fptr, "global_boxl:%lf\n", global_boxl);
fprintf(fptr, "global_boxr:%lf\n", global_boxr);
fprintf(fptr, "global_npart:%d\n", global_npart);
fprintf(fptr, "global_minvol:%lf\n", global_minvol);
fprintf(fptr, "global_temperature:%lf\n", global_temperature);
fprintf(fptr, "global_scaling:%lf\n", global_scaling);
fprintf(fptr, "global_num_configs:%d\n", 50000);
fprintf(fptr, "global_num_steps:%d\n", 10);
fprintf(fptr, "global_init_step:%d\n", 1000);
fprintf(fptr, "run time(s) for one loop:%lf\n",-1);
fprintf(fptr, "global_para:%lf\n",global_LJpara);
fclose(fptr);



/*    typedef coreDoubleRMC rmc; 
    
    double rndstep = 10;
    if(argc>=2) rndstep=atof(argv[1]);
    
    vector<double> arrX;// = {0,10,20,30};arrX.push_back(100);
    arrX = linspace(0,30,31);arrX.push_back(1000);
    //coreDoubleRMC* rmc = new coreDoubleRMC(arrX);
    rmc::initPotentialX(arrX);
    //static void asignConfig(double talpha,int tnum_configs,int tnum_steps,
    //                    int tinit_step,int tnum_particles,
    //                    int tmove_scale,int tmove_cut)
    
    rmc::asignConfig(0.5,100,10,5000,200,2*rndstep,rndstep);
    rmc::readyToRun();
    rmc::run(1000);*/
    
//     for(int i=0;i<10;i++)
//     {
//         rmc::run(10000);
//         
//         char buffer[50]; // <- danger, only storage for 256 characters.
//         sprintf(buffer, "savedpotential%d_10000_5k.dat", i);
//         rmc::writePot(buffer);
//     }
    
//     arrX = linspace(0,30,10);arrX.push_back(100);
//     rmc::updatePotX(arrX);
//     rmc::run(4);
//     
//     arrX = linspace(0,30,31);arrX.push_back(100);
//     rmc::updatePotX(arrX);
//     rmc::run(2);
    
    
    
    return 0; 
} 




















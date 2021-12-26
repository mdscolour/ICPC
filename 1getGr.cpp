//#include "/home/li/bin/coreMolecularMC.h"
#include "coreMolecularMC.h"
//#include "/home/li/bin/coreMolecularMC.cpp"//only one instance
//#include "coreMolecularMC.cpp"//only one instance


double LJPotential(double r,double epsilon=1,double sigma=1) {
    return 4*epsilon*(pow(sigma/r,12)-pow(sigma/r,6));}
double QuarticBond(double r,double k0,double k1,double R){
    return 0.5*k0*pow((r-R),2)+0.25*k1*pow(r-R,4);}
double FeneBond(double r,double k,double drm,double r0){
    return -0.5*k*drm*drm*log(1-pow((r-r0)/drm,2));}
//double potential(double r){return LJPotential(r,1,1);}

double potential(double r){
    //return FeneBond(r,1,30,3)+QuarticBond(r,1,1,10)+FeneBond(r,-1,30,8);
    return LJPotential(r,3,3);
}
void getCanGr(int itarea,const char *nchr,const char *likebedname,const char *fprefix,int seclen,int inclen,int lgr=1000,double dx=5)
{
    double st = itarea*seclen;
    double ed = itarea*seclen+seclen+inclen;
    // function of potential, temperature, step size
	coreMolecularMC MMC=coreMolecularMC(&potential,1,0.45,100);
    char saveGrName[50];
    sprintf(saveGrName, "%s/chr%s-%d.midgrpre",fprefix,nchr,itarea);
    MMC.setBox(st,ed);

    FILE *fptr = fopen(likebedname, "r");
    fscanf(fptr, "%*[^\n]\n");
    fscanf(fptr, "%*[^\n]\n");
    int x;int y;int z;
    char s[50];
    double midxy=0;
    int npartcount = 0;
    double lenpartcount = 0;
    while(fscanf(fptr, "%s %d %d %d\n",&s,&x,&y,&z)==4)
    {
        //if(y<(itarea*150000))continue;
        //if(x>(itarea+1)*150000)break;
        //npartcount++;lenpartcount+=(y-x);
        //for(int pos=x;pos<y;pos++)if(pos>=(itarea*150000) && pos<=(itarea+1)*150000)MMC.part.push_back(pos);
        
        midxy = double(x+y)/2.0;
        if(midxy<st)continue;
        if(midxy>ed)break;
        npartcount++;lenpartcount+=(y-x);
        if(midxy>=st && midxy<=ed)MMC.part.push_back(midxy);
    }
    MMC.readyEverything();
    lenpartcount /= npartcount;
    //cout<<npartcount<<endl;
    //cout<<lenpartcount<<endl;
    fclose(fptr);

    char bufferpre[50];
    sprintf(bufferpre,"%s/chr%s-%d.preconfig",fprefix,nchr,itarea);
    fptr = fopen(bufferpre, "w");
    fprintf(fptr, "global_boxl:%lf\n", 0.);
    fprintf(fptr, "global_boxr:%lf\n", ed-st);
    fprintf(fptr, "global_npart:%d\n", int(npartcount));
    fprintf(fptr, "global_minvol:%lf\n", lenpartcount);
    fclose(fptr);
//     FILE *fptr = fopen("chr2-0.bed", "r");
//     double x;double y;
//     while(fscanf(fptr, "%lf %lf\n", &x, &y)==2)
//     {
//         MMC.part.push_back(y);
//     }
//     printfVector(MMC.part);
//     MMC.readyEverything();
    
    //int lgr = 1000;///set as parameter
    //double dx=5;///set as parameter
    int ngr = int(double(lgr)/dx);
    vector<double> gr=MMC.getGr(ngr,dx);

    write2DVector(saveGrName,arange(0,lgr,dx),gr);
 
}

int main(int argc,char *argv[])
{    /**************  random seed  *******************************/
SeedByTime();

//int maxlen = atoi(argv[2]);
//int maxsection = int(maxlen/50000);
//if(maxsection*50000+25000 > maxlen) maxsection--;
int maxsection = atoi(argv[2]);
int argseclen = atoi(argv[5]);
int arginclen = atoi(argv[6]);

for(int i=0;i<maxsection;i++)
    getCanGr(i,argv[1],argv[3],argv[4],argseclen,arginclen);

    return 0; 
} 




















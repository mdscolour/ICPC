#include "MyLib.h"
#define PARA_FUNC getGr
#define ROLLING_SIZE 5000

typedef double (*type_pot)(double);

double LJlike(double r, double sigma, double epsilon, double a, double b)
{
    return 4 * epsilon * (pow(sigma / r, a) - pow(sigma / r, b));
}
////// class declaration
/// class of coreMolecularMC, which will map from potential function to radial distribution function by simulating the configurations
class coreMolecularMC
{
public:
    type_pot funcNam;
    vector<double> disFor;
    vector<double> disBck;
    double minvol;

    double alpha;
    vector<double> part;
    vector<int> ptype;
    int npart;
    double scaling;
    double box_l;
    double box_r;
    double box_s;

    coreMolecularMC(type_pot tfuncNam, double talpha, double tscaling = 0.45,
                    double tminvol = PRECISION)
    {
        funcNam = tfuncNam;
        alpha = talpha;
        minvol = tminvol;
        scaling = tscaling;
    }
    double getSpacing()
    {
        double s = part[npart - 1] - part[0];
        if (s < 0)
            return s + box_s;
        return s;
    }
    void setBox(double tboxl, double tboxr)
    {
        box_l = tboxl;
        box_r = tboxr;
        box_s = tboxr - tboxl;
    }
    bool passExcludedVolume(int ipart, double step)
    {
    if (step < disFor[ipart] - minvol)
        if (-step < disBck[ipart] - minvol)
            return true;
    return false;
    }
    void addPart(double pos) { part.push_back(pos); }
    int getNpart() { return npart; }
    double getPart(int i) { return part[i]; }
    vector<double> getPartAll() { return part; }
    void printfPart()
    {
        for (int i = 0; i < part.size(); i++)
        {
            std::cout << part[i] << " ";
        }
        std::cout << '\n';
    }

    vector<double> getGr(int ngr = 160, double dx = 1)
    {
    vector<double> gr(ngr, 0.0);
    //double dismax = fabs(part[0] - part[npart-1]);
    double distemp;
    int indtemp;
    for (int i = 0; i < npart; i++)
    {
        for (int j = 0; j < npart; j++)
        {
            distemp = fabs(part[i] - part[j]);
            indtemp = (int)(ceil(distemp / dx)); //rounding up
            if (indtemp < ngr)
                gr[indtemp] += 1.0; //i&j j&i twice
            //else break;//only count until ngr*dx

            /////periodic boundary condition
            distemp = fabs(part[i] - part[j] + box_s);
            indtemp = (int)(ceil(distemp / dx)); //rounding up
            if (indtemp < ngr)
                gr[indtemp] += 1.0;
            distemp = fabs(part[i] - part[j] - box_s);
            indtemp = (int)(ceil(distemp / dx)); //rounding up
            if (indtemp < ngr)
                gr[indtemp] += 1.0;
        }
    }
    VectorDividedByScalar(gr, ((double)npart));
    VectorDividedByScalar(gr, dx);
    VectorDividedByScalar(gr, 2.0 * (((double)npart) / box_s));
    //gr = gr/npart;
    //gr = gr/dx;
    //gr = gr/(npart/dismax);
    gr[0] = 0;
    gr[1] = 0;
    return gr;
    }

    double rollingCount(int pos, int rolsiz)
    {
    double count = 0;
    double test = 0;
    int halfvol = int(minvol) / 2.0;
    double left = box_l + pos;
    double right = box_l + pos + rolsiz;
    if (right > (box_l + box_s))
        throw range_error("rolling box error");
    //     if(right<=(box_l+box_s))
    //     {
    //         for(int i=0;i<npart;i++)for(double exv=-halfvol;exv<=halfvol;exv++)
    //         {if(part[i]+exv>=left && part[i]+exv<=right)count+=1;}
    //     }
    //     else
    //     {
    //         for(int i=0;i<npart;i++)for(double exv=-halfvol;exv<=halfvol;exv++)
    //         {
    //             if(part[i]+exv>=left && part[i]+exv<=right)count+=1;
    //             else if((part[i]+exv+box_s)>=left && (part[i]+exv+box_s)<=right)count+=1;
    //         }
    //     }

    for (int i = 0; i < npart; i++)
        for (double exv = -halfvol; exv <= halfvol; exv++)
        {
            if (part[i] + exv >= left && part[i] + exv <= right)
                count += 1;
        }

    //    for(int i=0;i<npart;i++)
    //    {if(part[i]>=left && part[i]<=right)count+=1;}

    return count;
    }

    vector<double> GetAutocorrelation(int nu, double t_intersteps = 1)
    {
    int intersteps = (int)t_intersteps;
    int n = int((box_s - ROLLING_SIZE + 1) / t_intersteps);
    double tmpval;
    std::vector<double> val = {};
    double valmiu = 0.;
    double C0 = 0.;
    double C = 0.;
    //double tao_int = 0.5;
    //unsigned M = 0;

    std::vector<double> autocor = {};
    autocor.push_back(1);

    // Read in the values for val
    for (int k = 0; k < n; k++)
    {

        tmpval = rollingCount(k * intersteps, ROLLING_SIZE);
        val.push_back(tmpval);
        valmiu += tmpval;
    }

    if (n != val.size())
    {
        std::cout << "Alert! " << n << " and " << val.size() << std::endl;
    }

    //std::cout << "val[0] = " << val[0] << std::endl;

    //unsigned nu = n - n/4;

    valmiu /= n;

    // Calculate C(0)
    for (int i = 0; i < n; ++i)
        C0 += (val[i] - valmiu) * (val[i] - valmiu);

    C0 /= n;

    bool tmp = false;

    // Calculate C(t) for t>0
    for (unsigned t = 1; t < nu; ++t)
    {

        C = 0;

        // Go through all the data
        for (unsigned i = 0; i < (n - t); ++i)
            C += (val[i] - valmiu) * (val[i + t] - valmiu);

        C /= (n - t) * C0;

        ///// Integrate w/ very simple rule
        //tao_int += intersteps*C;

        autocor.push_back(C);
        // 		fout << t*intersteps << "    " << C << "    " << tao_int << std::endl;
        // 		if ( (tmp==false) && ((t*intersteps) >= 10*tao_int) ) {
        // 			M = (unsigned int)tao_int;
        // 			tmp = true;
        // 			break;
        // 		}
    }
    //printfVector(val);
    //printfVector(autocor);
    return autocor;
    }
    
    double getPotMin(double maxX = 30, double dx = 0.01)
    {
    vector<double> x = arange(0, maxX, dx);
    int numRange = int(x.size());
    double potminy = 1e99;
    int potminx = -1;
    for (int i = 0; i < numRange; i++)
    {
        double tmpy = funcNam(x[i]);
        if (tmpy <= potminy)
        {
            potminy = tmpy;
            potminx = i;
        }
    }
    return x[potminx];
    }

    void addAllPartEvenly(int num)
    {
    double spacing = box_s / double(num);
    for (int i = 0; i < num; i++)
    {
        addPart(box_l + 0.5 * spacing + spacing * i);
    }
    if (part[num - 1] >= box_r)
        throw range_error("initializaiton in outside");
    if (part[0] <= box_l)
        throw range_error("initializaiton in outside");
    }

    bool illedPot(double maxX=5,double sep=0.1,double tol=0.1)
    {
    vector<double> x = arange(minvol, maxX+minvol,sep);
    int numRange = int(x.size());
    double potminy = 1e99;
    int potminx = -1;
    for (int i = 0; i < numRange; i++)
    {
        double tmpy = funcNam(x[i]);
        if (tmpy <= potminy)
        {
            potminy = tmpy;
            potminx = i;
        }
    }
    if(x[potminx]-minvol<=tol)return true;
    else return false;
    }


    bool readyEverything()
    {
        updateAllDis();
        return true;
    }

    void updateAllDis()
    {
    npart = (int)(part.size());

    for (int i = 0; i < npart; i++)
        if (part[i] < box_l || part[i] > box_r)
            throw range_error("outside box 1");

    double dis;
    dis = part[0] - part[npart - 1] + box_s;
    if (dis < 0)
        throw range_error("outside box 2");
    disBck.push_back(dis);

    for (int i = 0; i < npart - 1; i++)
    {
        dis = fabs(part[i] - part[i + 1]);
        disFor.push_back(dis);
        disBck.push_back(dis);
    }

    dis = part[0] - part[npart - 1] + box_s;
    if (dis < 0)
        throw range_error("outside box 3");
    disFor.push_back(dis);
    }

    double proposal(int ipart)
    {
    double step;
    ///////     double scaling = max(disFor[ipart],disBck[ipart])*0.9;

    //    double rnda;double rndb;rnda = -disBck[ipart]+minvol;rndb = disFor[ipart]-minvol;

    for (int count = 0; count < 1; count++) //only once now
    {
        step = RNG_NAME() * 2 * scaling - scaling;
        //        step = (rndb - rnda) * RNG_NAME() + rnda;
        return step;
    }
    //cout<<ipart<<endl<<disFor[ipart]<<endl<<disBck[ipart]<<endl;
    throw logic_error("error 1"); //should not reach here
    }

    inline bool acceptOrNot(double oldp, double newp)
    {
    if (newp <= oldp)
        return true;
    if (RNG_NAME() < exp(alpha * (oldp - newp)))
        return true;
    return false;
    }

    bool runOne()
    {
    bool resFlag = false;
    for (int count = 0; count < (npart * 0.8); count++)
    {
        int ipart = Random_integer_uniform(0, npart);

        //if no space to move
        //should return false normally but tolerate here to check system collapse or not,
        if ((disFor[ipart] + disBck[ipart]) <= 2 * minvol + PRECISION)
        {
            continue;
        } //tolerate here

        double prop = proposal(ipart);

        if (passExcludedVolume(ipart, prop) == false)
            return false;

        resFlag = acceptOrNot(getOriPot(ipart), getNewPot(ipart, prop));
        if (resFlag)
            updateOne(ipart, prop);
        return resFlag;
    }
    printfVector(part);
    throw logic_error("error 2: system collapse"); //should not reach here
    }

    int run(int num)
    {
    int resCount = 0;
    bool resFlag;
    for (int i = 0; i < num; i++)
    {
        for (int mcsi = 0; mcsi < npart; mcsi++)
        {
            resFlag = runOne();
            //cout<<resFlag<<'\t';
            if (resFlag)
                resCount++;
        }
    }
    return resCount;
    //printf("%d step    ",num);
    //printf("scaling: %f, success rate: %f\n",scaling,double(resCount)/double(num*npart));
    //cout<<scaling<<" "<<double(resCount)/double(num*npart)<<endl;
    }

    double getOriPot(int ipart)
    {
    double dis1 = disFor[ipart];
    double dis2 = disBck[ipart];
    //if(fabs(dis1-0)<PRECISION)return funcNam(dis2);
    //if(fabs(dis2-0)<PRECISION)return funcNam(dis1);
    return funcNam(dis1) + funcNam(dis2);
    }

    double getNewPot(int ipart, double prop)
    {
    double dis1 = disFor[ipart];
    double dis2 = disBck[ipart];
    //if(fabs(dis1-0)<PRECISION)return funcNam(dis2+prop);
    //if(fabs(dis2-0)<PRECISION)return funcNam(dis1-prop);
    return funcNam(dis1 - prop) + funcNam(dis2 + prop);
    }

    void updateOne(int i, double prop)
    {
    part[i] = part[i] + prop;
    if (part[i] > box_r)
        part[i] -= box_s;
    if (part[i] < box_l)
        part[i] += box_s;

    disFor[i] = disFor[i] - prop;
    disBck[i] = disBck[i] + prop;
    if (i == 0)
    {
        disFor[npart - 1] = disFor[npart - 1] + prop;
        disBck[i + 1] = disBck[i + 1] - prop;
    }
    else if (i == npart - 1)
    {
        disFor[i - 1] = disFor[i - 1] + prop;
        disBck[0] = disBck[0] - prop;
    }
    else
    {
        disFor[i - 1] = disFor[i - 1] + prop;
        disBck[i + 1] = disBck[i + 1] - prop;
    }
    }
};
class coreDirectRun
{
public:
    static vector<double> targetX;
    static vector<double> targetY;
    static vector<double> curpara;
    static coreMolecularMC *MMC;
    static double dx;
    static vector<double> savegr;///static so changing every time

    static double global_scaling;
    static double global_boxl;
    static double global_boxr;
    static int global_npart;
    static double global_minvol;
    static double global_temperature;
    static int global_num_configs;
    static int global_num_steps;
    static int global_init_step;

    static double modelPot(double r)
    {
        //return LJlike(r, 162,curpara[1],9,7);
        return LJlike(r, curpara[0], curpara[1], curpara[2], curpara[3]);
    }
    static double assignAndRun(double p1,double p2,double p3,double p4)
    {
        assignPara(p1,p2,p3,p4);
        return runOne();
    }

    static void assignPara(double p1,double p2,double p3,double p4)
    {
        //if(int(curpara.size())!=4)
        //{
            curpara.clear();
            curpara.push_back(p1);
            curpara.push_back(p2);
            curpara.push_back(p3);
            curpara.push_back(p4);
        //}
        //else
        //{
       // curpara[0]=p1;
       // curpara[1]=p2;
       // curpara[2]=p3;
       // curpara[3]=p4;
       // }
    }
    coreDirectRun()
    {
    }

    static void newMMC()
    {
        MMC = new coreMolecularMC(&modelPot, global_temperature, 
            global_scaling, global_minvol);
        MMC->setBox(global_boxl,global_boxr);
        MMC->addAllPartEvenly(global_npart);
        MMC->readyEverything();
    }

    static void loadTarget(const char *name)
    {
        FILE *fptr = fopen(name, "r");
        double x;
        double y;
        targetX.clear();
        targetY.clear();

        while (fscanf(fptr, "%lf %lf\n", &x, &y) == 2)
        {
            targetX.push_back(x);
            targetY.push_back(y);
        }
        dx = targetX[1] - targetX[0];
        //cout<<"load targetGr.gr ~success~"<<endl;
    }
    
    static double getDiff(vector<double> gr)
    {
        int numGr = int(gr.size());
        if (numGr != targetY.size())
            printf("Gr size not equal:%d %d\n", numGr, int(targetY.size()));

        double sum = 0;
        for (int i = 0; i < numGr; i++)
            sum += pow(gr[i] - targetY[i], 2);
        //for(int i=0;i<numGr;i++)sum+=(gr[i]-targetY[i])*(gr[i]-targetY[i]);
        return sum / (double(numGr));
        //return sum;
    }

    static double runOne()
    {
        double illedFlag = MMC->illedPot();
        if (illedFlag)
        {   
            //printf("illed potential\n");
            return 99999;
        }
        savegr.clear();
        savegr = runMMC();
        return getDiff(savegr);
    }
    static vector<double> runMMC()
    {
        deleteMMC();
        newMMC();

        MMC->run(global_init_step);
        int ngrpoint = int(targetX.size());
        vector<double> gr = MMC->getGr(ngrpoint, dx);
        vector<double> tgr;
        double countgr = 1;
        for (int i = 0; i < global_num_configs; i++)
        {
            MMC->run(global_num_steps);
            tgr = MMC->getGr(ngrpoint, dx);
            VectorAddedByVector(gr, tgr);
            countgr++;
        }
        VectorDividedByScalar(gr, countgr);
        return gr;
    }
   
    static void readyToRun()
    {
        newMMC();
    } 
    
    static void deleteMMC() { delete MMC; };
    static void readConfig(const char *name = "config.config")
    {
        FILE *fptr = fopen(name, "r");
        fscanf(fptr, "global_boxl:%lf\n", &global_boxl);
        fscanf(fptr, "global_boxr:%lf\n", &global_boxr);
        fscanf(fptr, "global_npart:%d\n", &global_npart);
        fscanf(fptr, "global_minvol:%lf\n", &global_minvol);
        fscanf(fptr, "global_temperature:%lf\n", &global_temperature);
        fscanf(fptr, "global_scaling:%lf\n", &global_scaling);
        fscanf(fptr, "global_num_configs:%d\n", &global_num_configs);
        fscanf(fptr, "global_num_steps:%d\n", &global_num_steps);
        fscanf(fptr, "global_init_step:%d\n", &global_init_step);
        fclose(fptr);
    }
};

vector<double> coreDirectRun::targetX={};
vector<double> coreDirectRun::targetY={};
vector<double> coreDirectRun::curpara={};
coreMolecularMC *coreDirectRun::MMC=0;
double coreDirectRun::dx=0;
vector<double> coreDirectRun::savegr={};

double coreDirectRun::global_scaling=0;
double coreDirectRun::global_boxl=0;
double coreDirectRun::global_boxr=0;
int coreDirectRun::global_npart=0;
double coreDirectRun::global_minvol=0;
double coreDirectRun::global_temperature=0;
int coreDirectRun::global_num_configs=0;
int coreDirectRun::global_num_steps=0;
int coreDirectRun::global_init_step=0;

class coreTableRun
{
public:
    coreTableRun() {}
    //// ersearchy ersearchx must be set outside.
    static vector<double> ersearchy;
    static vector<double> ersearchx;
    static vector<double> savegr;///static so changing every time
    //static coreMolecularMC* MMC;////change to local one

    static double erPotential(double r)
    {
        return interpolate(ersearchx, ersearchy, r, true);
    }

    static vector<double> runPotInput(double t, double scal, 
        double minvol, double tbox_l, double tbox_r, int npart, 
        int num_configs, int num_steps, int init_step, int ngr, 
        double dx)
    {
        //set ersearchx and ersearchy in advance

        coreMolecularMC MMC = coreMolecularMC(&erPotential, t, scal, minvol);

        MMC.setBox(tbox_l, tbox_r);
        MMC.addAllPartEvenly(npart);
        MMC.readyEverything();

        MMC.run(init_step);

        //int lgr = int(double(ngr)*dx);
        vector<double> gr = MMC.getGr(ngr, dx);

        for (int i = 0; i < num_configs; i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr, MMC.getGr(ngr, dx));
        }
        VectorDividedByScalar(gr, num_configs+1);
        savegr.clear();
        savegr = gr;
        return gr;
        //double newp = coreDoubleRMC::getDiff(gr);
        //cout<<newp<<endl;
    }
};
vector<double> coreTableRun::ersearchy = {};
vector<double> coreTableRun::ersearchx = {};
vector<double> coreTableRun::savegr = {};

//all inline and static except run()
class ClassErSearch
{
public:
    vector<double> ersearchymid;
    vector<double> ersearchymodel;
    vector<double> vecval;
    int vnum;
    FILE *fptr;

    ClassErSearch() {}

    void doErSearch(double outer_define,
    const char* tag,const char* cog,const char* savenam)
    {
        coreDirectRun::loadTarget(tag);
        coreDirectRun::readConfig(cog);
        coreDirectRun::readyToRun();
        fptr = fopen(savenam, "w");
        
        //double outer_define; 
        //if(argc>=2)outer_define = atof(argv[1]);
        //else printf("Error, no outer input!\n");
        ersearchymodel = {outer_define,-99, -99, -99}; //-99 refer to ersearchymid
        ersearchymid = {-99, -99, -99};
        vecval = arange(1, 19, 1);
        vnum = (int)vecval.size();

        setCombiAndRun(int(ersearchymid.size()) - 1);
        fclose(fptr);
    }
    
    void tolRunErSearch()
    {
        vector<double> para={};
        para.clear();
        int pcount = 0;
        for(int i=0;i<ersearchymodel.size();i++)
        {
            if(ersearchymodel[i]!=-99){para.push_back(ersearchymodel[i]);}
            else {para.push_back(ersearchymid[pcount]);pcount++;}
        }
        double newp = coreDirectRun::assignAndRun(para[0],para[1],para[2],para[3]);
        fprintf(fptr, "%.5f %.5f %.5f %.5f %.5f %.5f\n",
        0.0,para[0],para[1],para[2],para[3],newp);
    }

    //void setCombiAndRun(int lv);
    void setCombiAndRun(int lv)
    {
        if (lv == 0)
            for (int i = 0; i < vnum; i++)
            {
                ersearchymid[lv] = vecval[i];
                tolRunErSearch();
            }
        else
        {
            for (int i = 0; i < vnum; i++)
            {
                ersearchymid[lv] = vecval[i];
                setCombiAndRun(lv - 1);
            }
        }
    }
//    static void updateErSearchy()
//    {
//        ersearchy.clear();
//        int numx = ersearchx.size();
//        for (int j = 0; j < numx; j++)
//            ersearchy.push_back(LJlike(ersearchx[j],
//                                       ersearchymodel[0], ersearchymid[0], ersearchymid[1], ersearchymid[2]));
//        if (ersearchy.size() != numx)
//            printf("input of ersearchy error\n");
//    }
    double runSpecificPot(bool writeGr = false)
    {
        //set ersearchx and ersearchy in advance
        coreDirectRun::loadTarget("targetGr.gr");
        coreMolecularMC MMC = coreMolecularMC(&coreTableRun::erPotential, 1, 0.45);

        MMC.setBox(0, 800);
        int npart = 200;
        MMC.addAllPartEvenly(npart);
        MMC.readyEverything();

        int num_configs = 500;
        int num_steps = 1;
        int init_step = 500;
        MMC.run(init_step);

        //int lgr = 35;
        //double dx=0.01;
        //int ngr = int(double(lgr)/dx);
        int ngr = int(coreDirectRun::targetX.size());
        double dx = coreDirectRun::dx;
        int lgr = int(double(ngr) * dx);
        vector<double> gr = MMC.getGr(ngr, dx);

        FILE *fptr;
        if (writeGr)
            fptr = fopen("test.out", "w");
        for (int i = 0; i < num_configs; i++)
        {
            MMC.run(num_steps);
            VectorAddedByVector(gr, MMC.getGr(ngr, dx));

            if (writeGr)
            {
                fprintf(fptr, "ITEM: TIMESTEP\n%d\n", init_step + i * num_steps);
                fprintf(fptr, "ITEM: NUMBER OF ATOMS\n%d\n", npart);
                fprintf(fptr, "ITEM: BOX BOUNDS pp pp pp\n");
                fprintf(fptr, "0 1\n");
                fprintf(fptr, "0 1\n");
                fprintf(fptr, "0 1\n");
                fprintf(fptr, "ITEM: ATOMS id xs ys zs\n");
                for (int i = 0; i < npart; i++)
                {
                    fprintf(fptr, "%d %.5f %.5f %.5f\n", i, MMC.getPart(i), 0.0, 0.0);
                }
            }
        }
        VectorDividedByScalar(gr, num_configs);
        double newp = coreDirectRun::getDiff(gr);
        //cout<<newp<<endl;

        if (writeGr)
            fclose(fptr);
        if (writeGr)
            write2DVector("MC result.gr", arange(0, lgr, dx), gr);
        return newp;
    }
    void randomSearch()
    {
        //set ersearchx in advance outside
        coreDirectRun::loadTarget("targetGr.gr");
        int numx = int(coreTableRun::ersearchx.size());
        
        vector<double> para(4, 0.0);
        vector<double> tgr={},bestpara={};
        double newp,bestp;
        
        //for(int i=0;i<100000;i++)
        for (int i = 0; i < 10000; i++)
        {
            coreTableRun::ersearchy.clear();
            for (int k = 0; k < 4; k++)
                para[k] = 19 * RNG_NAME() + 1;
            for (int j = 0; j < numx; j++)
                coreTableRun::ersearchy.push_back(LJlike(coreTableRun::ersearchx[j],
                                           para[0], para[1], para[2], para[3]));
            //ersearchy.push_back(0);
            if (numx != coreTableRun::ersearchy.size())
                cout << "error 111" << endl;

            tgr = coreTableRun::runPotInput(1,40,100,0,15000,74,500,10,500,
                int(coreDirectRun::targetX.size()),coreDirectRun::dx);
            newp = coreDirectRun::getDiff(tgr);
            if (newp <= bestp)
            {
                bestpara = para;
                bestp = newp;
                cout << bestp << "   ";
                printfVector(bestpara);
            }
        }
        //ersearchy.clear();ersearchy=best;
        //runSpecificPot(true);
        //write2DVector("erPotential.rndpot",ersearchx,best);
    }
};

class classAnneal
{
    public:
        vector<double> oldpara;
        vector<double> newpara;
        vector<double> vRndSize;
        double energy;
        bool digitFlag=true;
        void readyToRun(const char* tag,const char* cog,vector<double> toldpara,vector<double> tvRndSize)
        {
            coreDirectRun::loadTarget(tag);
            coreDirectRun::readConfig(cog);
            coreDirectRun::readyToRun();
            oldpara.clear();oldpara.assign(toldpara.begin(), toldpara.end());
            newpara = oldpara;
            energy = 99999;
            vRndSize.clear();vRndSize.assign(tvRndSize.begin(), tvRndSize.end());
            //digitFlag = true;
        }

        vector<double> anneal(double T,double T_min,double lambda,int SAcount_max,
            const char* name="potparaAll.SApot")
        {
            //double T = 100;
            //double T_min = 0.0006;
            //double lambda = 0.992;
            
            int SAcount = 0;
            int i;
            double new_E;
            bool apFlag;
            FILE *fptr;
            
            fptr = fopen(name, "w");
            while(T > T_min && SAcount<SAcount_max)
            {
                i = 1;
                while(i <= 1)
                {
                    energy = coreDirectRun::assignAndRun(oldpara[0],oldpara[1],oldpara[2],oldpara[3]);
                    newpara = updateNew(oldpara);
                    new_E = coreDirectRun::assignAndRun(newpara[0],newpara[1],newpara[2],newpara[3]);
                    apFlag = acceptOrNot(energy, new_E, T);
                    if(apFlag)
                    {
                        oldpara = newpara;
                        energy = new_E;
                    }
                    i++;
                    SAcount++;
                }
                T *= lambda;
                fprintf(fptr, "%e %.5f %.5f %.5f %.5f %.5f\n",T,oldpara[0],oldpara[1],oldpara[2],oldpara[3],energy);
                //printf("%e %.5f %.5f %.5f %.5f,%.5f\n",T,oldpara[0],oldpara[1],oldpara[2],oldpara[3],energy);
            }
            fclose(fptr);
            return oldpara;
        }
        bool acceptOrNot(double oldp, double newp, double T)
        {
            if (newp <= oldp)return true;
            if (RNG_NAME() < exp((oldp - newp)/T))return true;
            return false;
        }
        vector<double> updateNew(vector<double> v1)
        {
            //// parameter domains set manually
            vector<double> v2={};
            if(digitFlag)
            {
                if(vRndSize[0]==0)v2.push_back(v1[0]);
                else v2.push_back( fPeriodRestrict(v1[0] + round(RNG_NAME() * (2*vRndSize[0])-vRndSize[0]),100,350));
                
                if(vRndSize[1]==0)v2.push_back(v1[1]);
                else v2.push_back( fPeriodRestrict(v1[1] + round(RNG_NAME() * (2*vRndSize[1])-vRndSize[1]),1,30));
                
                if(vRndSize[2]==0)v2.push_back(v1[2]);
                else v2.push_back( fPeriodRestrict(v1[2] + round(RNG_NAME() * (2*vRndSize[2])-vRndSize[2]),1,20));
                
                if(vRndSize[3]==0)v2.push_back(v1[3]);
                else v2.push_back( fPeriodRestrict(v1[3] + round(RNG_NAME() * (2*vRndSize[3])-vRndSize[3]),1,20));
            }
            else
            {
                vector<double> v2={};
                if(vRndSize[0]==0)v2.push_back(v1[0]);
                else v2.push_back( fPeriodRestrict(v1[0] + (RNG_NAME() * (2*vRndSize[0])-vRndSize[0]),100,350));
                
                if(vRndSize[1]==0)v2.push_back(v1[1]);
                else v2.push_back( fPeriodRestrict(v1[1] + (RNG_NAME() * (2*vRndSize[1])-vRndSize[1]),1,30));
                
                if(vRndSize[2]==0)v2.push_back(v1[2]);
                else v2.push_back( fPeriodRestrict(v1[2] + (RNG_NAME() * (2*vRndSize[2])-vRndSize[2]),1,20));
                
                if(vRndSize[3]==0)v2.push_back(v1[3]);
                else v2.push_back( fPeriodRestrict(v1[3] + (RNG_NAME() * (2*vRndSize[3])-vRndSize[3]),1,20));
            }
            // v2.push_back(12);
            // v2.push_back(6);
            return v2;
        }

};

class classHillClimb
{
    public:
        vector<double> oldpara;
        vector<double> newpara;
        vector<double> bestpara;
        vector<double> vRndSize;
        double energy;
        //bool digitFlag=true;
        void readyToRun(const char* tag,const char* cog,vector<double> toldpara,vector<double> tvRndSize)
        {
            coreDirectRun::loadTarget(tag);
            coreDirectRun::readConfig(cog);
            coreDirectRun::readyToRun();
            oldpara.clear();oldpara.assign(toldpara.begin(), toldpara.end());
            newpara.clear();newpara.assign(toldpara.begin(), toldpara.end());
            bestpara.clear();bestpara.assign(toldpara.begin(), toldpara.end());
            energy = 99999;
            vRndSize.clear();vRndSize.assign(tvRndSize.begin(), tvRndSize.end());
            //digitFlag = true;
        }

        vector<double> climb(const char* name="potparaAll.HCpot")
        {
            double new_E;
            FILE *fptr;
            fptr = fopen(name, "w");
            int runcount = 0;
            while(true)
            {
            energy = coreDirectRun::assignAndRun(oldpara[0],oldpara[1],oldpara[2],oldpara[3]);

            ////// block search
            // for(int i1=-vRndSize[0];i1<=vRndSize[0];i1++)
            // for(int i2=-vRndSize[1];i2<=vRndSize[1];i2++)
            // for(int i3=-vRndSize[2];i3<=vRndSize[2];i3++)
            // for(int i4=-vRndSize[3];i4<=vRndSize[3];i4++)
            // {
            //     //if(i1==0&&i2==0&&i3==0&&i4==0)continue;
            //     newpara.clear();
            //     newpara.push_back(oldpara[0]+i1);
            //     newpara.push_back(oldpara[1]+i2);
            //     newpara.push_back(oldpara[2]+i3);
            //     newpara.push_back(oldpara[3]+i4);
            ///// line search
            for(int ip=0;ip<4;ip++)
            for(int add=-vRndSize[ip];add<=vRndSize[ip];add++)
            {
                if(add==0)continue;
                newpara[ip] += add;

                new_E = coreDirectRun::assignAndRun(newpara[0],newpara[1],newpara[2],newpara[3]);
                if(new_E <= energy)
                {
                    energy = new_E;
                    bestpara = newpara;
                }
                newpara[ip] = oldpara[ip];
            }
                
            fprintf(fptr, "%f %.5f %.5f %.5f %.5f %.5f\n",float(runcount),bestpara[0],bestpara[1],bestpara[2],bestpara[3],energy);
            //printf("%e %.5f %.5f %.5f %.5f,%.5f\n",T,oldpara[0],oldpara[1],oldpara[2],oldpara[3],energy);
            if(bestpara == oldpara)break;
            oldpara = bestpara;
            runcount++;
            }//// end while
            fclose(fptr);
            return oldpara;
        }
};

//void readyToRun(const char* tag,const char* cog,vector<double> toldpara,vector<double> tvRndSize)
//vector<double> anneal(double T,double T_min,double lambda,int SAcount_max,const char* name="potparaAll.SApot")
void SAMethod1(double a,double b,const char* tagnam="targetGr.gr",
const char* lowname="config.midconfig",const char* savenam="potparaAll.SApot")
{
    classAnneal* ann;
    ann = new classAnneal();
    //char buffer[50];
    vector<double> told,tsiz;
    //vector<vector<double>> recPara = {};
    vector<double> tmpPara = {};

    //for(double a=2;a<19;a++)
    //for(double b=1;b<a;b++)
    if(true)
    {
        //sprintf(buffer, "potSAM1%d_%d.LJpot",int(a),int(b));
        delete ann;
        ann = new classAnneal();
        told = {180.,10,a,b};
        tsiz = {30.,3.,0.0,0.0};//// 0 will keep the value fixed
        
        ann->readyToRun(tagnam,lowname,told,tsiz);
        ann->digitFlag=true;
        tmpPara = ann->anneal(100,0.001,0.94,150,savenam);
        //fprintf(fptr, "%.5f %.5f %.5f %.5f %.5f %.5f\n",
        //    0.0,tmpPara[0],tmpPara[1],tmpPara[2],tmpPara[3],ann->energy);
        //recPara.push_back(tmpPara);
        //printf("a=%d b=%d finished.\n",a,b);
    }

// coreDirectRun::loadTarget(tagnam);
// coreDirectRun::readConfig(highname);
// coreDirectRun::readyToRun();

// int nrec = int(recPara.size());
// double oldenergy = 99999;
// double energy;
// vector<double> best={};
// double p1=0,p2=0,p3=0,p4=0;
// FILE *fptr;
// //char buffer[50];
// sprintf(buffer, "potparaRes.SAM1exppot"); 

// fptr = fopen(buffer, "w");
// for(double ii=0;ii<nrec;ii+=1)
// {
//     p1=recPara[ii][0];
//     p2=recPara[ii][1];
//     p3=recPara[ii][2];
//     p4=recPara[ii][3];
//     energy = coreDirectRun::assignAndRun(p1,p2,p3,p4);
//     fprintf(fptr, "%.5f %.5f %.5f %.5f %.5f\n",p1,p2,p3,p4,energy);
//     printf("%.5f %.5f %.5f %.5f,%.5f\n",p1,p2,p3,p4,energy);
//     if(energy<oldenergy){best=recPara[ii];oldenergy=energy;}
// }
// cout<<"best fit:";
// printfVector(best);
}

void SAMethod1B(double a,double b,const char* tagnam="targetGr.gr",
const char* lowname="config.midconfig",const char* savenam="potparaAll.SApot")
{
    classAnneal* ann;
    ann = new classAnneal();
    //char buffer[50];
    vector<double> told,tsiz;
    //vector<vector<double>> recPara = {};
    vector<double> tmpPara = {};

    if(true)
    {
        //sprintf(buffer, "potSAM1%d_%d.LJpot",int(a),int(b));
        delete ann;
        ann = new classAnneal();
        told = {180.,10,a,b};
        tsiz = {30.,3.,0.0,0.0};//// 0 will keep the value fixed
        
        ann->readyToRun(tagnam,lowname,told,tsiz);
        ann->digitFlag=false;
        tmpPara = ann->anneal(100,0.001,0.94,150,savenam);
    }
}

void SAMethod2(const char* tagnam="targetGr.gr",
const char* midname="config.midconfig",const char* savenam="potparaAll.SApot")
{
    classAnneal* ann;
    ann = new classAnneal();
    vector<double> told,tsiz;
    vector<double> tmpPara = {};

    //for(double nrep=0;nrep<20;nrep++)
    if(true)
    {
        delete ann;
        ann = new classAnneal();
        tsiz = {30.,3.,3.,3.};
        told.clear();
        told.push_back(int((230-130)*RNG_NAME()+130));
        told.push_back(int((20-1)*RNG_NAME()+1));
        told.push_back(int((20-1)*RNG_NAME()+1));
        told.push_back(int((20-1)*RNG_NAME()+1));
        
        ann->readyToRun(tagnam,midname,told,tsiz);
        ann->digitFlag=false;
        tmpPara = ann->anneal(100,0.001,0.98,500,savenam);
    }

// coreDirectRun::loadTarget("standardLJTest.gr");
// coreDirectRun::readConfig("high.config");
// coreDirectRun::readyToRun();

// int nrec = int(recPara.size());
// double oldenergy = 99999;
// double energy;
// vector<double> best={};
// double p1=0,p2=0,p3=0,p4=0;
// FILE *fptr;
// //char buffer[50];
// sprintf(buffer, "potparaRes.SAM2pot"); 

// fptr = fopen(buffer, "w");
// for(double ii=0;ii<nrec;ii+=1)
// {
//     p1=recPara[ii][0];
//     p2=recPara[ii][1];
//     p3=recPara[ii][2];
//     p4=recPara[ii][3];
//     energy = coreDirectRun::assignAndRun(p1,p2,p3,p4);
//     fprintf(fptr, "%.5f %.5f %.5f %.5f %.5f\n",p1,p2,p3,p4,energy);
//     printf("%.5f %.5f %.5f %.5f,%.5f\n",p1,p2,p3,p4,energy);
//     if(energy<oldenergy){best=recPara[ii];oldenergy=energy;}
// }
// cout<<"best fit:";
// printfVector(best);
}

void HCMethod1(double a,double b,const char* tagnam="targetGr.gr",
const char* lowname="config.midconfig",const char* savenam="potparaAll.SApot")
{
    classHillClimb* hill;
    hill = new classHillClimb();
    //char buffer[50];
    vector<double> told,tsiz;
    //vector<vector<double>> recPara = {};
    vector<double> tmpPara = {};

    //for(double a=2;a<19;a++)
    //for(double b=1;b<a;b++)
    if(true)
    {
        //sprintf(buffer, "potSAM1%d_%d.LJpot",int(a),int(b));
        delete hill;
        hill = new classHillClimb();
        told = {180.,10,a,b};
        tsiz = {15.,3.,0.0,0.0};//// 0 will keep the value fixed
        
        hill->readyToRun(tagnam,lowname,told,tsiz);
        //hill->digitFlag=true;
        tmpPara = hill->climb(savenam);
    }
}

void HCMethod2(const char* tagnam="targetGr.gr",
const char* midname="config.midconfig",const char* savenam="potparaAll.SApot")
{
    classHillClimb* hill;
    hill = new classHillClimb();
    vector<double> told,tsiz;
    vector<double> tmpPara = {};

    //for(double nrep=0;nrep<20;nrep++)
    if(true)
    {
        delete hill;
        hill = new classHillClimb();
        tsiz = {15.,3.,3.,3.};
        told.clear();
        told.push_back(int((230-130)*RNG_NAME()+130));
        told.push_back(int((20-1)*RNG_NAME()+1));
        told.push_back(int((20-1)*RNG_NAME()+1));
        told.push_back(int((20-1)*RNG_NAME()+1));
        
        hill->readyToRun(tagnam,midname,told,tsiz);
        //hill->digitFlag=false;
        tmpPara = hill->climb(savenam);
    }
}
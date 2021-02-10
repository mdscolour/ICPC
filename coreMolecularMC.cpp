#include "coreMolecularMC.h"

//////////////////////////////////////////////////////////////////////////// Real content starts here
///////////////////////////////////////////////////////////////////

///// class function for coreMolecularMC
void coreMolecularMC::printfDisAll()
{
    for (double n : disFor)
    {
        std::cout << n << " ";
    }
    std::cout << '\n';
    for (double n : disBck)
    {
        std::cout << n << " ";
    }
    std::cout << '\n';
}
bool coreMolecularMC::thermalOrNotByLastThree(vector<double> &v, double spacing, int sizV)
{
    v.push_back(spacing);
    int numV = int(v.size());
    if (numV < sizV)
        return false;
    if (numV > sizV)
        v.erase(v.begin());
    double minV = 0, maxV = 0;
    for (int i = 1; i < sizV; i++)
    {
        if (v[i] <= v[minV])
            minV = i;
        if (v[i] >= v[maxV])
            maxV = i;
    }
    if (minV > sizV - 3 || maxV > sizV - 3)
        return false;
    return true;
}
int coreMolecularMC::runTilThermalByLastThree(int each_step)
{ //not updated
    int nsteps = 0;
    double spacing;
    vector<double> v(0);
    while (true)
    {
        run(each_step);
        nsteps += each_step;
        //spacing = part[npart-1]-part[0];
        printf("initializaiton %.2f mil. steps....done, %.6f spacing\n", double(nsteps) / 1000000.0, spacing / npart);
        if (thermalOrNotByLastThree(v, spacing, 10))
            break;
    }
    return nsteps;
}
int coreMolecularMC::runTilThermalBySpacing(int each_step, bool upgoing, int sizV)
{ //not updated
    int nsteps = 0;
    if (upgoing)
    {
        double spacing = 0.0, sum_spacing = 0.0, sum_spacing_old = 0;
        while (true)
        {
            for (int i = 0; i < sizV; i++)
            {
                run(each_step);
                //spacing = part[npart-1]-part[0];
                sum_spacing += spacing;
                nsteps += each_step;
                printf("initializaiton %.2f mil. steps....done, %.6f spacing\n", double(nsteps) / 1000000.0, spacing / npart);
            }
            if (sum_spacing <= sum_spacing_old)
                break;
            else
            {
                sum_spacing_old = sum_spacing;
                sum_spacing = 0;
            }
        }
    }
    else
    {
        double spacing = 0.0, sum_spacing = 0.0, sum_spacing_old = 1e99;
        while (true)
        {
            for (int i = 0; i < sizV; i++)
            {
                run(each_step);
                //spacing = part[npart-1]-part[0];
                sum_spacing += spacing;
                nsteps += each_step;
                printf("initializaiton %.2f mil. steps....done, %.6f spacing\n", double(nsteps) / 1000000.0, spacing / npart);
            }
            if (sum_spacing >= sum_spacing_old)
                break;
            else
            {
                sum_spacing_old = sum_spacing;
                sum_spacing = 0;
            }
        }
    }
    return nsteps;
}
double coreMolecularMC::getPotMin(double maxX, double dx)
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
void coreMolecularMC::addAllPartByMinPot(int num, double potRange, double potPrecision)
{
    double spacing = getPotMin(potRange, potPrecision);
    for (int i = 1; i <= num; i++)
    {
        addPart(box_l + spacing * i);
    }
    if (part[num - 1] >= box_r)
        throw range_error("initializaiton in outside");
    if (part[0] <= box_l)
        throw range_error("initializaiton in outside");
}
void coreMolecularMC::addAllPartEvenly(int num)
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
vector<double> coreMolecularMC::getGr(int ngr, double dx)
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
double coreMolecularMC::rollingCount(int pos, int rolsiz)
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
vector<double> coreMolecularMC::GetAutocorrelation(int nu, double t_intersteps)
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

void coreMolecularMC::updateAllDis()
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
bool coreMolecularMC::passExcludedVolume(int ipart, double step)
{
    if (step < disFor[ipart] - minvol)
        if (-step < disBck[ipart] - minvol)
            return true;
    return false;
}
//essential function, generate a random step
double coreMolecularMC::proposal(int ipart)
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
inline bool coreMolecularMC::acceptOrNot(double oldp, double newp)
{
    if (newp <= oldp)
        return true;
    if (RNG_NAME() < exp(alpha * (oldp - newp)))
        return true;
    return false;
}
double coreMolecularMC::getOriPot(int ipart)
{
    double dis1 = disFor[ipart];
    double dis2 = disBck[ipart];
    //if(fabs(dis1-0)<PRECISION)return funcNam(dis2);
    //if(fabs(dis2-0)<PRECISION)return funcNam(dis1);
    return funcNam(dis1) + funcNam(dis2);
}
double coreMolecularMC::getNewPot(int ipart, double prop)
{
    double dis1 = disFor[ipart];
    double dis2 = disBck[ipart];
    //if(fabs(dis1-0)<PRECISION)return funcNam(dis2+prop);
    //if(fabs(dis2-0)<PRECISION)return funcNam(dis1-prop);
    return funcNam(dis1 - prop) + funcNam(dis2 + prop);
}
void coreMolecularMC::updateOne(int i, double prop)
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
//essential function, which goes one step
bool coreMolecularMC::runOne()
{
    bool resFlag = false;
    for (int count = 0; count < (npart * 0.8); count++)
    {
        int ipart = Random_integer_uniform(0, npart);

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
int coreMolecularMC::run(int num)
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

///// static member for coreDoubleRMC
vector<double> coreDoubleRMC::targetX = {};
vector<double> coreDoubleRMC::targetY = {};
vector<double> coreDoubleRMC::potX = {};
vector<double> coreDoubleRMC::curinpY = {};
vector<double> coreDoubleRMC::propinpY = {};
coreMolecularMC *coreDoubleRMC::MMC = 0;
double coreDoubleRMC::dx = 0;
bool coreDoubleRMC::potFlag = 0;
int coreDoubleRMC::numX = 0;

double coreDoubleRMC::oldp = 0;
double coreDoubleRMC::alpha = 0;
int coreDoubleRMC::num_configs = 0;
int coreDoubleRMC::num_steps = 0;
int coreDoubleRMC::init_step = 0;
int coreDoubleRMC::num_particles = 0;
double coreDoubleRMC::move_scale = 0;
double coreDoubleRMC::move_cut = 0;

///// class function for coreDoubleRMC
void coreDoubleRMC::newMMC()
{
    MMC = new coreMolecularMC(&MMCPotential, 1);
    MMC->setBox(0, 800);
    //MMC = new coreMolecularMC(&potential,1);
    //for(int i=0;i<num_particles;i++)MMC->addPart(3*i+0.5*RNG_NAME());
    MMC->addAllPartEvenly(num_particles);
    //MMC->addAllPartByMinPot(num_particles,30,0.01);
    MMC->readyEverything();
    //cout<<"Initialization success with #particles:"<<(MMC->npart)<<endl;
}
void coreDoubleRMC::loadTarget(const char *name)
{
    FILE *fptr = fopen(name, "r");
    double x;
    double y;

    while (fscanf(fptr, "%lf %lf\n", &x, &y) == 2)
    {
        targetX.push_back(x);
        targetY.push_back(y);
    }
    dx = targetX[1] - targetX[0];
    //cout<<"load targetGr.gr ~success~"<<endl;
}
void coreDoubleRMC::initGusZero()
{
    numX = int(potX.size());
    for (int i = 0; i < numX; i++)
    {
        curinpY.push_back(0.0);
        propinpY.push_back(0.0);
    }
    if (potX.size() != curinpY.size())
        throw range_error("Error in initGusTarget!");
    if (potX.size() != propinpY.size())
        throw range_error("Error in initGusTarget!");
}
void coreDoubleRMC::initGusTarget()
{
    int cut = 0;
    vector<double> tmprdf;
    vector<double> otmprdf;
    for (double n : targetY)
        tmprdf.push_back(-n);
    for (double n : targetY)
        otmprdf.push_back(n + 10);

    vector<int> peaks = argpeakdet(tmprdf,
                                   0.3 * (*max_element(tmprdf.begin(), tmprdf.end()) - *min_element(tmprdf.begin(), tmprdf.end())));
    if (int(peaks.size()) >= 2)
        cut = peaks[1];

    numX = int(potX.size());
    double tmpPotY = 0.0; //as infinite far away

    curinpY.push_back(1e10);
    propinpY.push_back(1e10);
    for (int i = 1; i < numX - 3; i++)
    {
        if (i != 1 && cut != 0)
        {
            //using i-1 to avoid staying at the top, safe since here i != 0
            if (potX[i] > targetX[cut])
            {
            } //do nothing, not update tmpPotY
            //if(potX[i-1]>targetX[cut]){}//do nothing, not update tmpPotY
            else
            {
                tmpPotY = interpolate(targetX, otmprdf, potX[i], true);
            }
        }
        else
        {
            tmpPotY = interpolate(targetX, otmprdf, potX[i], true);
        }

        curinpY.push_back(-log(tmpPotY));
        propinpY.push_back(-log(tmpPotY));
    }
    curinpY.push_back(0.0);
    propinpY.push_back(0.0); //as infinite far away
    curinpY.push_back(0.0);
    propinpY.push_back(0.0); //as infinite far away
    curinpY.push_back(0.0);
    propinpY.push_back(0.0); //as infinite far away
    if (potX.size() != curinpY.size())
        throw range_error("Error in initGusTarget!");
    if (potX.size() != propinpY.size())
        throw range_error("Error in initGusTarget!");
}
double coreDoubleRMC::MMCPotential(double r)
{
    if (potFlag)
        return interpolate(potX, curinpY, r, true);
    else
        return interpolate(potX, propinpY, r, true);
}
void coreDoubleRMC::updatePotX(vector<double> arrX)
{
    vector<double> tmpPotY;
    numX = int(arrX.size());
    for (int i = 0; i < numX - 1; i++)
    {
        tmpPotY.push_back(interpolate(potX, curinpY, arrX[i], true));
    }
    tmpPotY.push_back(0.0); //as infinite far away
    potX = arrX;
    curinpY = tmpPotY;
    propinpY = tmpPotY;
    if (potX.size() != arrX.size())
        throw range_error("Error in updatePotX!");
    if (potX.size() != curinpY.size())
        throw range_error("Error in updatePotX!");
    if (potX.size() != propinpY.size())
        throw range_error("Error in updatePotX!");
}
//essential function, generate a random step
inline int coreDoubleRMC::proposal()
{
    int rnd_yindex = Random_integer_uniform(0, numX);
    double rnd_yval = RNG_NAME() * move_scale - move_cut;
    propinpY[rnd_yindex] += rnd_yval;
    //cout<<"ind:"<<rnd_yindex<<" val:"<<rnd_yval;
    return rnd_yindex;
}
inline void coreDoubleRMC::proposalGaussian()
{
    int numX = int(potX.size());

    double miu = (potX[numX - 4] - potX[1]) * RNG_NAME() + potX[1];
    double sigma = 3; //test
    //double sigma = ((potX[numX-4]-potX[1])-0.01)*RNG_NAME()+0.01;
    double A = RNG_NAME() * move_scale - move_cut;

    for (int i = 1; i < numX - 3; i++)
    {
        propinpY[i] += A * GaussianPDF(potX[i], miu, sigma);
    }
}
inline int coreDoubleRMC::proposalGlobal()
{
    double rnd_yval;
    double rnd_yval_all;
    int numX = int(potX.size());
    for (int i = 0; i < numX; i++)
    {
        rnd_yval = RNG_NAME() * move_scale - move_cut;
        propinpY[i] += rnd_yval;
        rnd_yval_all += rnd_yval;
    }
    return rnd_yval_all;
}
double coreDoubleRMC::getDiff(vector<double> gr)
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
inline bool coreDoubleRMC::acceptOrNot(double oldp, double newp)
{
    if (newp <= oldp)
        return true;
    if (RNG_NAME() < exp(alpha * (oldp - newp)))
        return true;
    return false;
}
//essential function, goes one step
bool coreDoubleRMC::runOneGlobal()
{
    //potFlag=true;
    //double oldp = getDiff(runMMC());

    //int rnd_yval_all = proposalGlobal();
    proposalGaussian();

    potFlag = false;
    double spacing = MMC->getPotMin(30, 0.01);
    if (spacing < PRECISION)
    { /*printf("trival step\n");*/
        return false;
    }

    double newp = getDiff(runMMC());
    potFlag = true;
    //cout<<" old_SSD:"<<oldp<<" new_SSD:"<<newp;
    bool resFlag = acceptOrNot(oldp, newp);
    //printf("%f  %f  %d \n",oldp,newp,resFlag);
    if (resFlag)
    {
        curinpY.clear();
        curinpY.assign(propinpY.begin(), propinpY.end());
        oldp = newp;
    }
    else
    {
        propinpY.clear();
        propinpY.assign(curinpY.begin(), curinpY.end());
    }
    //cout<<"  "<<resFlag<<endl;
    return resFlag;
}
bool coreDoubleRMC::runOne()
{
    //potFlag=true;
    //double oldp = getDiff(runMMC());

    int rnd_yindex = proposal();

    potFlag = false;
    double spacing = MMC->getPotMin(30, 0.01);
    if (spacing < PRECISION)
    { /*printf("trival step\n");*/
        return false;
    }

    double newp = getDiff(runMMC());
    potFlag = true;
    //cout<<" old_SSD:"<<oldp<<" new_SSD:"<<newp;
    bool resFlag = acceptOrNot(oldp, newp);
    if (resFlag)
        updateOne(rnd_yindex, newp);
    else
        propinpY[rnd_yindex] = curinpY[rnd_yindex];
    //cout<<"  "<<resFlag<<endl;
    return resFlag;
}
vector<double> coreDoubleRMC::runMMC()
{
    deleteMMC();
    newMMC();

    MMC->run(init_step);
    int ngrpoint = int(targetX.size());
    vector<double> gr = MMC->PARA_FUNC(ngrpoint, dx);
    vector<double> tgr;
    double countgr = 1;
    for (int i = 0; i < num_configs; i++)
    {
        MMC->run(num_steps);
        tgr = MMC->PARA_FUNC(ngrpoint, dx);
        VectorAddedByVector(gr, tgr);
        countgr++;
    }
    VectorDividedByScalar(gr, countgr);
    return gr;
}

vector<double> ClassErSearch::ersearchy = {};
vector<double> ClassErSearch::ersearchx = {};
vector<double> ClassErSearch::ersearchymid = {};
vector<double> ClassErSearch::ersearchymodel = {};
vector<double> ClassErSearch::vecval = {};
vector<double> ClassErSearch::best = {};
vector<double> ClassErSearch::bestpara = {};
double ClassErSearch::bestp = 1e99;
int ClassErSearch::vnum = -1;
double ClassErSearch::global_scaling = 0;
double ClassErSearch::global_boxl = 0;
double ClassErSearch::global_boxr = 0;
int ClassErSearch::global_npart = 0;
double ClassErSearch::global_minvol = 0;
double ClassErSearch::global_temperature = 0;
int ClassErSearch::global_num_configs = 0;
int ClassErSearch::global_num_steps = 0;
int ClassErSearch::global_init_step = 0;
//coreMolecularMC* ClassErSearch::MMC=0;

vector<double> coreSampleMCMC::targetX = {};
vector<double> coreSampleMCMC::targetY = {};
vector<double> coreSampleMCMC::oldpara = {};
vector<double> coreSampleMCMC::curpara = {};
coreMolecularMC *coreSampleMCMC::MMC = 0;
double coreSampleMCMC::dx = 0;
double coreSampleMCMC::oldp = 0;
double coreSampleMCMC::rndSize = 0;
double coreSampleMCMC::alpha = 0;

double coreSampleMCMC::global_LJpara = 0;
double coreSampleMCMC::global_scaling = 0;
double coreSampleMCMC::global_boxl = 0;
double coreSampleMCMC::global_boxr = 0;
int coreSampleMCMC::global_npart = 0;
double coreSampleMCMC::global_minvol = 0;
double coreSampleMCMC::global_temperature = 0;
int coreSampleMCMC::global_num_configs = 0;
int coreSampleMCMC::global_num_steps = 0;
int coreSampleMCMC::global_init_step = 0;

vector<double> coreDirectRun::targetX={};
vector<double> coreDirectRun::targetY={};
vector<double> coreDirectRun::curpara={};
coreMolecularMC *coreDirectRun::MMC=0;
double coreDirectRun::dx=0;

double coreDirectRun::global_scaling=0;
double coreDirectRun::global_boxl=0;
double coreDirectRun::global_boxr=0;
int coreDirectRun::global_npart=0;
double coreDirectRun::global_minvol=0;
double coreDirectRun::global_temperature=0;
int coreDirectRun::global_num_configs=0;
int coreDirectRun::global_num_steps=0;
int coreDirectRun::global_init_step=0;


/*extern "C" { 
    //type_pot getPotential(){return &potential;}
    
    coreMolecularMC* coreMolecularMC_new(double talpha,double tminvol=0.001)
    { return new coreMolecularMC(&potential,talpha,tminvol); } 
    
    void coreMolecularMC_addPart(coreMolecularMC* MMC,double pos)
    {return MMC ->addPart(pos); }
    int coreMolecularMC_getNpart(coreMolecularMC* MMC)
    {return MMC ->getNpart();}
    double coreMolecularMC_getPart(coreMolecularMC* MMC,int i)
    {return MMC ->getPart(i);}
    vector<double> coreMolecularMC_getPartAll(coreMolecularMC* MMC)
    {return MMC ->getPartAll();}
    void coreMolecularMC_printfDisAll(coreMolecularMC* MMC)
    {return MMC ->printfDisAll();}
    void coreMolecularMC_printfPart(coreMolecularMC* MMC)
    {return MMC ->printfPart();}
    vector<double> coreMolecularMC_getGr(coreMolecularMC* MMC,int ngr=160,double dx=1)
    {return MMC ->getGr(ngr,dx);}
    bool coreMolecularMC_readyEverything(coreMolecularMC* MMC)
    {return MMC ->readyEverything();}
    void coreMolecularMC_run(coreMolecularMC* MMC,int num)
    {return MMC ->run(num);}
}*/

/*coreDoubleRMC* rmc;//tbd global
double potentialDoubleRMC(double r)
{
    return rmc->MMCPotential(r);
}*/

/**************  random seed  *******************************/

/*int main() 
{    
SeedByTime();

	coreMolecularMC MMC=coreMolecularMC(&potential,1);
    int num = 20;
    for(int i=0;i<num;i++)MMC.addPart(3*i+0.5*RNG_NAME());
    MMC.readyEverything();
    
    int num_configs=100;
    int num_steps=200;
    int init_step = 10000;
    
    int npart = MMC.getNpart();
    int nsteps = 0;
    
    std::cout << "particle\n";
    MMC.printfPart();
    std::cout << "distance\n";
    MMC.printfDisAll();
    
    MMC.run(init_step);
    nsteps+=init_step;
    
    vector<double> gr=MMC.getGr(160);
    printfVector(gr);
    
    FILE *fptr;
    fptr = fopen("test.out", "w");
    for(int i=0;i<num_configs;i++)
    {
        MMC.run(num_steps);
        nsteps+=num_steps;
        VectorAddedByVector(gr,MMC.getGr(160));
        
        fprintf(fptr, "ITEM: TIMESTEP\n%d\n",nsteps);
        fprintf(fptr, "ITEM: NUMBER OF ATOMS\n%d\n",npart);
        fprintf(fptr, "ITEM: BOX BOUNDS pp pp pp\n");
        fprintf(fptr, "0 1\n");
        fprintf(fptr, "0 1\n");
        fprintf(fptr, "0 1\n");
        fprintf(fptr, "ITEM: ATOMS id xs ys zs\n");
        for(int i=0;i<npart;i++){
            fprintf(fptr, "%d %.5f %.5f %.5f\n",i, MMC.getPart(i), 0.0, 0.0);}
    }
    fclose(fptr);
    VectorDividedByScalar(gr,num_configs);
    write2DVector("savegr",arange(0,160),gr);

    
    
    
    
    typedef coreDoubleRMC rmc; 
    
    vector<double> arrX = {0,10,20,30};arrX.push_back(100);
    //coreDoubleRMC* rmc = new coreDoubleRMC(arrX);
    rmc::initPotentialX(arrX);
    rmc::asignConfig(0.5,5,200,1000,20,20,10);
    rmc::readyToRun();
    
    rmc::run(400);
    
    arrX = linspace(0,30,10);arrX.push_back(100);
    rmc::updatePotX(arrX);
    rmc::run(400);
    
    arrX = linspace(0,30,31);arrX.push_back(100);
    rmc::updatePotX(arrX);
    rmc::run(200);
    
    rmc::writePot();
    
    return 0; 
} */

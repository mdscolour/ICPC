#PBS -l nodes=1:ppn=1:medium_buster
#PBS -q medium_buster
#PBS -e ncxxx.errout
#PBS -o ncxxx.resout
#PBS -l walltime=11:00:00
#PBS -l mem=500mb,vmem=500mb

echo "START_TIME           = `date +'%y-%m-%d %H:%M:%S %s'`"

tdr=SAM2B_program3_xxx_yyy  # at least one xxx remain
underWS=/remote/pi310b/li/MolecularMC/LJTestInt/reverseMC/

cptarget1=${underWS}/*.cpp
cptarget2=${underWS}/*.h
cptarget3=${underWS}/*config
cptarget4=${underWS}/*gr

cpback1=${underWS}/chr2-0LJSAM2B/
#tdr=xxx_yyy_prolong

# in case intermediate interruption
cd /tmp
if [ -d li_${tdr} ]
then
    rm -rf li_${tdr}
fi
mkdir li_${tdr}

if [ -f li_${tdr}.tar.gz ]
then
    rm li_${tdr}.tar.gz
fi
#cd ../
#cd /tmp/chu_${tdr}
cd li_${tdr}

##### start command
cp ${cptarget1} .
cp ${cptarget2} .
cp ${cptarget3} .
cp ${cptarget4} .

g++ runSAM2.cpp
./a.out 0 xxx

##### data transfer back
#cp *.like_bed $cpback1
#cp *_ctcfbin $cpback1
#cp *_ctcfbin47* $cpback1
cp *pot $cpback1
#cp alpha* $cpback1
#cp _chrxxx.like_wig $cpback1
#cp chrxxx $underWS
#cp ../*.out $cpback1

cd ../
#tar -czvf li_${tdr}.tar.gz li_${tdr}/*
rm -rf ./li_${tdr}
#cp ./li_${tdr}.tar.gz $underWS
#rm ./li_${tdr}.tar.gz
### cd out, important in some cases
cd ../


echo "END_TIME             = `date +'%y-%m-%d %H:%M:%S %s'`"

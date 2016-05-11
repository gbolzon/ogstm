#! /bin/bash
OGSTM_ARCH=x86_64
OGSTM_OS=LINUX
OGSTM_COMPILER=intel
DEBUG=
DEBUG=.dbg
#DEBUG=   # for production

export OPENMP_FLAG=     #-fopenmp
export MODULEFILE=machine_modules/pico.intel

OCEANVAR=true
DEBUG_OCEANVAR=

usage() {
echo "SYNOPSYS"
echo "Build BFM and ogstm model"
echo "ogstm_builder.sh [ BFMDIR ] [ OGSTMDIR ]"
echo ""
echo " Dirs have to be expressed as full paths "
echo "EXAMPLE"
echo " ./ogstm_builder.sh $PWD/bfm $PWD/ogstm "

}

if [ $# -lt 2 ] ; then
   usage
   exit 1
fi

BFMDIR=$1
OGSTMDIR=$2


############### MODULES AND ENVIRONMENT VARIABLES
source $MODULEFILE


# -------------- 3d_var _____
if [ $OCEANVAR == true ] ; then
   cd 3d_var
   INC_FILE=${OGSTM_ARCH}.${OGSTM_OS}.${OGSTM_COMPILER}${DEBUG_OCEANVAR}.inc
   cp $INC_FILE compiler.inc
   gmake
   if [ $? -ne 0 ] ; then  echo  ERROR; exit 1 ; fi
   export DA_INC=$PWD
   echo DA_INC= $DA_INC
fi

# ----------- BFM library ---------------------
cd $BFMDIR
A=`svn info 2>&1 `  # exit status 1 if bfmv5
if [ $? == 0 ] ; then 
   BFMversion=BFMv2
else
   BFMversion=bfmv5
fi

INC_FILE=${OGSTM_ARCH}.${OGSTM_OS}.${OGSTM_COMPILER}${DEBUG}.inc

if [ $BFMversion == BFMv2 ] ; then
   cd ${BFMDIR}/compilers
   cp $INC_FILE compiler.inc

   ###################################################### just because R1.3 does not have include/
   mkdir -p  ${BFMDIR}/include

   cd ${BFMDIR}/build
   ./config_BFM.sh -a ${OGSTM_ARCH} -c ogstm
   cd BLD_OGSTMBFM
   gmake

else
   # in-place replace the entire ARCH line
   sed -i "s/.*ARCH.*/        ARCH    = '$INC_FILE'  /"  build/configurations/OGS_PELAGIC/configuration
   cd $BFMDIR/build
   ./bfm_configure.sh -gc -o ../lib/libbfm.a -p OGS_PELAGIC
   if [ $? -ne 0 ] ; then  echo  ERROR; exit 1 ; fi
fi

export BFM_INC=${BFMDIR}/include
export BFM_LIB=${BFMDIR}/lib




# ------------ OGSTM builder -----------------

cd ${OGSTMDIR}/compilers
INC_FILE=${OGSTM_ARCH}.${OGSTM_OS}.${OGSTM_COMPILER}${DEBUG}.inc
cp $INC_FILE compiler.inc

cd ${OGSTMDIR}/build
./config_OGSTM.sh ${OGSTM_ARCH} 
cd BLD_OGSTM
make -f MakeLib
rm -f get_mem_mod.o
gmake

if [ $? -ne 0 ] ; then  echo  ERROR; exit 1 ; fi

### OGSTM NAMELIST GENERATION (also by Frequency Control )

mkdir -p ${OGSTMDIR}/ready_for_model_namelists/

if [ $BFMversion == bfmv5 ] ; then
   cp ${BFMDIR}/build/tmp/OGS_PELAGIC/namelist.passivetrc ${OGSTMDIR}/bfmv5/
   cd ${OGSTMDIR}/bfmv5/
   python ogstm_namelist_gen.py #generates namelist.passivetrc_new

   cp ${OGSTMDIR}/src/namelists/namelist*    ${OGSTMDIR}/ready_for_model_namelists/ 
   cp namelist.passivetrc_new                ${OGSTMDIR}/ready_for_model_namelists/namelist.passivetrc #overwriting
   cp ${BFMDIR}/build/tmp/OGS_PELAGIC/*.nml  ${OGSTMDIR}/ready_for_model_namelists/
else
   #V2
   cp ${OGSTMDIR}/src/namelists/namelist*    ${OGSTMDIR}/ready_for_model_namelists/
   cp ${BFMDIR}/src/namelist/*.nml           ${OGSTMDIR}/ready_for_model_namelists/
fi

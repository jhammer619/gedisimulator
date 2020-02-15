REM # Download/Install miniconda 
call conda create -n -y gedisimulator
call conda activate gedisimulator

pushd %CONDA_PREFIX%
mkdir src
cd src

call conda install -y git mercurial
call conda install -y geotiff hdf5 gdal 
call conda install -y -c statiskit clang 
call conda install -y -c conda-forge make gsl

REM install visual_studio_sdk

git clone https://bitbucket.org/caiohamamura/gedisimulator.git
hg clone https://bitbucket.org/StevenHancock/tools
hg clone https://bitbucket.org/StevenHancock/libclidar

call conda install -y -c menpo wget

wget --no-check-certificate https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz
tar -xf cmpfit-1.2.tar.gz


cd gedisimulator
SET HDF5_LIB=../../Library/lib
SET HANCOCKTOOLS_ROOT=../tools
SET CMPFIT_ROOT=../cmpfit-1.2
SET LIBCLIDAR_ROOT=../libclidar
SET GEDIRAT_ROOT=.
SET GSL_ROOT=../../Library/include
SET CC=clang
call "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
SET LIB=%LIB%;../../Library/lib


FOR /F "tokens=* USEBACKQ" %%g IN (`where link`) do (SET "linker_path=%%g")
mkdir %temp_path%\bin\HostX64\x64
copy "%linker_path%" .
set VCToolsInstallDir=%temp_path%
set PATH=%path%;%linker_path%

make SHELL=cmd CC=clang CFLAGS="-I../../Library/include"
make THIS=gediMetric SHELL=cmd CC=clang
make THIS=lvisBullseye SHELL=cmd CC=clang
make THIS=mapLidar SHELL=cmd CC=clang
make THIS=lasPoints SHELL=cmd CC=clang
make THIS=gediMetric SHELL=cmd CC=clang install
make THIS=lvisBullseye SHELL=cmd CC=clang install
make THIS=mapLidar SHELL=cmd CC=clang install
make THIS=lasPoints SHELL=cmd CC=clang install
# README #

### What is this repository for? ###

This is a set of programs to simulate large-footprint full-waveform lidar from airborne lidar and to process it and perform various other tasks. The parts are:


**gediRat**: simulates GEDI waveforms from ALS .las files and outputs ASCII or HDF5 waveforms.

**gediMetric**: processes full-waveform data (LVIS or simulated GEDI) and outputs metrics.

**mapLidar**: produces geotiffs from .las files of different properties. Also produces the bounding boxes of .las files for use in overlapLasFiles.csh.

**lasPoints**: outputs .pts files from .las files for selected areas.

**lvisBullseye**: produces the correlation bullseye plots of ALS to full-waveform data following Blair and Hofton (1999).

**addNoiseHDF**: Reads waveform data from HDF5 files and adds noise of a chosen level.


To find options, type the above command with "-help". These are explained in full detail later. 


There are some C and bash shells to control the above:

**gediRatList.csh**: batch processes gediRat.

**filtForR.csh**: converts the output .txt files into .csv files for reading in to R.

**overlapLasFiles.cs**h: determines which .las files are needed for a given simulation.

**orbitTracks.bash**: produces lists of footprints from GEDI orbital simulations.

**listALS.csh**: produces the lists of .las files needed to read multiple files.

The other .c files are either small test programs in the development of the above or are libraries called by the above. Important libraries are:

**gediIO.c**: contains the GEDI simulator functions.

### How do I get set up? ###

Clone the repository and point to its location with the environment variable:

    GEDIRAT_ROOT


All programs depend on these libraries:

* Gnu Scientific Library
* Geotiff
* HDF5
* GDAL

Minpack (Levenberg-Maruqardt):  https://www.physics.wisc.edu/~craigm/idl/cmpfit.html


Point to these with the environment variable:

    GSL_ROOT
    HDF5_LIB
    CMPFIT_ROOT

Once they're installed, it also requires two libraries without package managers:

[**C-tools**](https://bitbucket.org/StevenHancock/tools)

[**libClidar**](https://bitbucket.org/StevenHancock/libclidar)

Clone these from bitbucket and point to their locations with the environment variables:

    HANCOCKTOOLS_ROOT
    LIBCLIDAR_ROOT

To compile, type:

  **make THIS=gediRat**

Make sure that ~/bin/$ARCH exists and is in the path, where $ARCH is the result of `uname -m`.

  **make THIS=gediRat install**

Replace "**gediRat**" with each of the commands above to compile and install.


Make sure that all **.csh** and **.bash** files are also in your path.


# Function operation #

## gediRat ##

Program to create GEDI waveforms from ALS las or pts files. laz not yet supported. Data is output either as ASCII files or as a HDF5 file, both of which can be ready by gediMetric below.


##### Input output filenames and format
    -input name;     lasfile input filename
    -inList list;    input file list (ASCII file) for multiple files
    -output name;    output filename
    -ground;         record separate ground and canopy waveforms
    -hdf;            write output as HDF5. Best with gridded or list of coords
    -ascii;          write output as ASCII (default). Good for quick tests
    -waveID id;      supply a waveID to pass to the output (only for single footprints)

##### Single footprint, list of footprints, or grid of footprints
    -coord lon lat;  footprint coordinate in same system as lasfile
    -listCoord name; list of coordinates
    -gridBound minX maxX minY maxY;     make a grid of waveforms in this box
    -gridStep res;   grid step size

##### Lidar characteristics. Defaults are expected GEDI values.
    -pSigma sig;     set Gaussian pulse width as 1 sigma
    -pFWHM fhwm;     set Gaussian pulse width as FWHM in ns
    -readPulse file; read pulse shape and width from a file insteda of making Gaussian
    -fSigma sig;     set footprint width
    -wavefront file; read wavefront shape from file instead of setting Gaussian. Note that footprint width is still set by fSigma
    -res res;        range resolution of waveform digitisation to output, in units of ALS data
    -LVIS;           use LVIS pulse length, sigma=6.25m
    -topHat;         use a top hat wavefront
    -sideLobe;       use side lobes
    -lobeAng ang;    lobe axis azimuth


##### Input data quality filters
    -checkCover;     check that the footprint is covered by ALS data. Do not output if not
    -maxScanAng ang; maximum scan angle, degrees
    -decimate x;     probability of accepting an ALS beam


##### Computational speed options
    -pBuff s;        point reading buffer size in Gbytes
    -maxBins;        Optional: for HDF5, limit number of bins to save trimming.
    -countOnly;      only use count method
    -pulseAfter;     apply the pulse smoothing after binning for computational speed, at the risk of aliasing (default)
    -pulseBefore;    apply the pulse smoothing before binning to avoid the risk of aliasing, at the expense of computational speed
    -noNorm;         don't normalise for ALS density

##### Octree
    -noOctree;       do not use an octree
    -octLevels n;    number of octree levels to use
    -nOctPix n;      number of octree pixels along a side for the top level


##### Using full-waveform input data (not tested)
    -decon;          deconvolve
    -indDecon;       deconvolve individual beams
    -readWave;       read full-waveform where available


##### Miscellaneous
    -listFiles;      list files. Do not read them
    -keepOld;        do not overwrite old files, if they exist
    -useShadow;      account for shadowing in discrete return data through voxelisation
    -polyGround;     find mean ground elevation and slope through fitting a polynomial
    -nnGround;       find mean ground elevation and slope through nearest neighbour




## gediMetric ##

Program to process large-footprint lidar data (real or simulated) and produce standard waveform metrics. It can add noise to simulations and alter pulse shapes (increase length only). It reads either ASCII or HDF5 files created by gediRat, or can read LVIS data in either HDF5 or .lgw format. It will be updated to read GEDI data when that is available. Take care when reading ASCII data as some options are mutually exclusive (different gediRat options can change the column order). This outputs an ASCII file with the first row defining the contenst of each column. Output variable names are defined below.

##### Input output
    -input name;      waveform  input filename
    -outRoot name;    output filename root
    -inList list;     input file list for multiple files
    -writeFit;        write fitted waveform
    -writeGauss;      write Gaussian parameters
    -readBinLVIS;     input is an LVIS binary file
    -readHDFlvis;     read LVIS HDF5 input
    -readHDFgedi;     read GEDI simulator HDF5 input
    -level2 name;     level2 filename for LVIS ZG
    -bounds minX minY maxX maxY;    only analyse data within bounds

##### Switches
    -ground;          read true ground from file
    -useInt;          use discrete intensity instead of count
    -useFrac;         use fractional hits rather than counts
    -rhRes r;         percentage energy resolution of RH metrics
    -laiRes res;      lai profile resolution in metres
    -laiH h;          height to calculate LAI to
    -noRHgauss;       do not fit Gaussians
    -gTol tol;        ALS ground tolerance. Used to calculate slope.
    -fhdHistRes res;  waveform intensity resolution to use when calculating FHD from histograms
    -forcePsigma;     do not read pulse sigma from file
    -bayesGround;     use Bayseian ground finding
    -dontTrustGround; don't trust ground in waveforms, if included
    -noRoundCoord;    do not round up coords when outputting

##### Adding noise:
    -dcBias n;        mean noise level
    -nSig sig;        noise sigma
    -seed n;          random number seed
    -hNoise n;        hard threshold noise as a fraction of integral
    -linkNoise linkM cov;     apply Gaussian noise based on link margin at a cover
    -linkFsig sig;    footprint width to use when calculating and applying signal noise
    -linkPsig sig;    pulse width to use when calculating and applying signal noise
    -trueSig sig;     true sigma of background noise
    -bitRate n;       digitisation bit rate
    -maxDN max;       maximum DN
    -renoise;         remove noise from truth before applying new noise level
    -newPsig sig;     new value for pulse width, when lengthening pulse
    -oldPsig sig;     old value for pulse width if not defined in waveform file, when lengthening pulse
    -addDrift xi;     apply detector background drift
    -missGround;      assume ground is missed to assess RH metrics
    -minGap gap;      delete signal beneath min detectable gap fraction

##### Photon counting
    -photonCount;     output point cloud from photon counting
    -nPhotons n;      mean number of photons
    -photonWind x;    window length for photon counting search, metres
    -noiseMult x;     noise multiplier for photon-counting

##### Denoising:
    -meanN n;         mean noise level, if using a predefined mean level
    -thresh n;        noise threshold, if using a predefined noise threshold
    -varNoise;        use a variable noise threshold
    -varScale x;      variable noise threshold scale (multiple of stdev above mean to set threshold)
    -statsLen len;    length to calculate noise stats over for varNoise
    -noiseTrack;      use noise tracking
    -sWidth sig;      smoothing width, after densoising
    -psWidth sigma;   smoothing width, before denoising
    -msWidth sig;     smoothing width, after noise stats, before denoising
    -preMatchF;       matched filter before denoising
    -postMatchF;      matched filter after denoising
    -pFile file;      read pulse file, for deconvoltuion and matched filters
    -gWidth sig;      Gaussian parameter selection smoothing width
    -minGsig sig;     minimum Gaussian sigma to fit
    -minWidth n;      minimum feature width in bins
    -medNoise;        use median stats rather than mean
    -varDrift;        correct detector drift with variable factor
    -driftFac xi;     fix drift with constant drift factor
    -rhoG rho;        ground reflectance
    -rhoC rho;        canopy reflectance
    -pSigma sig;      pulse width to smooth by if using Gaussian pulse
    -gold;            deconvolve with Gold's method
    -deconTol;        deconvolution tolerance

### gediMetric variable names ###
Note that some metrics are "true" and will not be available to GEDI. They are included to assess errors and sensitivities.


##### Metrics available to GEDI
    gHeight - ground elevation (m) from Gaussian fitting
    maxGround - ground elevation (m) from lowest maximum
    inflGround - ground elevation (m) from inflection points.
    signal top - elevation of first point above noise (may include noise tracking).
    signal bottom - elevation of last return above noise (may include noise tracking).
    cover - canopy cover (fraction) from area of Gaussian fitted ground. Uses rho_v=0.57 and rho_g=0.4.
    leading edge ext - leading edge extent (m), from Lefksy et al (2007).
    trailing edge extent - trailing edge extent (m), from Lefksy et al (2007).
    rhGauss 0-100 - RH metrics, 0%-100%, using ground from Gaussian fitting (m).
    rhMax 0-100 - RH metrics, 0%-100%, using ground from lowest maximum (m).
    rhInfl 0-100 - RH metrics, 0%-100%, using ground from inflection points (m).
    gaussHalfCov - canopy cover (fraction) from double the energy beneath the Gaussian ground. Uses rho_v=0.57 and rho_g=0.4.
    maxHalfCov - canopy cover (fraction) from double the energy beneath the lowest maximum ground. Uses rho_v=0.57 and rho_g=0.4.
    infHalfCov - canopy cover (fraction) from double the energy beneath the inflection point ground. Uses rho_v=0.57 and rho_g=0.4.
    bayHalfCov - canopy cover (fraction) from double the energy beneath the experimental "Bayesian" ground. Uses rho_v=0.57 and rho_g=0.4.
    lon - footprint centre longitude in projection of ALS data (m).
    lat - footprint centre latitude in projection of ALS data (m).
    waveEnergy - total energy within waveform (will be 1 scaled by noise for simulations).
    blairSense - Blair's sensitivity metric. Canopy cover at which this SNR would have90% chance of detecting ground (does not account for rho_v/rho_g).
    FHD - Foliage height diversity
    niM2 - Wenge Ni's biomass metric, equal to the sum of the RH metrics to the power of 2 (unpublished)
    niM2.1 - Wenge Ni's biomass metric, equal to the sum of the RH metrics to the power of 2.1 (unpublished)

##### Metrics unavailable to GEDI
    wave ID - waveform label, relates to plot name and footprint number.
    true ground - ground elevation (m) from ALS. Centre of gravity of ground points within footprint
    true top - elevation of highest point of waveform (m), without noise. Includes pulse blurring.
    ground slope - effective ground slope (degrees), from width of ground return. Includes roughness.
    ALS cover - canopy cover (fraction) from ALS data. Uses rho_v=0.57 and rho_g=0.4.
    rhReal 0-100 - RH metrics, 0%-100%, using "true" ground from ALS data (m).
    groundOverlap - fraction of ground return overlapping with canopy return. A measure of understorey.
    groundMin - depth of minimum between ground and canopy return. A measure of understorey.
    groundInfl - d2y/dx2 of inflection point between ground and canopy return. A measure of understorey.
    pointDense - average ALS point density within GEDI footprint.
    beamDense - average ALS beam density within GEDI footprint.

##### System settings
    pSigma - GEDI system pulse width, sigma (m).
    fSigma - GEDI footprint width, sigma (m).
    linkM - link margin if noise is added (db).
    linkCov - canopy cover at which the above link margin is true (fraction).
    filename - name of input waveform file.




##### Signal processing description


##### Gaussian fitting
Used for "gHeight", "rhGauss" and "gaussHalfCov". The waveform is denoised (mean+5*sigma, noise tracking to avoid truncation), smoothed (pSigma*0.75) and Gaussians fitted with Levenberg-Marquardt optimisation. The center of the lowest Gaussian containing at least 0.5% of the waveform energy is selected as the ground.

##### Maximum
Used for "maxGround", "rhMax" and "maxHalfCov". The waveform is denoised (mean+5*sigma, noise tracking to avoid truncation), smoothed (pSigma*0.75). The lowest maximum is taken as the ground.

##### Inflection points
Used for "inflGround", "rhInfl" and "inflHalfCov". The waveform is denoised (mean+5*sigma, noise tracking to avoid truncation), smoothed (pSigma*0.75). The centre of gravity between the lowest two inflection points is taken as the ground.

##### Half covers
Used for "gaussHalfCov", "maxHalfCov" and "inflHalfCov". Sum energy beneath estimated ground position. Double that is the ground energy. Calculate canopy cover, correcting for rho_v and rho_g.

cover=Ecan/(Ecan+Eg*rho_v/rho_g)

Where Ecan is the canopy energy, Eg is the ground energy, rho_v is the vegetation reflectance and rho_g is the ground reflectance.


##### Edge extents
These are described in:

Lefsky, Michael A., Michael Keller, Yong Pang, Plinio B. De Camargo, and Maria O. Hunter. "Revised method for forest canopy height estimation from Geoscience Laser Altimeter System waveforms." Journal of Applied Remote Sensing 1, no. 1 (2007): 013537-013537.



## mapLidar ##
Generates a geotiff from las file properties, combining multiple files.

## lasPoints ##
Extracts a point cloud as a pts for a bounding box within a collection of las files.


### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

svenhancock@gmail.com

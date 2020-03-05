## MAF install notes ##

2020-03-04 WIC - I used these steps to successfully install MAF **from
source** on an HP Z840 running Ubuntu 18.0.4 with **bash** as the
default shell. The installation steps in the MAF documentation **do**
work - and these steps are lightly edited version of that material --
but are somewhat sensitive to environment variables. These notes are a
record of the steps that worked for me.

## Per-session Incantations to set up once installed ##

source /raid1/soft/lsst_stack/loadLSST.bash  
setup sims_maf -t $USER  

And, if using the sims_maf_contrib:

setup sims_maf_contrib -t $USER -t sims  

## Links to instructions ##

Pre-requisites: https://pipelines.lsst.io/install/newinstall.html

Installation steps: https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF

sims_maf github page: https://github.com/lsst/sims_maf

Including maf_contrib: https://github.com/LSST-nonproject/sims_maf_contrib 

## Steps ##

### 0. Set up a dedicated user account ###

I wasn't able to install from source on my regular user account
without conflicts, probably due to the various customizations and
environment variables I use. I also found that the anaconda python on
my regular account was causing conflicts with sims_maf. I created a
new user account for sims_maf and ran the steps logged in as that
account.  

Simply switching to that account by **su [my-lsst-account]** failed to
resolve the conflicts: I had to log out completely, log back in to the
fresh account, and perform the following steps. When logged in to the
fresh account, however, everything ran pretty much as advertised on
the installation documentation.

### 1. Install the pre-requisites ###

From here: https://pipelines.lsst.io/install/newinstall.html

sudo apt-get install \
    bison \
    ca-certificates \
    cmake \
    curl \
    default-jre \
    flex \
    gettext \
    git \
    libbz2-dev \
    libcurl4-openssl-dev \
    libfontconfig1 \
    libglib2.0-dev \
    libncurses5-dev \
    libreadline6-dev \
    libx11-dev \
    libxrender1 \
    libxt-dev \
    m4 \
    make \
    perl-modules \
    rsync \
    zlib1g-dev

### 2. Create a directory to hold the lsst stack and install it. ###

From here: https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF

On my system, this all lives on a large data disk (so that all the
necessary maps have plenty of space).

mkdir /raid1/soft/lsst_stack  
cd /raid1/soft/lsst_stack  

cd ~/lsst  
curl -OL https://raw.githubusercontent.com/lsst/lsst/master/scripts/newinstall.sh  
bash newinstall.sh -ct  

source ./loadLSST.bash  
eups distrib install lsst_sims -t **sims_w_2020_05**  
curl -sSL https://raw.githubusercontent.com/lsst/shebangtron/master/shebangtron | python  
 
* The "sims_w_2020_05" is the weekly tag I most recently installed (from 2020 Feb 1). When using a different weekly "tag", replace this with the version you actually installed. The tags can be found here: https://github.com/lsst/afw/releases

* As indicated on the confluence and maf page, this took about an hour to install from source, with about 150 packages installed. In my fresh lsst user account, this ran all the way through without problems. 

### 3. Set up to run sims_maf ###

setup sims_maf -t **sims_w_2020_05**

### 4. Bring down the new sims_maf from github and incorporate it ###

mkdir /raid1/soft/maf_github  
cd /raid1/soft/maf_github  
git clone git@github.com:lsst/sims_maf.git  

cd sims_maf  
eups declare -r . -t $USER  
setup sims_maf -t $USER  
scons  

### 5. Bring down maf_contrib and incorporate it ###

From here: https://github.com/LSST-nonproject/sims_maf_contrib

mkdir /raid1/soft/mafcontrib_github  
cd /raid1/soft/mafcontrib_github  
git clone  git@github.com:LSST-nonproject/sims_maf_contrib.git  
cd sims_maf_contrib  
eups declare sims_maf_contrib -r . -t $USER  
setup sims_maf_contrib -t $USER -t sims  


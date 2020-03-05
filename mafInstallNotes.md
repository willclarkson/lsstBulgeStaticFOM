## MAF install notes ##

2020-03-04 WIC - I used these steps to successfully install MAF **from
source** on an HP Z840 running Ubuntu 18.0.4 with **bash** as the
default shell. The installation steps in the MAF documentation **do**
work - and these steps are taken from the links below -- but are
somewhat sensitive to environment variables. These notes are a static
record of the steps that worked on my system.

## Incantations to set up sims_maf once installed ##

The steps below must be done for each new session (after the one in which the software was installed).
```
source /raid1/soft/lsst_stack/loadLSST.bash  
setup sims_maf -t $USER  
```

And, if using the sims_maf_contrib:

```
setup sims_maf_contrib -t $USER -t sims  
```

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
my regular account was causing conflicts with sims_maf. 

As suggested by the documentation, I created a fresh user account for
use with sims_maf, and ran the steps below logged in to that account.

Simply switching in the terminal from my personal account to that
account in the terminal (by **su [my-lsst-account]**) failed to
resolve the conflicts: I had to log out completely, log back in to the
fresh account, and perform the steps below (this also applies to the
per-session install steps above). When logged in to the fresh account,
however, everything ran pretty much as advertised on the installation
documentation.

### 1. Install the pre-requisites ###

From here: https://pipelines.lsst.io/install/newinstall.html

```
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
```

### 2. Create a directory to hold the lsst stack and install it. ###

From here: https://confluence.lsstcorp.org/display/SIM/Catalogs+and+MAF

On my system, this all resides on a large data disk (so that all the
necessary maps have plenty of space). To configure the installation
for my system:

```
mkdir /raid1/soft/lsst_stack  
cd /raid1/soft/lsst_stack  

cd ~/lsst  
curl -OL https://raw.githubusercontent.com/lsst/lsst/master/scripts/newinstall.sh  
bash newinstall.sh -ct  
```

And then to install:

```
source ./loadLSST.bash  
eups distrib install lsst_sims -t sims_w_2020_06  
curl -sSL https://raw.githubusercontent.com/lsst/shebangtron/master/shebangtron | python  
```

Notes:   
* The "sims_w_2020_06" is the weekly tag I most recently installed (from 2020 Feb 07). When using a different weekly "tag", replace this with the version actually installed. The tags can be found at the following link, but the install tag seems to differ slightly from the tag at the link. For example, entry ''w.2020.06 '' would be ''sims_w_2020_06'' (i.e. prepend with 'sims_' and replace '.' with '_'): https://github.com/lsst/afw/releases

* As suggested by the documentation, this took about an hour to install from source, with more than 120 packages installed. In my fresh lsst user account, this ran all the way through without problems. 

To set up to run:
```
setup sims_maf -t sims_w_2020_06
```

### 3. Bring down the new sims_maf from github and incorporate it ###

```
mkdir /raid1/soft/maf_github  
cd /raid1/soft/maf_github  
git clone git@github.com:lsst/sims_maf.git  

cd sims_maf  
eups declare -r . -t $USER  
setup sims_maf -t $USER  
scons  
```

### 4. Bring down maf_contrib and incorporate it ###

From here: https://github.com/LSST-nonproject/sims_maf_contrib

```
mkdir /raid1/soft/mafcontrib_github  
cd /raid1/soft/mafcontrib_github  
git clone  git@github.com:LSST-nonproject/sims_maf_contrib.git  

cd sims_maf_contrib  
eups declare sims_maf_contrib -r . -t $USER  
setup sims_maf_contrib -t $USER -t sims  
```

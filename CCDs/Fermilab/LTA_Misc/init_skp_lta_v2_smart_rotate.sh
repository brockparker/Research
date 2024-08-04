#!/bin/bash
source setup_lta.sh

seqRead=/home/astro/soft/ltaDaemon-javier/sequencers/sequencer_rotate.xml
volReadDefault=/home/astro/soft/ltaDaemon-javier/voltage_skp_lta_v2.sh


lta nomulti
lta sseq $seqRead

lta set sinit 30
lta set pinit 0
lta set ssamp 200
lta set psamp 200
lta set packSource 9
lta set cdsout 2 #sig-ped

# Physical CCD dimensions
lta CCDNROW 1246
lta CCDNCOL 742

# CCD number of rows
lta NROW 650
lta NCOL 450
lta NSAMP 1

if [ -n "$volRead" ]
then
    source $volRead
else
    source $volReadDefault
fi


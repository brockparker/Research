#!/bin/bash

daemonDir=/home/astro/soft/ltaDaemon-javier
seqClean=${daemonDir}/sequencers/sequencer_manu_fastclear_v4.xml
seqRead=${daemonDir}/sequencers/KWL_test_sequencer_rotate.xml #sequencer_rotate.xml

clean() {
# NROW should be a bit larger than for a normal readout (around +50/100)
# NCOL should be at least 1.5 times the usual
# check variables in sequencer file or set here!

        ./lta.sh 8888 sseq $seqClean
        source clean_erase_purge.sh
        ./lta.sh 8888 runseq $seqClean
        ./lta.sh 8888 sseq $seqRead
}
clean

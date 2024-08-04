#!/bin/bash


seqClean=/home/astro/soft/ltaDaemon-javier/sequencers/sequencer_manu_fastclear_v4.xml
seqRead=/home/astro/soft/ltaDaemon-javier/sequencers/sequencer_rotate.xml

clean() {

	local eap=${1:-yes}
	lta sseq $seqClean
	#if [ $eap == "yes" ]; then
		#echo "erase and purge..."
		#source eraseANDepurge.sh
	#fi
	echo "cleaning..."
	lta runseq $seqClean
	lta sseq $seqRead

	lta NCOL 450
	lta NROW 650
	lta NSAMP 1
}

clean "$1"

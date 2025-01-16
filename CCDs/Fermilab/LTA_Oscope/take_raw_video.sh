#!/bin/bash
# Takes raw video data 

amp=(0 1 2 3) # channels to take raw video data from

for vid in ${amp[@]}

do
    lta startseq
    sleep 0.35
    lta name ./raw_video/HDU_${vid}_
    lta set bufSel $vid
    lta set packSource 4
    lta set packStart 1
    lta getraw
    lta set packSource 9
    lta set packStart 0
done

python3 plot_raw_video.py 

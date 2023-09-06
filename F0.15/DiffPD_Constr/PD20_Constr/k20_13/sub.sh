#! /bin/bash

bsub  -J Conk50_9  -q kometshort -app Reserve1500M -W 300 -o outfile_Conk50_9 -e err_Conk50_9  ./out >&log_Conk50_9


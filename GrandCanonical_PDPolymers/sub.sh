#! /bin/bash

bsub  -J Test12  -q kometlong -app Reserve1500M -W 4000 -o outfile_Test12 -e err_Test12  ./out >&log_Test12


#!/bin/bash
g++ -w -O3 -o convert txt2bin.cc

# # ------------------------------------------------------------------------------
# #  Mnist 
# # ------------------------------------------------------------------------------
# dname=Mnist
# n=60000
# d=50
# qn=100
# dtype=uint8
# min=0
# max=255
# ipath=txt/${dname}/${dname}
# opath=${dname}/${dname}

# ./convert ${n} ${d} ${qn} ${dtype} ${ipath} ${opath} ${min} ${max}
# cp ${ipath}.stat ${opath}.stat 

# # ------------------------------------------------------------------------------
# #  Sift 
# # ------------------------------------------------------------------------------
# dname=Sift
# n=1000000
# d=128
# qn=100
# min=0
# max=255
# dtype=uint8
# ipath=txt/${dname}/${dname}
# opath=${dname}/${dname}

# ./convert ${n} ${d} ${qn} ${dtype} ${ipath} ${opath} ${min} ${max}
# cp ${ipath}.stat ${opath}.stat 

# # ------------------------------------------------------------------------------
# #  P53
# # ------------------------------------------------------------------------------
# dname=P53
# n=31059
# d=5408
# qn=100
# min=0
# max=10000
# dtype=uint16
# ipath=txt/${dname}/${dname}
# opath=${dname}/${dname}

# ./convert ${n} ${d} ${qn} ${dtype} ${ipath} ${opath} ${min} ${max}
# cp ${ipath}.stat ${opath}.stat 

# # ------------------------------------------------------------------------------
# #  Trevi
# # ------------------------------------------------------------------------------
# dname=Trevi
# n=100800
# d=4096
# qn=100
# min=0
# max=255
# dtype=uint8
# ipath=txt/${dname}/${dname}
# opath=${dname}/${dname}

# ./convert ${n} ${d} ${qn} ${dtype} ${ipath} ${opath} ${min} ${max}
# cp ${ipath}.stat ${opath}.stat 

# # ------------------------------------------------------------------------------
# #  Gist
# # ------------------------------------------------------------------------------
# dname=Gist
# n=1000000
# d=960
# qn=100
# min=0
# max=15000
# dtype=uint16
# ipath=txt/${dname}/${dname}
# opath=${dname}/${dname}

# ./convert ${n} ${d} ${qn} ${dtype} ${ipath} ${opath} ${min} ${max}
# cp ${ipath}.stat ${opath}.stat 

# ------------------------------------------------------------------------------
#  Sift10M
# ------------------------------------------------------------------------------
dname=Sift10M
n=10000000
d=128
qn=100
min=0
max=255
dtype=uint8
ipath=txt/${dname}/${dname}
opath=${dname}/${dname}

./convert ${n} ${d} ${qn} ${dtype} ${ipath} ${opath} ${min} ${max}
cp ${ipath}.stat ${opath}.stat

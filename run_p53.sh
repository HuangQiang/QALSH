#!/bin/bash
make
rm *.o

# ------------------------------------------------------------------------------
#  Parameters
# ------------------------------------------------------------------------------
dname=P53
n=31059
qn=100
d=5408
B=65536
c=2.0

dPath=./data/${dname}/${dname}
dFolder=./data/${dname}/

# ------------------------------------------------------------------------------
#  Running Scripts
# ------------------------------------------------------------------------------
p_list=(0.5 1.0 2.0 0.8 1.2 1.5) 
z_list=(1.0 0.0 0.0 0.0 0.0 0.0)
length=`expr ${#p_list[*]} - 1`

for j in $(seq 0 ${length})
do 
  p=${p_list[j]}
  z=${z_list[j]}
  oFolder=./results${c}/${dname}/L${p}/

  # ----------------------------------------------------------------------------
  #  Ground Truth
  # ----------------------------------------------------------------------------
  ./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds \
    -qs ${dPath}.q -ts ${dPath}.gt${p}

  # ----------------------------------------------------------------------------
  #  QALSH
  # ----------------------------------------------------------------------------
  ./qalsh -alg 3 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds \
    -df ${dFolder} -of ${oFolder}

  ./qalsh -alg 4 -qn ${qn} -d ${d} -qs ${dPath}.q -ts ${dPath}.gt${p} \
    -df ${dFolder} -of ${oFolder}

  # ----------------------------------------------------------------------------
  #  Linear Scan
  # ----------------------------------------------------------------------------
  ./qalsh -alg 5 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q \
    -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}
done
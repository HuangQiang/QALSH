make
rm *.o

# Parameters of QALSH:
#    -alg  (integer)   options of algorithms (0 - 3)
#    -n    (integer)   cardinality of the dataset
#    -qn   (integer)   number of queries
#    -d    (integer)   dimensionality of the dataset
#    -B    (integer)   page size
#    -p    (real)      lp-norm, where p in (0, 2]
#    -z    (real)      symmetric factor of p-stable distr. zeta in [-1, 1]
#    -c    (real)      approximation ratio (c > 1)
#    -ds   (string)    file path of the dataset
#    -qs   (string)    file path of the query set
#    -ts   (string)    file path of the ground truth set
#    -df   (string)    data folder to store new format of data
#    -of   (string)    output folder to store info of qalsh
#
# The options of algorithms (-alg) are:
#    0 - Ground-Truth
#        Parameters: -alg 0 -n -qn -d -p -ds -qs -ts
#
#    1 - Indexing
#        Parameters: -alg 1 -n -d -B -p -z -c -ds -df -of
#
#    2 - QALSH
#        Parameters: -alg 2 -qn -d -p -qs -ts -df -of
#
#    3 - Linear Scan
#        Parameters: -alg 3 -n -qn -d -B -p -qs -ts -df -of
#
# NOTE: Each parameter is required to be separated by one space

n=60000
qn=100
d=50
B=4096
c=2.0
dname=Mnist

dPath=./data/${dname}/${dname}
dFolder=./data/${dname}/
oFolder=./results/${dname}/


### p = 0.5 ###
p=0.5
z=1.0

./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.gt${p}

./qalsh -alg 1 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./qalsh -alg 2 -qn ${qn} -d ${d} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}

./qalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}


### p = 1.0 ###
p=1.0
z=0.0

./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.gt${p}

./qalsh -alg 1 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./qalsh -alg 2 -qn ${qn} -d ${d} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}

./qalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}


### p = 2.0 ###
p=2.0
z=0.0

./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.gt${p}

./qalsh -alg 1 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./qalsh -alg 2 -qn ${qn} -d ${d} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}

./qalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}


### p = 0.8 ###
p=0.8
z=0.0

./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.gt${p}

./qalsh -alg 1 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./qalsh -alg 2 -qn ${qn} -d ${d} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}

./qalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}


### p = 1.2 ###
p=1.2
z=0.0

./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.gt${p}

./qalsh -alg 1 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./qalsh -alg 2 -qn ${qn} -d ${d} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}

./qalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}


### p = 1.5 ###
p=1.5
z=0.0

./qalsh -alg 0 -n ${n} -qn ${qn} -d ${d} -p ${p} -ds ${dPath}.ds -qs ${dPath}.q -ts ${dPath}.gt${p}

./qalsh -alg 1 -n ${n} -d ${d} -B ${B} -p ${p} -z ${z} -c ${c} -ds ${dPath}.ds -df ${dFolder} -of ${oFolder}

./qalsh -alg 2 -qn ${qn} -d ${d} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}

./qalsh -alg 3 -n ${n} -qn ${qn} -d ${d} -B ${B} -p ${p} -qs ${dPath}.q -ts ${dPath}.gt${p} -df ${dFolder} -of ${oFolder}












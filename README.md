# QALSH: Query-Aware Locality-Sensitive Hashing

## Introduction

This package provides two external LSH schemes QALSH and QALSH<sup>+</sup> for high-dimensional ```c-Approximate Nearest Neighbor (c-ANN)``` search under l<sub>p</sub> norm from the following two papers, where 0 < p <= 2.

```bash
Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong Fang, Wilfred Ng. Query-Aware Locality-Sensitive 
Hashing for Approximate Nearest Neighbor Search. Proceedings of the VLDB Endowment (PVLDB), 
9(1): 1-12, 2015.

Qiang Huang, Jianlin Feng, Qiong Fang, Wilfred Ng, Wei Wang. Query-Aware Locality-Sensitive 
Hashing Scheme for l_p Norm. The VLDB Journal, 26(5): 683–708, 2017.
```

## Compilation

The package requires ```g++``` with ```c++11``` support. To download and compile the code, type:

```bash
$ git clone https://github.com/HuangQiang/QALSH.git
$ cd QALSH
$ make
```

## Datasets

We use four real-life datasets [Sift](https://drive.google.com/open?id=1Q3_dnblolD9GVis7OakP2mrqmBApytEL), [Gist](https://drive.google.com/open?id=1M3hJl5slY_pu50IQ7ie-t9E6RvzMizYT), [Trevi](https://drive.google.com/open?id=1RF1FJKWHv3y7W7aBrewnOMrWR15dNbJ3), and [P53](https://drive.google.com/open?id=15mzraPmxNRzcfhXsd_KWBgKclUFUZQEj) for comparison. The statistics of the datasets are summarized in the following table:

| Datasets | #Data Objects  | #Queries | Dimensionality | Page Size | Domain Size | Data Size |
| -------- | --------- | -------- | -------------- | --------- | ----------- | --------- |
| Sift     | 1,000,000 | 100      | 128            | 4 KB      | [0, 218]    | 337.8 MB  |
| Gist     | 1,000,000 | 100      | 960            | 16 KB     | [0, 14,772] | 4.0 GB    |
| Trevi    | 100,900   | 100      | 4,096          | 64 KB     | [0, 255]    | 1.5 GB    |
| P53      | 31,159    | 100      | 5,408          | 64 KB     | [0, 10,000] | 833.7 MB  |

## Run Experiments

```bash
Usage: qalsh [OPTIONS]

This package supports 6 options to evaluate the performance of QALSH, QALSH^+,
and Linear_Scan for c-ANN search. The parameters are introduced as follows.

  -alg    integer    options of algorithms (0 - 5)
  -n      integer    cardinality of dataset
  -d      integer    dimensionality of dataset and query set
  -qn     integer    number of queries
  -B      integer    page size
  -leaf   integer    leaf size of kd_tree
  -L      integer    number of projections for drusilla_select
  -M      integer    number of candidates  for drusilla_select
  -p      float      l_{p} norm, where 0 < p <= 2
  -z      float      symmetric factor of p-stable distribution (-1 <= z <= 1)
  -c      float      approximation ratio for c-ANN search (c > 1)
  -ds     string     address of data  set
  -qs     string     address of query set
  -ts     string     address of truth set
  -df     string     data folder to store new format of data
  -of     string     output folder to store output results
```

We provide the scripts to repeat experiments reported in VLDBJ 2017. A quick example is shown as follows (run QALSH<sup>+</sup> and QALSH on ```Mnist``` with ```Euclidean distance```, where ```p = 2.0``` and ```z = 0.0```):

```bash
# QALSH^+
./qalsh -alg 1 -n 60000 -d 50 -B 4096 -leaf 4000 -L 30 -M 10 -p 2.0 -z 0.0 -c 2.0 -ds data/Mnist/Mnist.ds -df data/Mnist/ -of results2.0/Mnist/L2.0/
./qalsh -alg 2 -qn 100 -d 50 -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.gt2.0 -df data/Mnist/ -of results/Mnist/L2.0/

# QALSH
./qalsh -alg 3 -n 60000 -d 50 -B 4096 -p 2.0 -z 0.0 -c 2.0 -ds data/Mnist/Mnist.ds -df data/Mnist/ -of results2.0/Mnist/L2.0/
./qalsh -alg 4 -qn 100 -d 50 -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.gt2.0 -df data/Mnist/ -of results/Mnist/L2.0/
```

If you would like to get more information to run other algorithms, please check the scripts in the package. When you run the package, please ensure that the path for the dataset, query set, and truth set is correct. Since the package will automatically create folder for the output path, please keep the path as short as possible.

## Parameter Settings

Finally, we introduce some tricks to set up parameters, i.e., ```B```, ```leaf```, ```L```, ```M```, ```p```, ```z```, and ```c```.

### The settings of ```B```

```B``` is the page size in bytes. It is determined by the dimensionality ```d``` of datasets. We use the following rules to set up ```B```:

- (1) ```d < 256```: we set ```B = 4096```.
- (2) ```256 ⩽ d < 512```: we set ```B = 8192```.
- (3) ```512 ⩽ d < 1024```: we set ```B = 16384```.
- (4) ```1024 ⩽ d < 4096```: we set ```B = 32768```.
- (5) ```4096 ⩽ d < 8192```: we set ```B = 65536```.

for the case ```d ⩾ 8192```, we can set up a corresponding larger ```B``` value following the rules above. 

### The settings of ```leaf```, ```L```, and ```M```

```leaf``` is the maximum leaf size of kd-tree. Thus, it should be smaller than the cardinality of dataset, i.e., ```leaf < n```. Let ```K``` be the number of blocks after kd-tree partitioning. Since we use kd-tree to divide the whole datasets into blocks, ```K``` is 2<sup>i</sup>, where ```i = ceil(log_2 (n/leaf))```. 

```L```, and ```M``` are two parameters introduced by Drusilla Select. Once ```leaf``` is determined, the actual leaf size n<sub>0</sub> (also known as the number of objects in each block) can be estimated as floor(n / 2<sup>i</sup>) or ceil(n / 2<sup>i</sup>). There are two conditions when we set up ```L``` and ```M```: 

- (1) L * M < n<sub>0</sub>. Since we run drusilla select for each block to select the representative objects, it is a natural condition to restrict its size (L * M) less than n<sub>0</sub>.
- (2) K * L * M ≈ n<sub>0</sub>. If the sample size (i.e., K*L*M) is large, we can well estimate which blocks are closer to the query, but it will a lot of extra time for estimation. If the sample size is small, the time to determine close blocks can be reduced, but these blocks may not be closer to the query than others. This condition is based on our observation. According to our experiments, we find that creating a sample set with cardinality similar to n<sub>0</sub> can achieve a good trade-off.

### The settings of ```p```, ```z```, and ```c```

```p``` and ```z``` determine the distance metric and the corresponding p-stable distribution. There are three common settings: 

- (1) Euclidean distance (l<sub>2</sub> distance): we set ```p=2.0```, ```z=0.0``` and apply standard Gaussian distribution.
- (2) Manhattan distance (l<sub>1</sub> distance): we set ```p=1.0```, ```z=0.0``` and apply standard Cauchy distribution.
- (3) l<sub>0.5</sub> distance: we set ```p=0.5```, ```z=1.0``` and apply standard Levy distribution. 

In addition, for other l<sub>p</sub> distance, users can set ```0 < p ⩽ 2``` and ```-1 ⩽ z ⩽ 1```.

```c``` is the approximation ratio for c-k-ANN search. We often set ```c=2.0```. But if the dataset is easy, it is also satisfied to set ```c=3.0```.

## Related Publications

If you use this package for publications, please cite the papers as follows.

```bib
@article{huang2017query,
    title={Query-aware locality-sensitive hashing scheme for $$ l\_p $$ norm}
    author={Huang, Qiang and Feng, Jianlin and Fang, Qiong and Ng, Wilfred and Wang, Wei},
    booktitle={The VLDB Journal},
    volumn={26},
    number={5},
    pages={683--708},
    year={2017},
    organization={Springer}
}

@article{huang2015query,
    title={Query-aware locality-sensitive hashing for approximate nearest neighbor search}
    author={Huang, Qiang and Feng, Jianlin and Zhang, Yikai and Fang, Qiong and Ng, Wilfred},
    booktitle={Proceedings of the VLDB Endowment},
    volumn={9},
    number={1},
    pages={1--12},
    year={2015},
    organization={VLDB Endowment}
}
```

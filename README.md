# QALSH: Query-Aware Locality-Sensitive Hashing

## Introduction

This package provides two external LSH schemes QALSH and QALSH<sup>+</sup> for high-dimensional ```c-Approximate Nearest Neighbor (c-ANN)``` search under l<sub>p</sub> norm from the following two papers, where 0 < p <= 2.

```bash
Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong Fang, Wilfred Ng. Query-Aware  Locality-Sensitive 
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

| Datasets | #Objects  | #Queries | Dimensionality | Page Size | Domain Size | Data Size |
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
  -nb     integer    number of blocks for c-ANN search
  -p      float      l_{p} norm, where 0 < p <= 2
  -z      float      symmetric factor of p-stable distribution (-1 <= z <= 1)
  -c      float      approximation ratio for c-ANN search (c > 1)
  -ds     string     address of data  set
  -qs     string     address of query set
  -ts     string     address of truth set
  -df     string     data folder to store new format of data
  -of     string     output folder to store output results
```

We provide the scripts to repeat experiments reported in VLDBJ 2017. A quick example is shown as follows (run QALSH<sup>+</sup> and QALSH on ```Mnist```, where ```p = 2.0```):

```bash
# QALSH^+
./qalsh -alg 1 -n 60000 -d 50 -B 4096 -leaf 4000 -L 30 -M 10 -p 2.0 -z 0.0 -c 2.0 -ds data/Mnist/Mnist.ds -df data/Mnist/ -of results2.0/Mnist/L2.0/
./qalsh -alg 2 -qn 100 -d 50 -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.gt2.0 -df data/Mnist/ -of results/Mnist/L2.0/

# QALSH
./qalsh -alg 3 -n 60000 -d 50 -B 4096 -p 2.0 -z 0.0 -c 2.0 -ds data/Mnist/Mnist.ds -df data/Mnist/ -of results2.0/Mnist/L2.0/
./qalsh -alg 4 -qn 100 -d 50 -qs data/Mnist/Mnist.q -ts data/Mnist/Mnist.gt2.0 -df data/Mnist/ -of results/Mnist/L2.0/
```

If you would like to get more information to run other algorithms, please check the scripts in the package. When you run the package, please ensure that the path for the dataset, query set, and truth set is correct. Since the package will automatically create folder for the output path, please keep the path as short as possible.

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

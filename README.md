# QALSH: Query-Aware Locality-Sensitive Hashing

Welcome to the **QALSH** GitHub!

**QALSH** is a package for the problem of Nearest Neighbor Search (NNS) over high-dimensional Euclidean spaces. Given a set of data points and a query, the problem of NNS aims to find the nearest data point to the query. It is a very fundamental probelm and has wide applications in many data mining and machine learning tasks.

This package provides the external memory implementations (disk-based) of QALSH and QALSH<sup>+</sup> for *c*-Approximate Nearest Neighbor Search (c-ANNS) under *l<sub>p</sub>* norm, where *0 < p ⩽ 2*. The internel memory version can be found [here](https://github.com/HuangQiang/QALSH_Mem).

If you want to get more details of QALSH and QALSH<sup>+</sup>, please refer to our works [Query-Aware Locality-Sensitive Hashing for Approximate Nearest Neighbor Search](https://dl.acm.org/doi/abs/10.14778/2850469.2850470) and [Query-Aware Locality-Sensitive Hashing Scheme for *l<sub>p</sub>* Norm](https://link.springer.com/article/10.1007/s00778-017-0472-7), which have been published in PVLDB 2015 and VLDBJ 2017, respectively.

## Datasets and Queries

We study the performance of QALSH and QALSH<sup>+</sup> over 4 real-life high-dimensional [datasets](https://drive.google.com/drive/folders/1tKMl0_iLSEeuT1ZJ7s4x1BbLHyX0D5OJ), i.e., **Sift**, **Gist**, **Trevi**, and **P53**. We also include a toy dataset **Mnist** for illustration and a large-scale dataset **Sift10M** for further validation. For each dataset, we provide 100 queries (randomly select from its test set or extract from the dataset itself) for evaluations. The statistics of datasets and queries are summarized as follows.

| Datasets | #Data Points (*n*) | #Queries | Dimensionality (*d*) | Range of Coordinates | Data Type |
| -------- | ------------ | ----- | ------- | ----------- | ------- |
| Sift     | 1,000,000    | 100   | 128     | [0, 255]    | uint8   |
| Gist     | 1,000,000    | 100   | 960     | [0, 15,000] | uint16  |
| Trevi    | 100,800      | 100   | 4,096   | [0, 255]    | uint8   |
| P53      | 31,059       | 100   | 5,408   | [0, 10,000] | uint16  |
| Mnist    | 60,000       | 100   | 50      | [0, 255]    | uint8   |
| Sift10M  | 10,000,000   | 100   | 128     | [0, 255]    | uint8   |

Note that all the datasets and queries are in a binary format, which can be considered as an array of `n·d` coordinates, where each coordinate is specified by the data type, e.g., `uint8` or `uint16`. We currently support four data types: `uint8`, `uint16`, `int32`, and `float32`. One can determine the data type of the dataset based on the range of its coordinates. If you want to support more data types, you can update the interface in the `main.cc` and re-compile the package.

## Compilation

The package requires `g++` with `c++11` support. To download and compile the c++ source codes, please run the commands as follows:

```bash
git clone git@github.com/HuangQiang/QALSH.git
cd QALSH/methods/
make -j
```

## Usages

Suppose you have cloned the project and you are in the folder `QALSH/`. We provide bash scripts to run experiments for the six real-life datasets.

### Step 1: Get the Datasets and Queries

Please download the [datasets](https://drive.google.com/drive/folders/1tKMl0_iLSEeuT1ZJ7s4x1BbLHyX0D5OJ) and copy them to the directory `data/`.

For example, when you get `Sift.ds` and `Sift.q`, please move them to the paths `data/Sift/Sift.ds` and `data/Sift/Sift.q`, respectively. We provide Mnist at this package, you can follow the same pattern to move the datasets to the right place.

### Step 2: Run Experiments

When you run the package, please ensure the paths for the dataset, query set, and truth set are correct. The package will automatically create folder for the output path, so please keep the output path *as short as possible*. All of the experiments can be run with the following commands:

```bash
cd methods/
bash run_all.sh
```

A gentle reminder is that when running QALSH and QALSH<sup>+</sup>, since they need the ground truth results for evaluation, please run `-alg 0` to get the ground truth results first.

### Step 3. Parameter Settings

Finally, if you would like to use this package for *c*-ANNS over other datasets, you may want to get more information about the parameters and know how to set them effectively.
Based on our experience when we conducted the experiments, we now share some tricks on setting up the parameters, i.e., `B`, `lf`, `L`, `M`, `p`, `z`, and `c`. The illustration of the parameters are as follows.

```bash
Usage: qalsh [OPTIONS]

This package supports 6 options to evaluate the performance of QALSH, QALSH^+,
and Linear_Scan for c-k-ANNS. The parameters are introduced as follows.

  -alg    integer    options of algorithms (0 - 5)
  -n      integer    cardinality of dataset
  -d      integer    dimensionality of dataset and query set
  -qn     integer    number of queries
  -B      integer    page size
  -lf     integer    leaf size of kd_tree
  -L      integer    number of projections for drusilla_select
  -M      integer    number of candidates  for drusilla_select
  -p      float      l_{p} norm, where 0 < p ⩽ 2
  -z      float      symmetric factor of p-stable distribution (-1 ⩽ z ⩽ 1)
  -c      float      approximation ratio for c-k-ANNS (c > 1)
  -dt     string     data type (i.e., uint8, uint16, int32, float32)
  -pf     string     the prefix of dataset, query set, and truth set
  -of     string     output folder to store output results
```

#### The settings of `B`

`B` is the page size in bytes. It is determined by the data dimensiona `d`. Based on our experience, we use the following rules to set up `B`:

- (1) 1    ⩽ d < 256: we set **B = 4096**;
- (2) 256  ⩽ d < 512: we set **B = 8192**;
- (3) 512  ⩽ d < 1024: we set **B = 16384**;
- (4) 1024 ⩽ d < 4096: we set **B = 32768**;
- (5) 4096 ⩽ d < 8192: we set **B = 65536**;

for the case d ⩾ 8192, we can set up a corresponding larger `B` value following the rules above.

#### The settings of `lf`, `L`, and `M`

`lf` is the maximum leaf size of kd-tree. `L` and `M` are two parameters used for Drusilla_Select, where `L` is the number of random projections; `M` is the number of representative data points we select on each random projection.

Let `K` be the number of blocks, and let n<sub>0</sub> be the actual leaf size after kd-tree partitioning. Once `lf` is determined, `K` and **n<sub>0</sub>** can be computed as follows:

- **K = 2<sup>h</sup>**, where `h` is the height of the kd-tree, i.e., **h = ceil(log_2 (n / lf))**;
- **n<sub>0</sub> = floor(n / K)** or **n<sub>0</sub> = ceil(n / K)** (Note: if `n` is not divisible, these two cases can happen.)

When we set up these three parameters `lf`, `L` and `M`, it might be better to satisfy the following three conditions:

- **lf < n**: It is a natural condition that the maximum leaf size should be smaller than the cardinality of dataset;
- **L · M < n<sub>0</sub>**: It is a natural condition to restrict its size **(L · M)** less than **n<sub>0</sub>** when we run Drusilla_Select to select the representative data points on each block;
- **K · L · M ≈ n<sub>0</sub>**: This condition is the main trade-off between efficiency and accuracy to set up these three parameters.
  - On the one hand, if the total number of representative data points **(K · L · M)** is large, we can accurately identify the close blocks to the query, but it may introduce much time for the first-level close block search.
  - On the other hand, if **(K · L · M)** is small, the time for the first-level close block search can be reduced. However, since these blocks may not be really close to the query, the accuracy for the second-level c-ANNS may also be reduced.
  - In our experiments, we find that **selecting the representative data points with cardinality approximate to n<sub>0</sub>** can achieve a good trade-off between accuracy and efficiency.

#### The settings of `p`, `z`, and `c`

`p` and `z` determine the distance metric and the corresponding p-stable distribution. There are four cases as follows.

- *l<sub>2</sub>* distance: set up **p = 2.0** and **z = 0.0**, and apply standard Gaussian distribution.
- *l<sub>1</sub>* distance: set up **p = 1.0** and **z = 0.0** and apply standard Cauchy distribution.
- *l<sub>0.5</sub>* distance: set up **p = 0.5** and **z = 1.0** and apply standard Levy distribution.
- General *l<sub>p</sub>* distance: set up **0 < p ⩽ 2** and **-1 ⩽ z ⩽ 1**.

```c``` is the approximation ratio for *c*-ANNS. We set **c = 2.0** by default. If the dataset is easy, it is also satisfied to set **c = 3.0** or am even larger value.

## Reference

Please use the following BibTex to cite this work if you use **QALSH** for publications.

```tex
@article{huang2017query,
    title={Query-aware locality-sensitive hashing scheme for $$ l\_p $$ norm}
    author={Huang, Qiang and Feng, Jianlin and Fang, Qiong and Ng, Wilfred and Wang Wei},
    booktitle={The VLDB Journal},
    volumn={26},
    number={5},
    pages={683--708},
    year={2017}
}

@article{huang2015query,
    title={Query-aware locality-sensitive hashing for approximate nearest neighbor search}
    author={Huang, Qiang and Feng, Jianlin and Zhang, Yikai and Fang, Qiong and Ng, Wilfred},
    booktitle={Proceedings of the VLDB Endowment},
    volumn={9},
    number={1},
    pages={1--12},
    year={2015}
}
```

It is welcome to contact me (<huangq@comp.nus.edu.sg>) if you meet any issue. Thank you.

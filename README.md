# QALSH: Query-Aware Locality-Sensitive Hashing

Version: 1.3.0

Release date: 3 October 2017


Introduction
--------

This package is written in the C++ programming language. It provides an 
external LSH scheme (QALSH) for c-Approximate Nearest Neighbor (or simply 
c-ANN) search under l_p norm, where p \in (0, 2].

How to use this Package?
--------

We provide a Makefile and a script (i.e., run_mnist.sh) as a running example 
for comipling and running this package. Before start running this package, 
please ensure the input format of the dataset and query set is correct. We 
provide a sample dataset and query set (i.e., Mnist) for your reference.

We also provide the scripts (i.e., run_sift.sh, run_gist.sh, run_trevi.sh, 
and run_p53.sh) for the users who would like to reproduce our results presented 
in PVLDB 2015 and VLDBJ 2017. The datasets Sift, Gist, Trevi, and P53 we used 
can be downloaded from the following links:

* Sift: https://drive.google.com/open?id=1Q3_dnblolD9GVis7OakP2mrqmBApytEL

* Gist: https://drive.google.com/open?id=1M3hJl5slY_pu50IQ7ie-t9E6RvzMizYT

* Trevi: 

* P53: 


Authors
--------

* **Qiang Huang**

  Smart Systems Institute, National University of Singapore (NUS),
  
  Singapore, 119613 
  
  huangq2011@gmail.com, huangq25@mail2.sysu.edu.cn
  
  https://sites.google.com/site/qianghuang2017/
  

* **Jianlin Feng**

  School of Data and Computer Science, Sun Yat-Sen University (SYSU),
  
  Guangzhou, China, 510006
  
  fengjlin@mail.sysu.edu.cn, fengjl9@gmail.com
  
  http://ss.sysu.edu.cn/~fjl/


Relevant Papers
--------

The paper for the package of QALSH has been published in VLDB 2016 and The VLDB 
Journal, which are displayed as follows:

* **Qiang Huang, Jianlin Feng, Yikai Zhang, Qiong Fang, and Wilfred Ng. Query-Aware
Locality-Sensitive Hashing for Approximate Nearest Neighbor Search. In 
Proceedings of the VLDB Endowment, 9(1): 1 - 12, 2015.**

* **Qiang Huang, Jianlin Feng, Qiong Fang, Wilfred Ng, and Wei Wang. Query-Aware 
Locality-Sensitive Hashing Scheme for l_p Norm. The VLDB Journal, 26(5): 683 – 
708, 2017.**

If you use the package for publications, please cite the papers above.


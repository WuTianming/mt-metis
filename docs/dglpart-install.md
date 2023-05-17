## ParMETIS Installation

### Pre-requisites

The following software are required to compile ParMETIS:

```
sudo apt-get install build-essential
sudo apt-get install cmake
```

You will also need MPI compilers in your PATH.

### Install ParMETIS

Let `prefix` be a single directory where all of the following software will be installed to. For example,

```
prefix=/home/user/parmetis
```

Install GKlib:

```shell
git clone https://github.com/KarypisLab/GKlib.git && cd GKlib
make config prefix=$prefix
make -j$(nproc)
make install
```

Install METIS:

```shell
git clone https://github.com/KarypisLab/METIS.git && cd METIS
make config i64=1 prefix=$prefix
make -j$(nproc)
make install
```

We can install ParMETIS after installing METIS:

```shell
git clone https://github.com/KarypisLab/ParMETIS.git && cd ParMETIS
make config prefix=$prefix
make -j$(nproc)
make install
```

The command above installs the ParMETIS software suite into `$prefix`. The executables can be found under `$prefix/bin`. The executables include `pm_dglpart`, which is what we are going to use as our baseline.


## Preparing the dataset

`ParMETIS` uses a COO-based input format that can be generated using `testbench/dglpart_prepare.py`:

```shell
python testbench/dglpart_prepare.py \
    --dataset ogbn-products \
    --root /data1 \
    --output /data1/dglpart_data/
```



## Running

You can launch `pm_dglpart` with the following command:

```
export LD_LIBRARY_PATH=$prefix/lib:$LD_LIBRARY_PATH

mpirun -np <#procs> \
    /path/to/pm_dglpart \
    <dataset> \
    <partitions per process>
```

For example, to partition a graph into 8 partitions using 4 processes, the `<partitions per process>` argument should be 8/4 = 2.
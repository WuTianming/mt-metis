This is a fork of mt-metis (by Dominique LaSalle) focusing on reducing memory footprint.



## Installation

```shell
./configure --vertices64bit --edges64bit --weights64bit
make -j$(nproc)
```

The binary can be found under `/path/to/source/build/Linux-x86_64/bin`.



## Running

```shell
OMP_NUM_THREADS=<#nthreads> \
    <path/to/mtmetis> \
    -t \
    <path/to/graph data> \
    <#parts> \
    <path/to/output/partfile>
```

explanation:

- `#nthreads`: number of threads to use
- `-t`: print timing information
- `graph data`: can be METIS format (`.graph` file) or binary file (see `testbench/mtmetis_prepare.py`).
- `#parts`: total number of partitions
- `partfile`: the file to output partition assignments



Note:

- `#nthreads` and `#parts` do **not** require to divide each other.
- The program generates temporary files under the current work directory that can be very big (about several times the size of the original input graph). Please consider running `mt-metis` under a dedicated temporary directory that has sufficient disk space. You may want to manually delete the temporary directory after the program finishes.



## mt-metis testbench

Run the preparation script to generate binary-format input file (requires dgl python module):

```shell
python mtmetis_prepare.py --dataset ogbn-papers400M --root /data1 --balance-train --balance-edges --output /data1/binary/
```

Edit `testbench/launch_mtmetis.sh` to change `TMPDIR`,`LOGDIR` and `MTMETIS` path. Edit `testbench/experiment_mtmetis.sh` to change the dataset `prefix`. Run `testbench/experiment_mtmetis.sh` to run the testbench.

The testbench uses GNU `time` to measure the elapsed time and peak resident memory.

The testbench uses `nmon` to generate a runtime memory consumption profile which can be visualized later.

This is a fork of the **m**ulti-**t**hreaded **METIS** partitioning algorithm `mt-metis` by Dominique LaSalle, focused on reducing memory footprint.



## Installation

### mt-metis

In order to prepare the dataset, you need to have the following python modules installed:

```shell
pip install dgl
pip install ogb
```

`mt-metis` itself does not have any dependencies. To compile `mt-metis`:

```shell
./configure --vertices64bit --edges64bit --weights64bit
make -j$(nproc)
```

The executable can be found under `/path/to/source/build/Linux-x86_64/bin`.



## Prepare Dataset

The script `mtmetis_prepare.py` under `testbench/` converts a graph dataset into a CSR-based binary format that can be used directly by `mt-metis`. To run the preparation script:

```shell
python mtmetis_prepare.py \
    --dataset ogbn-papers400M \
    --root /data1 \
    [--balance-train] \
    [--balance-edges] \
    --output /path/to/data/binary/
```

where the dataset can be `ogbn-arxiv`, `ogbn-products`, `ogbn-papers100M` or `ogbn-papers400M`. The flag `--balance-train` and `--balance-edges` are optional. The flag `--balance-train` adds a constraint to balance the number of training nodes in each partition. The flag `--balance-edges` adds a constraint to balance the number of edges in each partition.

The script will download the dataset from the OGB website and convert it into binary format. The binary format is a collection of the following files:

- `[graphname]_meta.txt`: the metadata of the dataset, including the number of nodes, edges, and the number of constraints.
- `[graphname]_vwgt.bin`: the node weights. Used as constraints to balance the nodes in each partition.
- `[graphname]_indptr.bin` and `xxx_indices.bin`: the graph structure in CSR format.

You can then pass `[graphname]` as the input to `mt-metis`.




## Running

Use the following command to run `mt-metis` on a prepared binary format dataset:

```
OMP_NUM_THREADS=<#nthreads> \
    <path/to/mtmetis> \
    [-t] \
    <path/to/graph data> \
    <#parts> \
    [<path/to/output/partfile>]
```

Arguments:

- `#nthreads`: number of threads to use
- `-t`: print timing information (optional)
- `graph data`: can be METIS format (`.graph` file) or binary format
- `#parts`: total number of partitions
- `partfile`: the file to output partition assignments (optional)



Note:

- `#nthreads` and `#parts` do **not** require to divide each other.
- The program generates temporary files under the current work directory that can be very big (about several times the size of the original input graph). Please consider running `mt-metis` under a dedicated temporary directory that has sufficient disk space. You may want to manually delete the temporary directory after the program finishes.




## Benchmark

The benchmarks use GNU `time` to measure the elapsed time and peak resident memory.

### DGL Baseline

DGL provides parallel METIS partitioning functionality via `ParMETIS`. Refer to [dglpart-install](docs/dglpart-install.md) to install the required programs and prepare the dataset.

Before running the benchmark, first fill out the "CONFIGURATION" section in `testbench/launch_parmetis.sh`.

Make sure to put the input data in the current working directory, then run the benchmark:

```shell
./testbench/experiment_baseline.sh
```



### Ours

Edit `testbench/launch_mtmetis.sh` to change `TMPDIR`,`LOGDIR` and `MTMETIS` path. Edit `testbench/experiment_mtmetis.sh` to set the location of the datasets.

Then run `testbench/experiment_mtmetis.sh` to launch the benchmark.



## Results

(TODO: table)
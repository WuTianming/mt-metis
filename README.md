## mt-metis for refinement only

编译方法：

```
./configure ./configure --vertices64bit --edges64bit --weights64bit
make -j8
```





输入格式：兼容 xtrapulp 使用的 `.graph` 格式。需要额外附加一个 xtrapulp 切分好的 `.txt` 格式的 partition 文件，并命名为 "`xxxx.graph_where.txt`"，其中 `xxxx.graph` 是输入图文件的名称。这个 partition 文件需要与后面运行命令中的 "partition 数" 一致。





运行方法：

```
OMP_NUM_THREADS=16 ./build/Linux-x86_64/bin/mtmetis -t -b 1.05 \
    [XXXX.graph] \
    [partition 数] \
    [输出文件名.txt]
```

例如

```
OMP_NUM_THREADS=16 ./build/Linux-x86_64/bin/mtmetis -t -b 1.05 \
    ogbn-products_xtrapulp.graph \
    4 \
    输出文件名.txt
```

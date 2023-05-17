# generates **dglpart** input files

import code
import dgl
from tqdm import tqdm
import torch as th
import numpy as np
import os

# the load_* function returns (graph, train_mask, isdirected)

# note: dglpart 版本的 graph 中，双向边只出现一次，因此不需要手动增加反向边
def load_ogb(name, root):
    from ogb import nodeproppred
    print('load', name)
    data = nodeproppred.DglNodePropPredDataset(name=name, root=root)
    print('finish loading', name)

    graph, _ = data[0]

    print('generating train-mask')
    train_list = data.get_idx_split()['train']
    train_mask = th.zeros(graph.num_nodes(), dtype=th.bool)
    train_mask[train_list] = True
    print('finish generating train-mask')

    print(graph)
    # code.interact(local=locals())

    print('removing self loop')
    graph = dgl.remove_self_loop(graph)

    DirectedGraphs = ['ogbn-arxiv', 'ogbn-papers100M', 'ogbn-mag']

    isdirected = (name in DirectedGraphs)

    if not isdirected:
        print("will skip reverse edges when dumping edges")
        # not skipping is also ok, but slows down dglpart by 2x

    return graph, train_mask, isdirected

if __name__ == '__main__':
    # dataset = 'ogbn-arxiv'
    dataset = 'ogbn-products'
    # dataset = 'ogbn-papers100M'

    g, train_mask, isdirected = load_ogb(dataset, root="/data1")

    srcs, dsts = g.all_edges()

    # print(srcs, dsts)

    prefix = '/data1/mt-metis-final-exp/datasets/dglpart/' + dataset

    nn = g.num_nodes()
    ne = g.num_edges()
    ncon = 2            # first = trainable, second = degree

    with open(prefix + '_stats.txt', 'w') as f:
        if isdirected:
            print(nn, ne, ncon, file=f)
        else:
            print(nn, ne // 2, ncon, file=f)

    print('dump nodes...')
    with open(prefix + '_nodes.txt', 'w') as f:
        # for i in tqdm(range(nn)):
        #     print(0, 1, i, file=f)

        qwq = np.zeros((4, nn), dtype=np.int32)

        # first vertex weight
        qwq[1] = train_mask.numpy().astype(np.int32)

        # second vertex weight
        if not isdirected:
            qwq[2] = g.in_degrees().numpy()
        else:
            qwq[2] = g.in_degrees().numpy() + g.out_degrees().numpy()

        # type wise node id
        qwq[3] = np.arange(nn)

        # now dump the file using numpy
        # -- I expect this to be way faster than python loops
        np.savetxt(f, qwq.T, fmt='%d')

    print('dump edges...')
    with open(prefix + '_edges.txt', 'w') as f:
        # for i in tqdm(range(len(srcs))):
        #     print(srcs[i].item(), dsts[i].item(), i, 0, file=f)

        # pack srcs and dsts into numpy array
        # if not isdirected, skip the edges with even id
        # else, skip nothing
        qwq = np.zeros((2, ne), dtype=np.int32)
        qwq[0] = srcs.numpy()
        qwq[1] = dsts.numpy()
        if not isdirected:
            qwq = qwq[:, 1::2]
        np.savetxt(f, qwq.T, fmt='%d')
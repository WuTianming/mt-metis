import sys
import os
import code
import argparse
import time
import torch as th
import dgl
from utils.load_graph import load_papers400m_sparse, load_ogb

'''
usage:

python mtmetis_prepare.py --dataset ogbn-papers400M --root /data1 --balance-train --balance-edges --output /data1/binary/
'''

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('--root',
                           default='dataset/',
                           help='Path of the dataset.')
    argparser.add_argument(
        "--dataset",
        default="ogbn-papers400M",
        choices=["ogbn-arxiv", "ogbn-products", "ogbn-papers100M", "ogbn-papers400M"])
    argparser.add_argument("--balance-train",
                           action="store_true",
                           help="balance the training size in each partition.")
    argparser.add_argument("--balance-edges",
                           action="store_true",
                           help="balance the number of edges in each partition.")
    argparser.add_argument("--output",
                           type=str,
                           default="data",
                           help="Output directory to dump binary graph.")
    args = argparser.parse_args()

    start = time.time()
    if args.dataset == "ogbn-products":
        g, _ = load_ogb("ogbn-products", root=args.root)
        train_mask = g.ndata["train_mask"]
        g_coo = g
        del g
    elif args.dataset == "ogbn-arxiv":
        g, _ = load_ogb("ogbn-arxiv", root=args.root)
        train_mask = g.ndata["train_mask"]
        srcs, dsts = g.all_edges()
        del g
        # create dgl graph from srcs and dsts
        print('creating clean-slate graph from coo')
        g_coo = dgl.graph((srcs, dsts))
        print('adding reverse edges to coo')
        g_coo = dgl.add_edges(g_coo, dsts, srcs)
        del srcs, dsts
    elif args.dataset == "ogbn-papers100M":
        g, _ = load_ogb("ogbn-papers100M", root=args.root)
        train_mask = g.ndata["train_mask"]
        srcs, dsts = g.all_edges()
        del g
        # create dgl graph from srcs and dsts
        print('creating clean-slate graph from coo')
        g_coo = dgl.graph((srcs, dsts))
        print('adding reverse edges to coo')
        g_coo = dgl.add_edges(g_coo, dsts, srcs)
        del srcs, dsts
    elif args.dataset == "ogbn-papers400M":
        g, _ = load_papers400m_sparse(root='/data1', omit_features=True)
        train_mask = g.ndata["train_mask"]
        srcs, dsts = g.all_edges()
        del g
        # create dgl graph from srcs and dsts
        print('creating clean-slate graph from coo')
        g_coo = dgl.graph((srcs, dsts))
        print('adding reverse edges to coo')
        g_coo = dgl.add_edges(g_coo, dsts, srcs)
        del srcs, dsts

    print("load {} takes {:.3f} seconds".format(args.dataset,
                                                time.time() - start))
    print("|V|={}, |E|={}".format(g_coo.number_of_nodes(), g_coo.number_of_edges()))

    # instead of evoking mt-metis here (we don't yet have a python wrapper),
    # just dump the loaded graph into binary vectors: xadj,adjncy,vwgt

    # code.interact(local=locals())

    print('generating partitioning vertex weights')
    vwgt = []
    balance_cnt = 0
    if args.balance_train:
        print('retrieving vertex weight: train_mask')
        vwgt.append(train_mask)
        balance_cnt += 1
    if args.balance_edges:
        print('retrieving vertex weight: node degrees')
        vwgt.append(g_coo.out_degrees())
        balance_cnt += 1

    if len(vwgt) > 0:
        print('stacking vertex weights')
        vwgt = th.stack(vwgt).T # arrange the weights into vector
    else:
        vwgt = None

    # retrieve sparse vectors
    print('generating csr from coo')
    xadj, adjncy, _ = g_coo.adj_sparse('csr')
    del _

    # generate metadata file
    print(g_coo)

    outpath = os.path.join(args.output, args.dataset)
    print("generate metadata...")
    with open(outpath + "_meta.txt", "w") as f:
        f.write("{}\n{}\n{}\n{}\n".format(
            g_coo.number_of_nodes(),
            g_coo.number_of_edges()//2,
            g_coo.number_of_edges(),    # actual edge number (=xadj[n+1])
            balance_cnt))
    # then dump arrays
    if not vwgt is None:
        print("dump vwgt...")
        vwgt.numpy().astype('int64').tofile(outpath + "_vwgt.bin")
    print("dump xadj...")
    xadj.numpy().astype('int64').tofile(outpath + "_indptr.bin")
    print("dump adjncy...")
    adjncy.numpy().astype('int64').tofile(outpath + "_indices.bin")

    print("binary graph data dump done!")
    print("Please run mt-metis on the dumped files to generate partitioning results using the following command:")
    print("OMP_NUM_THREADS=<#threads> /path/to/mtmetis /path/to/{} <total #parts> /path/to/output".format(args.dataset))

if __name__ == "__main__":
    main()

import dgl
import torch as th
import numpy as np
import os


def load_papers400m_sparse(root="dataset", load_true_features=True, omit_features=False):
    if os.path.exists(os.path.join(root, 'papers400M_sparse.dgl')):
        print('load papers400M_sparse.dgl')
        graph = dgl.load_graphs(os.path.join(root,
                                             'papers400M_sparse.dgl'))[0][0]
        print('finish loading papers400M_sparse.dgl')

        if not omit_features:
            if load_true_features:
                print("load features from papers400M_features.npy")
                original_features = th.from_numpy(
                    np.load(os.path.join(root, 'papers400M_features.npy')))
            else:
                print("generate fake features")
                original_features = th.rand([111059956, 128]).float()

            print("copy the features quadruple")
            graph.ndata['features'] = th.cat([
                original_features, original_features, original_features,
                original_features
            ], dim=0)
            del original_features

            print('finish constructing papers400M_sparse features')

        return graph, 172
    else:
        from ogb import nodeproppred
        print('load papers100m.')
        data = nodeproppred.DglNodePropPredDataset(name="ogbn-papers100M",
                                                   root=root)
        print('finish loading papers100m.')
        splitted_idx = data.get_idx_split()
        original_graph, labels = data[0]
        original_labels = labels[:, 0]
        original_features = original_graph.ndata.pop('feat')

        original_src, original_dst = original_graph.edges()

        n_nodes = original_graph.number_of_nodes()
        n_edges = original_graph.number_of_edges()
        del original_graph
        # repeat original_graph for 4 times
        intra_src = th.cat([
            th.arange(n_nodes, dtype=th.int64).repeat(3).flatten(),
            th.arange(n_nodes, 2 * n_nodes,
                      dtype=th.int64).repeat(3).flatten(),
            th.arange(2 * n_nodes, 3 * n_nodes,
                      dtype=th.int64).repeat(3).flatten(),
            th.arange(3 * n_nodes, 4 * n_nodes,
                      dtype=th.int64).repeat(3).flatten()
        ])
        intra_dst = th.cat([
            th.arange(n_nodes, 4 * n_nodes, dtype=th.int64),
            th.cat([
                th.arange(0 * n_nodes, 1 * n_nodes, dtype=th.int64),
                th.arange(2 * n_nodes, 4 * n_nodes, dtype=th.int64)
            ]),
            th.cat([
                th.arange(0 * n_nodes, 2 * n_nodes, dtype=th.int64),
                th.arange(3 * n_nodes, 4 * n_nodes, dtype=th.int64)
            ]),
            th.arange(3 * n_nodes, dtype=th.int64)
        ])
        print(intra_src.shape)
        sm = th.randint(0, 4, (2 * n_edges, )).to(th.int64)
        dm = th.randint(0, 4, (2 * n_edges, )).to(th.int64)
        src = th.cat([
            original_src + sm[:n_edges].mul(n_nodes),
            original_dst + sm[n_edges:].mul(n_nodes), intra_src
        ],
                     dim=0)
        dst = th.cat([
            original_dst + dm[:n_edges].mul(n_nodes),
            original_src + dm[n_edges:].mul(n_nodes), intra_dst
        ],
                     dim=0)
        print(src.shape, dst.shape)
        graph_coo = dgl.graph((src, dst))
        del src, dst, intra_src, intra_dst, sm, dm
        out_degrees = graph_coo.out_degrees()
        print("Generating csc formats")
        indptr = graph_coo.adj_sparse('csc')[0]
        indices = graph_coo.adj_sparse('csc')[1]
        graph = dgl.graph(('csc', (indptr, indices, [])))
        graph.ndata['out_degrees'] = out_degrees
        del graph_coo
        np.save(os.path.join(root, 'papers400M_features.npy'),
                original_features)
        graph.ndata['labels'] = th.cat([
            original_labels, original_labels, original_labels, original_labels
        ],
                                       dim=0)

        num_labels = len(th.unique(labels[th.logical_not(th.isnan(labels))]))

        # Find the node IDs in the training, validation, and test set.
        train_nid, val_nid, test_nid = splitted_idx['train'], splitted_idx[
            'valid'], splitted_idx['test']
        train_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
        train_mask[train_nid] = True
        train_mask[train_nid + n_nodes] = True
        train_mask[train_nid + 2 * n_nodes] = True
        train_mask[train_nid + 3 * n_nodes] = True
        val_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
        val_mask[val_nid] = True
        val_mask[val_nid + n_nodes] = True
        val_mask[val_nid + 2 * n_nodes] = True
        val_mask[val_nid + 3 * n_nodes] = True
        test_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
        test_mask[test_nid] = True
        test_mask[test_nid + n_nodes] = True
        test_mask[test_nid + 2 * n_nodes] = True
        test_mask[test_nid + 3 * n_nodes] = True
        graph.ndata['train_mask'] = train_mask
        graph.ndata['val_mask'] = val_mask
        graph.ndata['test_mask'] = test_mask

        dgl.save_graphs(os.path.join(root, 'papers400M_sparse.dgl'), [graph])
        graph.ndata['features'] = th.cat([
            original_features, original_features, original_features,
            original_features
        ],
                                         dim=0)
        del original_features, original_labels

        print('finish constructing papers400m_sparse', graph)

    return graph, num_labels


def load_papers400m(root="dataset", load_true_features=True):
    if os.path.exists(os.path.join(root, 'papers400M.dgl')):
        print('load papers400M.dgl')
        graph = dgl.load_graphs(os.path.join(root, 'papers400M.dgl'))[0][0]
        print('finish loading papers400M.dgl')
        if load_true_features:
            original_features = th.from_numpy(
                np.load(os.path.join(root, 'papers400M_features.npy')))
        else:
            original_features = th.rand([111059956, 128]).float()
        graph.ndata['features'] = th.cat([
            original_features, original_features, original_features,
            original_features
        ],
                                         dim=0)
        print('finish constructing papers400M features')
        del original_features
        return graph, 172
    else:
        from ogb import nodeproppred
        print('load papers100m.')
        data = nodeproppred.DglNodePropPredDataset(name="ogbn-papers100M",
                                                   root=root)
        print('finish loading papers100m.')
        splitted_idx = data.get_idx_split()
        original_graph, labels = data[0]
        original_labels = labels[:, 0]
        original_features = original_graph.ndata.pop('feat')

        original_src, original_dst = original_graph.edges()

        n_nodes = original_graph.number_of_nodes()
        del original_graph
        # repeat original_graph for 4 times
        intra_src = th.cat([
            th.arange(n_nodes, dtype=th.int64).repeat(3).flatten(),
            th.arange(n_nodes, 2 * n_nodes,
                      dtype=th.int64).repeat(3).flatten(),
            th.arange(2 * n_nodes, 3 * n_nodes,
                      dtype=th.int64).repeat(3).flatten(),
            th.arange(3 * n_nodes, 4 * n_nodes,
                      dtype=th.int64).repeat(3).flatten()
        ])
        intra_dst = th.cat([
            th.arange(n_nodes, 4 * n_nodes, dtype=th.int64),
            th.cat([
                th.arange(0 * n_nodes, 1 * n_nodes, dtype=th.int64),
                th.arange(2 * n_nodes, 4 * n_nodes, dtype=th.int64)
            ]),
            th.cat([
                th.arange(0 * n_nodes, 2 * n_nodes, dtype=th.int64),
                th.arange(3 * n_nodes, 4 * n_nodes, dtype=th.int64)
            ]),
            th.arange(3 * n_nodes, dtype=th.int64)
        ])
        print(intra_src.shape)
        src = th.cat([
            original_src, original_src + n_nodes, original_src + 2 * n_nodes,
            original_src + 3 * n_nodes, intra_src
        ])
        dst = th.cat([
            original_dst, original_dst + n_nodes, original_dst + 2 * n_nodes,
            original_dst + 3 * n_nodes, intra_dst
        ])
        print(src.shape, dst.shape)
        graph_coo = dgl.graph((src, dst))
        del src, dst, intra_src, intra_dst
        out_degrees = graph_coo.out_degrees()
        print("Generating csc formats")
        indptr = graph_coo.adj_sparse('csc')[0]
        indices = graph_coo.adj_sparse('csc')[1]
        graph = dgl.graph(('csc', (indptr, indices, [])))
        graph.ndata['out_degrees'] = out_degrees
        del graph_coo
        np.save(os.path.join(root, 'papers400M_features.npy'),
                original_features)
        graph.ndata['labels'] = th.cat([
            original_labels, original_labels, original_labels, original_labels
        ],
                                       dim=0)

        num_labels = len(th.unique(labels[th.logical_not(th.isnan(labels))]))

        # Find the node IDs in the training, validation, and test set.
        train_nid, val_nid, test_nid = splitted_idx['train'], splitted_idx[
            'valid'], splitted_idx['test']
        train_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
        train_mask[train_nid] = True
        train_mask[train_nid + n_nodes] = True
        train_mask[train_nid + 2 * n_nodes] = True
        train_mask[train_nid + 3 * n_nodes] = True
        val_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
        val_mask[val_nid] = True
        val_mask[val_nid + n_nodes] = True
        val_mask[val_nid + 2 * n_nodes] = True
        val_mask[val_nid + 3 * n_nodes] = True
        test_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
        test_mask[test_nid] = True
        test_mask[test_nid + n_nodes] = True
        test_mask[test_nid + 2 * n_nodes] = True
        test_mask[test_nid + 3 * n_nodes] = True
        graph.ndata['train_mask'] = train_mask
        graph.ndata['val_mask'] = val_mask
        graph.ndata['test_mask'] = test_mask

        dgl.save_graphs(os.path.join(root, 'papers400M.dgl'), [graph])
        graph.ndata['features'] = th.cat([
            original_features, original_features, original_features,
            original_features
        ],
                                         dim=0)
        del original_features, original_labels

        print('finish constructing papers400m', graph)

    return graph, num_labels


def load_ogb(name, root="dataset"):
    from ogb.nodeproppred import DglNodePropPredDataset

    print('load', name)
    data = DglNodePropPredDataset(name=name, root=root)
    print('finish loading', name)
    splitted_idx = data.get_idx_split()
    graph, labels = data[0]
    labels = labels[:, 0]

    graph.ndata['features'] = graph.ndata.pop('feat')
    graph.ndata['labels'] = labels

    graph.ndata['out_degrees'] = graph.out_degrees()

    num_labels = len(th.unique(labels[th.logical_not(th.isnan(labels))]))

    # Find the node IDs in the training, validation, and test set.
    train_nid, val_nid, test_nid = splitted_idx['train'], splitted_idx[
        'valid'], splitted_idx['test']
    train_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
    train_mask[train_nid] = True
    val_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
    val_mask[val_nid] = True
    test_mask = th.zeros((graph.number_of_nodes(), ), dtype=th.bool)
    test_mask[test_nid] = True
    graph.ndata['train_mask'] = train_mask
    graph.ndata['val_mask'] = val_mask
    graph.ndata['test_mask'] = test_mask
    print('finish constructing', name)
    return graph, num_labels


def inductive_split(g):
    """Split the graph into training graph, validation graph, and test graph by training
    and validation masks.  Suitable for inductive models."""
    train_g = g.subgraph(g.ndata['train_mask'])
    val_g = g.subgraph(g.ndata['train_mask'] | g.ndata['val_mask'])
    test_g = g
    return train_g, val_g, test_g

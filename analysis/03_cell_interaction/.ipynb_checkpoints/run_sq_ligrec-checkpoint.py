import scanpy as sc
import squidpy as sq
import pickle

adata = sc.read("./data/adata_for_sqligrec.h5ad")
res= sq.gr.ligrec(
    adata,
    n_perms=10000,
    cluster_key="celltype",
    copy=True,
    use_raw=False,
    corr_method = "fdr_bh",
    threshold = 0.3, 
    alpha = 0.05,
    transmitter_params={"categories": "ligand"},
    receiver_params={"categories": "receptor"},
)


with open("./out/sq_ligrec_output.pickle", 'wb') as out:
    pickle.dump(res, out, pickle.HIGHEST_PROTOCOL)
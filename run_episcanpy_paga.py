import os
import os.path
import random
import sys
import tempfile
import functools as f
from contextlib import contextmanager, redirect_stderr, redirect_stdout

import altair as alt
import numpy as np
import pandas as pd
import scipy
import scipy.io
import scipy.sparse
import episcanpy as epi
import scanpy as sc
import anndata
import gzip

import umap

sys.path.append("/home/gw/code/utility/")

if True:  # In order to bypass isort when saving
    import altairThemes


# From https://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
# Ignore std output from a function.
@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)


random.seed(0)
np.random.seed(0)

labelFile = sys.argv[1]
output = sys.argv[2]
paths = sys.argv[3:]

# register the custom theme under a chosen name
alt.themes.register("publishTheme", altairThemes.publishTheme)

# enable the newly registered theme
alt.themes.enable("publishTheme")

def decompressFile(path1, path2):
    """Decompress a file from path1 to path2."""

    with gzip.open(path1, "rb") as f1:
        with open(path2, "wb") as f2:
            contents = f1.read()
            f2.write(contents)

    return()

def joinMats(dfs):
    df = f.reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True, how="outer"), dfs)

    item = list(df.columns)
    feature = list(df.index.values)

    mat = scipy.sparse.csr_matrix(df.values)

    annData = anndata.AnnData(mat.transpose(), obs=pd.DataFrame(index=item), var=pd.DataFrame(index=feature))
    annData.uns["omic"] = "ATAC"

    return(annData)

# Load matrix
def read10X(path):
    """Read a 10X matrix from a folder."""

    mat = os.path.join(path, "matrix.mtx")
    items = os.path.join(path, "barcodes.tsv")
    features = os.path.join(path, "genes.tsv")

    if os.path.isfile(os.path.join(path, "matrix.mtx.gz")):
        matTemp = tempfile.NamedTemporaryFile()
        itemsTemp = tempfile.NamedTemporaryFile()
        featuresTemp = tempfile.NamedTemporaryFile()

        decompressFile(mat + ".gz", matTemp.name)
        decompressFile(items + ".gz", itemsTemp.name)
        decompressFile(os.path.join(path, "features.tsv.gz"), featuresTemp.name)

        mat = matTemp.name
        items = itemsTemp.name
        features = featuresTemp.name

    # annData = epi.pp.read_ATAC_10x(mat, cell_names=items, var_names=features) # Does not work with large matrices

    # return(annData)

    # Instead:

    m = scipy.io.mmread(mat).tocsr()

    with open(items) as fi:
        i = fi.readlines()
        i = [x[:-1] for x in i]
    with open(features) as fi:
        f = fi.readlines()
        f = [x.split('\t')[0] for x in f]

    p = pd.DataFrame.sparse.from_spmatrix(m)
    p.columns = i
    p.index = f

    return(p)

# annData = anndata.AnnData.concatenate(*map(read10X, paths), join="inner")
annData = joinMats(list(map(read10X, paths)))
annData.obs_names_make_unique()

# Preprocessing
epi.pp.filter_cells(annData, min_features=10)
epi.pp.filter_features(annData, min_cells=10)

# Feature selection
epi.pl.cal_var(annData)
epi.pp.select_var_feature(annData)

# Reduction and neighbors
epi.pp.pca(annData)
epi.pp.neighbors(annData)
# If the above fails due to nearest neighbor issue:
# epi.pp.neighbors(annData, knn=False)
# annData.obsp["connectivities"]  = scipy.sparse.csr_matrix(annData.obsp["connectivities"])
# annData.obsp["distances"]  = scipy.sparse.csr_matrix(annData.obsp["distances"])
# If too small clusters due to homogeneity from features:
# sc.tl.diffmap(annData)
# sc.pp.neighbors(annData, use_rep='X_diffmap')
epi.tl.pca(annData)
epi.tl.umap(annData)
sc.tl.louvain(annData)
sc.tl.paga(annData, groups="louvain")

# Label matrix
almost = np.column_stack([annData.obs["louvain"], annData.obs_names])
colnames = np.array(["cluster", "item"])
final = np.row_stack([colnames, almost])

# Output results
np.savetxt(sys.stdout.buffer, final, delimiter=",", fmt="%s")

# Label matrix
labelsDf = pd.read_csv(labelFile)
labelsDf = labelsDf.astype("category")
itemsDf = pd.DataFrame({"item": annData.obs_names})
itemsDf = pd.merge(itemsDf, labelsDf, on = "item", how = "left", sort = False)
itemsDf = itemsDf.set_index("item")

# Sort labels by obs (https://stackoverflow.com/questions/23482668/sorting-by-a-custom-list-in-pandas)
def sorter(column):
    """Sort function"""
    items = annData.obs_names
    correspondence = {item: order for order, item in enumerate(items)}
    return column.map(correspondence)

itemsDf.sort_values(by="item", key=sorter)

annData.obs["label"] = itemsDf["label"]

# Set colors
colors = ["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"]
colorLabels = sorted(set(itemsDf["label"]))
if len(colorLabels) > 9:
    newColors = [x.hexcode for x in spectra.range(colors, len(colorLabels))]
else:
    newColors = colors[0:len(colorLabels)]
annData.uns["label_colors"] = newColors

# Plot results

#Create output directory
if not os.path.exists("./episcanpy_img"):
    os.makedirs("./episcanpy_img")
if not os.path.exists("./episcanpy_paga_img"):
    os.makedirs("./episcanpy_paga_img")
if not os.path.exists("figures/paga/"):
    os.makedirs("figures/paga/")
if not os.path.exists("figures/draw_graph_fr/"):
    os.makedirs("figures/draw_graph_fr/")

sc.pl.paga(annData, color=["louvain", "label"], cmap=newColors, save = "/" + output + "_episcanpy_cluster.svg")
sc.tl.umap(annData, init_pos="paga")
sc.tl.draw_graph(annData, init_pos="paga", layout="fr")
sc.pl.draw_graph(annData, color=["louvain", "label"], save = "/" + output + "_episcanpy_cluster.svg")

df = pd.DataFrame({
    "item": annData.obs_names,
    "PAGA 1": annData.obsm["X_draw_graph_fr"][:, 0],
    "PAGA 2": annData.obsm["X_draw_graph_fr"][:, 1],
    "cluster": annData.obs["louvain"]
})
df = pd.merge(itemsDf, df, on = "item", how = "inner", sort = False)

# Plot results
clusterChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="PAGA 1", y="PAGA 2", color=alt.Color("cluster:N",
                                            scale=alt.Scale(scheme="rainbow"),
                                            sort=list(map(lambda x: str(x),
                                                          sorted(map(lambda x:
                                                                     int(x),
                                                                     list(set(df["cluster"])))))),
                                            legend=alt.Legend(columns=2, symbolLimit=len(df)))).properties(
        width=180, height=180)
clusterChart.save("./episcanpy_img/" + output + "_paga_cluster.svg", webdriver="firefox")

pagaLabelChart = alt.Chart(df).mark_circle(opacity=1).encode(
    x="PAGA 1",
    y="PAGA 2",
    color=alt.Color("label",
        scale=alt.Scale(interpolate="rgb", range=["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"])
                    , legend=alt.Legend(columns=2, symbolLimit=len(df))
    )).properties(
                                                 width=180, height=180)
pagaLabelChart.save("./episcanpy_img/" + output + "_paga_label.svg", webdriver="firefox")
pagaLabelChart.save("./episcanpy_img/" + output + "_paga_label.html")

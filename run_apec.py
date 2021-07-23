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

from APEC import clustering,plot,generate,convert
import tempfile as tmp
from distutils.dir_util import copy_tree
import subprocess
import glob
import gzip

import umap

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

def decompressFile(path1, path2):
    """Decompress a file from path1 to path2."""

    with gzip.open(path1, "rb") as f1:
        with open(path2, "wb") as f2:
            contents = f1.read()
            f2.write(contents)

    return()

if __name__ == "__main__":

    random.seed(0)
    np.random.seed(0)

    labelFile = sys.argv[1]
    output = sys.argv[2]
    paths = sys.argv[3:]

    # register the custom theme under a chosen name
    alt.themes.register("publishTheme", altairThemes.publishTheme)

    # enable the newly registered theme
    alt.themes.enable("publishTheme")

    # Make copy
    tmpdir = tmp.mktemp("apec")
    os.mkdir(tmpdir)
    projectDir = os.path.join(tmpdir, "project")
    os.mkdir(projectDir)

    # Matrix join for copy
    pathArgs = list(pd.core.common.flatten([["-m", i] for i in paths]))
    args = ["too-many-cells", "matrix-output", "--mat-output", tmpdir,
            "--normalization", "NoneNorm", "--filter-thresholds", "(1,1)"] + pathArgs
    subprocess.run(args)
    list(map(lambda x: decompressFile(x, os.path.splitext(x)[0]), glob.glob(os.path.join(tmpdir, "*.gz"))))

    # Bed join for copy
    peakFileH = open(os.path.join(tmpdir, "peaks.bed"), "w")
    p1 = subprocess.Popen(["cat", os.path.join(tmpdir, "features.tsv")], stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["sed", "-e", "s/:/\t/g", "-e", "s/-/\t/g"], stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()
    p3 = subprocess.Popen(["cut", "-f1,2,3"], stdin=p2.stdout, stdout=peakFileH)
    p2.stdout.close()
    peakFileH.close()

    with suppress_stdout_stderr():
      # Load into APEC
      convert.convert_10X(tmpdir, projectDir)

      # Cluster
      clustering.build_accesson(projectDir)
      clustering.cluster_byAccesson(projectDir)
      plot.plot_tsne(projectDir, rs=0)

    # Input cluster results
    clusterDf = pd.read_csv(os.path.join(projectDir, "result",
                                         "cluster_by_APEC.csv"), sep = "\t")

    # Input TSNE results
    tsneDf = pd.read_csv(os.path.join(projectDir, "result",
                                      "TSNE_by_APEC.csv"), sep = "\t")

    # Output cluster results
    clusterDf = clusterDf.rename(columns={clusterDf.columns[0]: "item"})
    clusterDf.to_csv(sys.stdout, index = False)

    # Label matrix
    labelsDf = pd.read_csv(labelFile)
    itemsDf = pd.DataFrame({"item": clusterDf["item"]})
    itemsDf = pd.merge(itemsDf, labelsDf, on = "item", how = "left", sort = False)

    # Plot results
    df = pd.DataFrame({
        "item": clusterDf.iloc[:,0],
        "TSNE 1": tsneDf["TSNE1"],
        "TSNE 2": tsneDf["TSNE2"],
        "cluster": clusterDf["cluster"]
    })
    df = pd.merge(itemsDf, df, on = "item", how = "inner", sort = False)

    # # Create output directory
    if not os.path.exists("./apec_img"):
        os.makedirs("./apec_img")

    # Plot results
    clusterChart = alt.Chart(df).mark_circle(opacity=1).encode(
        x="TSNE 1", y="TSNE 2", color=alt.Color("cluster:N", scale=alt.Scale(scheme="rainbow"))).properties(
            width=180, height=180)
    clusterChart.save("./apec_img/" + output + "_apec_cluster.svg", webdriver="firefox")

    umapLabelChart = alt.Chart(df).mark_circle(opacity=1).encode(
        x="TSNE 1",
        y="TSNE 2",
        color=alt.Color("label",
            scale=alt.Scale(interpolate="rgb", range=["#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999"])
                        , legend=alt.Legend(columns=2, symbolLimit=len(df))
        )).properties(
                                                    width=180, height=180)
    umapLabelChart.save("./apec_img/" + output + "_tsne_label.svg", webdriver="firefox")
    umapLabelChart.save("./apec_img/" + output + "_tsne_label.html")

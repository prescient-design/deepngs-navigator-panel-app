import os
import pandas as pd
import numpy as np
from collections import Counter
from sklearn.decomposition import PCA
import matplotlib as mpl
from io import StringIO
from bokeh.models import HoverTool
import holoviews as hv
from holoviews import opts

# Load Spectral palette
from bokeh.models import LinearColorMapper
from bokeh.palettes import all_palettes

def get_color_schemes() -> dict[str, dict[str, str]]:
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    COLOR_SCHEMES_FILE = os.path.join(BASE_DIR, "utils", "COLOR_SCHEMES.tsv")
    df = pd.read_csv(COLOR_SCHEMES_FILE, delimiter="\t")

    name2color_scheme = {}
    for _, row in df.iterrows():
        name = row.iloc[0]
        color_scheme = {letter: color for letter, color in zip(df.columns[1:], row[1:])}
        name2color_scheme[name] = color_scheme

    return name2color_scheme

COLOR_SCHEMES = get_color_schemes()
colors = list(set(COLOR_SCHEMES["Clustal"].values()))
cmap1 = mpl.colors.ListedColormap(colors)
levels = [i for i, x in enumerate(colors)]
levels.append(len(colors))
colors_dict = dict(zip(colors, levels))




def hide_hook(plot, element):
    plot.handles["xaxis"].visible = False
    plot.handles["yaxis"].visible = False
    plot.handles["plot"].border_fill_color = None
    plot.handles["plot"].background_fill_color = None
    plot.handles["plot"].outline_line_color = None


def plot_consensus_hv2(seq_list, id_list, consensus, alignment_length, msa_count):
    
    data3 = []
    data4 = []
    for j, letter in enumerate(consensus):
        data3.append((j, 0, COLOR_SCHEMES["Clustal"][letter]))
        data4.append((j, 0, letter))
    data3 = pd.DataFrame(data3, columns=["x", "y", "val"])
    rects_data = [
        {"x0": x - 0.5, "y0": -5, "x1": x + 0.5, "y1": 5, "val": val}
        for i, (x, y, val) in data3.iterrows()
    ]

    heatmap2 = hv.Rectangles(rects_data, vdims="val").opts(color="val")

    labels2 = hv.Labels(data4)
    layout = heatmap2 * labels2
    layout.opts(
        opts.Labels(
            hooks=[hide_hook],
            xlabel="",
            ylabel="",
            text_color="black",
            fontsize=5,
            height=50,
            width=alignment_length * 30,
            xaxis=None,
            toolbar="disable",
        ),
        opts.HeatMap(
            cmap=cmap1,
            alpha=0.5,
            hooks=[hide_hook],
            xlabel="",
            ylabel="",
            xlim=(-0.5, alignment_length + 0.5),
            ylim=(-5, 5),
            height=50,
            width=alignment_length * 30,
            xaxis=None,
            toolbar="disable",
        ),
    )
    return layout


def plot_msa_hv(seq_list, id_list, consensus, alignment_length, msa_count):

    # Highlight differences in the viewer
    sequences = [
        "".join([("." if aa == cons_aa else aa) for aa, cons_aa in zip(seq, consensus)])
        for seq in seq_list
    ]
    
    # Create data for HeatMap
    data = []
    data2 = []

    for i, sequence in enumerate(sequences):
        for j, letter in enumerate(sequence):
            color = colors_dict[COLOR_SCHEMES["Clustal"][letter]]
            data.append((j, i, color))
            data2.append((j, i, letter))
    data = pd.DataFrame(data, columns=["x", "y", "val"])
    hm_opts = dict(kdims=["x", "y"], vdims=["val"])
    heatmap = hv.HeatMap(data, **hm_opts)
    if alignment_length<30:
        l=100
    else:
        l=10

    heatmap = heatmap.options(
        cmap=cmap1,
        alpha=0.5,
        xlabel="Positions",
        ylabel="",
        xlim=(-1, len(sequence)),
        ylim=(-2, len(sequences)),
        width=alignment_length * 30,
        height=msa_count * 14,
    )
    labels = hv.Labels(data2).opts(
        text_color="black",
        fontsize=7,
        ylabel="",
        xlim=(-1, alignment_length),
        ylim=(-2, msa_count),
    )
    layout2 = heatmap * labels
    layout2.opts(
        opts.HeatMap(
            height=int(msa_count * 14),
            width=alignment_length * 30,
            xaxis=None,
            toolbar="disable",
        ),
        opts.Labels(
            height=msa_count * 14,
            width=alignment_length * 30,
            tools=["hover"],
            toolbar="disable",
        ),
    )
    return layout2


def find_consensus(data2):
    strings=data2.AA.values
    transposed = zip(*strings)
    consensus = "".join(Counter(column).most_common(1)[0][0] for column in transposed)
    return consensus


def plot_msa(data2):

    data2 = data2.reset_index(drop=True)
    # Combine 'e1' and 'e2' into a single array for PCA
    dff = data2[["e1", "e2"]].values
    if dff.shape[0] > 10:
        # Perform PCA to find the principal axis (line of best fit)
        pca = PCA(n_components=2)
        pca.fit(dff)

        # Get the principal components (directions of maximum variance)
        principal_components = pca.components_

        # Sort the DataFrame based on the dot product with the principal axis
        data2["dot_product"] = data2[["e1", "e2"]].dot(principal_components.T[0])
        data3 = data2.sample(min(100, data2.shape[0])).sort_values(by="dot_product")
        mv = plot_msa_hv(
            data3.AA.values,
            data3.seq_id.values,
            find_consensus(data2),
            len(data3.AA.iloc[0]),
            data3.shape[0],
        )
    else:
        mv = plot_msa_hv(
            data2.AA.values,
            data2.seq_id.values,
            find_consensus(data2),
            len(data2.AA.iloc[0]),
            data2.shape[0],
        )
    return mv

def read_processed_files():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    INPUT_FILELIST = os.path.join(BASE_DIR, "utils", "processed_files_report.csv")
    processed_files_report=pd.read_csv(INPUT_FILELIST)
    processed_files_report.columns=['name','path']
    return processed_files_report

def plot_consensus(data2):

    data2 = data2.reset_index(drop=True)
    mv = plot_consensus_hv2(
        data2.AA.values,
        data2.seq_id.values,
        find_consensus(data2),
        len(data2.AA.iloc[0]),
        data2.shape[0],
    )
    return mv


# function to get data
def get_data_(data=None):
    global df
    processed_files_report=read_processed_files()

    if data:
        path_file = processed_files_report[processed_files_report['name']==data]['path'].values[0]
    else:
        path_file = processed_files_report[processed_files_report['name']=='default:All']['path'].values[0]
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    path_file = os.path.join(BASE_DIR, path_file)
    df = pd.read_csv(path_file, sep="\t")

    df = df.reset_index(drop=True)
    df = df.reset_index()

    df["seq_id"] = df.index.values
    
    if "picked_clone" not in df.columns:
        df["picked_clone"] = " "
    else:
        df["picked_clone"] = df["picked_clone"].fillna(" ")
        
    if "AA" not in df.columns:
        if "seq" in df.columns:
            df["AA"] = df["seq"]  # Assign 'seq' column to 'AA' if it exists
        else:
            df["AA"] = " "  # Default value if neither

        
    try:
        df = df.drop(["Unnamed: 0"], axis=1)
    except Exception:
        pass
    df['default size']=1

    
    df['default color']='1'


    df.reset_index(drop=True, inplace=True)
    df = df.sort_values("picked_clone", ascending=False)
    if df.shape[0]>50000:
        df1=df[~df['picked_clone'].isna()]
        df2=df.sample(50000)
        df=pd.concat([df1,df2])
    return df




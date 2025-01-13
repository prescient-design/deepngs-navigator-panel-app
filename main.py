import panel as pn
import holoviews as hv
from functions import *


# 'Spectral' palette with 256 colors
Spectral256 = all_palettes['Spectral'][11]

# Load Bokeh plotting extension
pn.extension(
    "bokeh",
    "matplotlib",
    design="material",
    template="material",
    loading_indicator=True,
)
hv.extension("bokeh", "matplotlib")

color_options_global = [
    "default color",
    "HCDR3",
    "target",
    "total_count",
    "animal",
    "charge_pH7_cdr3",
    "charge_pH7_fullSeq",
    "hydrophobicity_cdr3",
    "hydrophobicity_fullSeq",
    "v-shm",
    "log10_read_count"
]
size_options_global = ["default size","log10_num_neighbors", "total_count","None"]

# ____________define panel placeholders____________
blank_hv = hv.Scatter((1, 1))
graph_pane = pn.pane.HoloViews(blank_hv, visible=False)
graph_pane2 = pn.pane.HoloViews(blank_hv, visible=False)
selection_pane1 = pn.pane.HoloViews(blank_hv, visible=False)
selection_pane2 = pn.pane.HoloViews(blank_hv, visible=False)
selection_pane3 = pn.pane.HoloViews(blank_hv, visible=False)
download_pane = pn.widgets.FileDownload(
    button_type="default",
    auto=False,
    embed=False,
    name="double-click to download using 'Save as' dialog",
    filename="selected_seqeunces.csv",
)
selector_c = pn.widgets.Select(name="Select coloring property")
selector_s = pn.widgets.Select(name="Select sizing property")
df = pd.DataFrame()


processed_files_report=read_processed_files()
projects=list(processed_files_report.name.str.split(':').str[0].drop_duplicates().values)


# Define the project selector
project_selector = pn.widgets.Select(name='Project Browser', options=projects, value='')
# Initially, the sequence selector is not applicable, so it's empty or set to a default value
sequence_selector = pn.widgets.Select(name='group:Vgene-Length-cdr3Length', options=[''], value='')

# Function to update the sequence selector based on the project selected
def update_sequence_options(event):
    project = event.obj.value
    sequence_selector.options=list(processed_files_report[processed_files_report.name.str.split(':').str[0]==project].name.str.split(':').str[1].values)
    sequence_selector.value = sequence_selector.options[0]

# Watch for changes in the project selector to update the sequence selector
project_selector.param.watch(update_sequence_options, 'value')

# Function to handle selection, triggering when the sequence is selected
def handle_selection(event=None):
    # Construct the identifier for the selected project and sequence
    selected_project = project_selector.value
    selected_sequence = sequence_selector.value
    selected_item = f"{selected_project}:{selected_sequence}"
    print(f"Processing data for: {selected_item}")
    global df
    df = get_data_(data=selected_item)
    options=df.loc[:, df.nunique() > 0].columns
    selector_c.options = [c for c in color_options_global if c in options]
    selector_c.options.extend([x for x in options if 'enrichment' in x])
    selector_s.options = [c for c in size_options_global if c in options]
    if event:
        return f'Your selected project: "{selected_project}", number of seqeunces: {df.shape[0]}'

    else:

        return f"Currently no project is selected, test, is loaded (as default).\nnumber of seqeunces: {df.shape[0]}"


# Watch for changes
sequence_selector.param.watch(handle_selection, 'value')

# Organize the widgets into a layout
file_input = pn.Column(
    project_selector,
    sequence_selector,
    height=150,
)



text = """
#  Project selection
Click on the Projects button to select a ready project
"""

explanation = "Select properties for point size and point color"
use_markers = "experimentally confirmed + selected for screening"

explanation2 = "Lasso/box select in the interactive scatter map to see MSA of selected sequences and download selection"


selector_c = pn.widgets.Select(
    name="Select coloring property",
    options=[c for c in color_options_global if c in df.columns],
)
selector_s = pn.widgets.Select(
    name="Select sizing property",
    options=[c for c in size_options_global if c in df.columns ],
)

def display_scatter_plot(event=None):
    global df
    global download_pane
    # print(selector_c.value)
    if event:
        df["color_"] = df[selector_c.value]
        if pd.api.types.is_numeric_dtype(df["color_"]):
            df["color_"].fillna(df['color_'].min()-1,inplace=True)
            df['normalized_index'] = (df['color_'] - df['color_'].min()) / (df['color_'].max() - df['color_'].min()) * (len(Spectral256) - 1)
            # Assign color codes based on normalized index
            df['color'] = df['normalized_index'].apply(lambda x: Spectral256[int(x)])
            color_mapper = LinearColorMapper(palette=Spectral256, low=df['color_'].min(), high=df['color_'].max())

        else:
            top_colors = df["color_"].value_counts().head(29).index.tolist()
            # Categorize the 'color' column and set the rest to gray
            colors = [
                "#1c8c99",
                "#e4c3ab",
                "#707a3f",
                "#a68436",
                "#1a3531",
                "#f7a4e9",
                "#6e1108",
                "#18bbcc",
                "#6788d5",
                "#460ad8",
                "#9db8ec",
                "#232689",
                "#59100c",
                "#835902",
                "#4e8832",
                "#d76979",
                "#cd9ca1",
                "#22fdb3",
                "#faf5b6",
                "#e19da2",
                "#0cd702",
                "#9f20fd",
                "#99dcb7",
                "#7c0e36",
                "#07417f",
                "#24a9bc",
                "#a50b19",
                "#5ef73d",
                "#878414",
            ]
            color_mapping = {color: colors[i] for i, color in enumerate(top_colors)}
            color_mapping["gray"] = "#bababa"  # Set the color for 'gray' category
            df["color"] = df["color_"].apply(lambda x: x if x in top_colors else "gray")

        # normalize size property
        size_ = list(
            ((df[selector_s.value].values + 0.1) * 0.7 / (df[selector_s.value].max() + 0.1) + 0.5)
        )

        df["size_prop"] = [x * max(1.0, 10 - 2 * np.log10(df.shape[0] + 1)) for x in size_]

        # lable points
        df["symbol"] = "circle"
        labeled_df = df[(df["picked_clone"] != " ")]
        labeled_df = labeled_df[labeled_df["picked_clone"].groupby(labeled_df["picked_clone"]).transform("count")<= 350]
        # Separate labeled points based on category count
        if labeled_df.shape[0] > 0:
            large_category_df = labeled_df[
                labeled_df["picked_clone"].groupby(labeled_df["picked_clone"]).transform("count")
                > 150
            ][["e1", "e2", "picked_clone"]]
            small_category_df = labeled_df[
                labeled_df["picked_clone"].groupby(labeled_df["picked_clone"]).transform("count")
                <= 150
            ][["e1", "e2", "picked_clone"]]
            large_category_df["symbol"] = "triangle"
            symbols = [ 'triangle', 'square', 'cross', 'star','diamond', 'x', 'inverted_triangle', 'hexagon',  'circle_cross']

            # Create a mapping from categories to symbols
            category_to_symbol = {cat: sym for cat, sym in zip(small_category_df['picked_clone'].unique(), symbols)}
            print('category_to_symbol',category_to_symbol)
            # Assign symbols to a new column in the DataFrame based on the picked_clone category
            small_category_df['symbol'] = small_category_df['picked_clone'].map(category_to_symbol)

            merged_df = pd.merge(
                df,
                large_category_df,
                on=["e1", "e2"],
                how="outer",
                suffixes=("_df1", "_df2"),
            )
            merged_df = pd.merge(merged_df, small_category_df, on=["e1", "e2"], how="outer")
            # Use value from df3 if available, else use value from df2, else use value from df1
            merged_df["symbol"] = (
                merged_df["symbol"]
                .combine_first(merged_df["symbol_df2"])
                .combine_first(merged_df["symbol_df1"])
            )

            # Drop unnecessary columns
            merged_df.drop(["symbol_df2", "symbol_df1"], axis=1, inplace=True)
            merged_df["symbol_c"] = merged_df.groupby(["symbol"])["symbol"].transform("count")

        else:
            merged_df = df
            merged_df["symbol_c"] = merged_df.groupby(["symbol"])["symbol"].transform("count")

        merged_df.sort_values(["symbol_c"], ascending=False, inplace=True)
        merged_df["symbol_c"] = [
            0.01 if c == "circle" else 3 if c == "triangle" else 4
            for c in merged_df["symbol"].values
        ]
        merged_df["size_prop"] = merged_df["size_prop"] * 2 + 2 * merged_df["symbol_c"]
        merged_df["e1"] = merged_df["e1"].apply(lambda x: round(x, 5))
        merged_df["e2"] = merged_df["e2"].apply(lambda x: round(x, 5))

        # ___________plot scatter and collected selections____________
        def structure_df_plot(dff):
            sorter = list(dff["color"].value_counts().index.values)
            dff.sort_values(
                by="color", key=lambda column: column.map(lambda e: sorter.index(e)), inplace=True
            )
            d_s = {
                "e1": dff.e1.values,
                "e2": dff.e2.values,
                "color": dff["color"].values,
                "color_": dff["color_"].values,
                "size_prop": dff["size_prop"].values,
                "symbol": dff["symbol"].values,
                "symbol_c": dff["symbol_c"].values,
            }

            return d_s

        hover = HoverTool(tooltips=[("color_", "@color_")], mode="mouse")
        if pd.api.types.is_numeric_dtype(merged_df['color_']):
            merged_df=merged_df.sort_values(by='color_',ascending=False)
            points = hv.Points(
                structure_df_plot(merged_df),
                kdims=["e1", "e2"],
                vdims=["color", "color_", "size_prop", "symbol", "symbol_c"],
                ).opts(
                    color=hv.dim('color'),  # Applying coloring based on 'color_' value
                    cmap=Spectral256,  # Continuous color map
                    size=hv.dim('size_prop'),
                    marker='symbol',
                    line_width=hv.dim('symbol_c'),
                    line_color='black',
                    tools=[hover, 'box_select', 'lasso_select'],
                    width=800,
                    height=600,
                    legend_position='bottom'
                )

        else:
            points = hv.Points(
                structure_df_plot(merged_df),
                kdims=["e1", "e2"],
                vdims=["color", "color_", "size_prop", "symbol", "symbol_c"],
            ).opts(
                color=hv.dim("color").categorize(color_mapping),
                size="size_prop",
                marker="symbol",
                line_width="symbol_c",
                line_color="black",
                tools=[hover, "box_select", "lasso_select"],
                width=800,
                height=600,
                legend_position="bottom",
            )
        marker_df = merged_df[["picked_clone", "symbol"]].value_counts().reset_index().reset_index()

        marker_df = marker_df[marker_df["symbol"] != "circle"]
        points2 = hv.Points(
            marker_df,
            kdims=["count", "picked_clone"],
            vdims=["picked_clone", "symbol"],
            label="marker legend",
        ).opts(
            fontsize={"title": 10, "labels": 8, "xticks": 5, "yticks": 6},
            marker="symbol",
            line_color="black",
            tools=["hover"],
            size=10,
            width=300,
            height=600,
            legend_position="bottom",
            toolbar="disable",
        )

        selection = hv.streams.Selection1D(source=points)

        @pn.depends(selection.param.index)
        def update_scatter(selected_indices):
            if selected_indices:

                selected_scatter = plot_msa(
                    merged_df.iloc[selected_indices]
                )  # hv.DynamicMap(selected_info, streams=[selection]).opts(shared_axes=False,width=600)
                return selected_scatter
            else:
                # print(selected_indices)
                return hv.Div("No sequence selected")

        @pn.depends(selection.param.index)
        def update_scatter2(selected_indices):
            if selected_indices:

                selected_scatter = plot_consensus(
                    merged_df.iloc[selected_indices]
                )  # hv.DynamicMap(selected_info, streams=[selection]).opts(shared_axes=False,width=600)

                return selected_scatter
            else:
                # print(selected_indices)
                return hv.Div("No sequence selected")

        @pn.depends(selection.param.index)
        def update_scatter3(selected_indices):
            global download_pane
            if selected_indices:

                table = hv.Table(
                    merged_df.iloc[selected_indices],
                )
                table.opts(width=1200)

                sio = StringIO()
                merged_df.iloc[selected_indices].to_csv(sio)
                sio.seek(0)

                download_pane.file = sio

                return table
            else:
                download_pane.file = None
                return hv.Div("No sequence selected")

        selection_pane3.object = update_scatter
        selection_pane3.visible = True
        selection_pane2.object = update_scatter2
        selection_pane2.visible = True
        selection_pane1.object = update_scatter3
        selection_pane1.visible = True
        graph_pane.object = points
        graph_pane.visible = True
        graph_pane2.object = points2
        graph_pane2.visible = True
        
        
display_graph_but = pn.widgets.Button(name="display map", button_type="success")
display_graph_but.on_click(display_scatter_plot)

sidebar = pn.layout.WidgetBox(
    pn.pane.Markdown(text, margin=(0, 10)),
    file_input,
    explanation,
    selector_c,
    selector_s,
    display_graph_but,
    explanation2,
    max_width=350,
    sizing_mode="stretch_width",
)
static_text = pn.widgets.StaticText(
    name="2D projection", value="select some sequences to explore their MSA"
)


pn.Column(
    pn.Row(sidebar, graph_pane, graph_pane2),
    pn.Column(selection_pane2, selection_pane3, download_pane, selection_pane1),
).servable(area="main")

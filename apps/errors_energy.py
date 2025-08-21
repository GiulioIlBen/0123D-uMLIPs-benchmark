# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "ase==3.25.0",
#     "jinja2>=3.1.6",
#     "marimo>=0.13.8",
#     "matplotlib==3.10.3",
#     "pandas==2.3.0",
#     "plotly==6.2.0",
#     "seaborn==0.13.2",
# ]
# ///

import marimo

__generated_with = "0.14.8"
app = marimo.App(width="full")

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import matplotlib.pyplot as plt
    import pandas as pd
    import plotly.express as px
    import seaborn as sns
    from ase.data import atomic_numbers
    from explore_0123D_test_set import get_table_view, make_rotation_string, plot_agms


@app.cell
def _():
    df_true = get_table_view(lambda row: {"energy [eV/atom]": row.uncorrected_energy/ row.natoms, "formula":row.formula})
    return (df_true,)


@app.cell
def _():
    # df_pred = get_table_view(lambda row: {"energy [eV/atom]": row.uncorrected_energy/ row.natoms, "uMLIPs":row.id_theory},
    #                         relative_path="public/0123D_predictions.db")
    df_pred = pd.read_csv(mo.notebook_dir()/"public/df_pred_energy.csv")
    return (df_pred,)


@app.cell
def _(df_pred, df_true):
    df_deltas = (df_pred
        .merge(df_true, on=["mat_id", "id_dimension"], validate="m:1", suffixes=("", "_PBE"))
        .assign(delta_EnMLIP_minus_EnPBE_eVatoms = lambda f: f["energy [eV/atom]"] - f["energy [eV/atom]_PBE"])
    )
    return (df_deltas,)


@app.cell
def _():
    mo.md(r"""# Overview""")
    return


@app.cell
def _(df_deltas):
    flierprops = dict(marker=".", markerfacecolor="black", markersize=1, markeredgecolor="none")
    whiskerprops = dict(linestyle="-", linewidth=2)
    boxprops = dict(alpha=0.5, edgecolor=None)
    font_size_all = 20
    plt.rcParams.update({
        'figure.figsize':(20,20),
        "legend.title_fontsize": font_size_all,
        "font.size": font_size_all,
        "xtick.labelsize": font_size_all,
        "ytick.labelsize": font_size_all,
        "axes.labelsize": font_size_all,
    })
    ax = sns.boxplot(df_deltas.assign(delta_EnMLIP_minus_EnPBE_meVatoms=lambda f: f["delta_EnMLIP_minus_EnPBE_eVatoms"] * 1000), x="delta_EnMLIP_minus_EnPBE_meVatoms", y="uMLIPs", hue="id_dimension", **dict(
        flierprops=flierprops,
        showfliers=True,
        whiskerprops=whiskerprops,
        saturation=1.0,
        showcaps=True,
    ))
    ax.get_legend().set(loc=1)
    ax.set(xlim=(-50, 50))
    ax
    return


@app.cell
def _():
    mo.md(r"""# uMLIP single model error""")
    return


@app.cell
def _():
    dropdow_yval = mo.ui.dropdown(
       options=["energy [eV/atom]", "delta_EnMLIP_minus_EnPBE_eVatoms"], value="energy [eV/atom]", label="y axis",
    )
    return (dropdow_yval,)


@app.cell
def _(df_pred):
    dropdow_models = mo.ui.dropdown(
       options=df_pred.uMLIPs.unique(), value='esen-oam', label="uMLIPs",
    )
    return (dropdow_models,)


@app.cell
def _():
    dropdow_elements = mo.ui.text(
       value="All", label="element",
    )
    return (dropdow_elements,)


@app.cell
def _(df_deltas, dropdow_elements, dropdow_models, dropdow_yval):
    if dropdow_elements.value not in atomic_numbers.keys():
        data = df_deltas.query(f"uMLIPs=='{dropdow_models.value}'") 
        print(dropdow_elements.value + " not in elements allowed")
    else:
        data = df_deltas.query(f"uMLIPs=='{dropdow_models.value}'") .query(f"formula.str.match(r'(?![a-z]){dropdow_elements.value}(?![a-z])')")
    scatter = mo.ui.plotly(px.scatter(
        data,
        x="energy [eV/atom]_PBE", 
        y=dropdow_yval.value, 
        color="id_dimension",hover_data=["mat_id", "id_dimension", "formula"],
    ))
    return data, scatter


@app.cell
def _(scatter):
    max_sample = mo.ui.number(
        start=1, value=10, label=f"N Selected {len(scatter.points)}, MaxSamples (to keep low)"
    )
    ncols = mo.ui.number(start=1, value=8, label="NCols")
    ax_size = mo.ui.number(start=1, value=5, label="AxSize")
    control_subplots = mo.md(f"{max_sample}, {ncols}, {ax_size}")
    return ax_size, control_subplots, max_sample, ncols


@app.cell
def _():
    n_rotations = mo.ui.slider(1, 3, value=1, label="Add Rotation:", show_value=True)
    return (n_rotations,)


@app.cell
def _(n_rotations):
    rotation_dict = mo.ui.dictionary(
        {
            "Angles": mo.ui.array(
                [
                    mo.ui.slider(
                        0, 360, value=0, step=30, label="Angle", show_value=True
                    )
                    for _ in range(n_rotations.value)
                ]
            ),
            "Axis": mo.ui.array(
                [
                    mo.ui.dropdown(["x", "y", "z"], value="x", label="Axis")
                    for _ in range(n_rotations.value)
                ]
            ),
        }
    )
    return (rotation_dict,)


@app.cell
def _(n_rotations, rotation_dict):
    control_rotations = mo.md(
        f"""
        Rotate by: {make_rotation_string(rotation_dict)}\n
        {n_rotations}\n\n
        """
        + "\n\n".join(
            # Iterate over the elements and embed them in markdown
            [
                f"{angles} {ax_rot}"
                for angles, ax_rot in zip(
                    rotation_dict["Angles"], rotation_dict["Axis"]
                )
            ]
        )
    )
    return (control_rotations,)


@app.cell
def _(
    ax_size,
    control_rotations,
    control_subplots,
    data,
    dropdow_elements,
    dropdow_models,
    dropdow_yval,
    max_sample,
    ncols,
    rotation_dict,
    scatter,
):
    mae_energy = data.delta_EnMLIP_minus_EnPBE_eVatoms.abs().mean()
    rmse_energy = ((data.delta_EnMLIP_minus_EnPBE_eVatoms.abs()**2).sum() / len(data)) **0.5
    mo.vstack(
        [
            dropdow_models,
            dropdow_yval,
            dropdow_elements,
            mo.md(f"MAE: {mae_energy:0.3f}, RMSE: {rmse_energy:0.3f} eV/atoms"),
            scatter,
            control_subplots,
            plot_agms(
                scatter.points,
                rotation=make_rotation_string(rotation_dict),
                ncols=ncols.value,
                max_sample=max_sample.value,
                ax_size=ax_size.value,
            ),
            control_rotations,
        ]
    )
    return


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()

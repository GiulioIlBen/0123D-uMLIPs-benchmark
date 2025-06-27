# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "ase>=3.25.0",
#     "jinja2>=3.1.6",
#     "marimo>=0.13.8",
#     "plotly>=6.0.1",
#     "seaborn>=0.13.2",
# ]
# ///
import marimo

__generated_with = "0.13.8"
app = marimo.App(width="full", app_title="0123D dataset explore")

with app.setup:
    # Initialization code that runs before all other cells
    import random

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import plotly.express as px
    from ase.db import connect
    from ase.visualize.plot import plot_atoms


@app.function
def plot_agms(data, max_sample=10, ncols=5, rotation="", ax_size=5):
    if len(data) == 0:
        return mo.md("#No Data Selected")
    agms = [x["mat_id"] for x in data]
    if max_sample is not None:
        if len(agms) > max_sample:
            agms = agms[:max_sample]
    with connect(mo.notebook_location() / "public/0123D_pbe.db") as con:
        atoms_list = [
            con.get(f"mat_id={mat_id}").toatoms(
                add_additional_information=True
            )
            for mat_id in agms
        ]

    nplots = len(atoms_list)
    nrows = -(-nplots // ncols)
    print(nrows)
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(ax_size * ncols, ax_size * nrows),
        layout="tight",
    )
    for atoms, ax in zip(atoms_list, axes.ravel()):
        plot_atoms(atoms, ax=ax, rotation=rotation)
        dim = atoms.info["key_value_pairs"]["id_dimension"]
        agm_i = atoms.info["key_value_pairs"]["mat_id"]
        formula = atoms.get_chemical_formula()
        ax.set_title(f"{dim} - {agm_i}\n{formula}")
    if max_sample is not None:
        fig.suptitle(f"View only {max_sample} samples of {len(data)}")
    return fig


@app.cell
def _():
    df_TSNE = pd.read_csv(mo.notebook_location() / "public/df_TSNE.csv")
    return (df_TSNE,)


@app.cell
def _():
    marker_size = mo.ui.slider(
        start=1, stop=10, step=0.5, value=2, label="MarkerSize", show_value=True
    )
    return (marker_size,)


@app.cell
def _(df_TSNE, marker_size):
    scatter = mo.ui.plotly(
        px.scatter(
            df_TSNE,
            x="t-SNE 1",
            y="t-SNE 2",
            hover_data=["mat_id", "id_dimension"],
            color="id_dimension",
        ).update_traces(marker=dict(size=marker_size.value))
    )
    return (scatter,)


@app.function
def make_rotation_string(controls):
    return ",".join([f"{angles.value}{ax_rot.value}"
            for angles, ax_rot in zip(
                controls["Angles"], controls["Axis"])
                ])


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
    marker_size,
    max_sample,
    ncols,
    rotation_dict,
    scatter,
):
    mo.vstack(
        [
            marker_size,
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


@app.function
def get_table_view(extractor, relative_path="public/0123D_pbe.db"):
    data = []
    with connect(mo.notebook_location() / relative_path) as con:
        tot = con.count()
        for row in mo.status.progress_bar(con.select(), total=tot):
            data.append({
                "id_dimension":row.id_dimension,
                "mat_id":row.mat_id,
                **extractor(row),
                })
    return pd.DataFrame(data)


@app.function
def plot_natoms():
    df = get_table_view(lambda row: {"natoms": row.natoms})
    fig = px.histogram(
        df,
        x="natoms",
        # hover_data=["mat_id", "id_dimension"],
        facet_row="id_dimension",
        color="id_dimension",
    )
    return fig


@app.function
def plot_energy():
    df = get_table_view(lambda row: {"energy [eV/atom]": row.uncorrected_energy/ row.natoms})
    fig = px.histogram(
        df,
        x="energy [eV/atom]",
        # hover_data=["mat_id", "id_dimension"],
        facet_row="id_dimension",
        color="id_dimension",
    )
    return fig


@app.cell
def _():
    def get_dim_natms(row):
        atoms = row.toatoms()
        lengts = atoms.cell.lengths()
        if row.id_dimension == "3D":
            ret = row.volume **(1/3) / (row.natoms)
        elif row.id_dimension == "2D":
            lengts_l = lengts.tolist()
            to_pop = np.argmax(lengts)
            lengts_l.pop(to_pop)
            ret = (lengts_l[0] * lengts_l[1])**(1/2) / (row.natoms)
        elif row.id_dimension == "1D":
            ret = np.min(lengts) / (row.natoms)
        elif row.id_dimension == "0D":
            xyz_coords = atoms.get_positions()
            ret =  (np.linalg.norm(xyz_coords - xyz_coords.mean(axis=0), axis=1) ** 2).sum() ** 0.5 / (row.natoms)
        else:
            raise ValueError(f"{row.id_dimension}")    
        return {"[Vol**1/3|Area**1/2|L]_natoms|GyR [\\A/atom]": ret}

    def plot_volume_like_prop():
        df = get_table_view(get_dim_natms)
        fig = px.histogram(
            df,
            x="[Vol**1/3|Area**1/2|L]_natoms|GyR [\\A/atom]",
            # hover_data=["mat_id", "id_dimension"],
            facet_row="id_dimension",
            color="id_dimension",
        )
        return fig
    return (plot_volume_like_prop,)


@app.cell
def _(plot_volume_like_prop):
    description_plot = {
        "natoms":plot_natoms,
        "volume_like_prop":plot_volume_like_prop,
        "energy_natoms":plot_energy,
        None:lambda: "None"
    }
    dropdow_choice = mo.ui.dropdown(options=description_plot.keys(), label="Choose plot option:")
    return description_plot, dropdow_choice


@app.cell
def _(description_plot, dropdow_choice):
    mo.vstack(
        [
            dropdow_choice,
            description_plot[dropdow_choice.value]()
        ]
    )
    return


if __name__ == "__main__":
    app.run()

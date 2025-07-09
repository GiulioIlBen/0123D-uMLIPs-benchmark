import marimo

__generated_with = "0.14.8"
app = marimo.App(width="full")


@app.cell
def _():
    # Initialization code that runs before all other cells
    import random

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import plotly.express as px
    from ase.db import connect
    from ase.visualize.plot import plot_atoms
    from explore_0123D_test_set import get_table_view, make_rotation_string, plot_agms
    return connect, make_rotation_string, mo, np, plot_agms, plot_atoms, plt


@app.cell
def _(mo):
    mo.md(r"""# In search for nanotubes""")
    return


@app.cell
def _(connect, mo, np, plot_atoms, ransac_circle):
    def find_empty_square(positions, cell=None, min_number_atoms=0):
        if cell is None:
            cell = np.array([[1,0,0], [0,1,0]]).T
        pos_xy = np.dot(positions, cell / np.linalg.norm(cell, axis=0))
        c = pos_xy.mean(axis=0)
        pos_centred = pos_xy  - c
        rs = np.arange(0.1, 5, 0.1)
        p = np.prod(np.abs(pos_centred[:, :, np.newaxis]) < rs, axis=1)
        n_in = p.sum(axis=0)
        options = rs[np.where(n_in<=min_number_atoms)]
        if len(options) == 0:
            return 0
        return options.max()

    def test():
        with connect(mo.notebook_location() / "public/0123D_pbe.db") as con_2:
            row = list(con_2.select("id_dimension=1D,natoms=15"))[1]
        best_circle, best_inliers = ransac_circle(row.positions)
        print(best_circle, best_inliers)
        print(len(best_inliers), row.natoms, row.id_dimension)
        plot_atoms(row.toatoms(), rotation="0x")
    
        find_empty_square(row.positions, row.cell[:, [0,1]])

    return (find_empty_square,)


@app.cell
def _(mo):
    filter_square_edge = mo.ui.slider(0, 10, step=0.1, value=4, label="Square edge length / 2", show_value=True)
    return (filter_square_edge,)


@app.cell
def _(filter_square_edge, find_empty_square, mo):
    from ase.data.vdw import vdw_radii

    info_vdw_for_noobs = mo.md(f"TBN! lowest wdw radii: {vdw_radii[vdw_radii == vdw_radii].min()} | vdw of C: {vdw_radii[6]}")

    def filter_row(row, min_number_atoms=0, min_val=filter_square_edge.value):
        val = find_empty_square(row.positions, None, min_number_atoms)
        # row["key_value_pairs"]["l"] = val
        return val >= min_val
    return filter_row, info_vdw_for_noobs


@app.cell
def _(connect, filter_row, find_empty_square, info_vdw_for_noobs, mo, plt):
    with connect(mo.notebook_location() / "public/0123D_pbe.db") as connnn:
        # print(con.get(1))
        # data = [{"mat_id":r.mat_id} for r in con.select(selection="id_dimension=1D, C")]
        data = [{"mat_id":r.mat_id, "l":find_empty_square(r.positions)} for r in connnn.select(selection="id_dimension=1D", filter=filter_row)]

    ax = plt.gca()
    _ = ax.hist(list(map(lambda x: x["l"], data)))
    summary_hist = mo.hstack([info_vdw_for_noobs, ax])
    return data, summary_hist


@app.cell
def _(data, mo):
    max_sample = mo.ui.number(
        start=1, value=10, label=f"N Selected {len(data)}, MaxSamples (to keep low)"
    )
    return (max_sample,)


@app.cell
def _(max_sample, mo):
    ncols = mo.ui.number(start=1, value=8, label="NCols")
    ax_size = mo.ui.number(start=1, value=5, label="AxSize")
    control_subplots = mo.md(f"{max_sample}, {ncols}, {ax_size}")
    return ax_size, control_subplots, ncols


@app.cell
def _(mo):
    n_rotations = mo.ui.slider(1, 3, value=1, label="Add Rotation:", show_value=True)
    return (n_rotations,)


@app.cell
def _(mo, n_rotations):
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
def _(make_rotation_string, mo, n_rotations, rotation_dict):
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
    filter_square_edge,
    make_rotation_string,
    max_sample,
    mo,
    ncols,
    plot_agms,
    rotation_dict,
    summary_hist,
):
    mo.vstack(
        [
            filter_square_edge,
            summary_hist,
            control_subplots,
            plot_agms(
                data,
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
def _(mo):
    mo.md(r"""# In search for Nanorods""")
    return


@app.cell
def _(mo):
    filter_nano_thikness = mo.ui.slider(0.1, 5, 0.1, value=2, show_value=True, label="thinkness_threshold")
    return (filter_nano_thikness,)


@app.cell
def _(filter_nano_thikness):
    def filter_nanorods(row, thick=filter_nano_thikness.value):
        pos = row.positions
        if max(pos[:,0]) - min(pos[:,0]) > thick:
          return False
        if max(pos[:,1]) - min(pos[:,1]) > thick:
          return False
        if max(pos[:,2]) - min(pos[:,2]) > thick:
          return False
        return True
    return (filter_nanorods,)


@app.cell
def _(connect, filter_nanorods, mo):
    with connect(mo.notebook_location() / "public/0123D_pbe.db") as con:
        data_nanorods = [{"mat_id":r.mat_id} for r in con.select(selection="id_dimension=1D", filter=filter_nanorods)]
    return (data_nanorods,)


@app.cell
def _(
    ax_size,
    data_nanorods,
    filter_nano_thikness,
    make_rotation_string,
    max_sample,
    mo,
    ncols,
    plot_agms,
    rotation_dict,
):
    mo.vstack(
        [
            filter_nano_thikness,
            plot_agms(
                data_nanorods,
                rotation=make_rotation_string(rotation_dict),
                ncols=ncols.value,
                max_sample=max_sample.value,
                ax_size=ax_size.value,
            ),
        ]
    )
    return


if __name__ == "__main__":
    app.run()

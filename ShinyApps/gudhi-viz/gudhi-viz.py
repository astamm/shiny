#! /usr/bin/env python3

import numpy as np
import gudhi
import plotly
from plotly.graph_objs import graph_objs as go
import json

def compute_filtration(file, alpha):
    ac = gudhi.AlphaComplex(off_file = file)
    st = ac.create_simplex_tree()
    points = np.array([ac.get_point(i) for i in range(st.num_vertices())])
    triangles = np.array([s[0] for s in st.get_skeleton(2) if len(s[0]) == 3 and s[1] <= alpha])
    return {"points": points, "triangles": triangles}

def compute_figure(data):
    points = data["points"]
    triangles = data["triangles"]
    mesh = go.Mesh3d(
        x = points[:, 0], 
        y = points[:, 1], 
        z = points[:, 2], 
        i = triangles[:, 0], 
        j = triangles[:, 1], 
        k = triangles[:, 2]
    )
    fig = go.Figure(
        data = mesh, 
        layout = go.Layout(
            title = dict(
                text = 'Alpha Complex Representation of the 2-Torus'
            ), 
            scene = dict(
                xaxis = dict(nticks = 4, range = [-1.5, 1.5]), 
                yaxis = dict(nticks = 4, range = [-1.5, 1.5]), 
                zaxis = dict(nticks = 4, range = [-1.5, 1.5])
            )
        )
    )
    fig.data[0].i = triangles[:, 0]
    fig.data[0].j = triangles[:, 1]
    fig.data[0].k = triangles[:, 2]
    fig = json.dumps(fig, cls = plotly.utils.PlotlyJSONEncoder)
    return fig

#! /usr/bin/env python3

import numpy as np
import gudhi as gd
from plotly.graph_objs import graph_objs as go
from plotly.subplots import make_subplots
from plotly.offline import plot

def load_off_file(file):
    points = gd.read_points_from_off_file(off_file = file)
    points = np.array(points)
    points = points / np.max(points)
    return points

def get_diameter_lower_bound(complexes, dimension = 2):
    rips_skeleton = complexes["rips"]["simplex_tree"].get_skeleton(dimension)
    rips_min_filtration = min([s[1] for s in rips_skeleton if len(s[0]) == (dimension + 1)])
    alpha_skeleton = complexes["alpha"]["simplex_tree"].get_skeleton(dimension)
    alpha_min_filtration = min([s[1] for s in alpha_skeleton if len(s[0]) == (dimension + 1)])
    return max(rips_min_filtration, np.sqrt(3 * alpha_min_filtration))

def build_complexes(pts, max_diameter, dimension = 2):
    return dict(
        rips = rips_complex(pts, max_diameter, dimension),
        alpha = alpha_complex(pts)
    )

def rips_complex(pts, max_diameter, dimension = 2):
    sc = gd.RipsComplex(points = pts, max_edge_length = max_diameter)
    st = sc.create_simplex_tree(max_dimension = dimension)
    return dict(
        points = pts, 
        simplex_tree = st
    )

def alpha_complex(pts):
    sc = gd.AlphaComplex(points = pts)
    st = sc.create_simplex_tree()
    points = np.array([sc.get_point(i) for i in range(st.num_vertices())])
    return dict(
        points = points, 
        simplex_tree = st 
    )

def get_triangles(st, alpha):
    triangles = np.array([s[0] for s in st.get_skeleton(2) if len(s[0]) == 3 and s[1] <= alpha])
    return triangles

def compute_figure(points, triangles, complex_type, scene):
    mesh = go.Mesh3d(
        x = points[:, 0], 
        y = points[:, 1], 
        z = points[:, 2], 
        i = triangles[:, 0], 
        j = triangles[:, 1], 
        k = triangles[:, 2], 
        name = complex_type + " Complex", 
        scene = scene
    )
    return mesh

def sync_figures(complexes, alpha, fps = 0, file = "temp-plot.html"):
    js = '''
    <script>
        var gd = document.getElementById('{div_id}');
        var fpsInterval, startTime, now, then, elapsed;
        // initialize the timer variables and start the animation
        function startAnimating(fps) {{
            fpsInterval = 1000 / fps;
            then = Date.now();
            startTime = then;
            animate();
        }}
        function animate() {{
            // request another frame
            requestAnimationFrame(animate);
            // calc elapsed time since last loop
            now = Date.now();
            elapsed = now - then;
            // if enough time has elapsed, draw the next frame
            if (elapsed > fpsInterval) {{
                // Get ready for next frame by setting then=now, but also adjust for your
                // specified fpsInterval not being a multiple of RAF's interval (16.7ms)
                then = now - (elapsed % fpsInterval);
                // Put your drawing code here
                rotate('scene' , Math.PI / 180);
                rotate('scene2', Math.PI / 180);
            }}
        }}
        function rotate(id, angle) {{
            var scene = gd._fullLayout[id]._scene;
            var camera = scene.getCamera();
            var rtz = xyz2rtz(camera.eye);
            rtz.t += angle;
            camera.eye = rtz2xyz(rtz);
            Plotly.relayout(gd, id + '.camera', camera);
        }}
        function xyz2rtz(xyz) {{
            return {{
                r: Math.sqrt(xyz.x * xyz.x + xyz.y * xyz.y),
                t: Math.atan2(xyz.y, xyz.x),
                z: xyz.z
            }};
        }}
        function rtz2xyz(rtz) {{
            return {{
                x: rtz.r * Math.cos(rtz.t),
                y: rtz.r * Math.sin(rtz.t),
                z: rtz.z
            }};
        }}
        var isUnderRelayout = false;
        gd.on('plotly_relayout', () => {{
            console.log('relayout', isUnderRelayout)
            if (!isUnderRelayout) {{
                Plotly.relayout(gd, 'scene2.camera', gd._fullLayout['scene']._scene.getCamera()).then(() => {{ isUnderRelayout = false }}  )
            }}
            isUnderRelayout = true;
        }})
        startAnimating({fps_input});
    </script>
    '''
    fig = make_subplots(
        rows = 1, 
        cols = 2, 
        specs = [[{'is_3d': True}, {'is_3d': True}]], 
        print_grid = False, 
        subplot_titles = ['Rips Complex', 'Alpha Complex'], 
        vertical_spacing = 0.05, 
        horizontal_spacing = 0.05
    )
    max_rips  = np.max(np.abs(complexes["rips"]["points"]))
    max_alpha = np.max(np.abs(complexes["alpha"]["points"]))
    fig_max = max([max_rips, max_alpha])
    fig1 = compute_figure(
        points = complexes["rips"]["points"],
        triangles = get_triangles(
            st = complexes["rips"]["simplex_tree"],
            alpha = alpha
        ),
        complex_type = "Rips", 
        scene = "scene1"
    )
    fig.add_trace(fig1, row = 1, col = 1)
    fig2 = compute_figure(
        points = complexes["alpha"]["points"],
        triangles = get_triangles(
            st = complexes["alpha"]["simplex_tree"],
            alpha = alpha**2 / 3
        ),
        complex_type = "Alpha", 
        scene = "scene2"
    )
    fig.add_trace(fig2, row = 1, col = 2)
    fig = go.Figure(fig)
    axis_layout = dict(
        showbackground = False, 
        showticklabels = False, 
        title_text = "", 
        range = [-fig_max, fig_max]
    )
    scene_layout = dict(
        xaxis = axis_layout,
        yaxis = axis_layout, 
        zaxis = axis_layout
    )
    margin_value = 10
    fig.update_layout(
        scene = scene_layout,
        scene2 = scene_layout,
        margin = dict(
            r = margin_value, 
            l = margin_value, 
            b = margin_value, 
            t = 40
        ), 
        showlegend = True
    )
    
    # get the a div
    div = plot(fig, include_plotlyjs = False, output_type = 'div')
    # retrieve the div id (you probably want to do something smarter here with beautifulsoup)
    div_id = div.split('=')[1].split()[0].replace("'", "").replace('"', '')
    # your custom JS code
    js = js.format(div_id = div_id, fps_input = fps)
    # # merge everything
    div = '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>' + div + js
    with open(file, 'w') as f:
        f.write(div)


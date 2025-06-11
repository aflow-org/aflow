# Imports
from aflow_chull import CHull
import json
import re
import numpy as np
from random import sample
import plotly.graph_objs as go
from plotly.offline import init_notebook_mode, iplot


# Base class for plotting
class Plotter:
    def __init__(self, aflow_executable = 'aflow'):
        # Set-up notebook for offline mode
        init_notebook_mode(connected=True)

        # Load CHull in case a new json must be created
        self.chull = CHull(aflow_executable=aflow_executable)

    def plot(self, alloy):
        # seperate alloy into elements
        ems = re.findall('[A-Z][^A-Z]*', alloy)
        em1, em2, em3 = ems[0], ems[1], ems[2] if len(ems) == 3 else None
        # ----- Load CHull json -----
        hull = self.chull.get_hull(alloy)
        if not em3:
            # 2D Cartesian plot
            _Two_Species_Plotter().plot_CHull(hull, em1, em2)
        else:
            # Ternary plot
            _Three_Species_Plotter().plot_CHull(hull, em1, em2, em3)
        
        
# For plotting 2 species alloys (2D cartesian)
class _Two_Species_Plotter:
    def plot_CHull(self, hull, em1, em2):
        # Set colors, randomly picked from {red, blue, green}
        color_choices = [[1,0,0], [0,1,0], [0,0,1]]
        colors = sample(color_choices, 2)
        left_color, right_color = colors[0], colors[1]

        # ----- Create scatter for edges -----
        # Grab edges
        edges = [dic['vertices_position'] for dic in list(hull['facets_data'].values())[0]]
        # Create scatter for edges
        hull_scatter = self._create_edges_scatter(edges)

        # ----- Create scatter for points -----
        points_scatter = self._create_points_scatter(hull['points_data'], em1, left_color, right_color)

        # ----- Plot -----
        fig = self._create_figure(hull_scatter, points_scatter, em1, em2)
        iplot(fig)

    def _create_edges_scatter(self, edges):
        edge_xs = []
        edge_ys = []
        for edge in edges:
            point1 = edge[0]
            edge_xs.append(1.0 - point1[0])  # Flip x values because they are reversed in the json output for some reason
            edge_ys.append(point1[1] * 1000.0)  # Multiply by 1000 to match with the points

        hull_scatter = go.Scatter(
            x=np.array(edge_xs),
            y=np.array(edge_ys),
            mode='lines',
            name='hull',
            hoverinfo='none'
        )
        return hull_scatter

    def _create_points_scatter(self, points_data, em1, left_color, right_color):
        colors, auids, urls, c_name, xs, ys = [], [], [], [], [], []

        # First grab points from the points data
        for compound in points_data:
            # Add AUID and url
            urls.append(compound['url_entry_page'])
            auids.append(compound['auid'])
            # Add compound name
            c_name.append(compound['compound'])
            # Set y: enthalpy_formation_atom
            y = compound['enthalpy_formation_atom']
            # Set x
            frac = compound['fractional_compound']
            # If 1 element then figure out if first or second element
            if compound['nspecies'] == 1:
                x = 0.0 if frac == em1 else 1.0
            # o.w. mix of the 2 elements
            else:
                # Example search: 'Pd0.4Pt0.6' -> '0.6'
                x = float(re.search(r'(0.\d*)$', frac).group(0))
            # Add color
            new_color = [(x * right_color[i]) + ((1-x) * left_color[i]) for i in range(3)]
            colors.append(new_color)
            # Add point
            xs.append(x)
            ys.append(y)

        # Create scatter
        points_scatter = go.Scatter(
        x=np.array(xs),
        y=np.array(ys),
        mode='markers',
        name='enthalpy',
        hoverinfo='text',
        text=["""<a href="{}">{}</a>""".format(urls[i], c_name[i]) for i in range(len(urls))],
        hoverlabel=dict(bgcolor='black'),
        marker=dict(
            color=['rgb({}, {}, {})'.format(c[0], c[1], c[2]) for c in colors],
            line=dict(
                color='rgb(0, 0, 0)',
                width=1
                )
            ),
        )

        return points_scatter

    def _create_figure(self, hull_scatter, points_scatter, em1, em2):
        data = [hull_scatter, points_scatter]
        layout = go.Layout(showlegend=False,
                           hovermode='closest',
                           xaxis=dict(
                                title='{}-{}'.format(em1, em2),
                                autorange=True,
                                showgrid=True,
                                zeroline=False,
                                showline=True,
                                showticklabels=True,
                                titlefont=dict(
                                    size=18,
                                )
                            ),
                            yaxis=dict(
                                title='eV/Atom',
                                autorange=True,
                                showgrid=True,
                                zeroline=False,
                                showline=True,
                                showticklabels=True,
                                titlefont=dict(
                                    size=18,
                                )
                            ))
        fig = go.Figure(data=data, layout=layout)
        return fig
    
    
# For plotting 3 species alloys (Ternary)
class _Three_Species_Plotter:
    def plot_CHull(self, hull, em1, em2, em3):
        self.em1, self.em2, self.em3 = em1, em2, em3
        # ----- Create scatter for edges -----
        # Grab compounds that make up points for edges
        compounds = [dic['vertices_position'] for dic in hull['facets_data']['3-nary:{}-{}-{}'.format(em1, em2, em3)]]
        # Create scatter for edges
        hull_scatters = self._create_edges_scatter(compounds)

        # ----- Create scatter for points -----
        points_scatter = self._create_points_scatter(hull['points_data'])

        # ----- Plot -----
        fig = self._create_figure(hull_scatters, points_scatter)
        iplot(fig)

    def _create_edges_scatter(self, compounds):
        hull_scatters = []
        points_a, points_b, points_c = [], [], []
        for vertices in compounds:
            point1, point2, point3 = vertices[0], vertices[1], vertices[2]
            point1[2] = 1.0 - (point1[0] + point1[1])
            point2[2] = 1.0 - (point2[0] + point2[1])
            point3[2] = 1.0 - (point3[0] + point3[1])
            # Important note: As and Bs are switched in the json compared to ternary plot points in Plotly
            hull_scatter = {
                'type': 'scatterternary',
                'mode': 'lines',
                'name': 'hull',
                'hoverinfo': 'none',
                'marker': dict(color='rgba(150, 150, 150, 0.4)', size=0.5),
                'a': np.array([point1[1], point2[1], point3[1]]),
                'b': np.array([point1[0], point2[0], point3[0]]),
                'c': np.array([point1[2], point2[2], point3[2]])
            }
            hull_scatters.append(hull_scatter)
        
        return hull_scatters
    
    def _create_points_scatter(self, points_data):
        colors, auids, urls, c_name, As, Bs, Cs = [], [], [], [], [], [], []

        # First grab points from the points data
        for compound in points_data:
            # Add AUID and url
            urls.append(compound['url_entry_page'])
            auids.append(compound['auid'])
            # Add compound name
            c_name.append(compound['compound'])
            comp = compound['fractional_compound']
            if compound['nspecies'] == 1:
                if comp == self.em1: vals = [1, 0, 0]
                elif comp ==  self.em2: vals = [0, 1, 0]
                else: vals = [0, 0, 1]
            else:
                em_se = [re.search('{}(0.\d*)'.format(em), comp) for em in [self.em1, self.em2, self.em3]]
                vals = [float(v.group(1)) if v else 0.0 for v in em_se]
            As.append(vals[1])
            Bs.append(vals[0])
            Cs.append(vals[2])
            # Add color
            r, g, b = [1,0,0], [0,1,0], [0,0,1]
            new_color = [(As[-1] * b[i]) + (Bs[-1] * r[i]) + (Cs[-1] * g[i]) for i in range(3)]
            colors.append(new_color)
        
        points_scatter = {
                'type': 'scatterternary',
                'mode': 'markers',
                'name': 'points',
                'hoverinfo': 'text',
                'a': np.array(As),
                'b': np.array(Bs),
                'c': np.array(Cs),
                'text': ["""<a href="{}">{}</a>""".format(urls[i], c_name[i]) for i in range(len(urls))],
                'hoverlabel': dict(bgcolor='black'),
                'marker': dict(
                size = 7,
                color=['rgb({}, {}, {})'.format(c[0], c[1], c[2]) for c in colors],
                line=dict(
                    color='rgb(0, 0, 0)',
                    width=1
                    )
                )
            }

        return points_scatter

    def _create_figure(self, hull_scatters, points_scatter):
        data = hull_scatters
        data.append(points_scatter)
        
        layout = {
            'ternary':{
                'sum': 1.0,
                'aaxis': dict(title=self.em2, showgrid=False, showticklabels=False, titlefont=dict(color='blue')),
                'baxis': dict(title=self.em1, showgrid=False, showticklabels=False, titlefont=dict(color='red')),
                'caxis': dict(title=self.em3, showgrid=False, showticklabels=False, titlefont=dict(color='green'))
            },
            'showlegend': False,
            'hovermode': 'closest',
        }
        
        fig = {'data': data, 'layout': layout}
        return fig

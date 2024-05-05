import math
import numpy as np


'''
Module providing functions for marker in diagrams.
'''

def circle(subplots, ax, pos_x, pos_y, zorder, \
           segments_n, segments_filled, \
           size, color_name, size_offset_border=1):
    '''
    Create complete circle marker.

    Args:
        segments_n: number of total segments
        segments_filled: number of filled segments
        size_offset_border: width of white border around circles
    '''
    marker_style = 1 # 0: crystallographic symbols, 1: circles

    # white background circle used as white border
    subplots.append(ax.scatter(pos_x, pos_y, \
                s=(math.sqrt(size)+size_offset_border)**2, \
                marker=marker_component(segments_n, segments_filled, \
                                        0, marker_style), \
                linewidths=0, edgecolors='white', facecolors='white', \
                alpha=1, \
                zorder=zorder, label=''))

    # full circle
    subplots.append(ax.scatter(pos_x, pos_y, \
                s=size, \
                marker=marker_component(segments_n, segments_filled, \
                                        0, marker_style), \
                linewidths=0, edgecolors='white', facecolors=color_name, \
                alpha=1, \
                zorder=zorder, label=''))                  

    # empty sector
    subplots.append(ax.scatter(pos_x, pos_y, \
                s=size, \
                marker=marker_component(segments_n, segments_filled, \
                                        1, marker_style), \
                linewidths=0, edgecolors='white', facecolors='white', \
                alpha=0.7, \
                zorder=zorder, label=''))

    return


def marker_component(segments_n, segments_filled, marker_component, \
                     marker_style):
    '''
    Create marker component.

    Args:
        marker style: (0: crystallographic symbols, 1: circles)
        segments_n: number of total segments
        segments_filled: number of filled segments
        marker_component:
            True: return empty sector
            False: return fully filled symbol for use as background
    '''

    # style: crystallographic symbols
    if marker_style == 0:
        if segments_n == 2:
            w = 0.45

            path_sections = [[(-w, 0), (-w, 1), ( w, 1), ( w, 0)], \
                             [( w, 0), ( w,-1), (-w,-1), (-w, 0)]]

        elif segments_n == 3:
            h = math.sqrt(3)
            r = 1/math.sqrt(3)

            path_sections = [[(-0.5,(h/2)-r), ( 0, h-r), ( 0.5,(h/2)-r)], \
                             [( 0.5,(h/2)-r), ( 1,-r),   (   0,-r)     ], \
                             [(  0,-r),       (-1,-r),   (-0.5,(h/2)-r)]]

        elif segments_n == 4:
            path_sections = [[( 0, 1), ( 1, 0)], \
                             [( 1, 0), ( 0,-1)], \
                             [( 0,-1), (-1,0)], \
                             [(-1, 0), ( 0, 1)]]

        elif segments_n == 6:
            h = math.sqrt(3)/2

            path_sections = [[( 0,     h  ), ( 0.5, h), ( 0.75,  h/2)], \
                             [( 0.75,  h/2), (   1, 0), ( 0.75, -h/2)], \
                             [( 0.75, -h/2), ( 0.5,-h), ( 0,    -h)  ], \
                             [( 0,    -h  ), (-0.5,-h), (-0.75, -h/2)], \
                             [(-0.75, -h/2), (-1,   0), (-0.75,  h/2)], \
                             [(-0.75,  h/2), (-0.5, h), ( 0,     h)]]


    # style: circles
    elif marker_style == 1:
        s = round(180/segments_n)
        alpha = np.linspace(2.5*np.pi, 0.5*np.pi, segments_n*s+1)
        x = np.cos(alpha)
        y = np.sin(alpha)
        c = list(zip(x, y))

        path_sections = [ c[i*s:(i+1)*s+1] for i in range(0, segments_n) ]


    path = []

    # fully filled symbol
    if marker_component == 0:
        for i in range(0, segments_n):

            for s in path_sections[i]:
                if s not in path:
                    path.append(s)

    # empty sector of symbol
    elif marker_component == 1:
        path = [(0, 0)]

        for i in range(0, segments_n):
            if i not in range(0, segments_filled):

                for s in path_sections[i]:
                    if s not in path:
                        path.append(s)

    return path

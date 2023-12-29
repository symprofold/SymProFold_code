import matplotlib.colors as mcolors


'''
Module providing functions for handling of colors.
'''

def group_color(colordict, group_name):
    '''
    Assign a unique color to a given "group_name" of a element group in a
    diagram.
    '''
    if group_name in colordict:
        color_name = colordict[group_name]
    else:
        color_names = dict(mcolors.TABLEAU_COLORS, **mcolors.CSS4_COLORS)

        colors_used = []
        for c in colordict:
            colors_used.append(colordict[c])

        for cn in color_names:
            if cn not in colors_used:
                color_name = cn
                colordict[group_name] = cn
                break

    return colordict, color_name

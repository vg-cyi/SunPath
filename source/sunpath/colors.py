import matplotlib.pyplot as plt
import matplotlib.cm as cm  # todo
import matplotlib.colors

colorMain = '#567CB5'
colorGrid = '#B0B0B0FF'
colorAxes = '#303030'

# https://matplotlib.org/stable/tutorials/colors/colormaps.html
# cMap = plt.cm.turbo #todo

# uniform rainbow from ParaView
colorList = [
    '#0561ff', '#056cf7', '#0577f0', '#0582e8', '#058bdf',
    '#0595d5', '#059ecb', '#05a6bf', '#05afb3', '#05b8a8',
    '#05c19a', '#05ca8c', '#05d47f', '#05dc6d', '#06e55c',
    '#04ed4b', '#46f327', '#7ef51c', '#a4f90c', '#c2fb09',
    '#e1fd06', '#ffff03', '#fff414', '#ffe826', '#ffdd37',
    '#ffd137', '#ffc437', '#ffb837', '#ffac37', '#ff9f37',
    '#ff9337', '#ff8537', '#ff7737', '#ff6837', '#fe5536',
    '#fc4230', '#fd2636', '#f21e40', '#e6144a', '#da0a54',
    '#cc0b5b', '#bd0c62', '#ae0d6a'
]


def makeColorMap(colors):
    ans = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    ans.set_under(colors[0])
    ans.set_over(colors[-1])
    return ans


colorMapL = makeColorMap(colorList)
colorMapL_r = makeColorMap(colorList[::-1])  # reversed

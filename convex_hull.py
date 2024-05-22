###############################################################
# Prof. Luiz T. F. Eleno
# luizeleno@usp.br
# 2024.04.03
# code to plot ternary convex hull with tooltips
# CSV file must have the structure found in energy_data.csv
# Code generate_example_data.py generates an example
##############################################################

import numpy as np
import pandas
import scipy.spatial
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
import distance
import gibbs

# element labels
A, B, C = 'La', 'Si' , 'H'

# reading data
data = pandas.read_table('extended_convex_hull-20GPA.dat', usecols=[0, 1, 2, 3, 4], names=['ID', A, B, C, 'E'], delimiter='\\s+', comment="#")
data.dropna(inplace=True)  # drop out invalid rows

# finding reference states from data
idEA = data.loc[(data[B] == 0) & (data[C] == 0), 'E'].idxmin()
idEB = data.loc[(data[A] == 0) & (data[C] == 0), 'E'].idxmin()
idEC = data.loc[(data[A] == 0) & (data[B] == 0), 'E'].idxmin()
print('REFERENCE STATES:')
print(data.iloc[[idEA, idEB, idEC]])

# getting number of atoms of each type
xA, xB, xC = data[A].values, data[B].values, data[C].values
nat = xA + xB + xC

# getting IDs
IDs = data['ID'].values

# calculating formation energies
energy = data['E'].values  # energies already in eV/atom
EA = energy[idEA] #/ nat[idEA]
EB = energy[idEB] #/ nat[idEB]
EC = energy[idEC] #/ nat[idEC]
energy  = (energy - (xA * EA + xB * EB + xC * EC) / nat ) * 1000  # meV

# calculating atomic fractions
xA = xA / nat
xB = xB / nat
xC = xC / nat

# preparing for convex hull calculation
points = np.column_stack((xB, xC, energy))
points = points[energy <= 0]
IDs = IDs[energy <= 0]  # getting point IDs

# convex hull
hull = scipy.spatial.ConvexHull(points)

# findind distances
distances = distance.dist(hull)

# creating array with all data
points = np.column_stack((points, distances, IDs, np.zeros_like(distances)))  # adding column to indicate if point belongs to hull
points[hull.vertices, -1] = 1  # 1 if in hull, 0 otherwise

# saving excel file
df = pandas.DataFrame(points)
df.to_excel(excel_writer = "20GPa.xlsx")

# filtering out larger distances
dmax =  100  # np.inf
points = points[distances <= dmax]

# unpacking filtered data
x, y, energy, distances, IDs, inhull = points.T

# renumbering hull vertices:
vertices = np.where(inhull)
print("POINTS IN HULL:")
print(IDs[vertices])
data_hull = points[vertices]

# transform compositions to Gibbs triangle
x_tri, y_tri = gibbs.transform_coords(x, y)
xtri_hull, ytri_hull = gibbs.transform_coords(data_hull[:, 0], data_hull[:, 1])

# CREATING GRAPHICS

fig = plt.figure(figsize=(8, 7))
ax = fig.add_subplot(aspect='equal')

rc('font',**{'family':'Palatino','serif':['Palatino'], 'size':10})
rc('text', usetex=True)

mpl.rcParams['axes.linewidth'] = 3.5
mpl.rcParams['xtick.major.width'] = 3.5
mpl.rcParams['ytick.major.width'] = 3.5

# DRAWING GIBBS TRIANGLE
gibbs.triangle(ax, A, B, C)

# PLOTTING DATA

allpoints = ax.scatter(x_tri, y_tri, c=distances, marker='s', s=60, cmap='coolwarm', edgecolors='black', alpha=.4,  zorder=5)
hullpoints = ax.plot(xtri_hull, ytri_hull, 'ok', zorder=6, markersize=10)

# TRIANGULATION
for trio in hull.simplices:
    x1, x2 = hull.points[trio][:, :2].T
    x1_tri, x2_tri = gibbs.transform_coords(x1, x2)
    x1_tri = np.append(x1_tri, x1_tri[0])
    x2_tri = np.append(x2_tri, x2_tri[0])
    plt.plot(x1_tri, x2_tri, 'k-')

# creating tooltips
gibbs.annotate(data_hull, hullpoints)
# gibbs.annotate(points, allpoints)

# SHOW RESULT
cbar = plt.colorbar(allpoints)
cbar.set_label(r'$\Delta E_\textrm{hull}$\,(meV/atom)', fontsize=16)
plt.tight_layout()
plt.savefig('convex_hull.png', bbox_inches='tight', dpi=300)
plt.show()

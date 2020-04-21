# packages
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

from simupy.systems.symbolic import DynamicalSystem, dynamicsymbols
from simupy.block_diagram import BlockDiagram
from simupy.array import Array, r_

def shape_figure():
    fig = plt.gcf()
    fig.set_size_inches(10,10, forward=True)
    leg_artist = plt.legend(loc='lower right', bbox_to_anchor=(1.1, 1.0))
    mplpub.horizontal_center(fig)
    mplpub.vertical_aspect(fig,
            mplpub.golden_ratio, overlapping_extra_artists=[leg_artist])

def plot_x(result, label=''):
    plt.plot(result.t, result.y[:,0]*180/np.pi, label=label)
    plt.xlabel('time, s')
    plt.ylabel('position, degrees')


## define systems
x, v, u = dynamicsymbols('x v u')
l, m = sp.symbols('l m')

parameters = {l: 1, m: 1}

inertia = DynamicalSystem(
        state_equation = r_[v, u/(m*l**2)],
        state = r_[x,v],
        input_ = u,
        constants_values = parameters
)

g = sp.symbols('g')
parameters[g] = 9.8

gravity = DynamicalSystem(
        output_equation = -g*m*l*sp.sin(x),
        input_ = x,
        constants_values = parameters
)

## put them together
BD = BlockDiagram(inertia, gravity)
BD.connect(gravity, inertia)
BD.connect(inertia, gravity, outputs=[0])

plt.figure()
plot_x(BD.simulate(8), 'first simulation!')
shape_figure()
plt.savefig("first_sim.pdf")

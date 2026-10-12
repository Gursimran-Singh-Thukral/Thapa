import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import test1 as t

def tris(shape):
    v, f = shape.tessellate(0.05)
    v = np.array([(p.x, p.y, p.z) for p in v])
    return v[np.array(f)]

def draw(ax, name, p):
    plate, fluid = t.create(**p)
    pl = Poly3DCollection(tris(plate), facecolor='#d9dde3', edgecolor='none', alpha=0.12)
    fl = Poly3DCollection(tris(fluid), facecolor='#1f77b4', edgecolor='#0b3d66', linewidths=0.1)
    ax.add_collection3d(pl); ax.add_collection3d(fl)
    ax.set_xlim(0, t.L); ax.set_ylim(-t.W/2, t.W/2); ax.set_zlim(0, t.H*3)
    ax.set_box_aspect((t.L, t.W, t.H*3))
    ax.view_init(elev=38, azim=-65)
    ax.set_title(f"{name}: n={p['n']}, α={p['alpha']}°, d={p['d']}, h={p['h']} mm")
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('z')

fig = plt.figure(figsize=(11, 9))
draw(fig.add_subplot(211, projection='3d'), 'ICP (initial)', t.ICP)
draw(fig.add_subplot(212, projection='3d'), 'OCP (optimized)', t.OCP)
fig.tight_layout(); fig.savefig('out/leafcool_3d.png', dpi=120)

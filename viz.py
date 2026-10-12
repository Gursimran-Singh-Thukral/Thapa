import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path
import test1 as t

def draw(ax, name, p):
    prof = t.port_and_lines(p['n'], p['alpha'], p['d'])
    ax.add_patch(Rectangle((0, -t.W/2), t.L, t.W, fc='#d9dde3', ec='#555'))
    rings = [prof.exterior] + list(prof.interiors)
    verts, codes = [], []
    for r in rings:
        c = list(r.coords)
        verts += c; codes += [Path.MOVETO] + [Path.LINETO]*(len(c)-2) + [Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), fc='#1f77b4', ec='#0b3d66', lw=0.5))
    ax.set_xlim(-5, t.L+5); ax.set_ylim(-t.W/2-5, t.W/2+5); ax.set_aspect('equal')
    ax.set_title(f"{name}: n={p['n']}, α={p['alpha']}°, d={p['d']} mm, h={p['h']} mm")
    ax.set_xlabel('x (mm)  — flow →'); ax.set_ylabel('y (mm)')

fig, axes = plt.subplots(2, 1, figsize=(10, 11))
draw(axes[0], 'ICP (initial)', t.ICP)
draw(axes[1], 'OCP (optimized)', t.OCP)
fig.tight_layout(); fig.savefig('out/leafcool.png', dpi=130)

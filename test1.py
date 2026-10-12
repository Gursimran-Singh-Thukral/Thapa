"""LeafCool parametric CAD models — geometry parameters from Zhan et al. (2024).

Recreates a leaf-vein interpretation from the source's conceptual Fig. 1,
not an exact reconstruction of the authors' unpublished CAD.
Only the sourced values (L,W,H,a,b,w,x,n,alpha,d,h) are used as dimensional inputs.
The branch positions and line-to-collector join geometry are constructed
from the schematic and therefore necessarily approximate.

Requires cadquery + shapely.  Units: mm.
"""
from __future__ import annotations
import math, argparse, json, pathlib
from dataclasses import dataclass
import cadquery as cq
from cadquery import Vector
from shapely.geometry import LineString, box, Polygon
from shapely.ops import unary_union
from shapely.geometry.polygon import orient

# Paper Table 1, page 3: fixed dimensions
L = 200.0; W = 120.0; H = 6.0; a = 2.0
b = 5.0; w = 10.0; x = 20.0
# Initial structure described in paper page 3
ICP = {'n': 10, 'alpha': 110, 'd': 3.0, 'h': 3.0}
# Optimized parameters specified in paper section 3.7 page 9
OCP = {'n': 13, 'alpha': 91, 'd': 3.8, 'h': 4.0}
# Paper Table 4, p7: orthogonal four-level parameter study
LEVELS = {'n':[8,10,12,14], 'alpha':[90,110,130,150], 'd':[2.5,3.0,3.5,4.0], 'h':[2.0,3.0,4.0,5.0]}
# Paper Table 5, p7: indices 1-4 in order (n,alpha,d,h)
ORTHO_LEVEL_INDICES = [
  [1,1,1,1],[1,2,2,2],[1,3,3,3],[1,4,4,4],
  [2,1,2,3],[2,2,1,4],[2,3,4,1],[2,4,3,2],
  [3,1,3,4],[3,2,4,3],[3,3,1,2],[3,4,2,1],
  [4,1,4,2],[4,2,3,1],[4,3,2,4],[4,4,1,3]]

def port_and_lines(n, alpha, d):
    # Midrib, rim veins, and cross-collector connections as interpreted from Fig. 1.
    # W-a-b/2 places the rim centered 2 mm inward from each side wall.
    margin = a+b/2
    rim_y = W/2-margin
    rim_left = margin + b/2 + 2   # branches and rim have space for the inlet
    rim_right = L-rim_left
    pieces = [
        LineString([(0,0),(L,0)]).buffer(b/2, cap_style=2, join_style=2),
        box(0,-w/2,10,w/2), box(L-10,-w/2,L,w/2),
    ]
    for s in (-1,1):
        rim_center = s*rim_y
        pieces += [
            LineString([(rim_left, rim_center),(rim_right,rim_center)]).buffer(b/2,cap_style=2,join_style=2),
            LineString([(rim_left,0),(rim_left,rim_center)]).buffer(b/2,cap_style=2,join_style=2),
            LineString([(rim_right,0),(rim_right,rim_center)]).buffer(b/2,cap_style=2,join_style=2),
        ]
        # Alpha is the included angle relative to the upstream-facing trunk ray.
        # At 110 deg, each vein is at 70 deg to +x (Fig 1 orientation).
        acute = 180.0-alpha
        offset = (rim_y-b/2)*math.tan(math.radians(90-acute)) if acute != 90 else 0.
        # ie delta_x=delta_y/tan(acute), with vertical for alpha=90.
        # Near-outlet densification within paper's x=20 mm zone.
        dense_count = min(3, max(n-2,0))
        regular_count = n-dense_count
        xmin = max(16.0 + offset, 25.0)
        xmax = L-x-5
        if xmin > xmax: xmin = xmax
        regular = [(xmin + (xmax-xmin)*i/max(1,regular_count-1)) for i in range(regular_count)]
        dense = [L-x + 2.0+(x-12)*i/max(1,dense_count-1) for i in range(dense_count)]
        for ex in regular+dense:
            sx = ex-offset
            p0 = (sx, s*b/2)
            p1 = (ex, rim_center)
            pieces.append(LineString([p0,p1]).buffer(d/2,cap_style=2,join_style=2))
    profile = unary_union(pieces).buffer(0)
    # Channel polygon deliberately extends to plate end faces (inlet/outlet openings).
    profile = profile.intersection(box(0,-W/2+a,L,W/2-a)).buffer(0)
    if profile.geom_type != 'Polygon':
        raise RuntimeError(f'Fluid network not connected: {profile.geom_type}')
    return profile

def ring_wire(coords, z):
    pts = [Vector(float(px),float(py),float(z)) for px,py in list(coords)[:-1]]
    return cq.Wire.assembleEdges([cq.Edge.makeLine(pts[i],pts[(i+1)%len(pts)]) for i in range(len(pts))])

def create(n,alpha,d,h):
    if not (8<=n<=14 and 90<=alpha<=150 and 2.5<=d<=4 and 2<=h<=5):
        raise ValueError('Outside the paper Table 4 parameter bounds')
    if h >= H: raise ValueError('Channel must be shallower than plate thickness')
    profile=orient(port_and_lines(n,alpha,d),sign=1.0)
    # Remove unattached sub-mm slivers trapped at diagonal-vein crossings.
    # This is a CAD repair only (5 mm^2 solid-island cutoff), not a sourced design variable.
    profile=orient(Polygon(profile.exterior.coords,[r.coords for r in profile.interiors
                    if Polygon(r).area >= 5.0]),sign=1.0)
    z=(H-h)/2
    outer = ring_wire(profile.exterior.coords,z)
    holes = [ring_wire(r.coords,z) for r in profile.interiors]
    fluid = cq.Solid.extrudeLinear(outer,holes,Vector(0,0,h))
    stock = cq.Workplane('XY').box(L,W,H,centered=(False,True,False)).val()
    plate = stock.cut(fluid)
    return plate, fluid

def ascii_map(profile, cols=100, rows=30):
    """Top-down terminal view: '#' = fluid channel, '.' = solid plate."""
    from shapely.geometry import Point
    lines = []
    for r in range(rows):
        y = W/2 - (r+0.5)*W/rows
        lines.append(''.join('#' if profile.contains(Point((c+0.5)*L/cols, y)) else '.'
                             for c in range(cols)))
    return '\n'.join(lines)

def export(name, params, outdir=None, show_map=True, quiet=False):
    plate, fluid = create(**params)
    if outdir:
        outdir = pathlib.Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
        cq.exporters.export(cq.Workplane('XY').add(plate), str(outdir/f'{name}_plate.step'))
        cq.exporters.export(cq.Workplane('XY').add(fluid), str(outdir/f'{name}_fluid.step'))
    info = {'name': name, **params, 'plate_volume_mm3': plate.Volume(),
            'fluid_volume_mm3': fluid.Volume(), 'valid': bool(plate.isValid() and fluid.isValid())}
    if quiet: return info
    bar = '=' * 100
    print(f"{bar}\n {name}: n={params['n']}  alpha={params['alpha']} deg  d={params['d']} mm  h={params['h']} mm")
    print(f"   plate volume : {info['plate_volume_mm3']:>12,.1f} mm^3")
    print(f"   fluid volume : {info['fluid_volume_mm3']:>12,.1f} mm^3  ({100*info['fluid_volume_mm3']/(L*W*H):.1f}% of stock)")
    print(f"   valid solids : {info['valid']}")
    if outdir:
        print(f"   files        : {outdir/(name+'_plate.step')}, {outdir/(name+'_fluid.step')}")
    if show_map:
        print(f"\n Top view (x: 0..{L:.0f} mm, flow left -> right; y: +{W/2:.0f}..-{W/2:.0f} mm)")
        print(ascii_map(port_and_lines(params['n'], params['alpha'], params['d'])))
    return info

CASES = [
    ('Case 1 - ICP (initial)',   ICP),
    ('Case 2 - OCP (optimized)', OCP),
    ('Case 3 - custom',          {'n': 12, 'alpha': 130, 'd': 3.5, 'h': 5.0}),
]

def make_figure(cases, title=None, ncols=1):
    """Top-view image of each case's channel network."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import PathPatch, Rectangle
    from matplotlib.path import Path
    nrows = -(-len(cases)//ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.5*ncols, 4.0*nrows), squeeze=False)
    for ax in axes.flat[len(cases):]: ax.axis('off')
    for ax, (name, p) in zip(axes.flat, cases):
        prof = port_and_lines(p['n'], p['alpha'], p['d'])
        ax.add_patch(Rectangle((0, -W/2), L, W, fc='#d9dde3', ec='#555'))
        verts, codes = [], []
        for r in [prof.exterior] + list(prof.interiors):
            c = list(r.coords)
            verts += c; codes += [Path.MOVETO] + [Path.LINETO]*(len(c)-2) + [Path.CLOSEPOLY]
        ax.add_patch(PathPatch(Path(verts, codes), fc='#1f77b4', ec='#0b3d66', lw=0.5))
        ax.set_xlim(-5, L+5); ax.set_ylim(-W/2-5, W/2+5); ax.set_aspect('equal')
        ax.set_title(f"{name}\nn={p['n']}, alpha={p['alpha']}, d={p['d']}, h={p['h']}", fontsize=9)
        ax.set_xlabel('x (mm) - flow ->', fontsize=8); ax.set_ylabel('y (mm)', fontsize=8)
    if title: fig.suptitle(title, fontweight='bold')
    fig.tight_layout()
    return fig

def sweep_cases(base, param):
    """One-at-a-time study: vary `param` over its Table 4 levels, others held at `base`."""
    return [(f'{param} = {v}', {**base, param: v}) for v in LEVELS[param]]

def result_table(infos):
    lines = [f"{'case':<28}{'n':>4}{'alpha':>7}{'d':>6}{'h':>6}{'fluid mm^3':>14}{'fluid %':>9}{'valid':>7}"]
    for r in infos:
        lines.append(f"{r['name']:<28}{r['n']:>4}{r['alpha']:>7}{r['d']:>6}{r['h']:>6}"
                     f"{r['fluid_volume_mm3']:>14,.0f}{100*r['fluid_volume_mm3']/(L*W*H):>9.1f}{str(r['valid']):>7}")
    return '\n'.join(lines)

def write_csv(path, rows):
    import csv
    for r in rows: r['fluid_pct_of_stock'] = round(100*r['fluid_volume_mm3']/(L*W*H), 2)
    cols = list(rows[0])
    with open(path, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=cols); wr.writeheader()
        wr.writerows({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()} for r in rows)

def main():
    ap = argparse.ArgumentParser(description="LeafCool parametric CAD models")
    ap.add_argument('--set', choices=['sweep','cases','icp','ocp','ortho'], default='sweep',
                    help='sweep: change one parameter at a time around --base (default)')
    ap.add_argument('--base', choices=['icp','ocp'], default='ocp', help='baseline for --set sweep')
    ap.add_argument('--out', default=None, help='optional: folder to also save STEP files')
    ap.add_argument('--results', default='results', help='folder for the image (PNG) and table (CSV) files')
    ap.add_argument('--no-files', action='store_true', help='terminal output only: write no PNG/CSV files')
    ap.add_argument('--no-images', action='store_true', help='skip the images')
    ap.add_argument('--show', action='store_true', help='also open the images in windows')
    args = ap.parse_args()
    save = not args.no_files
    resdir = pathlib.Path(args.results)
    if save: resdir.mkdir(parents=True, exist_ok=True)
    if save and not args.show:
        import matplotlib; matplotlib.use('Agg')
    all_infos = []
    figs = []   # (name, figure)

    if args.set == 'sweep':
        base = ICP if args.base == 'icp' else OCP
        print(f"One-at-a-time study around {args.base.upper()} baseline: {base}")
        for param in ('n', 'alpha', 'd', 'h'):
            cases = sweep_cases(base, param)
            print('\n' + '='*86 + f"\n SWEEP OF {param}  (others fixed: " +
                  ', '.join(f'{k}={v}' for k, v in base.items() if k != param) + ')\n' + '='*86)
            infos = [export(f'{param} = {p[param]}', p, args.out, quiet=True) for _, p in cases]
            print(result_table(infos))
            all_infos += [{'sweep': param, **r} for r in infos]
            if save: write_csv(resdir/f'table_sweep_{param}.csv', infos)
            if not args.no_images:
                figs.append((f'sweep_{param}', make_figure(cases, f'Sweep of {param} (others fixed)', ncols=2)))
    elif args.set in ('cases', 'icp', 'ocp'):
        cases = {'cases': CASES, 'icp': CASES[:1], 'ocp': CASES[1:2]}[args.set]
        infos = [export(name, p, args.out) for name, p in cases]
        if len(infos) > 1:
            print('\n' + '='*86); print(result_table(infos))
        if not args.no_images:
            figs.append(('cases', make_figure(cases)))
    else:
        keys = ['n','alpha','d','h']
        infos = []
        for i, idx in enumerate(ORTHO_LEVEL_INDICES, 1):
            p = {k: LEVELS[k][j-1] for k, j in zip(keys, idx)}
            infos.append(export(f'ortho_{i:02d}', p, args.out, quiet=True))
        print('\n' + '='*86); print(result_table(infos))

    if save:
        if all_infos: write_csv(resdir/'table_all_sweeps.csv', all_infos)
        elif infos: write_csv(resdir/f'table_{args.set}.csv', infos)
        for name, fig in figs:
            fig.savefig(resdir/f'{name}.png', dpi=110)
        written = sorted(p.name for p in resdir.iterdir())
        print(f"\nSaved {len(written)} files to {resdir.resolve()}:\n  " + '\n  '.join(written))
    if figs and args.show:
        import matplotlib.pyplot as plt
        plt.show()

if __name__ == '__main__':
    main()

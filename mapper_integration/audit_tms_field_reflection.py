#!/usr/bin/env python3
from pathlib import Path
import csv
import math
import numpy as np

MAPPER = Path('Mapper_global_DIAGNOSTIC.txt')
RAW_INNER = Path('Field_maps/Plate1_15.fld')
RAW_OUTER = Path('Field_maps/Plate2_15.fld')
Z_LOCAL_MM = -4000.0
Y_GLOBAL_OFFSET_MM = -1352.23
Z_GLOBAL_OFFSET_MM = 15202.5
OUTDIR = Path('TMS_REFLECTION_AUDIT')


def numeric_rows(path):
    with path.open('r', errors='ignore') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            try:
                vals = [float(v) for v in s.replace(',', ' ').split()]
            except ValueError:
                continue
            yield vals


def read_raw_fld(path):
    rows = []
    for vals in numeric_rows(path):
        if len(vals) >= 6:
            rows.append(tuple(vals[:6]))
    if not rows:
        raise RuntimeError(f'No field rows found in {path}')
    return np.asarray(rows, dtype=float)


def raw_stats(name, arr):
    out = {'name': name, 'n': len(arr)}
    for comp, col in (('Bx',3),('By',4),('Bz',5)):
        med = float(np.median(arr[:,col]))
        out[f'{comp}_median'] = med
        out[f'{comp}_sign'] = '+' if med > 0 else '-' if med < 0 else '0'
    out['x_min'] = float(np.min(arr[:,0])); out['x_max'] = float(np.max(arr[:,0]))
    out['y_min'] = float(np.min(arr[:,1])); out['y_max'] = float(np.max(arr[:,1]))
    return out


def read_mapper_slice(path, z_local):
    rows = numeric_rows(path)
    meta = next(rows)
    if len(meta) < 6:
        raise RuntimeError('Bad mapper header')
    dx,dy,dz = meta[3:6]
    z_global = z_local + Z_GLOBAL_OFFSET_MM
    tol = max(abs(dz)*1e-6, 1e-6)
    data = {}
    for vals in rows:
        if len(vals) < 6:
            continue
        xg,yg,zg,bx,by,bz = vals[:6]
        if abs(zg-z_global) > tol:
            continue
        xl = xg
        yl = yg - Y_GLOBAL_OFFSET_MM
        data[(round(xl,6), round(yl,6))] = np.array([bx,by,bz], dtype=float)
    if not data:
        raise RuntimeError('No mapper rows found at requested z slice')
    return data,(dx,dy,dz),z_global


def expected(b):
    bx,by,bz = b
    return {
        'x': np.array([ bx,-by,-bz]),
        'y': np.array([-bx, by,-bz]),
        'xy': np.array([-bx,-by, bz]),
    }


def audit_pairs(data):
    rows = []
    for (x,y),bq in data.items():
        if x <= 0 or y <= 0:
            continue
        keys = {
            'x': (round(-x,6), round(y,6)),
            'y': (round(x,6), round(-y,6)),
            'xy': (round(-x,6), round(-y,6)),
        }
        if not all(k in data for k in keys.values()):
            continue
        exp = expected(bq)
        for mode in ('x','y','xy'):
            act = data[keys[mode]]
            res = act - exp[mode]
            rows.append({
                'x_mm':x,'y_mm':y,'reflection':mode,
                'Bx_src':bq[0],'By_src':bq[1],'Bz_src':bq[2],
                'Bx_actual':act[0],'By_actual':act[1],'Bz_actual':act[2],
                'Bx_expected':exp[mode][0],'By_expected':exp[mode][1],'Bz_expected':exp[mode][2],
                'dBx':res[0],'dBy':res[1],'dBz':res[2],
            })
    return rows


def summarize(rows, mode):
    sub = [r for r in rows if r['reflection'] == mode]
    if not sub:
        return None
    arr = np.asarray([[r['dBx'],r['dBy'],r['dBz']] for r in sub], dtype=float)
    return {
        'n':len(sub),
        'median_abs':np.median(np.abs(arr),axis=0),
        'rms':np.sqrt(np.mean(arr**2,axis=0)),
        'max_abs':np.max(np.abs(arr),axis=0),
    }


def main():
    for p in (MAPPER,RAW_INNER,RAW_OUTER):
        if not p.exists():
            raise SystemExit(f'Missing {p}\nRun from /Users/h988h767/CODES/NEW_MAPPER_WITH__NEW_GEOMETRY')

    OUTDIR.mkdir(exist_ok=True)
    inner = raw_stats('Plate1_15 (inner)', read_raw_fld(RAW_INNER))
    outer = raw_stats('Plate2_15 (outer)', read_raw_fld(RAW_OUTER))
    mapper,spacing,zglobal = read_mapper_slice(MAPPER,Z_LOCAL_MM)
    rows = audit_pairs(mapper)

    if rows:
        csv_path = OUTDIR/'symmetric_pair_residuals.csv'
        with csv_path.open('w',newline='') as f:
            w = csv.DictWriter(f,fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    lines = []
    lines.append('TMS MAGNETIC-FIELD REFLECTION AUDIT\n')
    lines.append('='*72+'\n\n')
    lines.append('A) ORIGINAL SUPPLIED +X,+Y QUARTER FIELD\n')
    lines.append('-'*72+'\n')
    for s in (inner,outer):
        lines.append(f"{s['name']}\n")
        lines.append(f"  N raw points : {s['n']}\n")
        lines.append(f"  x range [m]  : {s['x_min']:.4f} to {s['x_max']:.4f}\n")
        lines.append(f"  y range [m]  : {s['y_min']:.4f} to {s['y_max']:.4f}\n")
        lines.append(f"  median Bx    : {s['Bx_median']:+.6f} T\n")
        lines.append(f"  median By    : {s['By_median']:+.6f} T\n")
        lines.append(f"  median Bz    : {s['Bz_median']:+.6f} T\n\n")

    lines.append('Raw-quarter dominant By signs:\n')
    lines.append(f"  INNER Plate1: {inner['By_sign']}Y\n")
    lines.append(f"  OUTER Plate2: {outer['By_sign']}Y\n\n")

    lines.append('B) FULL MAPPER PURE AXIAL-REFLECTION TEST\n')
    lines.append('-'*72+'\n')
    lines.append(f'local z = {Z_LOCAL_MM:.1f} mm\n')
    lines.append(f'global z = {zglobal:.1f} mm\n')
    lines.append(f'grid spacing = ({spacing[0]:g}, {spacing[1]:g}, {spacing[2]:g}) mm\n')
    lines.append(f'stored XY points = {len(mapper)}\n\n')
    lines.append('Expected pure axial-vector mirror rules:\n')
    lines.append('  x  : (Bx,By,Bz) -> ( Bx,-By,-Bz)\n')
    lines.append('  y  : (Bx,By,Bz) -> (-Bx, By,-Bz)\n')
    lines.append('  xy : (Bx,By,Bz) -> (-Bx,-By, Bz)\n\n')

    for mode in ('x','y','xy'):
        s = summarize(rows,mode)
        if s is None:
            lines.append(f'{mode}: no symmetric pairs found\n')
            continue
        med=s['median_abs']; rms=s['rms']; mx=s['max_abs']
        lines.append(f"{mode} reflection, N={s['n']}\n")
        lines.append(f"  median abs residual [T]: Bx={med[0]:.3e}, By={med[1]:.3e}, Bz={med[2]:.3e}\n")
        lines.append(f"  RMS residual [T]       : Bx={rms[0]:.3e}, By={rms[1]:.3e}, Bz={rms[2]:.3e}\n")
        lines.append(f"  max abs residual [T]   : Bx={mx[0]:.3e}, By={mx[1]:.3e}, Bz={mx[2]:.3e}\n\n")

    lines.append('C) INTERPRETATION\n')
    lines.append('-'*72+'\n')
    lines.append('Residuals near zero mean the generated mapper follows the pure axial-vector reflection rule internally.\n')
    lines.append('That does not prove that pure geometric reflection is the final physical relation between the two powered TMS halves.\n')
    lines.append('The final physical sign convention still needs confirmation from the intended coil/current polarity and agreed field baseline.\n')

    report = OUTDIR/'reflection_audit.txt'
    report.write_text(''.join(lines))
    print(report.read_text())
    print('Saved:',report)
    if rows:
        print('Saved:',OUTDIR/'symmetric_pair_residuals.csv')


if __name__ == '__main__':
    main()

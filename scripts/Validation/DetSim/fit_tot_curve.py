"""Fit of the fast-shaper pulse order and the bench threshold to the A5202 ToT-versus-charge curve.

Source: CAEN UM7945 (A5202/DT5202 user manual) rev. 4, sec. 9.3.2: ToT = 31.514 ln(Qin/pC) + 45.127 (0.5 ns channels), TD threshold
290, Qin 0.5-13.5 pC. Model: unipolar pulse (t/tp)^n exp(n(1-t/tp)), tp = 15 ns (Citiroc-1A fast shaper), first lobe only.
Result (2026-10-12): n = 2.1 (n = 2: rms 2.8 ns; n = 1: 11.8 ns), threshold 0.48 pC = 1.7 pe at 1.7e6 gain (0.272 pC/pe);
used for Sim.FrontEnd.FastShaperOrder / DiscriminatorThreshold in config/TMS_Readout_Default_Config.toml.
Run: PYTHONDONTWRITEBYTECODE=1 python3 fit_tot_curve.py   (needs numpy)
"""
import numpy as np
tp = 15.0
Q = np.exp(np.linspace(np.log(0.5), np.log(13.5), 40))
tot_data = 0.5 * (31.514 * np.log(Q) + 45.127)          # ns
t = np.arange(0, 600, 0.05)
def tot_model(n, qth):                                   # qth = charge (pC) whose pulse peak just reaches threshold
    x = t / tp; shape = x**n * np.exp(n * (1 - x))
    out = []
    for q in Q:
        above = (q / qth) * shape >= 1.0
        out.append(0.05 * above.sum() if above.any() else 0.0)   # first pulse only (unipolar, single lobe)
    return np.array(out)
best = None
for n in np.arange(0.5, 4.01, 0.1):
    for qth in np.exp(np.linspace(np.log(0.02), np.log(8.0), 160)):
        r = tot_model(n, qth) - tot_data
        c = np.sqrt(np.mean(r**2))
        if best is None or c < best[0]: best = (c, n, qth)
c, n, qth = best
print("best n=%.2f  qth=%.3f pC  rms residual=%.2f ns" % (n, qth, c))
for nn in (1, 2):
    b = min(((np.sqrt(np.mean((tot_model(nn, q) - tot_data)**2)), q) for q in np.exp(np.linspace(np.log(0.02), np.log(8.0), 300))))
    print("n=%d: best qth=%.3f pC rms=%.2f ns" % (nn, b[1], b[0]))
m = tot_model(n, qth)
for i in range(0, 40, 6): print("Q=%.2f pC  data %.1f ns  model %.1f ns" % (Q[i], tot_data[i], m[i]))
print("pC per PE (gain 1.7e6): %.3f ; qth in PE: %.2f" % (1.7e6 * 1.602e-7, qth / (1.7e6 * 1.602e-7)))

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-field O(4) bounce (SFV/dSB) with robust figure exports (v4h).

Fixes the action-density plot:
  - Evaluates on a geometric r-grid strictly > 0
  - Uses symlog y-scale (no artificial clipping)
  - Overlays KE and PE contributions for sanity
Saves:
  data/background_profile.csv
  figs/bounce_fields.png
  figs/bounce_action_density.png
  figs/bounce_potential_slice.png
  figs/bounce_residuals.png
  figs/bounce_mesh.png
"""

import os, sys, time, argparse, traceback
import numpy as np

# Non-interactive backend
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.integrate import solve_bvp, trapezoid

# ----------------- Utility -----------------
def ensure_dirs():
    for d in ("figs", "data"):
        os.makedirs(d, exist_ok=True)
        print(f"[init] ensured ./{d} exists")

def confirm(path):
    ok = os.path.exists(path) and os.path.getsize(path) > 0
    print(f"[write] {'OK ' if ok else 'ERR'} -> {path}")
    return ok

# ----------------- Potential -----------------
def potential_rescaled(Phi_t, phi_t, p):
    V_phi   = p['lam_phi']/4.0*(Phi_t**2 - p['v_phi_t']**2)**2 + p['bias_t']*Phi_t**2
    V_brane = p['lam_t']/4.0*(phi_t**2 - 1.0)**2 - p['mu2_t']*phi_t**2
    V_int   = p['g_portal_t']*Phi_t**2*phi_t**2
    return V_phi + V_brane + V_int

def derivatives_rescaled(Phi_t, phi_t, p):
    Phi_t = np.clip(Phi_t, -1e6, 1e6)
    phi_t = np.clip(phi_t, -1e6, 1e6)
    dVdPhi_t = p['lam_phi']*Phi_t*(Phi_t**2 - p['v_phi_t']**2) + 2*p['bias_t']*Phi_t + 2*p['g_portal_t']*Phi_t*phi_t**2
    dVdphi_t = p['lam_t'] *phi_t*(phi_t**2 - 1)             - 2*p['mu2_t']*phi_t + 2*p['g_portal_t']*Phi_t**2*phi_t
    return dVdPhi_t, dVdphi_t

# ----------------- ODEs, Jacobian, BCs -----------------
def bounce_odes_rescaled(r, y, p):
    Phi_t, Phi_p, phi_t, phi_p = y
    r_safe = np.maximum(r, 1e-9)
    dVdPhi_t, dVdphi_t = derivatives_rescaled(Phi_t, phi_t, p)
    dPhi_p_dr = dVdPhi_t - (3/r_safe)*Phi_p
    dphi_p_dr = dVdphi_t - (3/r_safe)*phi_p
    return np.vstack([Phi_p, dPhi_p_dr, phi_p, dphi_p_dr])

def bounce_odes_jac(r, y, p):
    Phi_t, _, phi_t, _ = y
    r_safe = np.maximum(r, 1e-9)
    d2VdPhi2    = p['lam_phi']*(3*Phi_t**2 - p['v_phi_t']**2) + 2*p['bias_t'] + 2*p['g_portal_t']*phi_t**2
    d2Vdphi2    = p['lam_t']  *(3*phi_t**2 - 1)              - 2*p['mu2_t'] + 2*p['g_portal_t']*Phi_t**2
    d2VdPhidphi = 4*p['g_portal_t']*Phi_t*phi_t
    J = np.zeros((4,4,r.shape[0]))
    J[0,1,:] = 1
    J[1,0,:] = d2VdPhi2
    J[1,1,:] = -3/r_safe
    J[1,2,:] = d2VdPhidphi
    J[2,3,:] = 1
    J[3,0,:] = d2VdPhidphi
    J[3,2,:] = d2Vdphi2
    J[3,3,:] = -3/r_safe
    return J

def boundary_conditions_rescaled(ya, yb, p):
    # Negative-phi false vacuum at infinity
    return np.array([
        ya[0]-p['Phi_TV_t'],
        ya[2]-p['phi_TV_t'],
        yb[0]-p['Phi_FV_t'],
        yb[2]+p['phi_FV_t']
    ])

# ----------------- Helpers -----------------
def initial_guess_tanh(mesh, p, R0, w):
    phi_FV_target = -p['phi_FV_t']
    Phi_guess = 0.5*(p['Phi_TV_t']-p['Phi_FV_t'])*(1-np.tanh((mesh-R0)/w)) + p['Phi_FV_t']
    phi_guess = 0.5*(p['phi_TV_t']-phi_FV_target)*(1-np.tanh((mesh-R0)/w)) + phi_FV_target
    dPhi_guess = np.gradient(Phi_guess, mesh)
    dphi_guess = np.gradient(phi_guess, mesh)
    dPhi_guess[0]=dphi_guess[0]=0.0
    dPhi_guess[-1]=dphi_guess[-1]=0.0
    Phi_guess[-1], phi_guess[-1] = p['Phi_FV_t'], phi_FV_target
    return np.vstack([Phi_guess, dPhi_guess, phi_guess, dphi_guess])

def wall_radius_from_solution(p, sol):
    r, y = sol.x, sol.y
    Phi, Phi_p, phi, phi_p = y
    V_false = potential_rescaled(p['Phi_FV_t'], p['phi_FV_t'], p)
    KE = 0.5*(Phi_p**2 + phi_p**2)
    PE = potential_rescaled(Phi, phi, p) - V_false
    dens4 = 2.0*np.pi**2 * r**3 * (KE + PE)
    return float(r[np.argmax(dens4)])

def calculate_action(p, sol):
    r, y = sol.x, sol.y
    Phi, Phi_p, phi, phi_p = y
    V_false = potential_rescaled(p['Phi_FV_t'], p['phi_FV_t'], p)
    KE = 0.5*(Phi_p**2 + phi_p**2)
    PE = potential_rescaled(Phi, phi, p) - V_false
    return 2*np.pi**2*trapezoid(r**3*(KE+PE), r)

def residual_profile(p, sol, r_eval=None):
    if r_eval is None:
        r_eval = sol.x
    Y = sol.sol(r_eval)
    dY = np.gradient(Y, r_eval, axis=1)
    f  = bounce_odes_rescaled(r_eval, Y, p)
    res = np.linalg.norm(dY - f, axis=0)
    return r_eval, res

# ----------------- Figures -----------------
def save_background_csv(sol, path="data/background_profile.csv"):
    r = sol.x
    Phi, _, phi, _ = sol.y
    import csv
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["r","Phi","phi"])
        for i in range(len(r)):
            w.writerow([f"{r[i]:.10g}", f"{Phi[i]:.10g}", f"{phi[i]:.10g}"])
    confirm(path)

def fig_fields(sol, p, out="figs/bounce_fields.png"):
    r = sol.x
    Phi, _, phi, _ = sol.y
    plt.figure(figsize=(7.0,4.6))
    plt.plot(r, Phi, label="Phi")
    plt.plot(r, phi, label="phi")
    plt.axhline(-p['phi_FV_t'], color='0.6', ls='--', lw=1, label='target -phi_FV')
    plt.xlabel(r"$\tilde r$")
    plt.ylabel("fields")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    confirm(out)

def fig_action_density(sol, p, out="figs/bounce_action_density.png"):
    """
    Corrected action-density plot:
      I(r) = 2π^2 r^3 [ KE(r) + PE(r) ], with KE=0.5(Phi'^2+phi'^2), PE=V-V_false.
      Evaluate on geom-spaced r>0; use symlog to show small negative regions cleanly.
    """
    # Geometric r-grid strictly > 0 to avoid log artifacts
    rmin = max(sol.x[1], 1e-8)
    rmax = sol.x[-1]
    r_eval = np.geomspace(rmin, rmax, 4000)

    Y = sol.sol(r_eval)
    Phi, phi = Y[0], Y[2]
    # smooth derivatives on dense grid
    dPhi = np.gradient(Phi, r_eval)
    dphi = np.gradient(phi, r_eval)

    V_false = potential_rescaled(p['Phi_FV_t'], p['phi_FV_t'], p)
    KE = 0.5*(dPhi**2 + dphi**2)
    PE = potential_rescaled(Phi, phi, p) - V_false

    I_tot = 2*np.pi**2 * (r_eval**3) * (KE + PE)
    I_ke  = 2*np.pi**2 * (r_eval**3) * KE
    I_pe  = 2*np.pi**2 * (r_eval**3) * PE

    # Print a quick consistency check
    S_num = 2*np.pi**2*trapezoid((r_eval**3)*(KE+PE), r_eval)
    print(f"[check] S_E from dense grid ≈ {S_num:.6g}")

    plt.figure(figsize=(7.2,4.8))
    plt.plot(r_eval, I_tot, label=r"$I_{\rm total}(r)$", lw=2)
    plt.plot(r_eval, I_ke,  label=r"$I_{\rm KE}(r)$", alpha=0.7)
    plt.plot(r_eval, I_pe,  label=r"$I_{\rm PE}(r)$", alpha=0.7)
    # Symmetric log handles tiny negative bands gracefully
    plt.yscale('symlog', linthresh=1e-20, linscale=1.0)
    plt.xlabel(r"$\tilde r$")
    plt.ylabel(r"$I(\tilde r)$")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    confirm(out)

def fig_potential_slice(sol, p, out="figs/bounce_potential_slice.png"):
    Phi_path = sol.y[0]; phi_path = sol.y[2]
    Phi_min, Phi_max = np.min(Phi_path), np.max(Phi_path)
    phi_min, phi_max = np.min(phi_path), np.max(phi_path)
    padP = 0.2*(Phi_max - Phi_min + 1e-9)
    padp = 0.2*(phi_max - phi_min + 1e-9)
    P = np.linspace(Phi_min-padP, Phi_max+padP, 260)
    Q = np.linspace(phi_min-padp, phi_max+padp, 260)
    PP, QQ = np.meshgrid(P, Q)
    V_false = potential_rescaled(p['Phi_FV_t'], p['phi_FV_t'], p)
    Vgrid = potential_rescaled(PP, QQ, p) - V_false

    plt.figure(figsize=(7.0,5.6))
    cs = plt.contourf(P, Q, Vgrid, levels=50)
    plt.colorbar(cs, label=r"$V(\Phi,\phi)-V_{\rm false}$")
    plt.plot(Phi_path, phi_path, 'k-', lw=2, label="bounce path")
    plt.xlabel(r"$\Phi$")
    plt.ylabel(r"$\phi$")
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    confirm(out)

def fig_residuals(sol, p, g_hist, out="figs/bounce_residuals.png"):
    r_prof = np.linspace(max(sol.x[1],1e-8), sol.x[-1], 2000)
    _, res_prof = residual_profile(p, sol, r_eval=r_prof)
    plt.figure(figsize=(7.6,5.0))
    ax1 = plt.gca()
    ax1.semilogy(r_prof, np.maximum(res_prof, 1e-20), label="|residual|(final g)")
    ax1.set_xlabel(r"$\tilde r$")
    ax1.set_ylabel("residual (semilogy)")
    if g_hist is not None and len(g_hist) > 1:
        ax2 = ax1.twinx()
        ax2.plot(g_hist, np.arange(len(g_hist)), 'o-', alpha=0.5)
        ax2.set_ylabel("continuation step index")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    confirm(out)

def fig_mesh(sol, out="figs/bounce_mesh.png"):
    r = sol.x
    if len(r) > 1:
        dr = np.diff(r); rc = 0.5*(r[1:] + r[:-1])
    else:
        dr = np.array([1.0]); rc = np.array([r[0]])
    plt.figure(figsize=(7.0,4.6))
    plt.semilogy(rc, np.maximum(dr, 1e-20))
    plt.xlabel(r"$\tilde r$")
    plt.ylabel(r"$\Delta \tilde r$ (node spacing)")
    plt.tight_layout()
    plt.savefig(out, dpi=200)
    plt.close()
    confirm(out)

# ----------------- Simulation -----------------
def run_simulation(v_phys, v_phi_phys, lam_phi_phys, lam_phys, g_portal_max, R0, w, r_max=100.0, tol=1e-5):
    print(f"\n{'='*62}\nRUN: v_phi={v_phi_phys:.2e}, lam_phi={lam_phi_phys:.3f}, lam={lam_phys:.3f}, g_max={g_portal_max:.2f}\n{'='*62}")
    p = {'v_phi_t': v_phi_phys/v_phys, 'lam_phi': lam_phi_phys, 'lam_t': lam_phys}
    p['mu2_t']   = lam_phys
    energy_diff  = (p['mu2_t']**2/p['lam_t']) + p['mu2_t'] + (p['lam_t']/4)
    p['bias_t']  = -1.01*energy_diff/p['v_phi_t']**2
    p['Phi_FV_t'] = 0.0
    p['phi_FV_t'] = np.sqrt(1 + 2*p['mu2_t']/p['lam_t'])
    p['Phi_TV_t'] = p['v_phi_t']
    p['phi_TV_t'] = 0.0
    p['g_portal_t'] = 0.0

    # wall-clustered composite mesh
    mesh = np.unique(np.concatenate([
        np.linspace(1e-9, max(1e-9,R0-6*w), 300),
        R0 + w*np.tanh(np.linspace(-4,4,2500)),
        np.linspace(R0+6*w, r_max, 300)
    ]))

    guess = initial_guess_tanh(mesh, p, R0, w)
    print("[bvp] solving at g=0 …")
    sol = solve_bvp(
        lambda r,y: bounce_odes_rescaled(r,y,p),
        lambda ya,yb: boundary_conditions_rescaled(ya,yb,p),
        mesh, guess, fun_jac=lambda r,y: bounce_odes_jac(r,y,p),
        verbose=1, max_nodes=500000, tol=tol
    )
    if not sol.success:
        print("!!! FAILURE: base solution (g=0) did not converge")
        return False, None, None, None

    R0_calc = wall_radius_from_solution(p, sol)
    S0_calc = calculate_action(p, sol)
    print(f"[cont] OK   g=0.000  R≈{R0_calc:.3f}  S_E≈{S0_calc:.6g}")

    g_hist = [0.0]; S_hist = [S0_calc]
    g = 0.0; step = max(0.04, g_portal_max/50)

    t0 = time.time()
    while g < g_portal_max - 1e-12:
        g_next = min(g + step, g_portal_max)
        p['g_portal_t'] = g_next
        guess_next = sol.y.copy()
        sol_next = solve_bvp(
            lambda r,y: bounce_odes_rescaled(r,y,p),
            lambda ya,yb: boundary_conditions_rescaled(ya,yb,p),
            sol.x, guess_next, fun_jac=lambda r,y: bounce_odes_jac(r,y,p),
            verbose=0, max_nodes=500000, tol=tol
        )
        if sol_next.success:
            sol = sol_next
            g = g_next
            S_now = calculate_action(p, sol)
            R_now = wall_radius_from_solution(p, sol)
            g_hist.append(g); S_hist.append(S_now)
            print(f" -> reached g={g:.4f}, R≈{R_now:.3f}, S_E={S_now:.4f}")
        else:
            step *= 0.5
            print(f" -> step reduced to {step:.5f}")
            if step < 1e-4:
                print(" -> continuation stalled; stopping.")
                break
    print(f"[done] max g={g:.4f}, S_E={S_hist[-1]:.4f}  (elapsed {time.time()-t0:.1f}s)")
    return True, sol, p.copy(), g_hist

# ----------------- Main -----------------
def main():
    parser = argparse.ArgumentParser(description="Two-field O(4) bounce with corrected action-density figure.")
    parser.add_argument("--v", type=float, default=4.2e-5, help="v (unit for tilde r)")
    parser.add_argument("--vphi", type=float, default=9.0e-5, help="v_phi (bulk)")
    parser.add_argument("--lamphi", type=float, default=0.1, help="lambda_Phi")
    parser.add_argument("--lam", type=float, default=1.3e-4, help="lambda (brane)")
    parser.add_argument("--gmax", type=float, default=2.0, help="target portal coupling")
    parser.add_argument("--R0", type=float, default=11.9, help="initial wall center")
    parser.add_argument("--w", type=float, default=3.0, help="initial wall width")
    parser.add_argument("--rmax", type=float, default=100.0, help="box size")
    parser.add_argument("--tol", type=float, default=1e-5, help="BVP tolerance")
    args = parser.parse_args()

    ensure_dirs()
    try:
        ok, sol, params, g_hist = run_simulation(
            v_phys=args.v, v_phi_phys=args.vphi, lam_phi_phys=args.lamphi,
            lam_phys=args.lam, g_portal_max=args.gmax, R0=args.R0, w=args.w,
            r_max=args.rmax, tol=args.tol
        )
        if not ok:
            sys.exit(2)
    except Exception:
        traceback.print_exc()
        sys.exit(1)

    # CSV
    try:
        rS = calculate_action(params, sol)
        print(f"[report] S_E (solver grid) ≈ {rS:.6g}")
        save_background_csv(sol, "data/background_profile.csv")
    except Exception:
        print("!!! Failed writing background_profile.csv")
        traceback.print_exc()

    # Figures
    for fn, f in [
        ("figs/bounce_fields.png",              lambda: fig_fields(sol, params)),
        ("figs/bounce_action_density.png",      lambda: fig_action_density(sol, params)),
        ("figs/bounce_potential_slice.png",     lambda: fig_potential_slice(sol, params)),
        ("figs/bounce_residuals.png",           lambda: fig_residuals(sol, params, g_hist)),
        ("figs/bounce_mesh.png",                lambda: fig_mesh(sol)),
    ]:
        try:
            f()
        except Exception:
            print(f"!!! Failed to create {fn}")
            traceback.print_exc()

    print("\n[summary]")
    print("  data/background_profile.csv  ->", "OK" if os.path.exists("data/background_profile.csv") else "MISSING")
    for img in [
        "figs/bounce_fields.png",
        "figs/bounce_action_density.png",
        "figs/bounce_potential_slice.png",
        "figs/bounce_residuals.png",
        "figs/bounce_mesh.png",
    ]:
        print(f"  {img:35s} -> {'OK' if os.path.exists(img) else 'MISSING'}")

if __name__ == "__main__":
    main()

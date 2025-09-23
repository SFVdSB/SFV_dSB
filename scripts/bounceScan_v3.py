import numpy as np
import matplotlib.pyplot as plt
import time
import csv
from datetime import datetime
from scipy.integrate import solve_bvp, trapezoid

# ==============================================================================
#  Bounce Solver with Parameter Scanning and Action Analysis
# ==============================================================================

def potential_rescaled(Phi_t, phi_t, p_t):
    V_phi = p_t['lam_phi']/4.0*(Phi_t**2 - p_t['v_phi_t']**2)**2 + p_t['bias_t']*Phi_t**2
    V_brane = p_t['lam_t']/4.0*(phi_t**2 - 1.0)**2 - p_t['mu2_t']*phi_t**2
    V_int = p_t['g_portal_t']*Phi_t**2*phi_t**2
    return V_phi + V_brane + V_int

def derivatives_rescaled(Phi_t, phi_t, p_t):
    # Clamp large values to prevent overflow
    Phi_t = np.clip(Phi_t, -1e6, 1e6)
    phi_t = np.clip(phi_t, -1e6, 1e6)
    dVdPhi_t = p_t['lam_phi']*Phi_t*(Phi_t**2 - p_t['v_phi_t']**2) \
               + 2*p_t['bias_t']*Phi_t + 2*p_t['g_portal_t']*Phi_t*phi_t**2
    dVdphi_t = p_t['lam_t']*phi_t*(phi_t**2 - 1) - 2*p_t['mu2_t']*phi_t \
               + 2*p_t['g_portal_t']*Phi_t**2*phi_t
    return dVdPhi_t, dVdphi_t

def bounce_odes_rescaled(r, y, p_t):
    Phi_t, Phi_p, phi_t, phi_p = y
    r_safe = np.maximum(r, 1e-9)
    dVdPhi_t, dVdphi_t = derivatives_rescaled(Phi_t, phi_t, p_t)
    dPhi_p_dr = dVdPhi_t - (3/r_safe)*Phi_p
    dphi_p_dr = dVdphi_t - (3/r_safe)*phi_p
    return np.vstack([Phi_p, dPhi_p_dr, phi_p, dphi_p_dr])

def bounce_odes_jac(r, y, p_t):
    Phi_t, _, phi_t, _ = y
    r_safe = np.maximum(r, 1e-9)
    d2VdPhi2 = p_t['lam_phi']*(3*Phi_t**2 - p_t['v_phi_t']**2) + 2*p_t['bias_t'] + 2*p_t['g_portal_t']*phi_t**2
    d2Vdphi2 = p_t['lam_t']*(3*phi_t**2 - 1) - 2*p_t['mu2_t'] + 2*p_t['g_portal_t']*Phi_t**2
    d2VdPhidphi = 4*p_t['g_portal_t']*Phi_t*phi_t
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

def boundary_conditions_rescaled(ya, yb, p_t):
    return np.array([
        ya[0]-p_t['Phi_TV_t'],
        ya[2]-p_t['phi_TV_t'],
        yb[0]-p_t['Phi_FV_t'],
        yb[2]-p_t['phi_FV_t']
    ])

def initial_guess_tanh(mesh, p_t, R0, w):
    Phi_guess = 0.5*(p_t['Phi_TV_t']-p_t['Phi_FV_t'])*(1-np.tanh((mesh-R0)/w)) + p_t['Phi_FV_t']
    phi_guess = 0.5*(p_t['phi_TV_t']-p_t['phi_FV_t'])*(1-np.tanh((mesh-R0)/w)) + p_t['phi_FV_t']
    dPhi_guess = np.gradient(Phi_guess, mesh)
    dphi_guess = np.gradient(phi_guess, mesh)
    dPhi_guess[0], dphi_guess[0] = 0, 0
    dPhi_guess[-1], dphi_guess[-1] = 0, 0
    Phi_guess[-1], phi_guess[-1] = p_t['Phi_FV_t'], p_t['phi_FV_t']
    return np.vstack([Phi_guess, dPhi_guess, phi_guess, dphi_guess])

def calculate_action(p_t, solution):
    r, y = solution.x, solution.y
    Phi, Phi_p, phi, phi_p = y
    V_false = potential_rescaled(p_t['Phi_FV_t'], p_t['phi_FV_t'], p_t)
    KE = 0.5*(Phi_p**2 + phi_p**2)
    PE = potential_rescaled(Phi, phi, p_t) - V_false
    return 2*np.pi**2*trapezoid(r**3*(KE+PE), r)

def analyze_solution(p_t, solution, g_history=None, action_history=None, title_prefix="", make_plots=False, save_prefix=None):
    """
    Print a brief analysis; only plot if make_plots=True.
    """
    r, y = solution.x, solution.y
    Phi, Phi_p, phi, phi_p = y
    S_E = calculate_action(p_t, solution)
    print(f"\n--- {title_prefix} Analysis ---")
    print(f"Euclidean Action: S_E={S_E:.4f}")

    if not make_plots:
        return

    V_false = potential_rescaled(p_t['Phi_FV_t'], p_t['phi_FV_t'], p_t)
    KE = 0.5*(Phi_p**2 + phi_p**2)
    PE = potential_rescaled(Phi, phi, p_t) - V_false

    plt.figure(figsize=(14, 10))

    plt.subplot(2,2,1)
    plt.plot(r, Phi, label='Phi')
    plt.plot(r, phi, label='phi')
    plt.legend(); plt.xlabel('r'); plt.ylabel('fields'); plt.grid(True)
    plt.title(f"{title_prefix} Field Profiles")

    plt.subplot(2,2,2)
    plt.plot(r, r**3*(KE+PE), label='Action Density')
    plt.grid(True); plt.xlabel('r'); plt.title(f"{title_prefix} r^3 (KE+ΔV)")
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(r, PE, label='PE = V - V_FV')
    plt.axhline(0, color='grey', linestyle='--')
    plt.grid(True); plt.xlabel('r'); plt.ylabel('ΔV'); plt.title(f"{title_prefix} Potential Energy")
    plt.legend()

    if g_history is not None and action_history is not None and len(g_history)==len(action_history):
        plt.subplot(2,2,4)
        plt.plot(g_history, action_history, 'o-')
        plt.xlabel('g_portal'); plt.ylabel('S_E'); plt.title('Action vs g')
        plt.grid(True)

    plt.tight_layout()
    if save_prefix:
        plt.savefig(f"{save_prefix}.png", dpi=150)
    else:
        plt.show()
def run_simulation(v_phys, v_phi_phys, lam_phi_phys, lam_phys, g_portal_max, R0, w):
    """
    Runs the entire bounce simulation for a given set of physical parameters.
    Returns success status, final solution, and action history.
    """
    start_time = time.time()
    print(f"\n{'='*60}")
    print(f"RUNNING: v_phi={v_phi_phys:.2e}, lam_phi={lam_phi_phys:.3f}, lam={lam_phys:.3f}, g_max={g_portal_max:.2f}")
    print(f"{'='*60}")

    # --- Setup dimensionless parameters ---
    params_tilde = {
        'v_phi_t': v_phi_phys/v_phys,
        'lam_phi': lam_phi_phys,
        'lam_t': lam_phys
    }
    params_tilde['mu2_t'] = lam_phys
    energy_diff = (params_tilde['mu2_t']**2/params_tilde['lam_t']) + params_tilde['mu2_t'] + (params_tilde['lam_t']/4)
    params_tilde['bias_t'] = -1.01*energy_diff/params_tilde['v_phi_t']**2
    params_tilde['Phi_FV_t'] = 0.0
    params_tilde['phi_FV_t'] = np.sqrt(1 + 2*params_tilde['mu2_t']/params_tilde['lam_t'])
    params_tilde['Phi_TV_t'] = params_tilde['v_phi_t']
    params_tilde['phi_TV_t'] = 0.0

    # --- Setup Mesh ---
    r_max = 100.0
    mesh = np.unique(np.concatenate([
        np.linspace(1e-9, max(1e-9,R0-6*w), 300),
        R0 + w*np.tanh(np.linspace(-4,4,2500)),
        np.linspace(R0+6*w, r_max, 300)
    ]))

    # --- Solve decoupled case g=0 ---
    params_decoupled = params_tilde.copy()
    params_decoupled['g_portal_t'] = 0.1
    guess = initial_guess_tanh(mesh, params_decoupled, R0, w)
   
    solution = solve_bvp(
        lambda r,y: bounce_odes_rescaled(r,y,params_decoupled),
        lambda ya,yb: boundary_conditions_rescaled(ya,yb,params_decoupled),
        mesh, guess,
        fun_jac=lambda r,y: bounce_odes_jac(r,y,params_decoupled),
        verbose=1, max_nodes=500000, tol=1e-5
    )
    if not solution.success:
        print("!!! FAILURE: Decoupled solution (g=0) not viable for these parameters !!!")
        return False, None, None, None, None
    print("Decoupled solution OK. Starting adaptive continuation...")

    # *** NEW: Analyze the initial g=0 solution to check for overshoot ***
    decoupled_params_full = params_decoupled.copy()
    decoupled_params_full['v_phys'] = v_phys
    decoupled_params_full['v_phi_phys'] = v_phi_phys
    decoupled_params_full['lam_phys'] = lam_phys
    analyze_solution(decoupled_params_full, solution, title_prefix="Initial g=0", make_plots=False)


    # --- Adaptive continuation in g_portal ---
    g_values = [0.0]
    action_values = [calculate_action(params_decoupled, solution)]
   
    final_params = params_decoupled.copy()
    g_step = g_portal_max/50
    g_val = 0.0
    while g_val < g_portal_max:
        g_val_next = min(g_val + g_step, g_portal_max)
        params_tilde['g_portal_t'] = g_val_next
        guess = solution.y.copy() # Use last successful solution as guess
        sol_next = solve_bvp(
            lambda r,y: bounce_odes_rescaled(r,y,params_tilde),
            lambda ya,yb: boundary_conditions_rescaled(ya,yb,params_tilde),
            solution.x, guess,
            fun_jac=lambda r,y: bounce_odes_jac(r,y,params_tilde),
            verbose=0, max_nodes=500000, tol=1e-5
        )
        if sol_next.success:
            solution = sol_next
            g_val = g_val_next
            final_params = params_tilde.copy()
           
            # --- Track action at each step ---
            current_action = calculate_action(params_tilde, solution)
            g_values.append(g_val)
            action_values.append(current_action)
            print(f" -> Successfully reached g_portal={g_val:.4f}, S_E={current_action:.4f}")
        else:
            g_step /= 2
            if g_step < 1e-4:
                print(f" -> Continuation failed at g ~ {g_val:.4f}. Stopping run for these parameters.")
                break
   
    end_time = time.time()
    final_action = action_values[-1]
    print(f"\nMax stable g_portal = {g_val:.4f}")
    print(f"Final Action S_E = {final_action:.4f}")
    print(f"Calculation took {end_time-start_time:.2f}s")
   
    return True, solution, final_params, g_values, action_values


# ==============================================================================
# Main Execution Block
# ==============================================================================
if __name__=="__main__":
   
    # --- Base Physics Parameters ---
    v_phys = 4.2e-5
   
    # --- Initial guess shape parameters ---
    R0 = 11.9
    w = 3.0

    # --- MODE 1: Run a single simulation (the "Golden Run") ---
    print("--- Starting Single 'Golden Run' Simulation ---")
    v_phi_phys_single = 5.4e-5
    lam_phi_phys_single = 0.1
    lam_phys_single = 0.01
    g_portal_max_single = 2.0
   
    success, final_solution, final_params, g_hist, action_hist = run_simulation(
        v_phys, v_phi_phys_single, lam_phi_phys_single, lam_phys_single,
        g_portal_max_single, R0, w
    )

    if success:
        # Add physical parameters needed for plotting title
        final_params['v_phys'] = v_phys
        final_params['v_phi_phys'] = v_phi_phys_single
        final_params['lam_phys'] = lam_phys_single
        analyze_solution(final_params, final_solution, g_history=g_hist, action_history=action_hist, title_prefix="Final g=2.0", make_plots=False)


    # --- MODE 2: Run a 2D parameter scan ---
    # Uncomment the block below to perform a detailed scan.
    # This will take a significant amount of time to run.

    print("\n\n--- Starting 2D Parameter Scan ---")
   
    # Define the parameter ranges to scan
    # *** NEW: High-resolution scan tightly focused around the golden run values ***
    v_phi_phys_scan = np.linspace(8.5e-5, 9.5e-5, 9) # 5 very small steps downwards
    lam_phys_scan = np.linspace(0.001, 0.0001, 30)      # 5 very small steps upwards

    results = []

    for v_phi_p in v_phi_phys_scan:
        for lam_p in lam_phys_scan:
           
            # Use fixed "golden run" values for other parameters during scan
            lam_phi_phys_scan_val = 0.1
            g_portal_max_scan_val = 2.0

            success, final_solution, final_params, g_hist, action_hist = run_simulation(
                v_phys, v_phi_p, lam_phi_phys_scan_val, lam_p,
                g_portal_max_scan_val, R0, w
            )
           
            if success:
                final_action = action_hist[-1]
                results.append({
                    'v_phys': v_phys,
                    'v_phi_phys': v_phi_p,
                    'lam_phi_phys': lam_phi_phys_scan_val,
                    'lam_phys': lam_p,
                    'g_portal_max': g_portal_max_scan_val,
                    'R0': R0,
                    'w': w,
                    'final_action': final_action,
                    'max_g': g_hist[-1]
                })

    print("\n\n--- Scan Complete ---")
    if results:
        print("Successful runs found:")
        sorted_results = sorted(results, key=lambda x: x['final_action'])
        for res in sorted_results:
            print(f"v_phi={res['v_phi_phys']:.2e}, lam={res['lam_phys']:.3f}, S_E={res['final_action']:.2f} (at g={res['max_g']:.2f})")
        # Write CSV of scan results
        csv_name = f"scan_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        fieldnames = ["v_phys","v_phi_phys","lam_phi_phys","lam_phys","g_portal_max","R0","w","SE","max_g"]
        with open(csv_name, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in results:
                writer.writerow({
                    "v_phys": row.get("v_phys", v_phys),
                    "v_phi_phys": row.get("v_phi_phys"),
                    "lam_phi_phys": row.get("lam_phi_phys", 0.1),
                    "lam_phys": row.get("lam_phys"),
                    "g_portal_max": row.get("g_portal_max", 2.0),
                    "R0": row.get("R0", R0),
                    "w": row.get("w", w),
                    "SE": row.get("final_action"),
                    "max_g": row.get("max_g")
                })
        print(f"[CSV] Wrote {len(results)} rows to {csv_name}")
    else:
        print("No successful bounces found in the scan.")

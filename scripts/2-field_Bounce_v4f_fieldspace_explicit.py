import numpy as np
# === Field-space contour with bounce path (auto-saves fig_fieldspace_path.pdf) ===
def plot_fieldspace_with_path(p_t, solution, fname="fig_fieldspace_path.pdf", levels=40):
    import numpy as np
    import matplotlib.pyplot as plt
    r = solution.x
    Phi, _, phi, _ = solution.y
    Phi_points = np.array([p_t['Phi_TV_t'], p_t['Phi_FV_t'], *Phi])
    phi_points = np.array([p_t['phi_TV_t'], -p_t['phi_FV_t'], *phi])  # negative branch
    padPhi = 0.15*(Phi_points.max() - Phi_points.min() + 1e-9)
    padphi = 0.15*(phi_points.max() - phi_points.min() + 1e-9)
    Phi_min, Phi_max = Phi_points.min()-padPhi, Phi_points.max()+padPhi
    phi_min, phi_max = phi_points.min()-padphi, phi_points.max()+padphi
    n = 400
    Phi_grid = np.linspace(Phi_min, Phi_max, n)
    phi_grid = np.linspace(phi_min, phi_max, n)
    PP, pp = np.meshgrid(Phi_grid, phi_grid, indexing='xy')
    V = potential_rescaled(PP, pp, p_t)
    V_false = potential_rescaled(p_t['Phi_FV_t'], p_t['phi_FV_t'], p_t)
    Vrel = V - V_false
    fig, ax = plt.subplots(figsize=(7.5, 6.0), constrained_layout=True)
    ax.contourf(PP, pp, Vrel, levels=levels)
    ax.contour(PP, pp, Vrel, colors='k', linewidths=0.5, alpha=0.5)
    ax.plot(Phi, phi, lw=2.0, label='bounce path')
    ax.scatter([p_t['Phi_TV_t']], [p_t['phi_TV_t']], marker='*', s=180, label='TV')
    ax.scatter([p_t['Phi_FV_t']], [-p_t['phi_FV_t']], marker='o', s=90, label='FV (neg branch)')
    ax.set_xlabel(r'$\Phi$'); ax.set_ylabel(r'$\phi$')
    ax.set_title('Field-space $V(\\Phi,\\phi)-V_{\\rm false}$ with bounce path')
    ax.legend(loc='best')
    fig.savefig(fname, bbox_inches='tight')
    print("[info] Saved field-space plot:", fname)
    plt.close(fig)

import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_bvp, trapezoid

# ==============================================================================
#  Bounce Solver with Corrected Boundary Condition
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
    """
    Boundary conditions for the bounce solution.
    Crucially, we now force the phi field to settle in the NEGATIVE vacuum.
    """
    return np.array([
        ya[0]-p_t['Phi_TV_t'],
        ya[2]-p_t['phi_TV_t'],
        yb[0]-p_t['Phi_FV_t'],
        # The key change is here: yb[2] should equal -phi_FV_t
        yb[2]+p_t['phi_FV_t']
    ])

def initial_guess_tanh(mesh, p_t, R0, w):
    # The initial guess should also point towards the negative vacuum
    phi_FV_target = -p_t['phi_FV_t']
    Phi_guess = 0.5*(p_t['Phi_TV_t']-p_t['Phi_FV_t'])*(1-np.tanh((mesh-R0)/w)) + p_t['Phi_FV_t']
    phi_guess = 0.5*(p_t['phi_TV_t']-phi_FV_target)*(1-np.tanh((mesh-R0)/w)) + phi_FV_target
    dPhi_guess = np.gradient(Phi_guess, mesh)
    dphi_guess = np.gradient(phi_guess, mesh)
    dPhi_guess[0], dphi_guess[0] = 0, 0
    dPhi_guess[-1], dphi_guess[-1] = 0, 0
    Phi_guess[-1], phi_guess[-1] = p_t['Phi_FV_t'], phi_FV_target
    return np.vstack([Phi_guess, dPhi_guess, phi_guess, dphi_guess])

def calculate_action(p_t, solution):
    r, y = solution.x, solution.y
    Phi, Phi_p, phi, phi_p = y
    # IMPORTANT: The false vacuum potential V_false is the same for both +/- phi_FV
    # because V depends on phi^2 and phi^4.
    V_false = potential_rescaled(p_t['Phi_FV_t'], p_t['phi_FV_t'], p_t)
    KE = 0.5*(Phi_p**2 + phi_p**2)
    PE = potential_rescaled(Phi, phi, p_t) - V_false
    return 2*np.pi**2*trapezoid(r**3*(KE+PE), r)

def analyze_solution(p_t, solution, g_history=None, action_history=None, title_prefix=""):
    """Generates detailed plots for a single successful solution."""
    r, y = solution.x, solution.y
    Phi, Phi_p, phi, phi_p = y
    S_E = calculate_action(p_t, solution)
    print(f"\n--- {title_prefix} Analysis ---")
    print(f"Euclidean Action: S_E={S_E:.4f}")

    plt.figure(figsize=(14, 10))
    plt.subplot(2,2,1)
    plt.plot(r, Phi, label='Phi')
    plt.plot(r, phi, label='phi')
    plt.axhline(-p_t['phi_FV_t'], color='grey', linestyle='--', label='Target phi_FV')
    plt.legend()
    plt.title(f'{title_prefix} Field Profiles')
    plt.grid(True)

    V_false = potential_rescaled(p_t['Phi_FV_t'], p_t['phi_FV_t'], p_t)
    KE = 0.5*(Phi_p**2+phi_p**2)
    PE = potential_rescaled(Phi, phi, p_t) - V_false
    plt.subplot(2,2,2)
    plt.plot(r, r**3*(KE+PE), label='Action Density')
    plt.grid(True)
    plt.title(f'{title_prefix} r^3 * (KE + PE)')
    plt.legend()
    
    plt.subplot(2,2,3)
    plt.plot(r, PE, label='PE = V(r) - V_false')
    plt.axhline(0, color='grey', linestyle='--')
    plt.grid(True)
    plt.title(f'{title_prefix} Potential Energy Diagnostic')
    plt.legend()

    if g_history and action_history:
        plt.subplot(2,2,4)
        plt.plot(g_history, action_history, 'o-')
        plt.title("Action vs. Portal Coupling")
        plt.grid(True)

    plt.suptitle(f"{title_prefix} g_portal={p_t['g_portal_t']:.2f}, v_phi={p_t['v_phi_phys']:.2e}, lam={p_t['lam_phys']:.3f} -> S_E={S_E:.2f}", fontsize=16)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()

def run_simulation(v_phys, v_phi_phys, lam_phi_phys, lam_phys, g_portal_max, R0, w):
    start_time = time.time()
    print(f"\n{'='*60}")
    print(f"RUNNING: v_phi={v_phi_phys:.2e}, lam_phi={lam_phi_phys:.3f}, lam={lam_phys:.3f}, g_max={g_portal_max:.2f}")
    print(f"{'='*60}")

    params_tilde = {
        'v_phi_t': v_phi_phys/v_phys, 'lam_phi': lam_phi_phys, 'lam_t': lam_phys
    }
    params_tilde['mu2_t'] = lam_phys
    energy_diff = (params_tilde['mu2_t']**2/params_tilde['lam_t']) + params_tilde['mu2_t'] + (params_tilde['lam_t']/4)
    params_tilde['bias_t'] = -1.01*energy_diff/params_tilde['v_phi_t']**2
    params_tilde['Phi_FV_t'] = 0.0
    params_tilde['phi_FV_t'] = np.sqrt(1 + 2*params_tilde['mu2_t']/params_tilde['lam_t'])
    params_tilde['Phi_TV_t'] = params_tilde['v_phi_t']
    params_tilde['phi_TV_t'] = 0.0

    r_max = 100.0
    mesh = np.unique(np.concatenate([
        np.linspace(1e-9, max(1e-9,R0-6*w), 300),
        R0 + w*np.tanh(np.linspace(-4,4,2500)),
        np.linspace(R0+6*w, r_max, 300)
    ]))

    params_decoupled = params_tilde.copy()
    params_decoupled['g_portal_t'] = 0.0
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

    decoupled_params_full = params_decoupled.copy()
    decoupled_params_full.update({'v_phys': v_phys, 'v_phi_phys': v_phi_phys, 'lam_phys': lam_phys})
    analyze_solution(decoupled_params_full, solution, title_prefix="Initial g=0")

    g_values = [0.0]
    action_values = [calculate_action(params_decoupled, solution)]
    
    final_params = params_decoupled.copy()
    g_step = g_portal_max/50
    g_val = 0.0
    while g_val < g_portal_max:
        g_val_next = min(g_val + g_step, g_portal_max)
        params_tilde['g_portal_t'] = g_val_next
        guess = solution.y.copy()
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
            current_action = calculate_action(params_tilde, solution)
            g_values.append(g_val)
            action_values.append(current_action)
            print(f" -> Successfully reached g_portal={g_val:.4f}, S_E={current_action:.4f}")
        else:
            g_step /= 2
            if g_step < 1e-4:
                print(f" -> Continuation failed at g ~ {g_val:.4f}.")
                break
    
    end_time = time.time()
    final_action = action_values[-1]
    print(f"\nMax stable g_portal = {g_val:.4f}")
    print(f"Final Action S_E = {final_action:.4f}")
    print(f"Calculation took {end_time-start_time:.2f}s")
    
    return True, solution, final_params, g_values, action_values

if __name__=="__main__":
    v_phys = 4.2e-5
    R0 = 11.9
    w = 3.0

    print("--- Starting Single 'Golden Run' Simulation ---")
    v_phi_phys_single = 9e-5
    lam_phi_phys_single = 0.1
    lam_phys_single = 0.00015
    g_portal_max_single = 2.0
    
    success, final_solution, final_params, g_hist, action_hist = run_simulation(
        v_phys, v_phi_phys_single, lam_phi_phys_single, lam_phys_single,
        g_portal_max_single, R0, w
    )

    if success:
        final_params.update({'v_phys': v_phys, 'v_phi_phys': v_phi_phys_single, 'lam_phys': lam_phys_single})
        analyze_solution(final_params, final_solution, g_hist, action_hist, title_prefix="Final g=2.0")
        # Generate field-space plot
        plot_fieldspace_with_path(final_params, final_solution)

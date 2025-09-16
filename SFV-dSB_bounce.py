import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import solve_bvp, trapezoid

# ============================================================================== 
#  Golden Run Script with Improved Continuation
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

def analyze_solution(p_t, solution):
    r, y = solution.x, solution.y
    Phi, Phi_p, phi, phi_p = y
    S_E = calculate_action(p_t, solution)
    print(f"Euclidean Action: S_E={S_E:.4f}")
    plt.figure(figsize=(12,5))
    plt.subplot(1,2,1)
    plt.plot(r, Phi, label='Phi'); plt.plot(r, phi, label='phi'); plt.legend()
    plt.title('Field Profiles'); plt.grid(True)
    plt.subplot(1,2,2)
    KE = 0.5*(Phi_p**2+phi_p**2)
    PE = potential_rescaled(Phi, phi, p_t)-potential_rescaled(p_t['Phi_FV_t'], p_t['phi_FV_t'], p_t)
    plt.plot(r, r**3*(KE+PE), label='Action Density'); plt.grid(True); plt.title('r^3*(KE+PE)'); plt.legend()
    plt.tight_layout(); plt.show()

# ============================================================================== 
# Main Execution
# ============================================================================== 
if __name__=="__main__":
    start_time = time.time()
    
    # --- Physics parameters ---
    v_phys = 4.2e-5
    v_phi_phys = 5.4e-5
    lam_phi_phys = 0.1
    lam_phys = 0.01
    g_portal_max = 2.0
    R0 = 11.9
    w = 3.0

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

    # --- Mesh ---
    r_max = 100.0
    mesh = np.unique(np.concatenate([
        np.linspace(1e-9, max(1e-9,R0-6*w),300),
        R0 + w*np.tanh(np.linspace(-4,4,2500)),
        np.linspace(R0+6*w, r_max,300)
    ]))

    # --- Solve decoupled case g=0 ---
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
        print("!!! FAILURE: Decoupled solution not viable !!!")
        exit()
    print("Decoupled solution OK. Starting adaptive continuation...")

    # --- Adaptive continuation in g_portal ---
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
            print(f" -> Successfully reached g_portal={g_val:.4f}")
        else:
            g_step /= 2
            if g_step < 1e-4:
                print(f" -> Continuation failed at g ~ {g_val:.4f}")
                break

    end_time = time.time()
    print(f"\nMax stable g_portal = {g_val:.4f} (calculation took {end_time-start_time:.2f}s)")
    analyze_solution(params_tilde, solution)

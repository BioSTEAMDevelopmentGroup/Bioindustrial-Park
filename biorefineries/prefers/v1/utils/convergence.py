"""
Shared convergence utility for iterative titer optimization in biorefineries.
Encapsulates the standard loop: Titer Adjustment -> NH3 Optimization -> Production Rate Scaling -> Simulation.
"""
import biosteam as bst

def run_titer_convergence(
        system, 
        target_production, 
        target_titer_func=None, # Function to get/set target titer, or we assume it's on R302 unit
        adjust_glucose_func=None,
        optimize_nh3_func=None,
        set_prod_func=None,
        max_iter=15, 
        tol_pct=0.01,
        verbose=True
    ):
    """
    Executes the iterative convergence loop for titer and production rate.
    
    Parameters
    ----------
    system : bst.System
        The system object to simulate.
    target_production : float
        Target production rate in kg/hr (or system-specific unit).
    target_titer_func : callable or None
        Optional. If provided, used to get the target titer and actual titer for checking.
        If None, tries to find 'R302' in system.flowsheet.unit and use its .titer and .actual_titer.
        Return tuple: (actual_titer, target_titer)
    adjust_glucose_func : callable
        Function to adjust glucose/yield. Signature: f(system, verbose=bool)
    optimize_nh3_func : callable
        Function to optimize ammonia loading. Signature: f(system, verbose=bool)
    set_prod_func : callable
        Function to set production rate. Signature: f(system, target, verbose=bool)
    max_iter : int
        Maximum number of iterations.
    tol_pct : float
        Tolerance percentage for convergence (e.g. 0.01 for 1%).
    verbose : bool
        Whether to print status logs.
        
    Returns
    -------
    bool
        True if converged, False otherwise.
    """
    
    log = print if verbose else (lambda *args, **kwargs: None)
    
    # Validation
    if not adjust_glucose_func or not optimize_nh3_func or not set_prod_func:
        raise ValueError("All adjustment functions (adjust_glucose, optimize_nh3, set_prod) must be provided.")

    # Helper to get titers if no func provided
    def default_get_titers():
        try:
            R302 = system.flowsheet.unit.R302
            return R302.actual_titer, R302.target_titer
        except AttributeError:
            # Fallback or strict error
            raise AttributeError("Could not default-resolve R302 in system. Provide target_titer_func.")
            
    get_titers = target_titer_func if target_titer_func else default_get_titers
    
    # Initial Titer Check for Log
    try:
        _, target_titer_val = get_titers()
        log(f"\n2. Starting Iterative Convergence (Target Titer: {target_titer_val} g/L)...")
    except:
        log("\n2. Starting Iterative Convergence...")
    
    converged = False
    
    for i in range(max_iter):
        log(f"\n   --- Titer Convergence Iteration {i+1}/{max_iter} ---")
        
        # 1. Adjust Titer (Yield)
        adjust_glucose_func(system, verbose=verbose)
        
        # 2. Optimize NH3
        optimize_nh3_func(system, verbose=False)
        
        # 3. Set Production Rate
        # Verbose only on first iteration to avoid spamming scaler logs
        # But if global verbose is False, force False.
        set_prod_func(system, target_production, verbose=(verbose and (i==0)))
        
        # 4. Optimize NH3 again (after scaling)
        optimize_nh3_func(system, verbose=False)
        
        # 5. Simulate
        try:
            system.simulate()
        except Exception as e:
            log(f"   [ERROR] Simulation failed in iteration {i+1}: {e}")
            raise e
        
        # 6. Check Convergence
        actual_titer, target_titer = get_titers()
        
        # Avoid division by zero or None
        if target_titer is None or target_titer == 0:
             log(f"   [WARN] Target titer is None or zero. Cannot check convergence.")
             break
             
        tolerance_val = target_titer * tol_pct
        diff = abs(actual_titer - target_titer)
        
        if diff <= tolerance_val:
            log(f"\n   [CONVERGED] Titer: {actual_titer:.4f} g/L (Target: {target_titer} g/L, Tol: {tolerance_val:.4f})")
            converged = True
            break
        else:
            log(f"   [CONTINUING] Titer gap: {diff:.4f} g/L (Tol: {tolerance_val:.4f})")
            
    if not converged:
        log(f"\n   [WARN] Max iterations ({max_iter}) reached. Titer may not be converged.")
        
    return converged

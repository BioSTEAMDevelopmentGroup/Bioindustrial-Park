# Convergence Optimization Brainstorming

## 1. Numba Suitability Analysis
**Verdict: Not Recommended for this specific use case.**

### Why?
Numba is an excellent tool for "Just-In-Time" (JIT) compilation of numerical Python code. It works best when you have a function that performs heavy mathematical operations on NumPy arrays (e.g., matrix multiplication, custom loop logic on raw numbers).

However, your `run_titer_convergence` loop and the BioSTEAM simulations (`system.simulate()`) behave very differently:
1.  **Object-Oriented Overhead**: The code relies heavily on `biosteam` objects (`Unit`, `Stream`, `System`). Numba cannot compile arbitrary Python classes or dynamic method calls. It only understands specific types (integers, floats, arrays).
2.  **Inner Complexity**: The bottleneck is likely inside `system.simulate()`, which executes thousands of thermodynamic calculations (VLE, LLE) and unit operation mass balances. Optimizing the *outer* loop with Numba won't speed up the *inner* simulation.
3.  **Multiprocessing Compatibility**: While Numba releases the GIL (Global Interpreter Lock), allowing for better threading, your current setup uses `multiprocess` (process-based parallelism). Since each process has its own memory space and Python interpreter, Numba wouldn't provide additional concurrency benefits over what you already have.

**Conclusion**: Converting the high-level control logic in `convergence.py` to Numba would be extremely difficult (requiring a rewrite of BioSTEAM internals) and likely yield negligible performance gains.

---

## 2. Recommended Optimization Strategies

Since Numba is off the table, here are practical ways to accelerate your simulation without breaking existing multiprocessing logic.

### A. Algorithmic Optimization (Highest Impact)
The current method in `convergence.py` likely uses a custom "stepping" logic (coarse steps -> fine steps) to find the target titer/production. This effectively mimics a linear search or a primitive root-finding algorithm.

**Recommendation:** Replace the manual stepping logic with a robust numerical root-finding algorithm.
*   **Tool**: `scipy.optimize.root_scalar` (specifically **Brent's Method** or **Secant Method**).
*   **How it helps**: Instead of taking fixed steps (e.g., `+1%`, then `+0.1%`), Brent's method analyzes the function's slope and "jumps" directly to the estimated solution. It typically converges in 3-5 iterations compared to 10-20 manual steps.

**Proposed Change:**
```python
from scipy import optimize

def objective_function(titer_guess):
    # Set titer
    target_titer_func(titer_guess) 
    # Adjust glucose/NH3
    adjust_glucose_func(system)
    optimize_nh3_func(system)
    # Run sim
    system.simulate()
    # Calculate error (Production - Target)
    current_prod = some_get_prod_func()
    return current_prod - target_production

# Solve
sol = optimize.root_scalar(objective_function, bracket=[min_titer, max_titer], method='brentq')
```

### B. "Warm Start" for Monte Carlo
In Monte Carlo simulations, samples are often independent, but if you sort your samples (e.g., based on Titer or Yield), the solution for `Sample N` is likely very close to the solution for `Sample N-1`.

**Recommendation:**
*   Sort the samples by key parameters (if possible) before distribution.
*   Use the *final converged values* of the previous run as the *initial guess* for the current run.
*   **Caution**: In `multiprocess`, sharing state between workers is tricky. However, each worker processes a "chunk" of samples sequentially. You can implement warm-starting *within* each worker's loop.

### C. Relaxing Tolerances
You explicitly mentioned maintaining accuracy, but "accuracy" is relative.
*   BioSTEAM default tolerance (`molar_tolerance`) might be tighter than required for high-level techno-economic analysis (TEA).
*   **Recommendation**: Check `system.molar_tolerance` and `system.relative_tolerance`.
    *   If `molar_tolerance` is `1e-6` kmol/hr, broadening it to `1e-4` or `1e-3` can drastically reduce recycle loop iterations without meaningfully changing MSP or GWP.

### D. Profiling (Know the Enemy)
Before making complex changes, identify exactly *what* is slow.
*   **Action**: Run a single simulation with `cProfile`.
    ```python
    import cProfile
    import pstats
    
    with cProfile.Profile() as pr:
        model(sample) # Run one sample
    
    ps = pstats.Stats(pr).sort_stats('cumtime')
    ps.print_stats(20) # See top 20 time-consumers
    ```
*   **Look for**:
    *   Specific Unit Operations (e.g., `Flash`, `Distillation`).
    *   Thermodynamic property calls (e.g., `activity_coefficients`).
    *   If `optimize_NH3` takes 80% of the time, focus on optimizing *that* specific function (maybe call it every 3rd iteration instead of every iteration).

### E. Multiprocessing Tweaks
*   **Batch Size**: You are using `gen_data_mc.py`. Ensure your `chunk_size` is optimal. Too small = overhead from spawning processes. Too large = memory hoarding.
*   **Avoid Serialize/Deserialize Overhead**: Ensure that the `model` object doesn't contain massive unnecessary data that gets pickled/unpickled for every worker.

## 3. Summary of Action Plan

1.  **Do NOT use Numba** for this.
2.  **Profile** one run to find the bottleneck.
3.  **Refactor `convergence.py`** to use `scipy.optimize.root_scalar` (Brent's method) instead of manual stepping. This is safe for multiprocessing and likely provides the biggest speedup.
4.  **Implement Warm Starts** inside the worker loop (reuse converged state from previous sample in the chunk).

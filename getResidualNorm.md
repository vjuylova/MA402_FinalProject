# `KSP.getResidualNorm()`

**Documentation for petsc4py KSP method**

---

## Function Signature

```python
ksp.getResidualNorm()
```

---

## Summary

This returns the L2 norm of the residual vector at the end of the most recent KSP solve. The residual measures how well the computed solution satisfies the linear system, and is a way to measure of solution accuracy.

---

## Parameters

This method takes no arguments.

---

## Returns

**rnorm** : float
    The L2 (Euclidean) norm of the final residual vector: ‖r‖₂ where r = b - A·x.
    Returns 0.0 if `solve()` has not yet been called.

---

## Mathematical Definition

For the linear system:

```
A·x = b
```
The residual vector is defined as

```
r = b - A·x
```

And it measures the error in satisfying the equation. If r = 0, then x is the exact solution. Usually in practice, iterative solvers just stop when ‖r‖ is "small enough."

The L2 norm (Euclidean norm) of the residual is:

```
‖r‖₂ = √(r₁² + r₂² + ... + rₙ²) = √(rᵀ·r)
```

This is the value returned by `getResidualNorm()`.

The KSP solver uses the residual norm to decide when to stop. You configure this via `setTolerances()`:

```python
ksp.setTolerances(rtol=1e-8, atol=1e-12, dtol=1e30, max_it=1000)
```

### Stopping criteria:

1. **Relative tolerance (rtol)**: Stops when ‖r‖ / ‖b‖ < rtol

2. **Absolute tolerance (atol)**: Stops when ‖r‖ < atol

3. **Divergence tolerance (dtol)**: Stops (diverges) if ‖r‖ / ‖r₀‖ > dtol

4. **Max iterations (max_it)**: Stops after this many iterations

---

## Debugging with Residual Norm

### Case 1: The residual is not decreasing

```python
ksp.solve(b, x)
if ksp.getResidualNorm() > 1e-3:
    print("Residual is large! Possible causes:")
    print("  1. Matrix is singular or nearly singular")
    print("  2. Boundary conditions not applied correctly")
    print("  3. Initial guess is far from solution")
    print("  4. Preconditioner is ineffective")
```

### Case 2: Residual is NaN or Inf

```python
rnorm = ksp.getResidualNorm()
if np.isnan(rnorm) or np.isinf(rnorm):
    print("Residual is NaN/Inf!")
    print("Check for:")
    print("  - Zero diagonal entries in matrix")
    print("  - Unbounded growth in iterations")
    print("  - Numerical overflow in matrix assembly")
```

### Case 3: Residual oscillating

Some Krylov methods (like BiCGStab) can have non-monotone residual behavior:

```
Iter 10: ‖r‖ = 1.2e-03
Iter 11: ‖r‖ = 3.4e-03   ← increased!
Iter 12: ‖r‖ = 8.9e-04   ← decreased again
```

If this happens, use `ksp.setMonitor()` to track the full history.


---

## Residual Check After Solve

It's good practice to **manually verify** the residual after solving:

```python
ksp.solve(b, x)
reported_rnorm = ksp.getResidualNorm()

---

## C Function Mapping

This Python method wraps the PETSc C function:

```c
PetscErrorCode KSPGetResidualNorm(KSP ksp, PetscReal *rnorm)
```

**GitLab source**: [src/ksp/ksp/interface/itfunc.c](https://petsc.org/release/src/ksp/ksp/interface/iterativ.c.html#KSPGetResidualNorm)


---

## The Original DOcumentaion

```c
/*@
 13:   KSPGetResidualNorm - Gets the last (possibly approximate and/or preconditioned) residual norm that has been computed.

 15:   Not Collective

 17:   Input Parameter:
 18: . ksp - the iterative context

 20:   Output Parameter:
 21: . rnorm - residual norm

 23:   Level: intermediate

 25:   Notes:
 26:   For some methods, such as `KSPGMRES`, the norm is not computed directly from the residual.

 28:   The type of norm used by the method can be controlled with `KSPSetNormType()`

 30:   Certain solvers, under certain conditions, may not compute the final residual norm in an iteration, in that case the previous norm is returned.

 32: .seealso: [](ch_ksp), `KSP`, `KSPSetNormType()`, `KSPBuildResidual()`, `KSPNormType`
 33: @*/
 34: PetscErrorCode KSPGetResidualNorm(KSP ksp, PetscReal *rnorm)
 35: {
 36:   PetscFunctionBegin;
 38:   PetscAssertPointer(rnorm, 2);
 39:   *rnorm = ksp->rnorm;
 40:   PetscFunctionReturn(PETSC_SUCCESS);
 41: }

```

---

## References

- **PETSc Manual**: [KSP Convergence](https://petsc.org/release/manual/ksp/#ksp-convergence)
- **C API**: [KSPGetResidualNorm](https://petsc.org/release/manualpages/KSP/KSPGetResidualNorm/)
- **petsc4py docs**: [KSP class reference](https://petsc.org/release/src/ksp/ksp/interface/iterativ.c.html#KSPGetResidualNorm)
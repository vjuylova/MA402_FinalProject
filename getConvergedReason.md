# `KSP.getConvergedReason()`

**Documentation for petsc4py KSP method**

---

## Function Signature

```python
ksp.getConvergedReason()
```

---

## Summary

This returns an integer code indicating why the PETSc KSP (Krylov Subspace) linear solver terminated. This function tells you whether the solver successfully converged to a solution, diverged, or stopped for other reasons.

---

## Parameters

This method takes no arguments.

---

## Returns

**reason** : int
    An integer that follows convergence code. Positive values mean a successful convergence, negative values mean divergence or failure, and zero indicates the solver hasnt run yet.

---

## Context behind this function

The KSP solver attempts to solve the linear system below

```
A·x = b
```

where **A** is the stiffness matrix, **b** is the load vector, and **x** is the unknown solution. The solver is iterative, and would hopefully converge to the true solution.

At each iteration k, the solver computes the **residual**:

```
rₖ = b - A·xₖ
```

And the residual measures how close xₖ is to satisfying the equation. The solver stops when one of the following occurs:

1. ‖rₖ‖ < rtol · ‖b‖ (residual is small compared to RHS)
2. ‖rₖ‖ < atol (residual is small in absolute terms)
3. k ≥ max_iterations
4. ‖rₖ‖ is diverging or NaN/Inf encountered


---

## C Function Mapping

This Python method wraps the PETSc C function:

```c
PetscErrorCode KSPGetConvergedReason(KSP ksp, KSPConvergedReason *reason)
```

**GitLab source**: [src/ksp/ksp/interface/itfunc.c](https://petsc.org/release/src/ksp/ksp/interface/iterativ.c.html#KSPGetConvergedReason)


## The Original Documentation

```c
/*@
1890:   KSPGetConvergedReason - Gets the reason the `KSP` iteration was stopped.

1892:   Not Collective

1894:   Input Parameter:
1895: . ksp - the `KSP` context

1897:   Output Parameter:
1898: . reason - negative value indicates diverged, positive value converged, see `KSPConvergedReason` for the possible values

1900:   Options Database Key:
1901: . -ksp_converged_reason - prints the reason to standard out when the solve ends

1903:   Level: intermediate

1905:   Note:
1906:   If this routine is called before or doing the `KSPSolve()` the value of `KSP_CONVERGED_ITERATING` is returned

1908: .seealso: [](ch_ksp), `KSPConvergedReason`, `KSP`, `KSPSetConvergenceTest()`, `KSPConvergedDefault()`, `KSPSetTolerances()`,
1909:           `KSPConvergedReasonView()`, `KSPGetConvergedReasonString()`
1910: @*/
1911: PetscErrorCode KSPGetConvergedReason(KSP ksp, KSPConvergedReason *reason)
1912: {
1913:   PetscFunctionBegin;
1915:   PetscAssertPointer(reason, 2);
1916:   *reason = ksp->reason;
1917:   PetscFunctionReturn(PETSC_SUCCESS);
1918: }
```


---

## References

**PETSc Manual**: [KSP Chapter](https://petsc.org/release/manual/ksp/)
**C API**: [KSPGetConvergedReason](https://gitlab.com/petsc/petsc/-/blob/main/src/binding/petsc4py/src/petsc4py/PETSc/KSP.pyx?ref_type=heads#L1890)
**petsc4py docs**: [KSP class reference](https://petsc.org/release/src/ksp/ksp/interface/iterativ.c.html#KSPGetConvergedReason)

# `KSP.getIterationNumber()`

**Documentation for petsc4py KSP method**

---

## Function Signature

```python
ksp.getIterationNumber()
```

---

## Summary

This function returns the number of iterations executed by the PETSc KSP linear solver during the most recent call to `solve()`. It would tell you how much work the iterative solver performed before converging (or diverging, so that you can 

---

## Parameters

This method takes no arguments.

---

## Returns

**n_iter** : int
    An integer for the number of iterations completed. It will return 0 if `solve()` has not yet been called.

---

## Mathematical Context

Iterative Krylov subspace solvers (like GMRES, CG, BiCGStab) solve the linear system **A·x = b** by generating a sequence of improving approximations:

```
x₀ → x₁ → x₂ → ... → xₖ
```

Each step is called an **iteration**. The solver will stop when the residual ‖b - A·xₖ‖ is really small (and converged), or the maximum iteration count is reached, or divergence is detected.

The iteration count **k** tells you how hard the problem was (because well-conditioned systems converge in less iterations), how effective the preconditioner was (because good preconditioners reduce iteration count by a lot), or whether you hit the iteration limit.


## Factors That Will Effect Iteration Count

### 1. Matrix Conditioning
A **condition number** κ(A) = ‖A‖ · ‖A⁻¹‖ measures how sensitive the solution is to perturbations.

**κ(A) ≈ 1–100**: Well-conditioned

However,
**κ(A) > 10⁸**: Severely ill-conditioned 

**In elasticity** the condition number grows as
- The mesh is refined (more elements)
- The material becomes stiffer (larger μ, λ)
- Or the aspect ratio of elements increases

### 2. Preconditioner Quality
A **preconditioner** P transforms the system:

```
A·x = b  →  P⁻¹A·x = P⁻¹b
```

to reduce the condition number of P⁻¹A, accelerating convergence. This could include like using a Jacobian and only focus on the eigenvalues. This will effect how many iterations you will need, usually making it less.

### 3. Mesh Resolution
For FEM problems, the finer meshes mean larger systems which will call for more iterations. Having less of a resolution will give you less iterations and reduce likeliness of a divergence. 


---

## C Function Mapping

This Python method wraps the PETSc C function:

```c
PetscErrorCode KSPGetIterationNumber(KSP ksp, PetscInt *its)
```

**GitLab source**: [src/ksp/ksp/interface/itfunc.c](https://petsc.org/main/src/ksp/ksp/interface/iterativ.c.html#KSPGetIterationNumber)


---

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

- **PETSc Manual**: [KSP Chapter](https://petsc.org/release/manual/ksp/)
- **C API**: [KSPGetIterationNumber](https://gitlab.com/petsc/petsc/-/blob/main/src/binding/petsc4py/src/petsc4py/PETSc/KSP.pyx?ref_type=heads#L1890)
- **petsc4py docs**: [KSP class reference](https://petsc.org/release/src/ksp/ksp/interface/iterativ.c.html#KSPGetIterationNumber)

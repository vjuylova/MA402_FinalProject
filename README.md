# PETSc Example 17 - Linear Elasticity (petsc4py Translation)

This directory contains Python translations of the PETSc C example `ex17.c`, which demonstrates solving linear elasticity problems using finite elements with the DMPLEX unstructured mesh framework.

## Mathematical Context 

In linear elasticity the governing equation is:
```
∇·σ + f = 0
```
Which tells you how to find the euqlibrium state of a material. 

The stress tensor is:
```
σ_ij = λ δ_ij ε_kk + 2μ ε_ij
ε_ij = 1/2 (∂u_i/∂x_j + ∂u_j/∂x_i)
```

And you can use these to understand how an elastic material would deform under a force. 

The example solves the linear elasticity problem where 
- There is a eectangular domain [0,1]^2
- Solves the PDE ∇·σ = f, where σ is the stress tensor and f is the body force
- Uses an unstructured tetrahedral mesh with DMPLEX
- And uses a nonlinear solver (SNES) 

The implementation supports 2d quadratic linearly elastatic problem.

## My Experience, and AI Usage

I chose the PETSc C example `ex17.c`, because it simulates a lot of the work I do as a material science and soft matter physicist. In a lot of my work with characterizing granular materials, I have to use soft materials as my sample containers when compressing my samples. Understanding how my sample container deforms under pressure (or a constant vertical force) is crucial for properly designing my expiremental set up. 

To translate the C code, I used Claude Sonnet with the prompt: "Translate this PETSc C tutorial into a working Python script using petsc4py. Keep the same problem structure and solver setup." 
Personally, I really stuggled with using an LLM to properly translate the code. There were many times where the LLM would "hallucinate" examples and code specific examples of a linear problem while ignoring the source code, or try to implement more features that don't exist in the original C code. 

When using ChatGPT with the same prompt, none of the code worked, and half of the code produced was structured in C without translation into Python because of the request to "Keep the same problem structure and do not implement new features". Additionally, when asking it to check if the Claude code worked with the prompt "Is this a good translation of the file for just the quadratic 2dlastic problem or did claude hallucinate ?", the LLM shut down and was no longer able to provide full python translations of the original file code because of it hyperfixating on how and where "Claude" keeps going wrong. After erasing memory and previous chats, the problem prevailed with no useful code provided other than examples of how claude misinterpreted the original file. 

In the end, as ChatGPT had no longer been helpful to use, so I used Claude Haiku. Only by translating a few lines of code at a time, I recieved (in bits and pieces) a working python file in which I then was able to code in a way to output a solution and use it in a Jupyter Notebook. 

Prompt engineering was a difficult and frustrating experience that I have not dealt with before as a physicist, and as insightful this project was, I hope to never again need to rely upon. 

The following below will now have the rest of the information needed to use the elastic2D_SOLVER.


## Inputs, Material Properties

Implementations allow for inputs to material properties and simulation specs:
- **Lamé parameter λ**
- **Shear modulus μ**
- **Body Foce f**
- **Mesh size, nx by ny**

## Outputs, Graphs:
Based off the above inputs, the outputs will show:

- **A contoured colormap of the x displacement, y displamenet and total magnitude**
- **A visual of the mesh before and after force is enacted upon the material**


## Usage Documentation

The program supports the following format after importing the eslastic2D_SOLVER (values can be changed of course):

```python 
solver = ElasticitySolver(
    nx= 32,
    ny=32,
    mu=1.0,
    lam=1.0,
    body_force=[0.0, -9.8]
)

u = solver.solve()
solver.plot_solution(u)
solver.plot_deformed(u, scale=10.0)
```


## Implementation Notes

### Differences from C Version

1. **Python NumPy Arrays**: Solution/RHS evaluated using NumPy instead of raw C arrays
   - More Pythonic and readable
   - Easier numerical debugging


2. **Simplified Assembly**:
   - Focuses on problem setup and structure
   - Full element assembly would require DMPlex local assembly loops
   - Can be extended with `DMPlexComputeResidualAndJacobian` callbacks


## Building and Running

### Requirements
```bash
pip install petsc4py numpy
```

## References

- PETSc Documentation: https://petsc.org
- petsc4py: https://petsc4py.readthedocs.io
- Linear Elasticity Theory: https://en.wikipedia.org/wiki/Linear_elasticity
- DMPLEX mesh framework: https://petsc.org/release/docs/manual/dmplex/


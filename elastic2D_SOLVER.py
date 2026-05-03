"""
Linear Elasticity in 2D using Finite Elements (petsc4py)
=========================================================
Reusable module version of PETSc ex17.c translation
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

import petsc4py
from petsc4py import PETSc


# ===========================================================================
# Helper
# ===========================================================================

def compute_stress(grad_u, mu=1.0, lam=1.0):
    epsilon   = 0.5 * (grad_u + grad_u.T)
    trace_eps = np.trace(epsilon)
    sigma     = lam * trace_eps * np.eye(2) + 2.0 * mu * epsilon
    return sigma


# ===========================================================================
# Solver
# ===========================================================================

class ElasticitySolver:

    def __init__(self, nx=32, ny=32, mu=1.0, lam=1.0, body_force=None):
        self.nx  = nx
        self.ny  = ny
        self.mu  = mu
        self.lam = lam

        if body_force is None:
            body_force = [0.0, -1.0]

        self.body_force = np.array(body_force, dtype=float)

        self._build_mesh()

    # ---------------- Mesh ----------------

    def _build_mesh(self):
        nx, ny = self.nx, self.ny

        xs = np.linspace(0.0, 1.0, nx + 1)
        ys = np.linspace(0.0, 1.0, ny + 1)
        XX, YY = np.meshgrid(xs, ys)

        self.nodes = np.column_stack([XX.ravel(), YY.ravel()])
        self.n_nodes = self.nodes.shape[0]

        tris = []
        for j in range(ny):
            for i in range(nx):
                n00 = j*(nx+1) + i
                n10 = n00 + 1
                n01 = n00 + (nx+1)
                n11 = n01 + 1
                tris.append([n00, n10, n11])
                tris.append([n00, n11, n01])
        self.triangles = np.array(tris)

        tol = 1e-12
        c = self.nodes
        self.boundary_mask = (
            (np.abs(c[:, 0]) < tol) |
            (np.abs(c[:, 0] - 1.0) < tol) |
            (np.abs(c[:, 1]) < tol) |
            (np.abs(c[:, 1] - 1.0) < tol)
        )

    # ---------------- Assembly ----------------

    def _assemble(self):
        N = self.n_nodes
        lam, mu = self.lam, self.mu
        f_body = self.body_force

        K_values = {}
        f_vec = np.zeros(2 * N)

        for tri in self.triangles:
            x = self.nodes[tri]

            T = np.column_stack([np.ones(3), x])
            try:
                Tinv = np.linalg.inv(T)
            except np.linalg.LinAlgError:
                continue

            grads = Tinv[1:, :].T
            area  = 0.5 * abs(np.linalg.det(T))

            for a in range(3):
                for b in range(3):
                    ga, gb = grads[a], grads[b]
                    Kab = np.zeros((2, 2))

                    for i in range(2):
                        for j in range(2):
                            Kab[i, j] += mu * (ga[j]*gb[i] + ga[i]*gb[j])
                        for j in range(2):
                            Kab[i, j] += lam * ga[i] * gb[j]

                    Kab *= area

                    for i in range(2):
                        for j in range(2):
                            r = 2*tri[a] + i
                            c = 2*tri[b] + j
                            K_values[(r, c)] = K_values.get((r, c), 0.0) + Kab[i, j]

            f_local = f_body * area
            for a in range(3):
                for i in range(2):
                    f_vec[2*tri[a] + i] += f_local[i] / 3.0

        return K_values, f_vec

    # ---------------- PETSc ----------------

    def _build_petsc_system(self, K_values, f_vec, ndof):
        A = PETSc.Mat().create()
        A.setSizes([ndof, ndof])
        A.setType('aij')
        A.setPreallocationNNZ(50)
        A.setUp()

        for (r, c), val in K_values.items():
            A.setValue(r, c, val, addv=True)
        A.assemblyBegin()
        A.assemblyEnd()

        b = PETSc.Vec().create()
        b.setSizes(ndof)
        b.setUp()

        for i, val in enumerate(f_vec):
            b.setValue(i, val)

        b.assemblyBegin()
        b.assemblyEnd()

        return A, b

    def _apply_dirichlet(self, A, b):
        bc_dofs = []
        for node_id in np.where(self.boundary_mask)[0]:
            bc_dofs += [2*node_id, 2*node_id + 1]

        A.zeroRows(bc_dofs, diag=1.0)

        for dof in bc_dofs:
            b.setValue(dof, 0.0)

        b.assemblyBegin()
        b.assemblyEnd()

    # ---------------- Solve ----------------

    def solve(self):
        ndof = 2 * self.n_nodes

        print(f"\nElasticity Solver | mesh {self.nx}x{self.ny} | DOFs={ndof}")

        K_values, f_vec = self._assemble()
        A, b = self._build_petsc_system(K_values, f_vec, ndof)
        self._apply_dirichlet(A, b)

        ksp = PETSc.KSP().create()
        ksp.setOperators(A)
        ksp.setType('gmres')
        ksp.getPC().setType('ilu')
        ksp.setTolerances(rtol=1e-10)

        x = A.createVecRight()
        ksp.solve(b, x)

        u = x.getArray().copy().reshape(self.n_nodes, 2)

        A.destroy(); b.destroy(); x.destroy(); ksp.destroy()

        return u

    # ---------------- Plot ----------------

    def plot_solution(self, u, save_path=None):
        triang = mtri.Triangulation(
            self.nodes[:, 0], self.nodes[:, 1], self.triangles)

        mag = np.linalg.norm(u, axis=1)

        fig, axes = plt.subplots(1, 3, figsize=(16, 5))

        for ax, data, title, cmap in zip(
            axes,
            [u[:, 0], u[:, 1], mag],
            ['x-displacement $u_x$',
             'y-displacement $u_y$',
             'Magnitude $|\\mathbf{u}|$'],
            ['RdBu_r', 'RdBu_r', 'viridis']
        ):
            tc = ax.tricontourf(triang, data, levels=20, cmap=cmap)
            plt.colorbar(tc, ax=ax)
            ax.set_title(title, fontsize=12)
            ax.set_xlabel('x'); ax.set_ylabel('y')
            ax.set_aspect('equal')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"Saved to {save_path}")

        plt.show()

    def plot_deformed(self, u, scale=10.0, save_path=None):
        deformed = self.nodes + scale * u

        triang_orig = mtri.Triangulation(
            self.nodes[:, 0], self.nodes[:, 1], self.triangles)
        triang_def  = mtri.Triangulation(
            deformed[:, 0], deformed[:, 1], self.triangles)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        axes[0].triplot(triang_orig, color='steelblue', lw=0.5)
        axes[0].set_title('Original mesh', fontsize=12)
        axes[0].set_xlabel('x'); axes[0].set_ylabel('y')
        axes[0].set_aspect('equal')

        axes[1].triplot(triang_def, color='tomato', lw=0.5)
        axes[1].set_title(f'Deformed mesh (scale ×{scale})', fontsize=12)
        axes[1].set_xlabel('x'); axes[1].set_ylabel('y')
        axes[1].set_aspect('equal')

        for ax in axes:
            ax.set_aspect('equal')

        if save_path:
            plt.savefig(save_path, dpi=150)
            print(f"Saved to {save_path}")

        plt.show()


# ===========================================================================
# Optional test run
# ===========================================================================

if __name__ == "__main__":
    solver = ElasticitySolver(nx=64, ny=64, body_force=[0, -9.8])
    u = solver.solve()
    solver.plot_solution(u)
    solver.plot_deformed(u)
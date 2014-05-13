package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.boussinesq.ComputationalArrays;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

public class Solver {

	long startCompute, endCompute, startSolver, endSolver;

	/**
	 * Newton iteration.
	 * 
	 * @desc this method compute the Newton iteration.
	 * 
	 * @param arrb
	 *            the array of known terms
	 * @param arrT
	 *            the array of T for non Dirichlet cells
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the array of eta at the previous time step
	 * @param cg
	 *            the conjugate gradient class
	 * 
	 * @return the array of eta at the following time step
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public double[] newtonIteration(double[] arrb, double[] arrT, double[] eta,
			RCConjugateGradient cg, double tolerance)
			throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		DoubleMatrix1D matrixr;

		double maxResidual = 10;

		ComputeJr cJr = new ComputeJr();
		ComputeR cR = new ComputeR();

		ComputeBEqNoDirichlet.timeCompute = 0;
		ComputeBEqNoDirichlet.timeSolver = 0;

		do {

			startCompute = System.nanoTime();

			// compute Jr
			double[] jr = cJr.computeJr(arrT, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixJr = new SparseRCDoubleMatrix2D(ComputationalArrays.numberOfPolygon,
					ComputationalArrays.numberOfPolygon, ComputationalArrays.arrayMp,
					ComputationalArrays.arrayMj, jr);

			// compute the residual function
			double[] r = cR.computeR(arrT, arrb, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixr = new SparseDoubleMatrix1D(r);

			endCompute = System.nanoTime();

			startSolver = System.nanoTime();

			cg.solverCG(matrixr, matrixJr);

			endSolver = System.nanoTime();

			// compute the new eta for every cell
			for (int i = 0; i < ComputationalArrays.numberOfPolygon; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
			}

			// compute the max residual
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));

			System.out.println(maxResidual);
//			System.out.println(tolerance * 1000);

			ComputeBEqNoDirichlet.timeCompute = ComputeBEqNoDirichlet.timeCompute
					+ (endCompute - startCompute);
			ComputeBEqNoDirichlet.timeSolver = ComputeBEqNoDirichlet.timeSolver
					+ (endSolver - startSolver);

		} while (maxResidual > tolerance * 1000);

		return eta;

	}

}

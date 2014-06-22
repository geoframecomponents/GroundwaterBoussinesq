package org.boussinesq.rowCompressedForm;

import org.interfacesPDE.nonLinearParabolicPDE.solver.RowCompressedConjugateGradient;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleDiagonal;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleICC;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleILU;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleILUT;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoublePreconditioner;
import cern.colt.matrix.tdouble.algo.solver.preconditioner.DoubleSSOR;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

/**
 * The Class RCConjugateGradient.
 */
public class RCConjugateGradient implements RowCompressedConjugateGradient {

	/** The matrix_x. */
	DenseDoubleMatrix1D denseArray_x;

	/** The mat sol. */
	public DoubleMatrix1D arraySolution;

	protected DoublePreconditioner preconditioner;
	
	public void cgDoubleDiagonal_preconditioner(int equationNumber,
			SparseRCDoubleMatrix2D A) {
		preconditioner = new DoubleDiagonal(equationNumber);

		preconditioner.setMatrix(A);
	}

	public void cgICC_preconditioner(int equationNumber,
			SparseRCDoubleMatrix2D A) {
		
		preconditioner = new DoubleICC(equationNumber);

		preconditioner.setMatrix(A);

	}

	public void cgILU_preconditioner(int equationNumber,
			SparseRCDoubleMatrix2D A) {
		preconditioner = new DoubleILU(equationNumber);

		preconditioner.setMatrix(A);
	}

	public void cgILUT_preconditioner(int equationNumber,
			SparseRCDoubleMatrix2D A) {
		preconditioner = new DoubleILUT(equationNumber);

		preconditioner.setMatrix(A);
	}

	public void cgSSOR_preconditioner(int equationNumber,
			SparseRCDoubleMatrix2D A) {
		preconditioner = new DoubleSSOR(equationNumber);

		preconditioner.setMatrix(A);
	}
	
	/**
	 * Solver cg.
	 * 
	 * @param matrix_b
	 *            the matrix_b
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public void cgSolver(DenseDoubleMatrix1D matrix_b,
			SparseRCDoubleMatrix2D matrix_A)
			throws IterativeSolverDoubleNotConvergedException {

		int equationNumber = (int) matrix_b.size();
		
		denseArray_x = new DenseDoubleMatrix1D(equationNumber);
		arraySolution = new DenseDoubleMatrix1D(equationNumber);

		cgICC_preconditioner(equationNumber, matrix_A);
		
		DoubleCG conjugateGradient = new DoubleCG(denseArray_x);
		conjugateGradient.setPreconditioner(preconditioner);
		arraySolution = conjugateGradient.solve(matrix_A, matrix_b, denseArray_x);

	}

	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 * @throws IterativeSolverDoubleNotConvergedException
	 */
	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException {

		double[] b = { 1, 2 };
		int[] Mp = { 0, 2, 4 };
		int[] Mi = { 0, 1, 0, 1 };
		double[] Ml = { 4, 1, 1, 3 };

		DenseDoubleMatrix1D matrix_b = new DenseDoubleMatrix1D(b);
		SparseRCDoubleMatrix2D matrix_A = new SparseRCDoubleMatrix2D(b.length,
				b.length, Mp, Mi, Ml);
		// Mesh grid1 = new Mesh("Song");
		RCConjugateGradient cg = new RCConjugateGradient();

		cg.cgSolver(matrix_b, matrix_A);

		System.out.println(cg.arraySolution.get(0));

	}

}

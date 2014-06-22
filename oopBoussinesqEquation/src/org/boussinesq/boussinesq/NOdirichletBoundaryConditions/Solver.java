package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.computationalDomain.ActualComputationalDomain;
import org.boussinesq.rowCompressedForm.RCConjugateGradient;
import org.interfacesPDE.nonLinearParabolicPDE.LinearSystemAssembler;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

public class Solver implements LinearSystemAssembler {

	long startCompute, endCompute, startSolver, endSolver;

	private double[] array_Jr;
	private double[] array_r;

	private SparseRCDoubleMatrix2D sparseMatrix_Jr;
	private DenseDoubleMatrix1D denseArray_r;
//	private DenseDoubleMatrix1D denseArraySolution;
	private RCConjugateGradient conjugateGradient;

	public void computeCoefficientsMatrix(double[] array_T, double[] eta) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		array_Jr = new double[array_T.length];

		System.arraycopy(array_T, 0, array_Jr, 0, array_T.length);

		int endForLoop = ActualComputationalDomain.indexDiagonal.length;

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < endForLoop; i++) {

			// equation (A6)
			array_Jr[ActualComputationalDomain.indexDiagonal[i]] = array_T[ActualComputationalDomain.indexDiagonal[i]]
					+ PolygonGeometricalWetProperties.computeWetArea(eta[i],
							ActualComputationalDomain.bedRockElevation[i],
							ActualComputationalDomain.porosity[i],
							ActualComputationalDomain.planarArea[i]);

		}

	}

	public void computeSolutionArray(double[] array_T, double[] array_b, double[] eta) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		array_r = new double[ActualComputationalDomain.numberOfPolygons];

		for (int i = 0; i < ActualComputationalDomain.numberOfPolygons; i++) {

			for (int j = ActualComputationalDomain.array_Mp[i]; j < ActualComputationalDomain.array_Mp[i + 1]; j++) {
				sum += array_T[j] * eta[ActualComputationalDomain.array_Mj[j]];
			}

			double waterVolume = PolygonGeometricalWetProperties
					.computeWaterVolume(eta[i],
							ActualComputationalDomain.bedRockElevation[i],
							ActualComputationalDomain.porosity[i],
							ActualComputationalDomain.planarArea[i]);
			// equation (A3)
			array_r[i] = waterVolume + sum - array_b[i];

			sum = 0;

		}

	}

	/**
	 * Newton iteration.
	 * 
	 * @desc this method compute the Newton iteration.
	 * 
	 * @param array_b
	 *            the array of known terms
	 * @param array_T
	 *            the array of T for non Dirichlet cells
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the array of eta at the previous time step
	 * @param conjugateGradient
	 *            the conjugate gradient class
	 * 
	 * @return the array of eta at the following time step
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public double[] convergenceMethod(double[] array_b, double[] array_T,
			double[] eta, double tolerance)
			throws IterativeSolverDoubleNotConvergedException {

		double maxResidual = 10;

		ComputeBEq.timeCompute = 0;
		ComputeBEq.timeSolver = 0;

		do {

			startCompute = System.nanoTime();

			// compute Jr
			computeCoefficientsMatrix(array_T, eta);

			// convert array in sparse matrix for DoubleCG class
			sparseMatrix_Jr = new SparseRCDoubleMatrix2D(
					ActualComputationalDomain.numberOfPolygons,
					ActualComputationalDomain.numberOfPolygons,
					ActualComputationalDomain.array_Mp,
					ActualComputationalDomain.array_Mj, array_Jr);

			// compute the residual function
			computeSolutionArray(array_T, array_b, eta);

			// convert array in sparse matrix for DoubleCG class
			denseArray_r = new DenseDoubleMatrix1D(array_r);

			endCompute = System.nanoTime();

			startSolver = System.nanoTime();

			conjugateGradient.cgSolver(denseArray_r, sparseMatrix_Jr);

			endSolver = System.nanoTime();

			// compute the new eta for every cell
			for (int i = 0; i < ActualComputationalDomain.numberOfPolygons; i++) {
				eta[i] = eta[i] - conjugateGradient.arraySolution.get(i);
			}

			// compute the max residual
			maxResidual = Math.max(
					Math.abs(conjugateGradient.arraySolution.getMaxLocation()[0]),
					Math.abs(conjugateGradient.arraySolution.getMinLocation()[0]));

			ComputeBEq.timeCompute = ComputeBEq.timeCompute
					+ (endCompute - startCompute);
			ComputeBEq.timeSolver = ComputeBEq.timeSolver
					+ (endSolver - startSolver);

		} while (maxResidual > tolerance * 1000);

		return eta;

	}

}

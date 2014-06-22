package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.computationalDomain.AbstractDomain;
import org.boussinesq.boussinesq.computationalDomain.ActualComputationalDomain;
import org.boussinesq.rowCompressedForm.RCConjugateGradient;
import org.interfacesPDE.nonLinearParabolicPDE.LinearSystemAssembler;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

public class Solver extends IsNoValue implements LinearSystemAssembler {

	long startCompute, endCompute, startSolver, endSolver;

	double[] array_Jr;
	double[] array_r;

	SparseRCDoubleMatrix2D sparseMatrix_Jr;
	DenseDoubleMatrix1D denseArray_r;
	RCConjugateGradient conjugateGradient;

	public void computeCoefficientsMatrix(double[] arrT, double[] eta) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		array_Jr = new double[arrT.length];

		System.arraycopy(arrT, 0, array_Jr, 0, arrT.length);

		int endForLoop = ActualComputationalDomain.indexDiagonal.length;

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < endForLoop; i++) {

			if (isNoValue(ActualComputationalDomain.etaDirichlet[i],
					AbstractDomain.NOVALUE)) {
				// non Dirichlet cells
				// equation (A6)
				array_Jr[ActualComputationalDomain.indexDiagonal[i]] = arrT[ActualComputationalDomain.indexDiagonal[i]]
						+ PolygonGeometricalWetProperties.computeWetArea(
								eta[i],
								ActualComputationalDomain.bedRockElevation[i],
								ActualComputationalDomain.porosity[i],
								ActualComputationalDomain.planarArea[i]);

			} else {
				// Dirichlet cells
				array_Jr[ActualComputationalDomain.indexDiagonal[i]] = arrT[ActualComputationalDomain.indexDiagonal[i]];
			}
		}

	}

	public void computeSolutionArray(double[] arrT, double[] arrb, double[] eta) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		array_r = new double[ActualComputationalDomain.numberOfPolygons];

		for (int i = 0; i < ActualComputationalDomain.numberOfPolygons; i++) {
			if (isNoValue(ActualComputationalDomain.etaDirichlet[i],
					AbstractDomain.NOVALUE)) {

				// non Dirichlet cells
				for (int j = ActualComputationalDomain.array_Mp[i]; j < ActualComputationalDomain.array_Mp[i + 1]; j++) {
					sum += arrT[j] * eta[ActualComputationalDomain.array_Mj[j]];
				}

				double waterVolume = PolygonGeometricalWetProperties
						.computeWaterVolume(eta[i],
								ActualComputationalDomain.bedRockElevation[i],
								ActualComputationalDomain.porosity[i],
								ActualComputationalDomain.planarArea[i]);
				// equation (A3)
				array_r[i] = waterVolume + sum - arrb[i];

				sum = 0;
			} else {

				// Dirichlet cells
				array_r[i] = 0;
			}
		}

	}

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
	 * @param conjugateGradient
	 *            the conjugate gradient class
	 * 
	 * @return the array of eta at the following time step
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public double[] convergenceMethod(double[] arrb, double[] arrT,
			double[] eta, double tolerance)
			throws IterativeSolverDoubleNotConvergedException {

		double maxResidual = 10;

		ComputeBEqDirichlet.timeCompute = 0;
		ComputeBEqDirichlet.timeSolver = 0;

		conjugateGradient = new RCConjugateGradient();
		
		do {

			startCompute = System.nanoTime();

			// compute Jr
			computeCoefficientsMatrix(arrT, eta);

			// convert array in sparse matrix for DoubleCG class
			sparseMatrix_Jr = new SparseRCDoubleMatrix2D(
					ActualComputationalDomain.numberOfPolygons,
					ActualComputationalDomain.numberOfPolygons,
					ActualComputationalDomain.array_Mp,
					ActualComputationalDomain.array_Mj, array_Jr);

			// compute the residual function
			computeSolutionArray(arrT, arrb, eta);

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
			maxResidual = Math.max(Math.abs(conjugateGradient.arraySolution.getMaxLocation()[0]),
					Math.abs(conjugateGradient.arraySolution.getMinLocation()[0]));

			ComputeBEqDirichlet.timeCompute = ComputeBEqDirichlet.timeCompute
					+ (endCompute - startCompute);
			ComputeBEqDirichlet.timeSolver = ComputeBEqDirichlet.timeSolver
					+ (endSolver - startSolver);

		} while (maxResidual > tolerance * 100);

		return eta;

	}

}

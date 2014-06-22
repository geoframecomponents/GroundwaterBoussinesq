package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.AbstractComputeBEq;
//import org.boussinesq.boussinesq.AbstractComputeT;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.computationalDomain.AbstractDomain;
import org.boussinesq.boussinesq.computationalDomain.ActualComputationalDomain;
import org.boussinesq.rowCompressedForm.DeleteRowColumnNullDiagonalEntry;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEq extends AbstractComputeBEq {

	private double[] completeArray_T;

	private Solver newtonIteration;

	private DeleteRowColumnNullDiagonalEntry deleteSuperDryCells;
	private ActualComputationalDomain computeActualDomainArrays;
	
	public ComputeBEq() {

		eta = new double[AbstractDomain.Np];
		newtonIteration = new Solver();

		computeActualDomainArrays = new ActualComputationalDomain();
		deleteSuperDryCells = new DeleteRowColumnNullDiagonalEntry();
		
	}

	public double[] computeMatrixTerm(double[] eta) {

		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */
		double rowSum = 0;

		/* to identify the diagonal entry of matrix T in row-compressed form */
		int index = 0;

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[AbstractDomain.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < AbstractDomain.Np; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			for (int j = AbstractDomain.Mp[i]; j < AbstractDomain.Mp[i + 1]; j++) {

				if (AbstractDomain.Mi[j] != i) {
					// equation (21)
					arrayT[j] = -BoussinesqEquation.TIMESTEP
							* (1 / AbstractDomain.euclideanDistance[(int) AbstractDomain.Ml[j] - 1])
							* AbstractDomain.hydrConductivity[(int) AbstractDomain.Ml[j] - 1]
							* AbstractDomain.lengthSides[(int) AbstractDomain.Ml[j] - 1]
							* Math.max(
									Math.max(
											0,
											eta[AbstractDomain.Mi[j]]
													- AbstractDomain.bedRockElevation[AbstractDomain.Mi[j]]),
									Math.max(
											0,
											eta[i]
													- AbstractDomain.bedRockElevation[i]));

					rowSum += -arrayT[j];

				} else {
					index = j;
				}

			}
			// equation (20)
			arrayT[index] = rowSum;

			if (BoussinesqEquation.boundaryConditions.equals("NoDirichlet")
					&& rowSum == 0 && !unlockDeleteRowColumn)
				unlockDeleteRowColumn = true;

			rowSum = 0;
		}

		return arrayT;

	}

	public double[] computeArrayTerm(double[] eta) {

		double[] arrB = new double[ActualComputationalDomain.numberOfPolygons];

		for (int i = 0; i < ActualComputationalDomain.numberOfPolygons; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], ActualComputationalDomain.bedRockElevation[i],
					ActualComputationalDomain.porosity[i],
					ActualComputationalDomain.planarArea[i]);

			// delta t deve essere minore di 1/c
			arrB[i] = volume
					+ BoussinesqEquation.TIMESTEP
					* ActualComputationalDomain.planarArea[i]
					* ActualComputationalDomain.source[i]
					- BoussinesqEquation.TIMESTEP
					* ActualComputationalDomain.planarArea[i]
					* ActualComputationalDomain.ratingCurveCoeff_c[i]
					* Math.pow(volume / ActualComputationalDomain.planarArea[i],
							ActualComputationalDomain.ratingCurveCoeff_m[i]);

			if (arrB[i] < 0) {

				TextIO.putln("WARNING!!!\nThe element " + i
						+ " of the array of known terms is NEGATIVE");

			}

		}

		return arrB;

	}

	public void assemblePDE(double[] eta) {

		if (AbstractComputeBEq.unlockDeleteRowColumn) {

			completeArray_T = computeMatrixTerm(eta);
			completeArray_b = computeArrayTerm(eta);

			deleteSuperDryCells.computationalDomain(AbstractDomain.Np, completeArray_T,
					indexDiag, AbstractDomain.Mp, AbstractDomain.Mi);

			actualArray_T = deleteSuperDryCells.computeNewArrayT(completeArray_T);

			actualArray_b = deleteSuperDryCells.reduceToActualArray(completeArray_b);

		} else {

			actualArray_T = computeMatrixTerm(eta);
			actualArray_b = computeArrayTerm(eta);

			computeActualDomainArrays.completeDomain(eta, indexDiag);

		}

	}

	public double[] callSolution(double[] etaOld, double[] matT,
			double[] arrb) throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newtonIteration.convergenceMethod(arrb, matT, etaOld, tolerance);

		double[] etaNew = new double[AbstractDomain.eta.length];

		if (AbstractComputeBEq.unlockDeleteRowColumn) {

			AbstractComputeBEq.unlockDeleteRowColumn = false;

			etaNew = deleteSuperDryCells.addRemovedCells(eta,
					AbstractDomain.bedRockElevation);

		} else {

			System.arraycopy(eta, 0, etaNew, 0, eta.length);

		}

		return etaNew;

	}

	public void temporalLoop() throws IOException,
			IterativeSolverDoubleNotConvergedException {
		// allocate the memory for eta array

		firstThings();

		computeInitialVolume();

		// initialize eta array
		System.arraycopy(AbstractDomain.eta, 0, eta, 0,
				AbstractDomain.eta.length);

		for (double t = 0; t < BoussinesqEquation.SIMULATIONTIME; t += BoussinesqEquation.TIMESTEP) {

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			assemblePDE(eta);

			eta = callSolution(AbstractDomain.eta, actualArray_T, actualArray_b);

			computeOutputVariables(eta);

			writeSolution(t, eta, BoussinesqEquation.boundaryConditions,
					BoussinesqEquation.simulationType);

			computeVolumeConservation();

		}

	}

}

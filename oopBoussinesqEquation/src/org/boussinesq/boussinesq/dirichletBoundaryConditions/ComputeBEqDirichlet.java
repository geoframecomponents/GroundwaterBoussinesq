package org.boussinesq.boussinesq.dirichletBoundaryConditions;

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

public class ComputeBEqDirichlet extends AbstractComputeBEq {

	private double[] completeArray_T;
	private double[] completeArray_Tdirichlet;
	private double[] actualArray_TNoDirichlet;
	private double[] completeArray_TNoDirichlet1;

	private Solver newtonIteration;
	private ComputeTDirichlet computeTDirichlet;
	private ComputeTNoDirichlet computeTNoDirichlet;
	private IsNoValue verifyDirichletCells;

	private DeleteRowColumnNullDiagonalEntry deleteSuperDryCells;
	private ActualComputationalDomain computeActualArrays;
	
	public ComputeBEqDirichlet() {

		eta = new double[AbstractDomain.Np];
		newtonIteration = new Solver();

		verifyDirichletCells = new IsNoValue();

		computeTDirichlet = new ComputeTDirichlet();
		computeTNoDirichlet = new ComputeTNoDirichlet();

		computeActualArrays = new ActualComputationalDomain();
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

		double[] arrB = new double[AbstractDomain.Np];

		for (int i = 0; i < AbstractDomain.Np; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], AbstractDomain.bedRockElevation[i],
					AbstractDomain.porosity[i], AbstractDomain.planArea[i]);
			// equation (19)
			double sum = 0;
			for (int j = AbstractDomain.Mp[i]; j < AbstractDomain.Mp[i + 1]; j++) {
				sum += completeArray_Tdirichlet[j]
						* AbstractDomain.etaDirichlet[AbstractDomain.Mi[j]];
			}

			// delta t deve essere minore di 1/c
			arrB[i] = volume
					+ BoussinesqEquation.TIMESTEP
					* AbstractDomain.planArea[i]
					* AbstractDomain.source[i]
					- sum
					- BoussinesqEquation.TIMESTEP
					* AbstractDomain.planArea[i]
					* AbstractDomain.c[i]
					* Math.pow(volume / AbstractDomain.planArea[i],
							AbstractDomain.m[i]);

//			if (arrB[i] < 0) {
//
//				TextIO.putln("WARNING!!!\nThe element " + i
//						+ " of the array of known terms is NEGATIVE");
//
//			}

		}

		return arrB;

	}

	public void assemblePDE(double[] eta) {

		completeArray_T = computeMatrixTerm(eta);

		completeArray_Tdirichlet = computeTDirichlet.computeTDirichlet(completeArray_T, verifyDirichletCells);

		completeArray_TNoDirichlet1 = computeTNoDirichlet.computeTNoDirichlet(completeArray_T, indexDiag, verifyDirichletCells);

		if (AbstractComputeBEq.unlockDeleteRowColumn) {

			completeArray_b = computeArrayTerm(eta);

			deleteSuperDryCells.computationalDomain(AbstractDomain.Np,
					completeArray_TNoDirichlet1, indexDiag, AbstractDomain.Mp,
					AbstractDomain.Mi);

			actualArray_TNoDirichlet = deleteSuperDryCells
					.computeNewArrayT(completeArray_TNoDirichlet1);
			computeActualArrays.wetDomain(eta, indexDiag, deleteSuperDryCells);

			actualArray_b = deleteSuperDryCells.reduceToActualArray(completeArray_b);

		} else {

			actualArray_TNoDirichlet = computeTNoDirichlet
					.computeTNoDirichlet(completeArray_T, indexDiag, verifyDirichletCells);

			actualArray_b = computeArrayTerm(eta);

			computeActualArrays.completeDomain(eta, indexDiag);

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

		int endForLoop = AbstractDomain.etaDirichlet.length;

		for (int i = 0; i < endForLoop; i++) {

			if (!verifyDirichletCells.isNoValue(AbstractDomain.etaDirichlet[i],
					AbstractDomain.NOVALUE)) {

				etaNew[i] = AbstractDomain.etaDirichlet[i];

			}

			if (Math.abs(etaNew[i]) < tolerance)
				etaNew[i] = AbstractDomain.bedRockElevation[i];

		}

		return etaNew;

	}

	public void temporalLoop() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		firstThings();

		EtaInitialization etaInit = new EtaInitialization();

		computeInitialVolume();

		System.arraycopy(AbstractDomain.eta, 0, eta, 0,
				AbstractDomain.eta.length);

		for (double t = 0; t < BoussinesqEquation.SIMULATIONTIME; t += BoussinesqEquation.TIMESTEP) {

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			eta = etaInit.etaInitialization(eta, verifyDirichletCells);

			assemblePDE(eta);

			eta = callSolution(ActualComputationalDomain.eta,
					actualArray_TNoDirichlet, actualArray_b);

			computeOutputVariables(eta);

			writeSolution(t, eta, BoussinesqEquation.boundaryConditions,
					BoussinesqEquation.simulationType);

			computeVolumeConservation();

		}

	}

}

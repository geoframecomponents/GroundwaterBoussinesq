package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.RowCompressedForm.DeleteRowColumnNullDiagonalEntry;
import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeB;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq {

	double[] eta;
	double[] matT1, matT;
	double[] arrb;
	double[] arrb1;

	Solver newton;
	RCConjugateGradient cg;

	DeleteRowColumnNullDiagonalEntry deleteRowColumn;
	ComputationalArrays computationalArrays;

	ComputeT computeT;
	ComputeB cB;

	public ComputeBEqNoDirichlet() {

		eta = new double[ComputationalDomain.Np];
		newton = new Solver();
		cg = new RCConjugateGradient();
		computeT = new ComputeT();
		cB = new ComputeB();

		deleteRowColumn = new DeleteRowColumnNullDiagonalEntry();

	}

	public void computeBEqArrays(double[] eta) {

		computationalArrays = new ComputationalArrays();

		matT1 = computeT.computeT(eta, indexDiag);
		arrb1 = cB.computeB(eta);

		deleteRowColumn.computationalDomain(ComputationalDomain.Np, matT1,
				indexDiag, ComputationalDomain.Mp, ComputationalDomain.Mi);

		matT = computationalArrays.DefineDomainProperties(matT1, eta,
				indexDiag, deleteRowColumn);

		arrb = deleteRowColumn.computeCellsArray(arrb1);

	}

	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, etaOld, cg, tolerance);

		double[] etaNew = new double[ComputationalDomain.eta.length];

		if (ComputationalArrays.matrixReduction) {

			etaNew = deleteRowColumn.addRemovedElements(eta,
					ComputationalDomain.bedRockElevation);

		} else {

			System.arraycopy(eta, 0, etaNew, 0, eta.length);

		}

		return etaNew;

	}

	public void computeBEq(String boundaryCondition, String simulationType)
			throws IOException, IterativeSolverDoubleNotConvergedException {
		// allocate the memory for eta array

		firstThings();

		computeInitialVolume();

		// initialize eta array
		System.arraycopy(ComputationalDomain.eta, 0, eta, 0,
				ComputationalDomain.eta.length);

		for (int t = 0; t < BoussinesqEquation.SIMULATIONTIME; t += BoussinesqEquation.TIMESTEP) {

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			computeBEqArrays(eta);

			eta = solutionMethod(ComputationalArrays.arrayEta, matT, arrb);

			computeOutputFeatures(eta);

			writeSolution(t, eta, boundaryCondition, simulationType);

			ComputationalArrays.matrixReduction = false;

			computeVolumeConservation();

		}

	}

}

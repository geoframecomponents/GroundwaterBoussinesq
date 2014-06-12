package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq {

	double[] matT1;

	Solver newton;
	ComputeB cB;

	public ComputeBEqNoDirichlet() {

		eta = new double[ComputationalDomain.Np];
		newton = new Solver();
		cB = new ComputeB();

	}

	public void computeBEqArrays(double[] eta) {

		matT1 = computeT.computeT(eta, indexDiag);
		// arrb1 = cB.computeB(eta);
		//
		// // deleteRowColumn.computationalDomain(ComputationalDomain.Np, matT1,
		// // indexDiag, ComputationalDomain.Mp, ComputationalDomain.Mi);
		//
		// matT = computationalArrays.DefineDomainProperties(matT1, eta,
		// indexDiag, deleteRowColumn);
		//
		// arrb = deleteRowColumn.computeCellsArray(arrb1);

		if (ComputeT.unlockDeleteRowColumn) {

			// matT1 = computeT.computeT(eta, indexDiag);
			arrb1 = cB.computeB(eta);

			deleteRowColumn.computationalDomain(ComputationalDomain.Np, matT1,
					indexDiag, ComputationalDomain.Mp, ComputationalDomain.Mi);

			matT = deleteRowColumn.computeNewArrayT(matT1);

			arrb = deleteRowColumn.computeCellsArray(arrb1);

		} else {

			matT = computeT.computeT(eta, indexDiag);
			arrb = cB.computeB(eta);

			computationalArrays.completeDomain(eta, indexDiag);

		}

	}

	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, etaOld, cg, tolerance);

		double[] etaNew = new double[ComputationalDomain.eta.length];

		if (ComputeT.unlockDeleteRowColumn) {

			ComputeT.unlockDeleteRowColumn = false;

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

		for (double t = 0; t < BoussinesqEquation.SIMULATIONTIME; t += BoussinesqEquation.TIMESTEP) {

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			computeBEqArrays(eta);

			eta = solutionMethod(ComputationalArrays.arrayEta, matT, arrb);

			computeOutputFeatures(eta);

			writeSolution(t, eta, boundaryCondition, simulationType);

			computeVolumeConservation();

		}

	}

}

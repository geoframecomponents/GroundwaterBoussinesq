package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqDirichlet extends ComputeBEq {

	double[] matTDirichlet;
	double[] matTDirichlet1;
	double[] matTNoDirichlet;
	double[] matTNoDirichlet1;

	double volumeDirichlet;

	Solver newton;
	ComputeTDirichlet cTDirichlet;
	ComputeTNoDirichlet cTNoDirichlet;
	IsNoValue verifyDirichlet;

	ComputeB cB;

	public ComputeBEqDirichlet() {

		eta = new double[ComputationalDomain.Np];
		newton = new Solver();
		cB = new ComputeB();

		verifyDirichlet = new IsNoValue();

		cTDirichlet = new ComputeTDirichlet();
		cTNoDirichlet = new ComputeTNoDirichlet();

	}

	public void computeBEqArrays(double[] eta, int[] indexDiag) {

		matT = computeT.computeT(eta, indexDiag);

		matTDirichlet = cTDirichlet.computeTDirichlet(matT);

		matTNoDirichlet1 = cTNoDirichlet.computeTNoDirichlet(matT,
				indexDiag);
		
		if (ComputeT.unlockDeleteRowColumn) {

			arrb1 = cB.computeB(eta, matTDirichlet);

			deleteRowColumn.computationalDomain(ComputationalDomain.Np,
					matTNoDirichlet1, indexDiag, ComputationalDomain.Mp,
					ComputationalDomain.Mi);

			matTNoDirichlet = deleteRowColumn
					.computeNewArrayT(matTNoDirichlet1);
			computationalArrays.wetDomain(eta, indexDiag,
					deleteRowColumn);

			arrb = deleteRowColumn.computeCellsArray(arrb1);

		} else {

			matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT,
					indexDiag);
			
			arrb = cB.computeB(eta, matTDirichlet);

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

		int endForLoop = ComputationalDomain.etaDirichlet.length;

		for (int i = 0; i < endForLoop; i++) {

			if (!verifyDirichlet.isNoValue(ComputationalDomain.etaDirichlet[i],
					ComputationalDomain.NOVALUE)) {

				etaNew[i] = ComputationalDomain.etaDirichlet[i];

			}

			if (Math.abs(etaNew[i]) < tolerance)
				etaNew[i] = ComputationalDomain.bedRockElevation[i];

		}

		return etaNew;

	}

	public void computeBEq(String boundaryCondition, String simulationType)
			throws IOException, IterativeSolverDoubleNotConvergedException {

		firstThings();

		EtaInitialization etaInit = new EtaInitialization();

		computeInitialVolume();

		System.arraycopy(ComputationalDomain.eta, 0, eta, 0,
				ComputationalDomain.eta.length);

		for (double t = 0; t < BoussinesqEquation.SIMULATIONTIME; t += BoussinesqEquation.TIMESTEP) {

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(eta, indexDiag);

			eta = solutionMethod(ComputationalArrays.arrayEta, matTNoDirichlet,
					arrb);

			computeOutputFeatures(eta);

			writeSolution(t, eta, boundaryCondition, simulationType);

			computeVolumeConservation();

		}

	}
}

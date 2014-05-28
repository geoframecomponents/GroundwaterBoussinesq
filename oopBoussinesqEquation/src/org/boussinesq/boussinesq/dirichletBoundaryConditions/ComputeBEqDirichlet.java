package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.RowCompressedForm.DeleteRowColumnNullDiagonalEntry;
import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqDirichlet extends ComputeBEq implements TimeSimulation {

	double[] eta;
	double[] matT;
	double[] arrb1;
	double[] arrb;

	double[] matTDirichlet;
	double[] matTDirichlet1;
	double[] matTNoDirichlet;
	double[] matTNoDirichlet1;

	double volumeDirichlet;

	Solver newton;
	RCConjugateGradient cg;
	ComputeT computeT;
	ComputeTDirichlet cTDirichlet;
	ComputeTNoDirichlet cTNoDirichlet;
	DeleteRowColumnNullDiagonalEntry deleteRowColumn;
	ComputationalArrays computeNewArray;
	IsNoValue verifyDirichlet;

	ComputeB cB;

	int[] indexDiag;

	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;
	double tolerance;

	public ComputeBEqDirichlet() {

		eta = new double[ComputationalDomain.Np];

		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
		newton = new Solver();
		cg = new RCConjugateGradient();

		deleteRowColumn = new DeleteRowColumnNullDiagonalEntry();
		computeNewArray = new ComputationalArrays();
		verifyDirichlet = new IsNoValue();

		computeT = new ComputeT();
		cTDirichlet = new ComputeTDirichlet();
		cTNoDirichlet = new ComputeTNoDirichlet();
		cB = new ComputeB();

	}

	public void computeBEqArrays(double[] eta, int[] indexDiag) {

		matT = computeT.computeT(eta, indexDiag);

		matTNoDirichlet1 = cTNoDirichlet.computeTNoDirichlet(matT, indexDiag);
		matTDirichlet = cTDirichlet.computeTDirichlet(matT);

		arrb1 = cB.computeB(eta, matTDirichlet);

		deleteRowColumn.computationalDomain(ComputationalDomain.Np,
				matTNoDirichlet1, indexDiag, ComputationalDomain.Mp,
				ComputationalDomain.Mi);

		matTNoDirichlet = computeNewArray.DefineDomainProperties(
				matTNoDirichlet1, eta, indexDiag, deleteRowColumn);

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

		int endForLoop = ComputationalDomain.etaDirichlet.length;
		
		for (int i = 0; i < endForLoop; i++) {

			if (!verifyDirichlet.isNoValue(ComputationalDomain.etaDirichlet[i],
					ComputationalDomain.NOVALUE)) {

				etaNew[i] = ComputationalDomain.etaDirichlet[i];

			}
			
			if(Math.abs(etaNew[i]) < tolerance) etaNew[i] = ComputationalDomain.bedRockElevation[i];

		}

		return etaNew;

	}

	public void firstThings() {

		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
				ComputationalDomain.Np, ComputationalDomain.Mp,
				ComputationalDomain.Mi);

		tolerance = cMEd.computeMachineEpsilonDouble();

	}

	public void computeBEq(String boundaryCondition, String simulationType)
			throws IOException, IterativeSolverDoubleNotConvergedException {

		firstThings();

		EtaInitialization etaInit = new EtaInitialization();

		computeInitialVolume();

		System.arraycopy(ComputationalDomain.eta, 0, eta, 0,
				ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(eta, indexDiag);

			eta = solutionMethod(ComputationalArrays.arrayEta, matTNoDirichlet,
					arrb);

			computeOutputFeatures(eta);

			writeSolution(t, eta, boundaryCondition, simulationType);

			ComputationalArrays.matrixReduction = false;

			computeVolumeConservation();

		}

	}
}

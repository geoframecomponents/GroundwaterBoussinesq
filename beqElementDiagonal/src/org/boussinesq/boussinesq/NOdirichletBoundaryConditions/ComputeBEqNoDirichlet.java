package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeB;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq implements TimeSimulation {

	double[] eta;
	double[] matT;
	double[] arrb;

	Solver newton;
	RCConjugateGradient cg;

	// int[] indexDiag;

	// RCIndexDiagonalElement rcIndexDiagonalElement;
	// MachineEpsilon cMEd;
	// double tolerance;

	public ComputeBEqNoDirichlet() {

		eta = new double[ComputationalDomain.Np];
		// rcIndexDiagonalElement = new RCIndexDiagonalElement();
		// cMEd = new MachineEpsilon();
		newton = new Solver();
		cg = new RCConjugateGradient(ComputationalDomain.Np);

	}

	public void computeBEqArrays(double[] eta) {

		ComputeT computeT = new ComputeT();
		ComputeB cB = new ComputeB();
		matT = computeT.computeT(eta);
		arrb = cB.computeB(eta);

	}

	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, cg,
				tolerance);

		// int endForLoop = ComputationalDomain.Np;
		//
		// for (int i = 0; i < endForLoop; i++) {
		//
		//
		// eta[i] = Math.max(a, b)
		// if (Math.abs(eta[i]) < ComputationalDomain.bedRockElevation[i])
		// eta[i] = ComputationalDomain.bedRockElevation[i];
		//
		// }

		return eta;

	}

	public void computeBEq(String boundaryCondition, String simulationType)
			throws IOException, IterativeSolverDoubleNotConvergedException {
		// allocate the memory for eta array

		firstThings();

		computeInitialVolume();

		// initialize eta array
		System.arraycopy(ComputationalDomain.eta, 0, eta, 0,
				ComputationalDomain.eta.length);

		int contatore = 0;
		
//		int limit = 900/TIMESTEP;
		
//		int countLoop = limit;
		
		ComputationalDomain.source = new double[ComputationalDomain.Np];
		ComputationalDomain.outflow = new double[ComputationalDomain.Np];

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {

//			if (countLoop == limit) {
//
				for (int i = 0; i < ComputationalDomain.Np; i++) {

					ComputationalDomain.source[i] = ComputationalDomain.rainHour[contatore];

				}
				
				contatore++;
//				countLoop = 1;
//
//			}

			TextIO.putln("Time step " + (double) t / 3600);

			openTxtFile(t);

			computeBEqArrays(eta);

			eta = solutionMethod(eta, matT, arrb);

			computeOutputFeatures(eta);

			writeSolution(t, eta, boundaryCondition, simulationType);

			computeVolumeConservation();

//			countLoop++;

		}

	}

}

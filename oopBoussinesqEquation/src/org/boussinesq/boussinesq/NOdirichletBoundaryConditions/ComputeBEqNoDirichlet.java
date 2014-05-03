package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeB;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq implements TimeSimulation {

	double[] eta;
	double[] matT;
	double[] arrb;
	
	
	public ComputeBEqNoDirichlet() {
		
		eta = new double[ComputationalDomain.Np];


	}
	
	public void computeBEqArrays(double[] eta) {

		ComputeT computeT = new ComputeT();
		ComputeB cB = new ComputeB();
		matT = computeT.computeT(eta);
		arrb = cB.computeB(eta);

	}
	
	public void computeBEq(String boundaryCondition, String simulationType) throws IOException,
	IterativeSolverDoubleNotConvergedException {
		// 	allocate the memory for eta array
		
			
		firstThings();
		
		computeInitialVolume();
		
		// 	initialize eta array
		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
	
			TextIO.putln("Time step " + (double) t/3600);

			openTxtFile(t);
	
			computeBEqArrays(eta);

			eta = solutionMethod(eta, matT, arrb);

			computeOutputFeatures(eta);
	
			writeSolution(t, eta, boundaryCondition, simulationType);
	
			computeVolumeConservation();
	
		}


	}	

}

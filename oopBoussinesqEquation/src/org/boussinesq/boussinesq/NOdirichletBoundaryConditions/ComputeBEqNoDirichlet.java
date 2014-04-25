package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq implements TimeSimulation {

	double[] eta;
	
	int[] indexDiag;
	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;
	double tolerance;
	
	public ComputeBEqNoDirichlet() {
		
		eta = new double[ComputationalDomain.Np];
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
	}
	
	public void computeBEq() throws IOException,
	IterativeSolverDoubleNotConvergedException {
		// 	allocate the memory for eta array
		
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(ComputationalDomain.Np,
				ComputationalDomain.Mp, ComputationalDomain.Mi);

		tolerance = cMEd.computeMachineEpsilonDouble();
				
		computeInitialVolume();
		
		// 	initialize eta array
		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
	
			TextIO.putln("Time step " + (double) t/3600);

			openTxtFile(t);
	
			computeBEqArrays();

			eta = newton.newtonIteration(arrb, matT, indexDiag, eta, cg,
			tolerance);

			computeOutputFeatures();
	
			writeSolution(t);
	
			computeVolumeConservation();
	
		}


	}	

}

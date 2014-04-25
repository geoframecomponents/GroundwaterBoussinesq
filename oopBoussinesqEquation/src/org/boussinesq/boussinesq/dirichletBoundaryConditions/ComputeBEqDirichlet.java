package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.IOException;
import java.text.DecimalFormat;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
//import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.francescoS.usefulClasses.FileWrite;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqDirichlet extends ComputeBEq implements TimeSimulation {

	double[] eta;
	
	double[] matTDirichlet;
	double[] matTNoDirichlet;

	double volumeDirichlet;

	
	
	int[] indexDiag;
	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;
	double tolerance;
	
	public ComputeBEqDirichlet() {
		
		eta = new double[ComputationalDomain.Np];
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
	}
	

	
	
	
	
	public void computeBEqArrays(ComputeTDirichlet cTDirichlet,
			ComputeTNoDirichlet cTNoDirichlet, ComputeB cB) {

		matT = computeT(eta);
		matTDirichlet = cTDirichlet.computeTDirichlet(matT);
		matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT);
		arrb = cB.computeB(eta, matTDirichlet);

	}

	
	
	
	
	
	
	
	
	
	
	public void computeBEq() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		
		EtaInitialization etaInit = new EtaInitialization();

		ComputeTDirichlet cTDirichlet = new ComputeTDirichlet();
		ComputeTNoDirichlet cTNoDirichlet = new ComputeTNoDirichlet();
		
		computeInitialVolume();

		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
			
			TextIO.putln("Time step " + (double) t/3600);
			
			openTxtFile(t);

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(cTDirichlet, cTNoDirichlet, cB);

			eta = newton.newtonIteration(arrb, matTNoDirichlet, indexDiag, eta,
					cg, tolerance);

			computeOutputFeatures();
			
			writeSolution(t);
			
			computeVolumeConservation();
			
		}

	}
}

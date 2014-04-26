package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqDirichlet extends ComputeBEq implements TimeSimulation {

	double[] eta;
	double[] matT;
	double[] arrb;
	
	double[] matTDirichlet;
	double[] matTNoDirichlet;

	double volumeDirichlet;

	Solver newton;
	RCConjugateGradient cg;
	ComputeT computeT;
	ComputeTDirichlet cTDirichlet;
	ComputeTNoDirichlet cTNoDirichlet;
	
	
	int[] indexDiag;

	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;
	double tolerance;
	
	public ComputeBEqDirichlet() {
		
		eta = new double[ComputationalDomain.Np];
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
		newton = new Solver();
		cg = new RCConjugateGradient(ComputationalDomain.Np);
		
		computeT = new ComputeT();
		cTDirichlet = new ComputeTDirichlet();
		cTNoDirichlet = new ComputeTNoDirichlet();
	

	}
	

	
	
	
	
	public void computeBEqArrays(double[] eta) {

		
		matT = computeT.computeT(eta);
		matTDirichlet = cTDirichlet.computeTDirichlet(matT);
		matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT);
		ComputeB cB = new ComputeB();
		arrb = cB.computeB(eta, matTDirichlet);

	}
	
	
	
	
	
	
	
	
	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb) throws IterativeSolverDoubleNotConvergedException{
		
		double[] eta = new double[etaOld.length];
		
		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, cg,
				tolerance);
		
		return eta;
		
	}

	
	
	
	
	
	
	
	
	
	
	public void computeBEq() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		firstThings();
		
		EtaInitialization etaInit = new EtaInitialization();

		computeInitialVolume();

		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
			
			TextIO.putln("Time step " + (double) t/3600);
			
			openTxtFile(t);

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(eta);

			eta = solutionMethod(eta, matTNoDirichlet, arrb);
					
			computeOutputFeatures(eta);
			
			writeSolution(t, eta);
			
			computeVolumeConservation();
			
		}

	}
}

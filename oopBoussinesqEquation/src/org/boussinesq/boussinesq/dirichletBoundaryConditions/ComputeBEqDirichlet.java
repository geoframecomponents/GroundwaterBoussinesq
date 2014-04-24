package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
//import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDoman.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqDirichlet extends ComputeT implements TimeSimulation {

	// allocate the memory for eta array
	double[] eta;

	double[] matT;
	double[] matTDirichlet;
	double[] matTNoDirichlet;
	double[] arrb;
	double[] volume;
	
	double volumeOld = 0;
	double volumeNew = 0;

	
	
	
	
	
	
	
	
	
	public ComputeBEqDirichlet() {
		super();

	}
	
	
	
	
	
	
	
	
	
	
	/**
	 * Compute volume.
	 *
	 * @param index the index
	 * @param eta the eta
	 * @return the double
	 */
	public double computeVolume(int index,double eta){
		
		double volume;
		
		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta, ComputationalDomain.bedRockElevation[index], ComputationalDomain.porosity[index], ComputationalDomain.planArea[index]);
				
		volume = volume	- TIMESTEP * ComputationalDomain.planArea[index] * ComputationalDomain.source[index]
						+ TIMESTEP * ComputationalDomain.planArea[index] * ComputationalDomain.c[index] * Math.pow(volume/ComputationalDomain.planArea[index], ComputationalDomain.m[index]);
		
		return volume;
		
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
		ComputeB cB = new ComputeB();
		Solver newton = new Solver();

		RCConjugateGradient cg = new RCConjugateGradient(ComputationalDomain.Np);
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();

		int[] indexDiag = rcIndexDiagonalElement.computeIndexDiag(ComputationalDomain.Np,
				ComputationalDomain.Mp, ComputationalDomain.Mi);

		MachineEpsilon cMEd = new MachineEpsilon();
		double tolerance = cMEd.computeMachineEpsilonDouble();

		// allocate the memory for eta array
		eta = new double[ComputationalDomain.Np];
		volume = new double[ComputationalDomain.Np];
		
		for (int i = 0; i < ComputationalDomain.Np; i++){
			
			volume[i] = PolygonGeometricalWetProperties.computeWaterVolume(ComputationalDomain.eta[i], ComputationalDomain.bedRockElevation[i], ComputationalDomain.porosity[i], ComputationalDomain.planArea[i]);
			volumeOld = volumeOld + volume[i];
			
		}

		FileWriter Rstatfile = new FileWriter(ComputationalDomain.outputPathBeqDirichlet);
		PrintWriter errestat = new PrintWriter(Rstatfile);

		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(cTDirichlet, cTNoDirichlet, cB);

			eta = newton.newtonIteration(arrb, matTNoDirichlet, indexDiag, eta,
					cg, tolerance);

			System.out.println("Simulation time: " + (double) t / 3600);

			for (int j = 0; j < eta.length; j++) {
				
				volume[j] = computeVolume(j,eta[j]);
				volumeNew = volumeNew + volume[j];

			}
			
			TextIO.putln("Volume " + volumeNew + "at time step " + t/3600);
			
			volumeNew = 0;
			
		}
		for (int j = 0; j < eta.length; j++) {

			errestat.println(eta[j]);

		}

		errestat.println();
		System.out.println();
		Rstatfile.close();

		System.out.println("Exit code");
	}
}

package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.Mesh;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.machineEpsilon.MachineEpsilon;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEq extends ComputeT implements TimeSimulation {

	// allocate the memory for eta array
	double[] eta;

	double[] matT;
	double[] arrb;
	double[] volume;
	
	double volumeOld = 0;
	double volumeNew = 0;

	
	
	
	
	
	
	
	
	
	public ComputeBEq() {
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
		
		volume = (eta-Mesh.bedRockElevation[index])*Mesh.porosity[index]*Mesh.planArea[index];
		
		return volume;
		
	}

	
	
	
	
	
	
	
	
	
	public void computeBEqArrays(ComputeB cB) {

		matT = computeT(eta);
		arrb = cB.computeB(eta);

	}

	
	
	
	
	
	
	
	
	
	public void computeBEq() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		ComputeB cB = new ComputeB();
		Solver newton = new Solver();
		RCConjugateGradient cg = new RCConjugateGradient(Mesh.Np);
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();
		MachineEpsilon cMEd = new MachineEpsilon();		
		
		FileWriter Rstatfile = new FileWriter(Mesh.outputPathBeqNoDirichlet);
		PrintWriter errestat = new PrintWriter(Rstatfile);		
		
		

		int[] indexDiag = rcIndexDiagonalElement.computeIndexDiag(Mesh.Np,
				Mesh.Mp, Mesh.Mi);

		double tolerance = cMEd.computeMachineEpsilonDouble();


		// allocate the memory for eta array
		eta = new double[Mesh.Np];
		volume = new double[Mesh.Np];

		for (int i = 0; i < Mesh.Np; i++){
			
			volume[i] = computeVolume(i,Mesh.eta[i]);
			volumeOld = volumeOld + volume[i];
			
		}
		
		// initialize eta array
		System.arraycopy(Mesh.eta, 0, eta, 0, Mesh.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {

			computeBEqArrays(cB);

			eta = newton.newtonIteration(arrb, matT, indexDiag, eta, cg,
					tolerance);

			System.out.println("Simulation time: " + t / 3600);

		}
		for (int j = 0; j < eta.length; j++) {

			errestat.println(eta[j]);
			
			
			volume[j] = computeVolume(j,eta[j]);
			volumeNew = volumeNew + volume[j];

		}

		errestat.println();
		System.out.println();
		Rstatfile.close();

		System.out.println("Exit code");

	}

}

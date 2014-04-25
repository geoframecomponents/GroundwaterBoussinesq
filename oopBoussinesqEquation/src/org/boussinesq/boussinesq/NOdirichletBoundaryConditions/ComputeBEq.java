package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

//import java.io.FileWriter;
import java.io.IOException;
//import java.io.PrintWriter;

import java.text.DecimalFormat;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.ComputeBEqDirichlet;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.francescoS.usefulClasses.FileWrite;
import org.francescoS.usefulClasses.TextIO;




//import cern.colt.Arrays;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEq extends ComputeT implements TimeSimulation {

	// allocate the memory for eta array
	double[] eta;
	double[] aquiferThickness;
	double[] volumeSource;

	double[] matT;
	double[] arrb;
	double[] volume;
	
	double volumeOld = 0;
	double volumeNew = 0;
	
	static long timeCompute;
	static long timeSolver;

	
	
	
	
	
	
	
	
	
	public ComputeBEq() {
		super();
	}
	
	
	
	
	
	
	
	
	
	
	public DecimalFormat computePattern(){
		
		int[] c = new int[String.valueOf(SIMULATIONTIME).length()];
		
		StringBuilder builder = new StringBuilder(c.length);
		
		for (int i:c){
			
			builder.append(c[i]);
			
		}
		
		String pattern = builder.toString();
		
		return new DecimalFormat(pattern);
		
	}
	
	
	
	
	
	
	
	
	
	
	public void writeSolution(int time) throws IOException{
		
		FileWrite.writeStringIntString("Iteration number", (int) time/TIMESTEP + 1, " [ ]");
		FileWrite.writeStringDoubleString("Timestep", TIMESTEP, " [s]");
		FileWrite.writeStringDoubleString("Time of simulation", SIMULATIONTIME, " [s]");
		FileWrite.writeStringDoubleString("Time of simulation", SIMULATIONTIME/60, " [min]");
		FileWrite.writeStringDoubleString("Initial volume", volumeOld, "[ m^3]");
		FileWrite.writeStringDoubleString("Total volume", volumeNew, " [m^3]");
		FileWrite.writeStringDoubleString("Time convergence of CG", ComputeBEq.timeSolver/1000000000.0, "[s]");
		FileWrite.writeFourStringColumn("PiezHead","AquifThick","WaterVolume","Source");
		FileWrite.writeFourStringColumn("[m]", "[m]", "[m^3]", "[m^3/s]");
		FileWrite.writeFourDoubleColumn(eta,aquiferThickness,volume,volumeSource);
		FileWrite.closeTxtFile();
		
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

	
	
	
	
	
	
	
	
	
	public void computeBEqArrays(ComputeB cB) {

		matT = computeT(eta);
		arrb = cB.computeB(eta);

	}

	
	
	
	
	
	
	
	
	
	public void computeBEq() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		ComputeB cB = new ComputeB();
		Solver newton = new Solver();
		RCConjugateGradient cg = new RCConjugateGradient(ComputationalDomain.Np);
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();
		MachineEpsilon cMEd = new MachineEpsilon();	
		
		
		
		String fileName = null;
		DecimalFormat myformatter = computePattern();
		
		int[] indexDiag = rcIndexDiagonalElement.computeIndexDiag(ComputationalDomain.Np,
				ComputationalDomain.Mp, ComputationalDomain.Mi);

		double tolerance = cMEd.computeMachineEpsilonDouble();


		// allocate the memory for eta array
		aquiferThickness = new double[ComputationalDomain.Np];
		volumeSource = new double[ComputationalDomain.Np];
		eta = new double[ComputationalDomain.Np];
		volume = new double[ComputationalDomain.Np];

		for (int i = 0; i < ComputationalDomain.Np; i++){
			
			volume[i] = PolygonGeometricalWetProperties.computeWaterVolume(ComputationalDomain.eta[i], ComputationalDomain.bedRockElevation[i], ComputationalDomain.porosity[i], ComputationalDomain.planArea[i]);
			volumeOld = volumeOld + volume[i];
			
		}
		
		TextIO.putln("Initial volume: " + volumeOld);
		
		// initialize eta array
		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
			
			TextIO.putln("Time step " + (double) t/3600);

			fileName = myformatter.format(t);
			FileWrite.openTxtFile(fileName.concat(".txt"), BoussinesqEquation.solutionDir, false);
			
			computeBEqArrays(cB);

			eta = newton.newtonIteration(arrb, matT, indexDiag, eta, cg,
					tolerance);

			for (int j = 0; j < eta.length; j++) {
				
				volumeSource[j] = ComputationalDomain.source[j] * ComputationalDomain.planArea[j];
				aquiferThickness[j] = Math.max(0, eta[j]-ComputationalDomain.bedRockElevation[j]);
				volume[j] = computeVolume(j,eta[j]);
				volumeNew = volumeNew + volume[j];

			}
			
			writeSolution(t);
			
			if (Math.abs(volumeNew-volumeOld) > Math.pow(10, -6)){
				
				TextIO.putln("WARNING!!! The system is losing mass");
				TextIO.putln("The difference between initial volume and compute volume is: " + Math.abs(volumeNew-volumeOld));
				
			}
			
			
			volumeNew = 0;
			
		}

//		System.out.println("Exit code");

	}

}

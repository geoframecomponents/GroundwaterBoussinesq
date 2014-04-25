package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.IOException;
import java.text.DecimalFormat;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
//import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDoman.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.francescoS.usefulClasses.FileWrite;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.jet.random.tdouble.Zeta;

public class ComputeBEqDirichlet extends ComputeT implements TimeSimulation {

	// allocate the memory for eta array
	double[] eta;
	double[] aquiferThickness;
	double[] volumeSource;

	double[] matT;
	double[] matTDirichlet;
	double[] matTNoDirichlet;
	double[] arrb;
	double[] volume;
	
	double volumeOld = 0;
	double volumeNew = 0;
	double volumeDirichlet;
	
	static long timeCompute;
	static long timeSolver;

	
	
	
	
	
	
	
	
	
	public ComputeBEqDirichlet() {
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
		
		FileWrite.writeStringIntString("Iteration number", (int) time/TIMESTEP + 1, "[ ]");
		FileWrite.writeStringDoubleString("Timestep", TIMESTEP, "[s]");
		FileWrite.writeStringDoubleString("Time of simulation", SIMULATIONTIME, "[s]");
		FileWrite.writeStringDoubleString("Time of simulation", SIMULATIONTIME/60, "[min]");
		FileWrite.writeStringDoubleString("Initial volume", volumeOld, "[ m^3]");
		FileWrite.writeStringDoubleString("Total volume", volumeNew, "[m^3]");
		FileWrite.writeStringDoubleString("Time convergence of CG", ComputeBEqDirichlet.timeSolver/1000000000.0, "[s]");
		FileWrite.writeFourStringColumn("PiezHead","AquifThick","WaterVolume","Source");
		FileWrite.writeFourStringColumn("[m]", "[m]", "[m^3]", "[m^3/s]");
		FileWrite.writeFourDoubleColumn(eta,aquiferThickness,volume,volumeSource);
		FileWrite.closeTxtFile();
		
	}
	
	
	
	
	
	
	
	
	
	public double computeVolumeDirichlet(double volume){
		
		for (int index = 0; index < ComputationalDomain.Np; index++){
			
			if (ComputationalDomain.etaDirichlet[index] != ComputationalDomain.NOVALUE){
				
				volume = volume + TIMESTEP * PolygonGeometricalWetProperties.computeWaterVolume(ComputationalDomain.etaDirichlet[index], ComputationalDomain.bedRockElevation[index], ComputationalDomain.porosity[index], ComputationalDomain.planArea[index]);
				
			}
			
		}
		
		return volume;
		
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

		String fileName = null;
		DecimalFormat myformatter = computePattern();
		
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
		


		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
			
			TextIO.putln("Time step " + (double) t/3600);
			
			volumeDirichlet = computeVolumeDirichlet(volumeDirichlet);
			
			fileName = myformatter.format(t);
			FileWrite.openTxtFile(fileName.concat(".txt"), BoussinesqEquation.solutionDir, false);

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(cTDirichlet, cTNoDirichlet, cB);

			eta = newton.newtonIteration(arrb, matTNoDirichlet, indexDiag, eta,
					cg, tolerance);

			for (int j = 0; j < eta.length; j++) {
				
				volumeSource[j] = ComputationalDomain.source[j] * ComputationalDomain.planArea[j];
				aquiferThickness[j] = Math.max(0, eta[j]-ComputationalDomain.bedRockElevation[j]);
				
				volume[j] = computeVolume(j,eta[j]);
				volumeNew = volumeNew + volume[j];

			}
			
			System.out.println(volumeDirichlet);
			System.out.println(volumeNew);
			
			volumeNew = volumeNew - volumeDirichlet;
			
			writeSolution(t);
			
			if (Math.abs(volumeNew-volumeOld) > Math.pow(10, -6)){
				
				TextIO.putln("WARNING!!! The system is losing mass");
				TextIO.putln("The difference between initial volume and compute volume is: " + Math.abs(volumeNew-volumeOld));
				
			}
			
			volumeNew = 0;
			
		}


		System.out.println("Exit code");
	}
}

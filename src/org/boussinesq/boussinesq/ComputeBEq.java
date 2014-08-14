package org.boussinesq.boussinesq;

import java.io.IOException;
import java.text.DecimalFormat;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.Solver;
<<<<<<< HEAD
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.wordpress.growworkinghard.usefulClasses.FileWrite;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEq implements TimeSimulation {
=======
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
import org.partialDifferentialEquation.nonLinearParabolicPDE.AbstractPde;
import org.wordpress.growworkinghard.usefulClasses.FileWrite;
//import org.wordpress.growworkinghard.usefulClasses.TextIO;
//
//import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public abstract class ComputeBEq extends AbstractPde {
>>>>>>> thesis_structure

	double[] aquiferThickness;
	double[] volumeSource;
	protected int[] indexDiag;
	protected double tolerance;

	double[] volume;
	
	double volumeOld = 0;
	double volumeNew = 0;
	
	String fileName = null;
	
	public static long timeCompute;
	public static long timeSolver;
	
	
	Solver newton;
	RCConjugateGradient cg;
	DecimalFormat myformatter;
	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;
	
	public ComputeBEq() {
		

<<<<<<< HEAD
//		cg = new RCConjugateGradient(ComputationalDomain.Np);
		newton = new Solver();
		
		aquiferThickness = new double[ComputationalDomain.Np];
		volumeSource = new double[ComputationalDomain.Np];
		
		volume = new double[ComputationalDomain.Np];
		
		myformatter = computePattern();
		
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
=======
//		cg = new RCConjugateGradient(mesh.Np);
//		newton = new Solver();
		

		myformatter = computePattern();
		
//		rcIndexDiagonalElement = new RCIndexDiagonalElement();
//		cMEd = new MachineEpsilon();
>>>>>>> thesis_structure
		
	}
	
	public DecimalFormat computePattern(){
		
<<<<<<< HEAD
		int[] c = new int[String.valueOf(SIMULATIONTIME).length()];
=======
		int[] c = new int[String.valueOf(TimeSimulation.SIMULATIONTIME).length()];
>>>>>>> thesis_structure
		
		StringBuilder builder = new StringBuilder(c.length);
		
		for (int i:c){
			
			builder.append(c[i]);
			
		}
		
		String pattern = builder.toString();
		
		return new DecimalFormat(pattern);
		
	}
	
<<<<<<< HEAD
	public void writeSolution(int time, double[] eta, String bc, String simulation) throws IOException{
=======
	public void writeSolution(int time, double[] eta, AbstractRCAdjacencyMatrixBased mesh) throws IOException{
>>>>>>> thesis_structure
		
//		FileWrite.writeStringString("Type of simulation", simulation);
//		FileWrite.writeStringString("Type of boundary conditions ", bc);
//		FileWrite.writeStringIntString("Iteration number", (int) time/TIMESTEP + 1, "");
//		FileWrite.writeStringDoubleString("Timestep", TIMESTEP, "[s]");
//		FileWrite.writeStringDoubleString("Time of simulation", SIMULATIONTIME, "[s]");
//		FileWrite.writeStringDoubleString("Time of simulation", SIMULATIONTIME/60, "[min]");
//		FileWrite.writeStringDoubleString("Initial volume", volumeOld, "[ m^3]");
//		FileWrite.writeStringDoubleString("Total volume", volumeNew, "[m^3]");
//		FileWrite.writeStringDoubleString("Time convergence of CG", timeSolver/1000000000.0, "[s]");
//		FileWrite.writeFourStringColumn("PiezHead","AquifThick","WaterVolume","Source");
//		FileWrite.writeFourStringColumn("[m]", "[m]", "[m^3]", "[m^3/s]");
//		FileWrite.writeFourDoubleColumn(eta,aquiferThickness,volume,volumeSource);
//		FileWrite.writeStringDoubleString("Time: ", time, "[s]");
<<<<<<< HEAD
		FileWrite.writeTwoDoubleColumn(eta,ComputationalDomain.outflow);
=======
		FileWrite.writeOneDoubleColumn(eta);
>>>>>>> thesis_structure
//		FileWrite.writeOneDoubleColumn(eta);
		FileWrite.closeTxtFile();
		
	}
	
<<<<<<< HEAD
	public double computeVolume(int index,double eta){
		
		double volume;
		
		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta, ComputationalDomain.bedRockElevation[index], ComputationalDomain.porosity[index], ComputationalDomain.planArea[index]);
				
		volume = volume	- TIMESTEP * ComputationalDomain.planArea[index] * ComputationalDomain.source[index]
						+ TIMESTEP * ComputationalDomain.planArea[index] * ComputationalDomain.c[index] * Math.pow(volume/ComputationalDomain.planArea[index], ComputationalDomain.m[index]);
=======
	public double computeVolume(int index,double eta, AbstractRCAdjacencyMatrixBased mesh){
		
		double volume;
		
		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta, mesh.bedRockElevation[index], mesh.porosity[index], mesh.planArea[index]);
				
		volume = volume	- TimeSimulation.TIMESTEP * mesh.planArea[index] * mesh.source[index]
						+ TimeSimulation.TIMESTEP * mesh.planArea[index] * mesh.c[index] * Math.pow(volume/mesh.planArea[index], mesh.m[index]);
>>>>>>> thesis_structure
		
		return volume;
		
	}
		
	public void openTxtFile(int time) throws IOException{
		
		fileName = myformatter.format(time);
		FileWrite.openTxtFile(fileName.concat(".txt"), BoussinesqEquation.solutionDir, false);
		
	}
	
<<<<<<< HEAD
	public void computeOutputFeatures(double[] eta){
		
		for (int j = 0; j < eta.length; j++) {
			
			volumeSource[j] = ComputationalDomain.source[j] * ComputationalDomain.planArea[j];
			aquiferThickness[j] = Math.max(0, eta[j]-ComputationalDomain.bedRockElevation[j]);
			volume[j] = computeVolume(j,eta[j]);
=======
	public void computeOutputFeatures(double[] eta, AbstractRCAdjacencyMatrixBased mesh){
		
		
		double[] aquiferThickness = new double[mesh.polygonsNumber];
				
		double[] volume = new double[mesh.polygonsNumber];
		
		double[] volumeSource = new double[mesh.polygonsNumber];
		
		for (int j = 0; j < eta.length; j++) {
			
			volumeSource[j] = mesh.source[j] * mesh.planArea[j];
			aquiferThickness[j] = Math.max(0, eta[j]-mesh.bedRockElevation[j]);
			volume[j] = computeVolume(j,eta[j], mesh);
>>>>>>> thesis_structure
			volumeNew = volumeNew + volume[j];

		}
		
	}
	
<<<<<<< HEAD
	public void computeInitialVolume(){
		
//		for (int i = 0; i < ComputationalDomain.Np; i++){
//			
//			volume[i] = PolygonGeometricalWetProperties.computeWaterVolume(ComputationalDomain.eta[i], ComputationalDomain.bedRockElevation[i], ComputationalDomain.porosity[i], ComputationalDomain.planArea[i]);
//			volumeOld = volumeOld + volume[i];
//	
//		}
//
//		TextIO.putln("Initial volume: " + volumeOld);
		
	}
	
	public void computeVolumeConservation(){
		
//		if (Math.abs(volumeNew-volumeOld) > Math.pow(10, -6)){
//			
//			TextIO.putln("WARNING!!! The system is losing mass");
//			TextIO.putln("The difference between initial volume and compute volume is: " + Math.abs(volumeNew-volumeOld));
//	
//		}


		volumeNew = 0;
		
	}
	
//	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb) throws IterativeSolverDoubleNotConvergedException{
//		
//		double[] eta = new double[etaOld.length];
//		
//		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, cg,
//				tolerance);
//		
//		return eta;
//		
//	}
	
	public void firstThings(){
		
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(ComputationalDomain.Np,
				ComputationalDomain.Mp, ComputationalDomain.Mi);
=======
	public void firstThings(AbstractRCAdjacencyMatrixBased mesh){
		
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(mesh.polygonsNumber,
				mesh.Mp, mesh.Mi);
>>>>>>> thesis_structure

		tolerance = cMEd.computeMachineEpsilonDouble();
		
	}
	
<<<<<<< HEAD
	public static void main(String[] args) {
		
		ComputeBEq test = new ComputeBEq();
		
		System.out.println(test.cMEd.computeMachineEpsilonDouble()*100);

	}

=======
>>>>>>> thesis_structure
}

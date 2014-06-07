package org.boussinesq.boussinesq;

import java.io.IOException;
//import java.text.DecimalFormat;

import org.boussinesq.RowCompressedForm.DeleteRowColumnNullDiagonalEntry;
import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;

import org.wordpress.growworkinghard.usefulClasses.FileWrite;
import org.wordpress.growworkinghard.usefulClasses.TextIO;


public class ComputeBEq {

	protected double[] eta;
	protected double[] matT;
	
	protected double[] arrb1;
	protected double[] arrb;
	
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

	// protected Solver newton;
	protected RCConjugateGradient cg;

	protected RCIndexDiagonalElement rcIndexDiagonalElement;
	protected MachineEpsilon cMEd;

	protected DeleteRowColumnNullDiagonalEntry deleteRowColumn;
	protected ComputationalArrays computationalArrays;

	protected ComputeT computeT;

	// protected ComputeB cB;

	public ComputeBEq() {

		computeT = new ComputeT();

		cg = new RCConjugateGradient();
		// newton = new Solver();

		computationalArrays = new ComputationalArrays();
		deleteRowColumn = new DeleteRowColumnNullDiagonalEntry();

		aquiferThickness = new double[ComputationalDomain.Np];
		volumeSource = new double[ComputationalDomain.Np];

		volume = new double[ComputationalDomain.Np];

		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();

	}

	public void writeSolution(int time, double[] eta, String bc,
			String simulation) throws IOException {

		FileWrite.writeStringString("Type of simulation", simulation);
		FileWrite.writeStringString("Type of boundary conditions ", bc);
		FileWrite.writeStringIntString("Iteration number", (int) time
				/ BoussinesqEquation.TIMESTEP + 1, "");
		FileWrite.writeStringDoubleString("Timestep",
				BoussinesqEquation.TIMESTEP, "[s]");
		FileWrite.writeStringDoubleString("Time of simulation",
				BoussinesqEquation.SIMULATIONTIME, "[s]");
		FileWrite.writeStringDoubleString("Time of simulation",
				BoussinesqEquation.SIMULATIONTIME / 60, "[min]");
		FileWrite
				.writeStringDoubleString("Initial volume", volumeOld, "[ m^3]");
		FileWrite.writeStringDoubleString("Total volume", volumeNew, "[m^3]");
		FileWrite.writeStringDoubleString("Time convergence of CG",
				timeSolver / 1000000000.0, "[s]");
		FileWrite.writeFourStringColumn("PiezHead", "AquifThick",
				"WaterVolume", "Source");
		FileWrite.writeFourStringColumn("[m]", "[m]", "[m^3]", "[m^3/s]");
		FileWrite.writeFourDoubleColumn(eta, aquiferThickness, volume,
				volumeSource);
		FileWrite.closeTxtFile();

	}

	public double computeVolume(int index, double eta) {

		double volume;

		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta,
				ComputationalDomain.bedRockElevation[index],
				ComputationalDomain.porosity[index],
				ComputationalDomain.planArea[index]);

		volume = volume
				- BoussinesqEquation.TIMESTEP
				* ComputationalDomain.planArea[index]
				* ComputationalDomain.source[index]
				+ BoussinesqEquation.TIMESTEP
				* ComputationalDomain.planArea[index]
				* ComputationalDomain.c[index]
				* Math.pow(volume / ComputationalDomain.planArea[index],
						ComputationalDomain.m[index]);

		return volume;

	}

	public void openTxtFile(int time) throws IOException {

		fileName = BoussinesqEquation.myformatter.format(time);
		FileWrite.openTxtFile(fileName.concat(".txt"),
				BoussinesqEquation.solutionDirectory, false);

	}

	public void computeOutputFeatures(double[] eta) {

		int endForLoop = eta.length;

		for (int j = 0; j < endForLoop; j++) {

			volumeSource[j] = ComputationalDomain.source[j]
					* ComputationalDomain.planArea[j];
			aquiferThickness[j] = Math.max(0, eta[j]
					- ComputationalDomain.bedRockElevation[j]);
			volume[j] = computeVolume(j, eta[j]);
			volumeNew = volumeNew + volume[j];

		}

	}

	public void computeInitialVolume() {

		for (int i = 0; i < ComputationalDomain.Np; i++) {

			volume[i] = PolygonGeometricalWetProperties.computeWaterVolume(
					ComputationalDomain.eta[i],
					ComputationalDomain.bedRockElevation[i],
					ComputationalDomain.porosity[i],
					ComputationalDomain.planArea[i]);
			volumeOld = volumeOld + volume[i];

		}

		TextIO.putln("Initial volume: " + volumeOld);

	}

	public void computeVolumeConservation() {

		if (Math.abs(volumeNew - volumeOld) > Math.pow(10, -6)) {

			TextIO.putln("WARNING!!! The system is losing mass");
			TextIO.putln("The difference between initial volume and compute volume is: "
					+ Math.abs(volumeNew - volumeOld));

		}

		volumeNew = 0;

	}

	public void firstThings() {

		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
				ComputationalDomain.Np, ComputationalDomain.Mp,
				ComputationalDomain.Mi);

		tolerance = cMEd.computeMachineEpsilonDouble();

	}

	public static void main(String[] args) {

		ComputeBEq test = new ComputeBEq();

		System.out.println(test.cMEd.computeMachineEpsilonDouble() * 100);

	}

}

package org.boussinesq.boussinesq;

import java.io.IOException;

import org.boussinesq.boussinesq.computationalDomain.AbstractDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.boussinesq.rowCompressedForm.RCIndexDiagonalElement;
import org.interfacesPDE.nonLinearParabolicPDE.ObjectAssembler;
import org.wordpress.growworkinghard.usefulClasses.FileWrite;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

public abstract class AbstractComputeBEq implements ObjectAssembler {

	protected double[] eta;
	protected double[] actualArray_T;

	protected double[] completeArray_b;
	protected double[] actualArray_b;

	protected double[] aquiferThickness;
	protected double[] volumeSource;
	protected int[] indexDiag;
	protected double tolerance;

	protected double[] volume;

	private double volumeOld = 0;
	private double volumeNew = 0;

	String fileName = null;

	public static long timeCompute;
	public static long timeSolver;
	public static boolean unlockDeleteRowColumn = false;

	protected RCIndexDiagonalElement rcCompleteIndecesDiagonalTerms;
	protected MachineEpsilon machineEpsilon;

	public AbstractComputeBEq() {

		aquiferThickness = new double[AbstractDomain.Np];
		volumeSource = new double[AbstractDomain.Np];

		volume = new double[AbstractDomain.Np];

		rcCompleteIndecesDiagonalTerms = new RCIndexDiagonalElement();
		machineEpsilon = new MachineEpsilon();

	}

	protected void writeSolution(double time, double[] eta, String bc,
			String simulation) throws IOException {

		// FileWrite.writeStringString("Type of simulation", simulation);
		// FileWrite.writeStringString("Type of boundary conditions ", bc);
		// FileWrite.writeStringIntString("Iteration number", (int) time
		// / BoussinesqEquation.TIMESTEP + 1, "");
		// FileWrite.writeStringDoubleString("Timestep",
		// BoussinesqEquation.TIMESTEP, "[s]");
		// FileWrite.writeStringDoubleString("Time of simulation",
		// BoussinesqEquation.SIMULATIONTIME, "[s]");
		// FileWrite.writeStringDoubleString("Time of simulation",
		// BoussinesqEquation.SIMULATIONTIME / 60, "[min]");
		// FileWrite
		// .writeStringDoubleString("Initial volume", volumeOld, "[ m^3]");
		// FileWrite.writeStringDoubleString("Total volume", volumeNew,
		// "[m^3]");
		// FileWrite.writeStringDoubleString("Time convergence of CG",
		// timeSolver / 1000000000.0, "[s]");
		// FileWrite.writeFourStringColumn("PiezHead", "AquifThick",
		// "WaterVolume", "Source");
		// FileWrite.writeFourStringColumn("[m]", "[m]", "[m^3]", "[m^3/s]");
		// FileWrite.writeFourDoubleColumn(eta, aquiferThickness, volume,
		// volumeSource);
		FileWrite.writeOneDoubleColumn(eta);
		FileWrite.closeTxtFile();

	}

	private double computeVolume(int index, double eta) {

		double volume;

		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta,
				AbstractDomain.bedRockElevation[index],
				AbstractDomain.porosity[index], AbstractDomain.planArea[index]);

		volume = volume
				- BoussinesqEquation.TIMESTEP
				* AbstractDomain.planArea[index]
				* AbstractDomain.source[index]
				+ BoussinesqEquation.TIMESTEP
				* AbstractDomain.planArea[index]
				* AbstractDomain.c[index]
				* Math.pow(volume / AbstractDomain.planArea[index],
						AbstractDomain.m[index]);

		return volume;

	}

	protected void openTxtFile(double time) throws IOException {

		fileName = BoussinesqEquation.myformatter.format(time);
		FileWrite.openTxtFile(fileName.concat(".txt"),
				BoussinesqEquation.solutionDirectory, true);

	}

	protected void computeOutputVariables(double[] eta) {

		int endForLoop = eta.length;

		for (int j = 0; j < endForLoop; j++) {

			volumeSource[j] = AbstractDomain.source[j]
					* AbstractDomain.planArea[j];
			aquiferThickness[j] = Math.max(0, eta[j]
					- AbstractDomain.bedRockElevation[j]);
			volume[j] = computeVolume(j, eta[j]);
			volumeNew = volumeNew + volume[j];

		}

	}

	protected void computeInitialVolume() {

		for (int i = 0; i < AbstractDomain.Np; i++) {

			volume[i] = PolygonGeometricalWetProperties.computeWaterVolume(
					AbstractDomain.eta[i], AbstractDomain.bedRockElevation[i],
					AbstractDomain.porosity[i], AbstractDomain.planArea[i]);
			volumeOld = volumeOld + volume[i];

		}

		TextIO.putln("Initial volume: " + volumeOld);

	}

	protected void computeVolumeConservation() {

		if (Math.abs(volumeNew - volumeOld) > Math.pow(10, -6)) {

			TextIO.putln("WARNING!!! The system is losing mass");
			TextIO.putln("The difference between initial volume and compute volume is: "
					+ Math.abs(volumeNew - volumeOld));

		}

		volumeNew = 0;

	}

	protected void firstThings() {

		indexDiag = rcCompleteIndecesDiagonalTerms.computeIndexDiag(AbstractDomain.Np,
				AbstractDomain.Mp, AbstractDomain.Mi);

		tolerance = machineEpsilon.computeMachineEpsilonDouble();

	}

}

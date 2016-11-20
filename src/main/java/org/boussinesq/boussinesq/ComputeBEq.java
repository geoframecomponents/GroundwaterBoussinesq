package org.boussinesq.boussinesq;

import java.io.IOException;
import java.text.DecimalFormat;

import com.blogspot.geoframe.differentialEquation.partialDifferentialEquation
		.timeDependent.parabolic.nonLinear.AbstractPde;
import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AbstractRCAdjacencyMatrixBasedMesh;
import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;
import org.wordpress.growworkinghard.usefulClasses.FileWrite;

public abstract class ComputeBEq extends AbstractPde {

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

	public ComputeBEq() {

		myformatter = computePattern();

	}

	public DecimalFormat computePattern() {

		int[] c = new int[String.valueOf(TimeSimulation.SIMULATIONTIME)
				.length()];

		StringBuilder builder = new StringBuilder(c.length);

		for (int i : c) {

			builder.append(c[i]);

		}

		String pattern = builder.toString();

		return new DecimalFormat(pattern);

	}

	public void writeSolution(int time, double[] eta,
			AbstractRCAdjacencyMatrixBasedMesh mesh) throws IOException {

		double[] etaPlot = new double[mesh.polygonsNumber];
		
		for (int i=0; i< mesh.polygonsNumber; i++){
			
			etaPlot[i] = eta[i] - ((CatchmentDomain) mesh).bedRockElevation[i];
			
		}
		
		// FileWrite.writeStringString("Type of simulation", simulation);
		// FileWrite.writeStringString("Type of boundary conditions ", bc);
		// FileWrite.writeStringIntString("Iteration number", (int)
		// time/TIMESTEP + 1, "");
		// FileWrite.writeStringDoubleString("Timestep", TIMESTEP, "[s]");
		// FileWrite.writeStringDoubleString("Time of simulation",
		// SIMULATIONTIME, "[s]");
		// FileWrite.writeStringDoubleString("Time of simulation",
		// SIMULATIONTIME/60, "[min]");
		// FileWrite.writeStringDoubleString("Initial volume", volumeOld,
		// "[ m^3]");
		// FileWrite.writeStringDoubleString("Total volume", volumeNew,
		// "[m^3]");
		// FileWrite.writeStringDoubleString("Time convergence of CG",
		// timeSolver/1000000000.0, "[s]");
		// FileWrite.writeFourStringColumn("PiezHead","AquifThick","WaterVolume","Source");
		// FileWrite.writeFourStringColumn("[m]", "[m]", "[m^3]", "[m^3/s]");
		// FileWrite.writeFourDoubleColumn(eta,aquiferThickness,volume,volumeSource);
		// FileWrite.writeStringDoubleString("Time: ", time, "[s]");
		FileWrite.writeOneDoubleColumn(etaPlot);
		// FileWrite.writeOneDoubleColumn(eta);
		FileWrite.closeTxtFile();

	}

//	public double computeVolume(int index, double eta,
//			AbstractRCAdjacencyMatrixBased mesh) {
//
//		double volume;
//
//		volume = PolygonGeometricalWetProperties.computeWaterVolume(eta,
//				mesh.bedRockElevation[index], mesh.porosity[index],
//				mesh.planArea[index]);
//
//		volume = volume - TimeSimulation.TIMESTEP * mesh.planArea[index]
//				* mesh.source[index] + TimeSimulation.TIMESTEP
//				* mesh.planArea[index] * mesh.c[index]
//				* Math.pow(volume / mesh.planArea[index], mesh.m[index]);
//
//		return volume;
//
//	}

	public void openTxtFile(int time) throws IOException {

		fileName = myformatter.format(time);
		FileWrite.openTxtFile(fileName.concat(".txt"),
				BoussinesqEquation.solutionDir, false);

	}

//	public void computeOutputFeatures(double[] eta,
//			AbstractRCAdjacencyMatrixBased mesh) {
//
//		double[] aquiferThickness = new double[mesh.polygonsNumber];
//
//		double[] volume = new double[mesh.polygonsNumber];
//
//		double[] volumeSource = new double[mesh.polygonsNumber];
//
//		for (int j = 0; j < eta.length; j++) {
//
//			volumeSource[j] = mesh.source[j] * mesh.planArea[j];
//			aquiferThickness[j] = Math
//					.max(0, eta[j] - mesh.bedRockElevation[j]);
//			volume[j] = computeVolume(j, eta[j], mesh);
//			volumeNew = volumeNew + volume[j];
//
//		}
//
//	}

//	public void firstThings(AbstractRCAdjacencyMatrixBased mesh) {
//
//		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
//				mesh.polygonsNumber, mesh.Mp, mesh.Mi);
//
//		tolerance = cMEd.computeMachineEpsilonDouble();
//
//	}

}

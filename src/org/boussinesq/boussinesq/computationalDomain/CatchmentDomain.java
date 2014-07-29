package org.boussinesq.boussinesq.computationalDomain;

import java.io.FileNotFoundException;

import org.wordpress.growworkinghard.usefulClasses.FileRead;
import org.wordpress.growworkinghard.usefulClasses.GUIpathFileRead;

public class CatchmentDomain {

	public CatchmentDomain() throws FileNotFoundException {

		GUIpathFileRead gui = new GUIpathFileRead();

		computeAdjacencyMatrixFeatures(gui);
		computeSidesProperties(gui);
		computePolygonsProperties(gui);

	}

	public void computeAdjacencyMatrixFeatures(GUIpathFileRead gui)
			throws FileNotFoundException {

		FileRead readMp = new FileRead();
		/** The Mp. */
		ComputationalDomain.Mp = readMp
				.readIntArray(gui.openDialog("Mp array"));

		FileRead readMj = new FileRead();
		/** The Mi. */
		ComputationalDomain.Mi = readMj
				.readIntArray(gui.openDialog("Mj array"));

		FileRead readMl = new FileRead();
		/** The Ml. */
		ComputationalDomain.Ml = readMl.readDoubleArray(gui
				.openDialog("Ml array"));

	}

	public void computeSidesProperties(GUIpathFileRead gui)
			throws FileNotFoundException {

		FileRead readLS = new FileRead();
		/** The length sides. */
		ComputationalDomain.lengthSides = readLS.readDoubleArray(gui
				.openDialog("LENGTH SIDES OF POLYGONS array"));

		FileRead readEucD = new FileRead();
		/** The euclidean distance. */
		ComputationalDomain.euclideanDistance = readEucD.readDoubleArray(gui
				.openDialog("EUCLIDEAN DISTANCE array"));

		FileRead readHydrC = new FileRead();
		/** The hydr conductivity. */
		ComputationalDomain.hydrConductivity = readHydrC.readDoubleArray(gui
				.openDialog("HYDRAULIC CONDUCTIVITY array"));

	}

	public void computePolygonsProperties(GUIpathFileRead gui)
			throws FileNotFoundException {

		FileRead readEta = new FileRead();
		/** The eta. */
		ComputationalDomain.eta = readEta.readDoubleArray(gui
				.openDialog("HYDRAULIC HEAD array"));

		FileRead readPlanArea = new FileRead();
		/** The plan area. */
		ComputationalDomain.planArea = readPlanArea.readDoubleArray(gui
				.openDialog("PLANIMETRIC POLYGONS AREA array"));

		ComputationalDomain.Np = ComputationalDomain.planArea.length;

		FileRead readSource = new FileRead();
		/** The source. per unit area of the polygon */
		ComputationalDomain.rainHour = readSource.readDoubleArray(gui
				.openDialog("SOURCE array"));

		FileRead readED = new FileRead();
		/** The eta. */
		ComputationalDomain.etaDirichlet = readED.readDoubleArray(gui
				.openDialog("HYDRAULIC HEAD DIRICHLET BC array"));

		FileRead readBRE = new FileRead();
		/** The bottom elevation. */
		ComputationalDomain.bedRockElevation = readBRE.readDoubleArray(gui
				.openDialog("BEDROCK ELEVATION array"));

		FileRead readPor = new FileRead();
		ComputationalDomain.porosity = readPor.readDoubleArray(gui
				.openDialog("POROSITY array"));

		FileRead readC = new FileRead();
		ComputationalDomain.c = readC.readDoubleArray(gui
				.openDialog("C array - coeff of flow rate"));

		FileRead readM = new FileRead();
		ComputationalDomain.m = readM.readDoubleArray(gui
				.openDialog("M array - coeff of flow rate"));

	}

	public static void main(String[] args) throws FileNotFoundException {

		CatchmentDomain test = new CatchmentDomain();
		
		GUIpathFileRead gui = new GUIpathFileRead();

		test.computeAdjacencyMatrixFeatures(gui);
		test.computeSidesProperties(gui);
		test.computePolygonsProperties(gui);

	}

}

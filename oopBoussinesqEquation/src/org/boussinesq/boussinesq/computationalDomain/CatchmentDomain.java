package org.boussinesq.boussinesq.computationalDomain;

import java.io.FileNotFoundException;

import org.wordpress.growworkinghard.usefulClasses.FileRead;
import org.wordpress.growworkinghard.usefulClasses.GUIpathFileRead;

public class CatchmentDomain extends AbstractDomain {

	private GUIpathFileRead gui;

	public CatchmentDomain() throws FileNotFoundException {

		gui = new GUIpathFileRead();

		computeAdjacencyMatrixFeatures();
		computeSidesProperties();
		computePolygonsProperties();

	}

	public void computeAdjacencyMatrixFeatures() throws FileNotFoundException {

		FileRead readMp = new FileRead();
		/** The Mp. */
		AbstractDomain.Mp = readMp.readIntArray(gui.openDialog("Mp array"));

		FileRead readMj = new FileRead();
		/** The Mi. */
		AbstractDomain.Mi = readMj.readIntArray(gui.openDialog("Mj array"));

		FileRead readMl = new FileRead();
		/** The Ml. */
		AbstractDomain.Ml = readMl.readDoubleArray(gui.openDialog("Ml array"));

	}

	public void computeSidesProperties() throws FileNotFoundException {

		FileRead readLS = new FileRead();
		/** The length sides. */
		AbstractDomain.lengthSides = readLS.readDoubleArray(gui
				.openDialog("LENGTH SIDES OF POLYGONS array"));

		FileRead readEucD = new FileRead();
		/** The euclidean distance. */
		AbstractDomain.euclideanDistance = readEucD.readDoubleArray(gui
				.openDialog("EUCLIDEAN DISTANCE array"));

		FileRead readHydrC = new FileRead();
		/** The hydr conductivity. */
		AbstractDomain.hydrConductivity = readHydrC.readDoubleArray(gui
				.openDialog("HYDRAULIC CONDUCTIVITY array"));

	}

	public void computePolygonsProperties() throws FileNotFoundException {

		FileRead readEta = new FileRead();
		/** The eta. */
		AbstractDomain.eta = readEta.readDoubleArray(gui.openDialog("HYDRAULIC HEAD array"));

		FileRead readPlanArea = new FileRead();
		/** The plan area. */
		AbstractDomain.planArea = readPlanArea.readDoubleArray(gui
				.openDialog("PLANIMETRIC POLYGONS AREA array"));

		AbstractDomain.Np = planArea.length;

		FileRead readSource = new FileRead();
		/** The source. per unit area of the polygon */
		AbstractDomain.source = readSource.readDoubleArray(gui.openDialog("SOURCE array"));

		FileRead readED = new FileRead();
		/** The eta. */
		AbstractDomain.etaDirichlet = readED.readDoubleArray(gui
				.openDialog("HYDRAULIC HEAD DIRICHLET BC array"));

		FileRead readBRE = new FileRead();
		/** The bottom elevation. */
		AbstractDomain.bedRockElevation = readBRE.readDoubleArray(gui
				.openDialog("BEDROCK ELEVATION array"));

		FileRead readPor = new FileRead();
		AbstractDomain.porosity = readPor.readDoubleArray(gui.openDialog("POROSITY array"));

		FileRead readC = new FileRead();
		AbstractDomain.c = readC.readDoubleArray(gui
				.openDialog("C array - coeff of flow rate"));

		FileRead readM = new FileRead();
		AbstractDomain.m = readM.readDoubleArray(gui
				.openDialog("M array - coeff of flow rate"));

	}

}

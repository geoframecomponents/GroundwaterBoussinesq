package org.boussinesq.boussinesq.computationalDomain;

import java.io.File;
import java.io.FileNotFoundException;

import org.wordpress.growworkinghard.usefulClasses.FileRead;
import org.wordpress.growworkinghard.usefulClasses.GUIpathFileRead;

public class CatchmentDomain extends ComputationalDomain {

	GUIpathFileRead gui = new GUIpathFileRead();
	
	public String dataFolder;
	public File dataPath;
	
	public CatchmentDomain() throws FileNotFoundException {
		
		NOVALUE = -9999;
		
		getAdjacencyMatrix();
		getGridProperties();
		getPolygonProperties();
		getSideProperties();
	}
	
	public void getAdjacencyMatrix() throws FileNotFoundException {
		
		FileRead readMp = new FileRead();
		/** The Mp. */
		Mp = readMp
				.readIntArray(gui.openDialog("Mp array"));

		FileRead readMj = new FileRead();
		/** The Mi. */
		Mi = readMj
				.readIntArray(gui.openDialog("Mj array"));

		FileRead readMl = new FileRead();
		/** The Ml. */
		Ml = readMl.readDoubleArray(gui
				.openDialog("Ml array"));
		
	}

	public void getGridProperties() throws FileNotFoundException {
		
		FileRead readLS = new FileRead();
		/** The length sides. */
		lengthSides = readLS.readDoubleArray(gui
				.openDialog("LENGTH SIDES OF POLYGONS array"));

		FileRead readEucD = new FileRead();
		/** The euclidean distance. */
		euclideanDistance = readEucD.readDoubleArray(gui
				.openDialog("EUCLIDEAN DISTANCE array"));
		
		FileRead readPlanArea = new FileRead();
		/** The plan area. */
		planArea = readPlanArea.readDoubleArray(gui
				.openDialog("PLANIMETRIC POLYGONS AREA array"));

		polygonsNumber = planArea.length;
		
	}

	public void getPolygonProperties() throws FileNotFoundException {
		
		FileRead readEta = new FileRead();
		/** The eta. */
		eta = readEta.readDoubleArray(gui
				.openDialog("INITIAL HYDRAULIC HEAD array"));

		FileRead readSource = new FileRead();
		/** The source. per unit area of the polygon */
		rainHour = readSource.readDoubleArray(gui
				.openDialog("SOURCE array"));

		FileRead readED = new FileRead();
		/** The eta. */
		etaDirichlet = readED.readDoubleArray(gui
				.openDialog("HYDRAULIC HEAD DIRICHLET BC array"));

		FileRead readBRE = new FileRead();
		/** The bottom elevation. */
		bedRockElevation = readBRE.readDoubleArray(gui
				.openDialog("BEDROCK ELEVATION array"));

		FileRead readPor = new FileRead();
		porosity = readPor.readDoubleArray(gui
				.openDialog("POROSITY array"));

		FileRead readC = new FileRead();
		c = readC.readDoubleArray(gui
				.openDialog("C array - coeff of flow rate"));

		FileRead readM = new FileRead();
		m = readM.readDoubleArray(gui
				.openDialog("M array - coeff of flow rate"));
		
	}

	public void getSideProperties() throws FileNotFoundException {
		
		FileRead readHydrC = new FileRead();
		/** The hydr conductivity. */
		hydrConductivity = readHydrC.readDoubleArray(gui
				.openDialog("HYDRAULIC CONDUCTIVITY array"));
		
	}

	public static void main(String[] args) throws FileNotFoundException{
		
		new CatchmentDomain();
		
	}
	
}

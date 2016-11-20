package org.boussinesq.boussinesq.computationalDomain;

import java.io.File;
import java.io.FileNotFoundException;

import org.wordpress.growworkinghard.usefulClasses.FileRead;
import org.wordpress.growworkinghard.usefulClasses.GUIpathFileRead;

public class CatchmentDomain extends ComputationalDomain {

	GUIpathFileRead gui = new GUIpathFileRead();
	
	public String dataFolder;
	public File dataPath;

    public double NOVALUE;

	public double[] eta;
	public double[] rainHour;
	public double[] etaDirichlet;
	public double[] bedRockElevation;
	public double[] porosity;
	public double[] c;
	public double[] m;

    public double[] hydrConductivity;
	public double[] source;
	public double[] outflow;

	public CatchmentDomain() throws FileNotFoundException {
		
		NOVALUE = -9999;
		
		getAdjacencyMatrix();
		getGridProperties();
		getPolygonProperties();
		getSideProperties();
	}
	
	public void getAdjacencyMatrix() {
		
		FileRead readMp = new FileRead();
		/** The Mp. */
		try {
			Mp = readMp
                    .readIntArray(gui.openDialog("Mp array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readMj = new FileRead();
		/** The Mi. */
		try {
			Mi = readMj
                    .readIntArray(gui.openDialog("Mj array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readMl = new FileRead();
		/** The Ml. */
		try {
			Ml = readMl.readDoubleArray(gui
                    .openDialog("Ml array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	public void getGridProperties() {
		
		FileRead readLS = new FileRead();
		/** The length sides. */
		try {
			lengthSides = readLS.readDoubleArray(gui
                    .openDialog("LENGTH SIDES OF POLYGONS array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readEucD = new FileRead();
		/** The euclidean distance. */
		try {
			euclideanDistance = readEucD.readDoubleArray(gui
                    .openDialog("EUCLIDEAN DISTANCE array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readPlanArea = new FileRead();
		/** The plan area. */
		try {
			planArea = readPlanArea.readDoubleArray(gui
                    .openDialog("PLANIMETRIC POLYGONS AREA array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		polygonsNumber = planArea.length;
		
	}

	public void getPolygonProperties() {
		
		FileRead readEta = new FileRead();
		/** The eta. */
		try {
			eta = readEta.readDoubleArray(gui
                    .openDialog("INITIAL HYDRAULIC HEAD array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readSource = new FileRead();
		/** The source. per unit area of the polygon */
		try {
			rainHour = readSource.readDoubleArray(gui
                    .openDialog("SOURCE array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readED = new FileRead();
		/** The eta. */
		try {
			etaDirichlet = readED.readDoubleArray(gui
                    .openDialog("HYDRAULIC HEAD DIRICHLET BC array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readBRE = new FileRead();
		/** The bottom elevation. */
		try {
			bedRockElevation = readBRE.readDoubleArray(gui
                    .openDialog("BEDROCK ELEVATION array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readPor = new FileRead();
		try {
			porosity = readPor.readDoubleArray(gui
                    .openDialog("POROSITY array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readC = new FileRead();
		try {
			c = readC.readDoubleArray(gui
                    .openDialog("C array - coeff of flow rate"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

		FileRead readM = new FileRead();
		try {
			m = readM.readDoubleArray(gui
                    .openDialog("M array - coeff of flow rate"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	public void getSideProperties() {
		
		FileRead readHydrC = new FileRead();
		/** The hydr conductivity. */
		try {
			hydrConductivity = readHydrC.readDoubleArray(gui
                    .openDialog("HYDRAULIC CONDUCTIVITY array"));
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}

	}

	public static void main(String[] args) throws FileNotFoundException{
		
		new CatchmentDomain();
		
	}
	
}

package org.boussinesq.boussinesq.computationalDomain;

<<<<<<< HEAD
=======
import java.io.File;
>>>>>>> thesis_structure
import java.io.FileNotFoundException;

import org.wordpress.growworkinghard.usefulClasses.FileRead;
import org.wordpress.growworkinghard.usefulClasses.GUIpathFileRead;

<<<<<<< HEAD
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
=======
public class CatchmentDomain extends ComputationalDomain {

	GUIpathFileRead gui = new GUIpathFileRead();
	
	public String dataFolder;
	public File dataPath;
	
	public CatchmentDomain() throws FileNotFoundException {
		
		NOVALUE = -999;
		
		getAdjacencyMatrix();
		getGridProperties();
		getPolygonProperties();
		getSideProperties();
	}
	
	public void getAdjacencyMatrix() throws FileNotFoundException {
		
		FileRead readMp = new FileRead();
		/** The Mp. */
		Mp = readMp
>>>>>>> thesis_structure
				.readIntArray(gui.openDialog("Mp array"));

		FileRead readMj = new FileRead();
		/** The Mi. */
<<<<<<< HEAD
		ComputationalDomain.Mi = readMj
=======
		Mi = readMj
>>>>>>> thesis_structure
				.readIntArray(gui.openDialog("Mj array"));

		FileRead readMl = new FileRead();
		/** The Ml. */
<<<<<<< HEAD
		ComputationalDomain.Ml = readMl.readDoubleArray(gui
				.openDialog("Ml array"));

	}

	public void computeSidesProperties(GUIpathFileRead gui)
			throws FileNotFoundException {

		FileRead readLS = new FileRead();
		/** The length sides. */
		ComputationalDomain.lengthSides = readLS.readDoubleArray(gui
=======
		Ml = readMl.readDoubleArray(gui
				.openDialog("Ml array"));
		
	}

	public void getGridProperties() throws FileNotFoundException {
		
		FileRead readLS = new FileRead();
		/** The length sides. */
		lengthSides = readLS.readDoubleArray(gui
>>>>>>> thesis_structure
				.openDialog("LENGTH SIDES OF POLYGONS array"));

		FileRead readEucD = new FileRead();
		/** The euclidean distance. */
<<<<<<< HEAD
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
=======
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
				.openDialog("HYDRAULIC HEAD array"));

		FileRead readSource = new FileRead();
		/** The source. per unit area of the polygon */
		rainHour = readSource.readDoubleArray(gui
>>>>>>> thesis_structure
				.openDialog("SOURCE array"));

		FileRead readED = new FileRead();
		/** The eta. */
<<<<<<< HEAD
		ComputationalDomain.etaDirichlet = readED.readDoubleArray(gui
=======
		etaDirichlet = readED.readDoubleArray(gui
>>>>>>> thesis_structure
				.openDialog("HYDRAULIC HEAD DIRICHLET BC array"));

		FileRead readBRE = new FileRead();
		/** The bottom elevation. */
<<<<<<< HEAD
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

=======
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
	
>>>>>>> thesis_structure
}

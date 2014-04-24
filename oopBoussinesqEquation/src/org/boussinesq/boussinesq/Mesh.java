package org.boussinesq.boussinesq;

import java.io.File;
import java.io.FileNotFoundException;

import org.francescoS.usefulClasses.*;

import cern.colt.Arrays;
// TODO: Auto-generated Javadoc
/**
 * The Class Grid.
 */
public class Mesh {

	public static String outputPathBeqDirichlet;
	public static String outputPathBeqNoDirichlet;
	
	public String inputPlanArea;
	public String inputSource;
	public String inputEta;
	public String inputEtaDirichlet;
	public String inputBedRockElevation;
	public String inputPorosity;
	public String inputC;
	public String inputm;
	public String inputLengthSides;
	public String inputEuclideanDistance;
	public String inputHydrConductivity;
	public String inputMp;
	public String inputMi;
	public String inputMl;
	public String dataFolder;
	public File dataPath;

	// POLYGONS PROPERTIES

	/** The number sides polygon. */
	static int[] numberSidesPolygon;

	/** The plan area. */
	public static double[] planArea;

	public static int Np;

	/** The source. per unit area of the polygon */
	public static double[] source;

	/** The eta. */
	public static double[] eta;

	/** The eta. */
	public static double[] etaDirichlet;

	public static double NOVALUE = -9999;

	/** The bottom elevation. */
	public static double[] bedRockElevation;

	public static double[] porosity;

	public static double[] c;

	public static double[] m;

	// SIDES PROPERTIES

	/** The length sides. */
	static double[] lengthSides;

	/** The euclidean distance. */
	static double[] euclideanDistance;

	/** The hydr conductivity. */
	static double[] hydrConductivity;

	// ADJACENCY MATRIX PROPERTIES

	/** The Mp. */
	public static int[] Mp;

	/** The Mi. */
	public static int[] Mi;

	/** The Ml. */
	public static double[] Ml;

	Mesh(String configuration) throws FileNotFoundException {

		if (configuration.equals("Song")) {
			Mesh.outputPathBeqDirichlet = "Dirichlet_20d_ts3600_ks001.txt";
			Mesh.outputPathBeqNoDirichlet = "NoDirchlet_20d_ts3600_ks001.txt";

			int dim = 150;
			// POLYGONS PROPERTIES
			//Mesh.numberSidesPolygon = new int[dim];
			Mesh.planArea = new double[dim];
			Mesh.source = new double[dim];
			Mesh.eta = new double[dim];
			Mesh.etaDirichlet = new double[dim];
			Mesh.bedRockElevation = new double[dim];
			Mesh.porosity = new double[dim];
			Mesh.c = new double[dim];
			Mesh.m = new double[dim];
			Mesh.lengthSides = new double[dim + 1];
			Mesh.euclideanDistance = new double[dim + 1];
			Mesh.hydrConductivity = new double[dim + 1];

			Mesh.NOVALUE = -999;

			Mesh.Np = planArea.length;

			for (int i = 0; i < dim; i++) {
//				Mesh.numberSidesPolygon[i] = 2;
				Mesh.planArea[i] = 1;
				Mesh.source[i] = 0;
				Mesh.eta[i] = 0;
				Mesh.etaDirichlet[i] = -999;
				Mesh.bedRockElevation[i] = 0;
				Mesh.porosity[i] = 0.4;
				Mesh.lengthSides[i] = 1;
				Mesh.euclideanDistance[i] = 1;
				Mesh.hydrConductivity[i] = 0.001;
				Mesh.c[i] = 0;
				Mesh.m[i] = 1;

			}

			Mesh.c[dim - 1] = 1;
			Mesh.m[dim - 1] = 1;
			Mesh.etaDirichlet[0] = 1;
			Mesh.lengthSides[dim] = 1;
			Mesh.euclideanDistance[dim] = 1;
			Mesh.hydrConductivity[dim] = 0.001;

			// ADJACENCY MATRIX PROPERTIES

			Mesh.Mp = new int[dim + 1];
			Mesh.Mi = new int[dim * 3 - 2];
			Mesh.Ml = new double[dim * 3 - 2];
			Mesh.Mp[0] = 0;
			Mesh.Mp[1] = 2;

			; // System.out.println("Son qui");

			for (int i = 2; i < (dim); i++) {
				Mesh.Mp[i] = Mesh.Mp[i - 1] + 3;
			}

			Mesh.Mp[dim] = Mesh.Mp[dim - 1] + 1;

			Mesh.Mi[0] = 0;
			Mesh.Mi[1] = 1;
			int index = 0;
			for (int i = 1; i < (dim); i++) {

				for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {

					Mesh.Mi[j] = index;
					index++;
				}
				index = index - 2;
			}

			// System.out.println(this.Mi.length);
			Mesh.Mi[Mesh.Mi.length - 1] = Mesh.Mi[Mesh.Mi.length - 2] + 1;

			Mesh.Ml[0] = -1;
			Mesh.Ml[1] = 2;

			int ind = 2;

			for (int i = 1; i < dim; i++) {
				/*
				 * nested for-loop to analyze diagonal entries, which are
				 * identified by a negative number
				 */

				for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {

					if (Mesh.Mi[j] == i) {
						Mesh.Ml[j] = -1;
					} else {
						Mesh.Ml[j] = ind;
						ind++;
					}

				}
				ind = (int) Mesh.Ml[Mesh.Mp[i + 1] - 1];
			}

			Mesh.Ml[Mesh.Ml.length - 1] = -1;

			// System.out.println(Arrays.toString(this.Mp));
			// System.out.println(Arrays.toString(this.Mi));
			// System.out.println(Arrays.toString(this.Ml));

		} else {
			
			
//			/** The Mp. */
//			int[] trMi;
//
//			/** The Mi. */
//			int[] trMj;
//
//			/** The Ml. */
//			double[] trMl;
			
//			outputPathBeqDirichlet = "/home/francesco/desktop/tesiTest/provaDirichlet.txt";
//			outputPathBeqNoDirichlet= "/home/francesco/desktop/tesiTest/provaNoDirichlet.txt";
			
//			dataFolder = ReadFromScreen.readText("Write the path of the data folder");
//			
//			File dataPath = new File(dataFolder);
//			
//			inputPlanArea = "vPlanarArea";
//			inputSource = "vSource";
//			inputEta = "vEtaInitialCondV";
//			inputEtaDirichlet = "vEtaDrichelet";
//			inputBedRockElevation = "vBedrock";
//			inputPorosity = "vPorosity";
//			inputC = "vC";
//			inputm = "vM";
//			inputLengthSides = "vLengthSides";
//			inputEuclideanDistance = "vEuclideanDistance";
//			inputHydrConductivity = "vHydrConductivity";
//			inputMp = "Mp";
//			inputMi = "Mj";
//			inputMl = "Ml";
			
			GUIpathFileRead gui = new GUIpathFileRead();
			// POLYGONS PROPERTIES

			FileRead readEta = new FileRead();
			/** The eta. */
			Mesh.eta = readEta.readDoubleArray(gui.openDialog("HYDRAULIC HEAD array"));
			
			FileRead readPlanArea = new FileRead();
			/** The plan area. */
			Mesh.planArea = readPlanArea.readDoubleArray(gui.openDialog("PLANIMETRIC POLYGONS AREA array"));

			Mesh.Np = planArea.length;

			FileRead readSource = new FileRead();
			/** The source. per unit area of the polygon */
			Mesh.source = readSource.readDoubleArray(gui.openDialog("SOURCE array"));

			FileRead readED = new FileRead();
			/** The eta. */
			Mesh.etaDirichlet = readED.readDoubleArray(gui.openDialog("HYDRAULIC HEAD DIRICHLET BC array"));

			FileRead readBRE = new FileRead();
			/** The bottom elevation. */
			Mesh.bedRockElevation = readBRE.readDoubleArray(gui.openDialog("BEDROCK ELEVATION array"));

			FileRead readPor = new FileRead();
			Mesh.porosity = readPor.readDoubleArray(gui.openDialog("POROSITY array"));

			FileRead readC = new FileRead();
			Mesh.c = readC.readDoubleArray(gui.openDialog("C array - coeff of flow rate"));

			FileRead readM = new FileRead();
			Mesh.m = readM.readDoubleArray(gui.openDialog("M array - coeff of flow rate"));

			// SIDES PROPERTIES

			FileRead readLS = new FileRead();
			/** The length sides. */
			Mesh.lengthSides = readLS.readDoubleArray(gui.openDialog("LENGTH SIDES OF POLYGONS array"));

			FileRead readEucD = new FileRead();
			/** The euclidean distance. */
			Mesh.euclideanDistance = readEucD.readDoubleArray(gui.openDialog("EUCLIDEAN DISTANCE array"));

			FileRead readHydrC = new FileRead();
			/** The hydr conductivity. */
			Mesh.hydrConductivity = readHydrC.readDoubleArray(gui.openDialog("HYDRAULIC CONDUCTIVITY array"));

			// ADJACENCY MATRIX PROPERTIES

			FileRead readMp = new FileRead();
			/** The Mp. */
			Mesh.Mp = readMp.readIntArray(gui.openDialog("Mp array"));

			FileRead readMj = new FileRead();
			/** The Mi. */
			Mesh.Mi = readMj.readIntArray(gui.openDialog("Mj array"));

			FileRead readMl = new FileRead();
			/** The Ml. */
			Mesh.Ml = readMl.readDoubleArray(gui.openDialog("Ml array"));
			
		}

	}

	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) throws FileNotFoundException {

		@SuppressWarnings("unused")
		Mesh gridTest = new Mesh("NoSong");
		
		System.out.println("\nEta initial condition");
		System.out.println(Arrays.toString(Mesh.eta));
		
		System.out.println("\nPlanar Area");
		System.out.println(Arrays.toString(Mesh.planArea));
		
		System.out.println("\nNumber of Polygons");
		System.out.println(Mesh.Np);
		
		System.out.println("\nSource Term");
		System.out.println(Arrays.toString(Mesh.source));
		
		System.out.println("\nPorosity");
		System.out.println(Arrays.toString(Mesh.porosity));

		System.out.println("\nMi");
		System.out.println(Arrays.toString(Mesh.Mi));
		
		System.out.println("\nMl");
		System.out.println(Arrays.toString(Mesh.Ml));
		
		System.out.println("\nMp");
		System.out.println(Arrays.toString(Mesh.Mp));

		System.exit(0);


	}

}
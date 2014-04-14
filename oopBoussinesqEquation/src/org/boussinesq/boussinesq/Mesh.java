package org.boussinesq.boussinesq;

import java.io.FileNotFoundException;

import org.francescoS.usefulClasses.*;
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

	public static double NOVALUE = -999;

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
				Mesh.numberSidesPolygon[i] = 2;
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
			
			FileRead reader = new FileRead();
			outputPathBeqDirichlet = "/home/francesco/desktop/tesiTest/provaDirichlet.txt";
			outputPathBeqNoDirichlet= "/home/francesco/desktop/tesiTest/provaNoDirichlet.txt";
			
			inputPlanArea = "/home/francesco/desktop/tesiTest/giuseppeTest/vPlanarArea";
			inputSource = "/home/francesco/desktop/tesiTest/giuseppeTest/vSource";
			inputEta = "/home/francesco/desktop/tesiTest/giuseppeTest/vEtaInitialCondV";
			inputEtaDirichlet = "/home/francesco/desktop/tesiTest/giuseppeTest/vEtaDrichelet";
			inputBedRockElevation = "/home/francesco/desktop/tesiTest/giuseppeTest/vBedrock";
			inputPorosity = "/home/francesco/desktop/tesiTest/giuseppeTest/vPorosity";
			inputC = "/home/francesco/desktop/tesiTest/giuseppeTest/vC";
			inputm = "/home/francesco/desktop/tesiTest/giuseppeTest/vM";
//			inputLengthSides = "/home/francesco/desktop/tesiTest/giuseppeTest/";
//			inputEuclideanDistance = "/home/francesco/desktop/tesiTest/giuseppeTest/";
//			inputHydrConductivity = "/home/francesco/desktop/tesiTest/giuseppeTest/";
			inputMp = "/home/francesco/desktop/tesiTest/giuseppeTest/Mp";
			inputMi = "/home/francesco/desktop/tesiTest/giuseppeTest/Mj";
			inputMl = "/home/francesco/desktop/tesiTest/giuseppeTest/Ml";
			// POLYGONS PROPERTIES

			/** The plan area. */
			Mesh.planArea = reader.readDoubleArray(inputPlanArea);

			Mesh.Np = planArea.length;

			/** The source. per unit area of the polygon */
			Mesh.source = reader.readDoubleArray(inputSource);

			/** The eta. */
			Mesh.eta = reader.readDoubleArray(inputEta);

			/** The eta. */
			Mesh.etaDirichlet = reader.readDoubleArray(inputEtaDirichlet);

			/** The bottom elevation. */
			Mesh.bedRockElevation = reader.readDoubleArray(inputBedRockElevation);

			Mesh.porosity = reader.readDoubleArray(inputPorosity);

			Mesh.c = reader.readDoubleArray(inputC);

			Mesh.m = reader.readDoubleArray(inputm);

			// SIDES PROPERTIES

//			/** The length sides. */
//			Mesh.lengthSides;
//
//			/** The euclidean distance. */
//			Mesh.euclideanDistance;
//
//			/** The hydr conductivity. */
//			Mesh.hydrConductivity;

			// ADJACENCY MATRIX PROPERTIES

			/** The Mp. */
			Mesh.Mp = reader.readIntArray(inputMp);

			/** The Mi. */
			Mesh.Mi = reader.readIntArray(inputMi);

			/** The Ml. */
			Mesh.Ml = reader.readDoubleArray(inputMl);
			
		}

	}

	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 */
	public static void main(String[] args) {

		// Grid gridTest = new Grid("Song");
		/*
		 * Grid gridTest = new Grid();
		 * 
		 * System.out.println("");
		 * System.out.println("Array POLYGONS properties length");
		 * System.out.println("");
		 * 
		 * System.out.println("numberSidesPolygon array length: "+
		 * gridTest.numberSidesPolygon.length);
		 * System.out.println("planArea array length: "+
		 * gridTest.planArea.length);
		 * System.out.println("sourceSink array length: "+
		 * gridTest.source.length); System.out.println("eta array length: "+
		 * gridTest.eta.length);
		 * System.out.println("topElevation array length: "+
		 * gridTest.topElevation.length);
		 * System.out.println("bottomElevation array length: "+
		 * gridTest.bottomElevation.length);
		 * 
		 * System.out.println("");
		 * System.out.println("Array SIDES properties length");
		 * System.out.println("");
		 * 
		 * System.out.println("lengthSides array length: "+
		 * gridTest.lengthSides.length);
		 * System.out.println("euclideanDistance array length: "+
		 * gridTest.euclideanDistance.length);
		 * System.out.println("hydrConductivity array length: "+
		 * gridTest.hydrConductivity.length);
		 * 
		 * System.out.println("");
		 * System.out.println("ADJACENCY MATRICES properties length");
		 * System.out.println("");
		 * 
		 * System.out.println("Mp array length: "+ gridTest.Mp.length);
		 * System.out.println("Mi array length: "+ gridTest.Mi.length);
		 * System.out.println("Ml array length: "+ gridTest.Ml.length);
		 */

	}

}
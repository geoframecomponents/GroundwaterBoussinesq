package org.boussinesq.boussinesq;

// TODO: Auto-generated Javadoc
/**
 * The Class Grid.
 */
public class Mesh {

	public static String outputPathBeqDirichlet;
	public static String outputPathBeqNoDirichlet;

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

	Mesh(String configuration) {

		if (configuration.equals("test0")) {
			// POLYGONS PROPERTIES

			/** The number sides polygon. */
			/*
			 * Mesh.numberSidesPolygon = new int[]{4,4,4,4,4,4,4,4,
			 * 4,4,4,4,4,4,4,4};
			 *//** The plan area. */
			/*
			 * Mesh.planArea = new double[]{6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
			 *//** The source. per unit area of the polygon */
			/*
			 * Mesh.source = new double[]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			 *//** The eta. */
			/*
			 * Mesh.eta = new
			 * double[]{12,12,12,12,12,14,12,12,23,12,12,12,12,12,12,12};
			 *//** The eta. */
			/*
			 * Mesh.etaDirichlet = new
			 * double[]{13,-999,-999,-999,-999,-999,-999,
			 * -999,-999,-999,-999,-999, -999,-999,-999,-999};
			 * 
			 * Mesh.NOVALUE = -999;
			 *//** The top elevation. */
			/*
			 * Mesh.topElevation = new
			 * double[]{16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
			 *//** The bottom elevation. */
			/*
			 * Mesh.bottomElevation = new
			 * double[]{10,10,10,10,10,10,10,10,10,10, 10,10,10,10,10,10};
			 * 
			 * Mesh.porosity = new
			 * double[]{0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,
			 * 0.10,0.10,0.10,0.10,0.10,0.10};
			 * 
			 * 
			 * //SIDES PROPERTIES
			 *//** The length sides. */
			/*
			 * Mesh.lengthSides = new
			 * double[]{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			 * 3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
			 *//** The euclidean distance. */
			/*
			 * Mesh.euclideanDistance = new double[]{1.5,1.5,1.5, 1.5,3,3,3,1.5,
			 * 3,3,3,3,3, 1.5,3,3,3,1.5, 1.5,1.5,1.5, 1,2,2,1, 1,2,2,2,2,1,
			 * 1,2,2,2,2,1, 1,2,2,1};
			 *//** The hydr conductivity. */
			/*
			 * Mesh.hydrConductivity = new
			 * double[]{10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
			 * 10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
			 * 10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
			 * 10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
			 * 10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
			 * 10^(-5),10^(-5),10^(-5),10^(-5),10^(-5)};
			 * 
			 * 
			 * //ADJACENCY MATRIX PROPERTIES
			 *//** The Mp. */
			/*
			 * Mesh.Mp = new
			 * int[]{0,3,7,10,13,18,23,28,31,34,39,44,49,52,55,59,62};
			 *//** The Mi. */
			/*
			 * Mesh.Mi = new int[]{0,1,4, 0,1,2,5, 1,2,6, 3,4,8, 0,3,4,5,9,
			 * 1,4,5,6,10, 2,5,6,7,11, 6,7,12, 3,8,9, 4,8,9,10,13, 5,9,10,11,14,
			 * 6,10,11,12,15, 7,11,12, 9,13,14, 10,13,14,15, 11,14,15};
			 *//** The Ml. */
			/*
			 * Mesh.Ml =new double[]{-1,23,5, 23,-1,24,6, 24,-1,7, -1,27,9,
			 * 5,27,-1,28,10, 6,28,-1,29,11, 7,29,-1,30,12, 30,-1,13, 9,-1,33,
			 * 10,33,-1,34,15, 11,34,-1,35,16, 12,35,-1,36,17, 13,36,-1,
			 * 15,-1,39, 16,39,-1,40, 17,40,-1};
			 */
		} else {

			Mesh.outputPathBeqDirichlet = "/home/francesco/desktop/tesiTest/OOPbeqDirichlet_1d_static.txt";
			Mesh.outputPathBeqNoDirichlet = "/home/francesco/desktop/tesiTest/beqNoDirchlet5d_static.txt";

			int dim = 1000;
			// POLYGONS PROPERTIES
			Mesh.numberSidesPolygon = new int[dim];
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
				Mesh.hydrConductivity[i] = 0.01;
				Mesh.c[i] = 0;
				Mesh.m[i] = 1;

			}

			Mesh.c[dim - 1] = 1;
			Mesh.m[dim - 1] = 1;
			Mesh.etaDirichlet[0] = 1;
			Mesh.lengthSides[dim] = 1;
			Mesh.euclideanDistance[dim] = 1;
			Mesh.hydrConductivity[dim] = 0.01;

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
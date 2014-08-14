package org.boussinesq.boussinesq.computationalDomain;

import java.io.File;
<<<<<<< HEAD
import java.io.FileNotFoundException;
=======
//import java.io.FileNotFoundException;

import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
>>>>>>> thesis_structure

/**
 * The Class Grid.
 */
<<<<<<< HEAD
public class ComputationalDomain {

	public String dataFolder;
	public File dataPath;

	// POLYGONS PROPERTIES

	/** The plan area. */
	public static double[] planArea;

	public static int Np;

	/** The source. per unit area of the polygon */
	
	public static double[] outflow;
	
	public static double[] rainHour;
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
	public static double[] lengthSides;

	/** The euclidean distance. */
	public static double[] euclideanDistance;

	/** The hydr conductivity. */
	public static double[] hydrConductivity;

	// ADJACENCY MATRIX PROPERTIES

	/** The Mp. */
	public static int[] Mp;

	/** The Mi. */
	public static int[] Mi;

	/** The Ml. */
	public static double[] Ml;

	public static void callSongDomain(){
		
		new SongDomain();
		
	}
	
	public static void callCatchmentDomain() throws FileNotFoundException{
		
		new CatchmentDomain();
		
	}

	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 * @throws FileNotFoundException 
	 */
	public static void main(String[] args) {

		System.exit(0);


	}
=======
public abstract class ComputationalDomain extends AbstractRCAdjacencyMatrixBased {



>>>>>>> thesis_structure

}
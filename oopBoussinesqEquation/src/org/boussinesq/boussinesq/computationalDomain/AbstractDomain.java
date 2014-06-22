package org.boussinesq.boussinesq.computationalDomain;

import org.interfacesPDE.nonLinearParabolicPDE.UnstructuredMeshDomain;

/**
 * The Class Grid.
 */
public abstract class AbstractDomain implements UnstructuredMeshDomain {

	// POLYGONS PROPERTIES

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

}
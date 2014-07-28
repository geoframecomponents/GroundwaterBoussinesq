package org.boussinesq.boussinesq.computationalDomain;

import java.io.FileNotFoundException;

import org.interfacesPDE.nonLinearParabolicPDE.BoundaryConditions;
import org.interfacesPDE.nonLinearParabolicPDE.InitialConditions;
import org.meshNumericalMethods.unstructuredMesh.UnstructuredMeshDomain;

/**
 * The Class Grid.
 */
public abstract class AbstractDomain implements UnstructuredMeshDomain,
		BoundaryConditions, InitialConditions {

	// POLYGONS PROPERTIES

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

	/** The hydr conductivity. */
	public static double[] hydrConductivity;

	// ADJACENCY MATRIX PROPERTIES

	/** The Mp. */
	public static int[] Mp;

	/** The Mi. */
	public static int[] Mi;

	/** The Ml. */
	public static double[] Ml;

	/** The length sides. */
	public static double[] lengthSides;

	/** The euclidean distance. */
	public static double[] euclideanDistance;
	
	/** The plan area. */
	public static double[] planArea;
	
	

	public abstract void getGridProperties() throws FileNotFoundException;

	public abstract void getSideProperties() throws FileNotFoundException;

	public abstract void getPolygonProperties() throws FileNotFoundException;

	
	
	public abstract void getBoundaryConditions() throws FileNotFoundException;
	
	
	
	public abstract void getInitialConditions() throws FileNotFoundException;
	
	
}
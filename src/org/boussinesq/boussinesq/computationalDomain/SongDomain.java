package org.boussinesq.boussinesq.computationalDomain;

import java.io.FileNotFoundException;

//import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;

public class SongDomain extends ComputationalDomain {
	
	
//	NOVALUE = -999;
	
//	double[] hydrConductivity = new double[polygonsNumber + 1];
	
	public SongDomain() throws FileNotFoundException{

		// SIDES PROPERTIES

		/** The length sides. */
//		public double[] lengthSides;

		/** The euclidean distance. */
//		public double[] euclideanDistance;

		/** The hydr conductivity. */
//		double[] hydrConductivity;

		
		polygonsNumber = 1200;
		
		getAdjacencyMatrix();
		getGridProperties();
		getPolygonProperties();
		getSideProperties();
		
	}
	
	public void getAdjacencyMatrix() throws FileNotFoundException {
		
		Mp = new int[polygonsNumber + 1];
		Mi = new int[polygonsNumber * 3 - 2];
		Ml = new double[polygonsNumber * 3 - 2];
		Mp[0] = 0;
		Mp[1] = 2;

		for (int i = 2; i < (polygonsNumber); i++) {
			Mp[i] = Mp[i - 1] + 3;
		}

		Mp[polygonsNumber] = Mp[polygonsNumber - 1] + 2;

		Mi[0] = 0;
		Mi[1] = 1;
		int index = 0;
		for (int i = 1; i < (polygonsNumber); i++) {

			for (int j = Mp[i]; j < Mp[i + 1]; j++) {

				Mi[j] = index;
				index++;
			}
			index = index - 2;
		}

		Mi[Mi.length - 1] = Mi[Mi.length - 2] + 1;

		Ml[0] = -1;
		Ml[1] = 2;

		int ind = 2;

		for (int i = 1; i < polygonsNumber; i++) {
			/*
			 * nested for-loop to analyze diagonal entries, which are
			 * identified by a negative number
			 */

			for (int j = Mp[i]; j < Mp[i + 1]; j++) {

				if (Mi[j] == i) {
					Ml[j] = -1;
				} else {
					Ml[j] = ind;
					ind++;
				}

			}
			ind = (int) Ml[Mp[i + 1] - 1];
		}

		Ml[Ml.length - 1] = -1;
		
	}

	public void getGridProperties() throws FileNotFoundException {
		
		lengthSides = new double[polygonsNumber + 1];
		euclideanDistance = new double[polygonsNumber + 1];
		planArea = new double[polygonsNumber];
		
		for (int i = 0; i < polygonsNumber; i++) {
			
			lengthSides[i] = 1;
			euclideanDistance[i] = 1;
			planArea[i] = 1;

		}
		
		lengthSides[polygonsNumber] = 1;
		euclideanDistance[polygonsNumber] = 1;
	
	}

	public void getPolygonProperties() throws FileNotFoundException {
		
		
		source = new double[polygonsNumber];
		eta = new double[polygonsNumber];
		etaDirichlet = new double[polygonsNumber];
		bedRockElevation = new double[polygonsNumber];
		porosity = new double[polygonsNumber];
		c = new double[polygonsNumber];
		m = new double[polygonsNumber];
		
		for (int i = 0; i < polygonsNumber; i++) {
			
			source[i] = 0;
			eta[i] = 0;
			etaDirichlet[i] = -999;
			bedRockElevation[i] = 0;
			porosity[i] = 1;
			c[i] = 0;
			m[i] = 1;

		}
		

		c[polygonsNumber - 1] = 1;
		m[polygonsNumber - 1] = 1;
		etaDirichlet[0] = 1;
		
	}

	public void getSideProperties() throws FileNotFoundException {
		
		hydrConductivity = new double[polygonsNumber + 1];
		
		for (int i = 0; i < (polygonsNumber +1); i++) {
			
			hydrConductivity[i] = 1;

		}
		
	}

	
	public static void main(String[] args) throws FileNotFoundException {

		new SongDomain();
		
	}


}

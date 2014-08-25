package org.boussinesq.boussinesq;

import org.boussinesq.boussinesq.TimeSimulation;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
import org.partialDifferentialEquation.AbstractPdeTerm;


public class PdeTermT extends AbstractPdeTerm {

	public PdeTermT() {

		matrix = true;

	}

	public double computeArrayTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex) {

		return 0;

	}

	public double computeMatrixTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex, int sideIndex) {

		double temp = 0;

	
		if (mesh.Mi[sideIndex] != polygonIndex) {

//			System.out.println("-----------------------------------------");
//			System.out.println("Polygon Index= " + polygonIndex);
//			System.out.println("Side Index= " + sideIndex);
			
			temp = -TimeSimulation.TIMESTEP
					* (1 / mesh.euclideanDistance[(int) mesh.Ml[sideIndex] - 1])
					* mesh.hydrConductivity[(int) mesh.Ml[sideIndex] - 1]
					* mesh.lengthSides[(int) mesh.Ml[sideIndex] - 1]
					* Math.max(
							Math.max(
									0,
									u[mesh.Mi[sideIndex]]
											- mesh.bedRockElevation[mesh.Mi[sideIndex]]),
							Math.max(
									0,
									u[polygonIndex]
											- mesh.bedRockElevation[polygonIndex]));
			
//			System.out.println("T = " + temp);

		}

		return temp;
	}

}

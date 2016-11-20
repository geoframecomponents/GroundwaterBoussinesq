package org.boussinesq.boussinesq;

import com.blogspot.geoframe.differentialEquation.partialDifferentialEquation
		.AbstractPdeTerm;
import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AdjacencyMatrixBased;
import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;

public class PdeTermT extends AbstractPdeTerm {

	public PdeTermT() {

		matrix = true;

	}

	public double computeArrayTerm(double[] u,
								   AdjacencyMatrixBased mesh,
								   int polygonIndex) {

		return 0;

	}

	public double computeMatrixTerm(double[] u,
			AdjacencyMatrixBased mesh, int polygonIndex, int
											sideIndex) {

		double temp = 0;

	
		if (((CatchmentDomain) mesh).Mi[sideIndex] != polygonIndex) {

//			System.out.println("-----------------------------------------");
//			System.out.println("Polygon Index= " + polygonIndex);
//			System.out.println("Side Index= " + sideIndex);
			
			temp = -TimeSimulation.TIMESTEP
					* (1 / mesh.euclideanDistance[(int) mesh.Ml[sideIndex] - 1])
					* ((CatchmentDomain) mesh).hydrConductivity[(int) mesh
					.Ml[sideIndex] - 1]
					* mesh.lengthSides[(int) mesh.Ml[sideIndex] - 1]
					* Math.max(
							Math.max(
									0,
									u[((CatchmentDomain) mesh).Mi[sideIndex]]
											- ((CatchmentDomain) mesh)
											.bedRockElevation[(
													(CatchmentDomain) mesh).Mi[sideIndex]]),
							Math.max(
									0,
									u[polygonIndex]
											- ((CatchmentDomain )mesh).bedRockElevation[polygonIndex]));
			
//			System.out.println("T = " + temp);

		}

		return temp;
	}

}

package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import com.blogspot.geoframe.differentialEquation.partialDifferentialEquation
		.AbstractPdeTerm;
import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AdjacencyMatrixBased;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

public class PdeTermB extends AbstractPdeTerm {

	PdeTermB() {

		matrix = false;

	}
	
	public double computeArrayTerm(double[] u,
								   AdjacencyMatrixBased mesh, int
										   polygonIndex) {
		
		double volume = PolygonGeometricalWetProperties.computeWaterVolume(
				u[polygonIndex], ((CatchmentDomain) mesh).bedRockElevation[polygonIndex],
				((CatchmentDomain) mesh).porosity[polygonIndex],
				mesh.planArea[polygonIndex]);

		// delta t deve essere minore di 1/c
		double temp = volume
				+ TimeSimulation.TIMESTEP
				* mesh.planArea[polygonIndex]
				* ((CatchmentDomain) mesh).source[polygonIndex]
				- TimeSimulation.TIMESTEP
				* mesh.planArea[polygonIndex]
				* ((CatchmentDomain) mesh).c[polygonIndex]
				* Math.pow(volume / mesh.planArea[polygonIndex],
				((CatchmentDomain) mesh).m[polygonIndex]);

		if (temp < 0) {

			TextIO.putln("WARNING!!!\nThe element " + polygonIndex
					+ " of the array of known terms is NEGATIVE");

		}

		return temp;
		
	}

	public double computeMatrixTerm(double[] u,
			AdjacencyMatrixBased mesh, int polygonIndex, int sideIndex) {
		
		return 0;
	}

}

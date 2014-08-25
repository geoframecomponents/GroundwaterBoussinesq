package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
import org.partialDifferentialEquation.AbstractPdeTerm;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

public class PdeTermB extends AbstractPdeTerm {

	PdeTermB() {

		matrix = false;

	}
	
	public double computeArrayTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex) {
		
		double volume = PolygonGeometricalWetProperties.computeWaterVolume(
				u[polygonIndex], mesh.bedRockElevation[polygonIndex],
				mesh.porosity[polygonIndex],
				mesh.planArea[polygonIndex]);

		// delta t deve essere minore di 1/c
		double temp = volume
				+ TimeSimulation.TIMESTEP
				* mesh.planArea[polygonIndex]
				* mesh.source[polygonIndex]
				- TimeSimulation.TIMESTEP
				* mesh.planArea[polygonIndex]
				* mesh.c[polygonIndex]
				* Math.pow(volume / mesh.planArea[polygonIndex],
						mesh.m[polygonIndex]);

		if (temp < 0) {

			TextIO.putln("WARNING!!!\nThe element " + polygonIndex
					+ " of the array of known terms is NEGATIVE");

		}

		return temp;
		
	}

	public double computeMatrixTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex, int sideIndex) {
		
		return 0;
	}

}

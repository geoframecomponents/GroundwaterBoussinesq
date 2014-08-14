package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
import org.partialDifferentialEquation.AbstractPdeTerm;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

public class PdeTermB extends AbstractPdeTerm {

	PdeTermB(){
		
		matrix = false;
		
	}
	
	/**
	 * Compute B.
	 * 
	 * @desc this method computes the elements of the array b of the known terms
	 *       of the linear system for every cell according to the equations (19)
	 *       of [Cordano & Rigon, 2012]. By considering head-based boundary
	 *       condi- tions (Dirichlet), it's necessary to subtract the product
	 *       between T matrix of Dirichlet cells and eta array of Dirichlet
	 *       cells. Thus the right side of equation (32) is solved. It's
	 *       considered the outflow scale of the closure point too, thus it is
	 *       subtract the outgoing flow from of last cell at the right side of
	 *       equation (32)
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the the piezometric head
	 * @param Tdrichelet
	 *            the matrix T computes for the Dirichlet cells
	 * @param etaDirichlet
	 *            the array of known eta in the Dirichlet cells
	 * 
	 * @return the array b
	 */
//	public double[] computeB(double[] eta, double[] Tdirichlet) {
//
//		// declaration of the array that holds the known terms of the linear
//		// system
//		double[] arrB = new double[ComputationalDomain.polygonsNumber];
//
//		for (int i = 0; i < ComputationalDomain.polygonsNumber; i++) {
//			// compute the water volume stored in the cell
//			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
//					eta[i], ComputationalDomain.bedRockElevation[i], ComputationalDomain.porosity[i],
//					ComputationalDomain.planArea[i]);
//			// equation (19)
//			double sum = 0;
//			for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {
//				sum += Tdirichlet[j] * ComputationalDomain.etaDirichlet[ComputationalDomain.Mi[j]];
//			}
//
//			// delta t deve essere minore di 1/c
//			arrB[i] = volume + TimeSimulation.TIMESTEP * ComputationalDomain.planArea[i]
//					* ComputationalDomain.source[i] - sum - TimeSimulation.TIMESTEP
//					* ComputationalDomain.planArea[i] * ComputationalDomain.c[i]
//					* Math.pow(volume / ComputationalDomain.planArea[i], ComputationalDomain.m[i]);
//			
//			if (arrB[i] < 0){
//				
//				TextIO.putln("WARNING!!!\nThe element " + i + " of the array of known terms is NEGATIVE");
//				
//			}
//
//		}
//
//		return arrB;
//	}

	public double computeArrayTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex) {

		double volume = PolygonGeometricalWetProperties.computeWaterVolume(
				u[polygonIndex], mesh.bedRockElevation[polygonIndex], mesh.porosity[polygonIndex],
				mesh.planArea[polygonIndex]);


		// delta t deve essere minore di 1/c
		double temp = volume + TimeSimulation.TIMESTEP * mesh.planArea[polygonIndex]
				* mesh.source[polygonIndex] - TimeSimulation.TIMESTEP
				* mesh.planArea[polygonIndex] * mesh.c[polygonIndex]
				* Math.pow(volume / mesh.planArea[polygonIndex], mesh.m[polygonIndex]);

		if (temp < 0){
			
			TextIO.putln("WARNING!!!\nThe element " + polygonIndex + " of the array of known terms is NEGATIVE");
			
		}
		
		return temp;
	}

	public double computeMatrixTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex, int sideIndex) {
		
		return 0;
	}
	
}

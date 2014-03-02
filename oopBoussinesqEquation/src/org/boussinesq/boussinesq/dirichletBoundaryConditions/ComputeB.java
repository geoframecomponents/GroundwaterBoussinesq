package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.Mesh;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;

public class ComputeB implements TimeSimulation {

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
	public double[] computeB(double[] eta, double[] Tdirichlet) {

		// declaration of the array that holds the known terms of the linear
		// system
		double[] arrB = new double[Mesh.Np];

		for (int i = 0; i < Mesh.Np; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], Mesh.bedRockElevation[i], Mesh.porosity[i],
					Mesh.planArea[i]);
			// equation (19)
			double sum = 0;
			for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {
				sum += Tdirichlet[j] * Mesh.etaDirichlet[Mesh.Mi[j]];
			}

			// delta t deve essere minore di 1/c
			arrB[i] = volume + TIMESTEP * Mesh.planArea[i]
					* Mesh.source[i] - sum - TIMESTEP
					* Mesh.planArea[i] * Mesh.c[i]
					* Math.pow(volume / Mesh.planArea[i], Mesh.m[i]);

		}

		return arrB;
	}
	
}

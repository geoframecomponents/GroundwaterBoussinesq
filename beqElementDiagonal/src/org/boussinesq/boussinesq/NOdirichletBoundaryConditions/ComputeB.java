package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

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
	 * @param etaDrichelet
	 *            the array of known eta in the Dirichlet cells
	 * 
	 * @return the array b
	 */
	public double[] computeB(double[] eta) {

		// declaration of the array that holds the known terms of the linear
		// system
		double[] arrB = new double[ComputationalDomain.Np];

		for (int i = 0; i < ComputationalDomain.Np; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], ComputationalDomain.bedRockElevation[i],
					ComputationalDomain.porosity[i],
					ComputationalDomain.planArea[i]);

			// delta t deve essere minore di 1/c
			arrB[i] = volume
					+ TIMESTEP
					* ComputationalDomain.planArea[i]
					* ComputationalDomain.source[i]
					- TIMESTEP
					* ComputationalDomain.planArea[i]
					* ComputationalDomain.c[i]
					* Math.pow(volume / ComputationalDomain.planArea[i],
							ComputationalDomain.m[i]);

			ComputationalDomain.outflow[i] = TIMESTEP
					* ComputationalDomain.planArea[i]
					* ComputationalDomain.c[i]
					* Math.pow(volume / ComputationalDomain.planArea[i],
							ComputationalDomain.m[i]);
			
//			TextIO.putln(ComputationalDomain.outflow[i]);

			if (arrB[i] < 0) {

				TextIO.putln("WARNING!!!\nThe element " + i
						+ " of the array of known terms is NEGATIVE");

			}

		}

		return arrB;
	}

}

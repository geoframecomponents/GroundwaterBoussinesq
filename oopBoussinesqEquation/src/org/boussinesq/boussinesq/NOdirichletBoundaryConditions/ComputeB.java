package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.boussinesq.BoussinesqEquation;
import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

public class ComputeB {

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
		double[] arrB = new double[ComputationalArrays.numberOfPolygon];

		for (int i = 0; i < ComputationalArrays.numberOfPolygon; i++) {
			// compute the water volume stored in the cell
			double volume = PolygonGeometricalWetProperties.computeWaterVolume(
					eta[i], ComputationalArrays.bedRockElevation[i], ComputationalArrays.porosity[i],
					ComputationalArrays.planarArea[i]);

			// delta t deve essere minore di 1/c
			arrB[i] = volume
					+ BoussinesqEquation.TIMESTEP
					* ComputationalArrays.planarArea[i]
					* ComputationalArrays.source[i]
					- BoussinesqEquation.TIMESTEP
					* ComputationalArrays.planarArea[i]
					* ComputationalArrays.coeffC[i]
					* Math.pow(volume / ComputationalArrays.planarArea[i],
							ComputationalArrays.coeffM[i]);

			if (arrB[i] < 0) {

				TextIO.putln("WARNING!!!\nThe element " + i
						+ " of the array of known terms is NEGATIVE");

			}

		}

		return arrB;
	}

}

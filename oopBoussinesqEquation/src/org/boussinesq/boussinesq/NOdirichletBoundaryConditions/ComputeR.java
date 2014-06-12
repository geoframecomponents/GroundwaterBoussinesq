package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;

public class ComputeR {

	/**
	 * Compute R.
	 * 
	 * @desc this method computes the values of the residual function at every
	 *       iteration of the Newton's method for every cell according to the
	 *       equation (A3) of [Cordano & Rigon, 2012]. By analyzing every cell,
	 *       the residual function is computed only for polygons that are not
	 *       Dirichlet cells, otherwise the residual is imposed equal to zero,
	 *       so the known piezometric head of Dirichlet cells remains constant
	 *       during the Newton's loop.
	 * 
	 * @param arrT
	 *            the array of T in Row Compressed Form
	 * @param arrb
	 *            the array of known terms
	 * @param zetaBedrock
	 *            the bedrock elevation
	 * @param porosity
	 *            the porosity
	 * @param Np
	 *            the number of polygons
	 * @param Mp
	 *            the array that holds the number of non-zero entries in
	 *            adjacency matrix
	 * @param Mi
	 *            the array that holds the column indices of non-zero entries
	 * @param eta
	 *            the piezometric head
	 * @param planimetricArea
	 *            the planimetric area of the cell
	 * 
	 * @return the array of the residual function
	 */
	public double[] computeR(double[] arrT, double[] arrb, double[] eta) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[ComputationalArrays.numberOfPolygon];

		for (int i = 0; i < ComputationalArrays.numberOfPolygon; i++) {

			for (int j = ComputationalArrays.arrayMp[i]; j < ComputationalArrays.arrayMp[i + 1]; j++) {
				sum += arrT[j] * eta[ComputationalArrays.arrayMj[j]];
			}

			double waterVolume = PolygonGeometricalWetProperties.computeWaterVolume(eta[i],
							ComputationalArrays.bedRockElevation[i],
							ComputationalArrays.porosity[i],
							ComputationalArrays.planarArea[i]);
			// equation (A3)
			arrR[i] = waterVolume + sum - arrb[i];

			sum = 0;

		}

		return arrR;
	}
	
}

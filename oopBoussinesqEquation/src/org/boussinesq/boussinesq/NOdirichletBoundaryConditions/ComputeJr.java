package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import org.boussinesq.boussinesq.ComputationalArrays;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;

public class ComputeJr {

	/**
	 * Compute Jr.
	 * 
	 * @desc this method computes the Jacobian matrix of the water volume stored
	 *       into every cell. In this case the array Jr, in Row Compressed Form,
	 *       is evaluated like sum between array T and the wet area, according
	 *       the equation (A6) and (A7) of [Cordano & Rigon, 2012]. The array Jr
	 *       is a copy of T where only diagonal entries are summed to P, because
	 *       P is a diagonal matrix in Row Compressed Form too. These operations
	 *       are made only in case the Jacobian is computed only in a non
	 *       Dirichlet cell. Otherwise the volume of water stored is constant
	 *       with eta and P is equal to zero.
	 * 
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param arrT
	 *            the array of T in Row Compressed Form
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the zeta bedrock
	 * @param porosity
	 *            the porosity
	 * @param planimetricArea
	 *            the planimetric area of a cell
	 * @param etaDirichlet
	 *            the eta of Dirichlet cells
	 * @param NOVALUE
	 *            the novalue
	 * 
	 * @return the Jacobian array of water volume stored in Row Compressed Form
	 */
	public double[] computeJr(double[] arrT, double[] eta) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		double[] arrJr = new double[arrT.length];

		System.arraycopy(arrT, 0, arrJr, 0, arrT.length);

		int endForLoop = ComputationalArrays.indexDiagonal.length;
		
		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < endForLoop; i++) {

			// equation (A6)
			arrJr[ComputationalArrays.indexDiagonal[i]] = arrT[ComputationalArrays.indexDiagonal[i]]
					+ PolygonGeometricalWetProperties.computeWetArea(eta[i],
							ComputationalArrays.bedRockElevation[i], ComputationalArrays.porosity[i],
							ComputationalArrays.planarArea[i]);

		}

		return arrJr;
	}

}

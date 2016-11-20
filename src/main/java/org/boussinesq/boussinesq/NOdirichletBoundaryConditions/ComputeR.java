package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AbstractRCAdjacencyMatrixBasedMesh;
import org.boussinesq.boussinesq.PolygonGeometricalWetProperties;
import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;

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
	 *
	 * @return the array of the residual function
	 */
	public double[] computeR(double[] arrT, double[] arrb, double[] eta,
			AbstractRCAdjacencyMatrixBasedMesh mesh) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[mesh.polygonsNumber];

		for (int i = 0; i < mesh.polygonsNumber; i++) {

			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
				sum += arrT[j] * eta[mesh.Mi[j]];
			}

			double waterVolume = PolygonGeometricalWetProperties
					.computeWaterVolume(eta[i], ((CatchmentDomain) mesh).bedRockElevation[i],
							((CatchmentDomain) mesh).porosity[i], mesh.planArea[i]);
			// equation (A3)
			arrR[i] = waterVolume + sum - arrb[i];

			sum = 0;

		}

		return arrR;
	}

}

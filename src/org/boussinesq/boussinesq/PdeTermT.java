package org.boussinesq.boussinesq;

import org.boussinesq.boussinesq.TimeSimulation;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
import org.partialDifferentialEquation.AbstractPdeTerm;


public class PdeTermT extends AbstractPdeTerm {

	public PdeTermT() {

		matrix = true;

	}

	/**
	 * Compute T.
	 * 
	 * @desc this method computes the elements of the matrix T in Row Compressed
	 *       Form according to the equations (20) e (21) of [Cordano & Rigon,
	 *       2012]
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the piezometric head
	 * 
	 * @return the matrix T like array in Row Compressed Form
	 */
	// public double[] computeT(double[] eta) {
	//
	// /*
	// * variable to which sum the terms of matrix T (T is an array because is
	// * in RC-F) that are outside the diagonal; after investigation of the
	// * row of the matrix the value is stored in the diagonal of matrix T
	// */
	// double rowSum = 0;
	//
	// /* to identify the diagonal entry of matrix T in row-compressed form */
	// int index = 0;
	//
	// /*
	// * the matrix T is an array because this code uses the Row Compressed
	// * Form to stored sparse matrix
	// */
	// double[] arrayT = new double[ComputationalDomain.Ml.length];
	//
	// /* for-loop to analyze the mesh cell by cell */
	// for (int i = 0; i < ComputationalDomain.polygonsNumber; i++) {
	// /*
	// * nested for-loop to analyze shared edges between the i-th cell and
	// * the Mi[j]-th cell
	// */
	// for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i +
	// 1]; j++) {
	//
	// if (ComputationalDomain.Mi[j] != i) {
	// // equation (21)
	// arrayT[j] = -TimeSimulation.TIMESTEP
	// * (1 / ComputationalDomain.euclideanDistance[(int)
	// ComputationalDomain.Ml[j] - 1])
	// * ComputationalDomain.hydrConductivity[(int) ComputationalDomain.Ml[j] -
	// 1]
	// * ComputationalDomain.lengthSides[(int) ComputationalDomain.Ml[j] - 1]
	// * Math.max(
	// Math.max(
	// 0,
	// eta[ComputationalDomain.Mi[j]]
	// - ComputationalDomain.bedRockElevation[ComputationalDomain.Mi[j]]),
	// Math.max(
	// 0,
	// eta[i]
	// - ComputationalDomain.bedRockElevation[i]));
	//
	// rowSum += -arrayT[j];
	//
	// } else {
	// index = j;
	// }
	//
	// }
	// // equation (20)
	// if (rowSum == 0) {
	//
	// arrayT[index] = 1;
	//
	// } else {
	//
	// arrayT[index] = rowSum;
	// rowSum = 0;
	//
	// }
	// }
	//
	// return arrayT;
	// }

	public double computeArrayTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex) {

		return 0;

	}

	public double computeMatrixTerm(double[] u,
			AbstractRCAdjacencyMatrixBased mesh, int polygonIndex, int sideIndex) {

		double temp = 0;

		if (mesh.Mi[sideIndex] != polygonIndex) {

			temp = -TimeSimulation.TIMESTEP
					* (1 / mesh.euclideanDistance[(int) mesh.Ml[sideIndex] - 1])
					* mesh.hydrConductivity[(int) mesh.Ml[sideIndex] - 1]
					* mesh.lengthSides[(int) mesh.Ml[sideIndex] - 1]
					* Math.max(
							Math.max(
									0,
									u[mesh.Mi[sideIndex]]
											- mesh.bedRockElevation[mesh.Mi[sideIndex]]),
							Math.max(
									0,
									u[polygonIndex]
											- mesh.bedRockElevation[polygonIndex]));

		}

		return temp;
	}

}

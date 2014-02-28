package org.boussinesq.boussinesq;

public class ComputeT {
	
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
	public static double[] computeT(double[] eta) {

		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */
		double rowSum = 0;

		/* to identify the diagonal entry of matrix T in row-compressed form */
		int index = 0;

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[Mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Mesh.Np; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {

				if (Mesh.Mi[j] != i) {
					// equation (21)
					arrayT[j] = -BoussinesqEquation.deltat
							* (1 / Mesh.euclideanDistance[(int) Mesh.Ml[j]])
							* Mesh.hydrConductivity[(int) Mesh.Ml[j]]
							* Mesh.lengthSides[(int) Mesh.Ml[j]]
							* Math.max(
									Math.max(0, eta[Mesh.Mi[j]]
											- Mesh.bottomElevation[Mesh.Mi[j]]),
									Math.max(0, eta[i]
											- Mesh.bottomElevation[i]));

					rowSum += -arrayT[j];

				} else {
					index = j;
				}

			}
			// equation (20)
			arrayT[index] = rowSum;
			rowSum = 0;
		}

		return arrayT;
	}

}

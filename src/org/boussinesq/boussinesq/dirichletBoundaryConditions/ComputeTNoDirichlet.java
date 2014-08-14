package org.boussinesq.boussinesq.dirichletBoundaryConditions;

<<<<<<< HEAD
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
=======
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
>>>>>>> thesis_structure

public class ComputeTNoDirichlet extends IsNoValue {

	/**
	 * Compute T for non Dirichlet cells.
	 * 
	 * @desc according the head-based Boundary Conditions (Dirichlet) at the
	 *       equation (29) of [Cordano & Rigon, 2012], the matrix T in row
	 *       compressed form for a non Dirichlet cell is computed only if
	 *       neither the cell that I'm observing nor the adjacency cell are a
	 *       Dirichlet cell. Other T is imposed equal to zero.
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param T
	 *            the array of T in Row Compressed Form
	 * 
	 * @return the array of T in RC-F for non Dirichlet cells
	 */
<<<<<<< HEAD
	public double[] computeTNoDirichlet(double[] T, int[] indexDiag) {
=======
	public double[] computeTNoDirichlet(double[] T, int[] indexDiag,
			AbstractRCAdjacencyMatrixBased mesh) {
>>>>>>> thesis_structure

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
<<<<<<< HEAD
		double[] arrayT = new double[ComputationalDomain.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < ComputationalDomain.Np; i++) {

			if (isNoValue(ComputationalDomain.etaDirichlet[i], ComputationalDomain.NOVALUE)) {
=======
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.polygonsNumber; i++) {

			if (isNoValue(mesh.etaDirichlet[i], mesh.NOVALUE)) {
>>>>>>> thesis_structure

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
<<<<<<< HEAD
				for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {

					if (isNoValue(ComputationalDomain.etaDirichlet[ComputationalDomain.Mi[j]], ComputationalDomain.NOVALUE)) {
=======
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (isNoValue(mesh.etaDirichlet[mesh.Mi[j]], mesh.NOVALUE)) {
>>>>>>> thesis_structure

						// adjacent non Dirichlet cells
						arrayT[j] = T[j];
					} else {

						// adjacent Dirichlet cells
						arrayT[j] = 0.0;
					}
<<<<<<< HEAD
					
					if (j == indexDiag[i] && arrayT[j] == 0.0){
						
						arrayT[j] = 1;
						
					}
					
				}

			} else {
				
				for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {

					if (j == indexDiag[i] && arrayT[j] == 0.0){
						
						arrayT[j] = 1;
						
					}
					
				}
				
			}			
=======

					if (j == indexDiag[i] && arrayT[j] == 0.0) {

						arrayT[j] = 1;

					}

				}

			} else {

				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (j == indexDiag[i] && arrayT[j] == 0.0) {

						arrayT[j] = 1;

					}

				}

			}
>>>>>>> thesis_structure

		}

		return arrayT;
	}
<<<<<<< HEAD
	
=======

>>>>>>> thesis_structure
}

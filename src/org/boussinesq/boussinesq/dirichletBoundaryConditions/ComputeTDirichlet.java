package org.boussinesq.boussinesq.dirichletBoundaryConditions;

<<<<<<< HEAD
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
=======
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
>>>>>>> thesis_structure

public class ComputeTDirichlet extends IsNoValue {

	/**
	 * Compute T for Dirichlet cells.
	 * 
	 * @desc according the head-based Boundary Conditions (Dirichlet) at the
	 *       equation (30) of [Cordano & Rigon, 2012], the matrix T in row
	 *       compressed form for a Dirichlet cell is computed only if the cell
	 *       that I'm observing or the adjacency cell are a Dirichlet cell.
	 *       Other T is imposed equal to zero.
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param T
	 *            the array of T in Row Compressed Form
	 * 
	 * @return the array of T in RC-F for Dirichlet cells
	 */
<<<<<<< HEAD
	public double[] computeTDirichlet(double[] T) {
=======
	public double[] computeTDirichlet(double[] T,
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

			if (!isNoValue(ComputationalDomain.etaDirichlet[i], ComputationalDomain.NOVALUE)) {

				// Dirichlet cells
				for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {
=======
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.polygonsNumber; i++) {

			if (!isNoValue(mesh.etaDirichlet[i], mesh.NOVALUE)) {

				// Dirichlet cells
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
>>>>>>> thesis_structure
					arrayT[j] = T[j];
				}
			} else {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
<<<<<<< HEAD
				for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {

					if (!isNoValue(ComputationalDomain.etaDirichlet[ComputationalDomain.Mi[j]], ComputationalDomain.NOVALUE)) {
=======
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (!isNoValue(mesh.etaDirichlet[mesh.Mi[j]], mesh.NOVALUE)) {
>>>>>>> thesis_structure

						// adjacent Dirichlet cell
						arrayT[j] = T[j];
					} else {

						// adjacent non Dirichlet cell
						arrayT[j] = 0.0;
					}
				}

			}

		}

		return arrayT;
	}
<<<<<<< HEAD
	
=======

>>>>>>> thesis_structure
}

package org.boussinesq.boussinesq.dirichletBoundaryConditions;

<<<<<<< HEAD
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;

public class EtaInitialization extends IsNoValue {

	public double[] etaInitialization(double[] eta){
		
		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(ComputationalDomain.etaDirichlet[i], ComputationalDomain.NOVALUE)) {
=======
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;

public class EtaInitialization extends IsNoValue {

	public double[] etaInitialization(double[] eta, AbstractRCAdjacencyMatrixBased mesh){
		
		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(mesh.etaDirichlet[i], mesh.NOVALUE)) {
>>>>>>> thesis_structure

				// not Dirichlet cells
				eta[i] = eta[i];
			} else {

				// Dirichlet cells
<<<<<<< HEAD
				eta[i] = ComputationalDomain.etaDirichlet[i];
=======
				eta[i] = mesh.etaDirichlet[i];
>>>>>>> thesis_structure
			}
		}
		
		return eta;
		
	}
	
}

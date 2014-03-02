package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.Mesh;

public class EtaInitialization extends IsNoValue {

	public double[] etaInitialization(double[] eta){
		
		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(Mesh.etaDirichlet[i], Mesh.NOVALUE)) {

				// not Dirichlet cells
				eta[i] = eta[i];
			} else {

				// Dirichlet cells
				eta[i] = Mesh.etaDirichlet[i];
			}
		}
		
		return eta;
		
	}
	
}

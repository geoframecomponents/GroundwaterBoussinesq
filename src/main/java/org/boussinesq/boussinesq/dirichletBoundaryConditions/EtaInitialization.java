package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;

public class EtaInitialization extends IsNoValue {

	public double[] etaInitialization(double[] eta,
									  CatchmentDomain mesh){
		
		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(mesh.etaDirichlet[i], mesh.NOVALUE)) {

				// not Dirichlet cells
				eta[i] = eta[i];
			} else {

				// Dirichlet cells
				eta[i] = mesh.etaDirichlet[i];
			}
		}
		
		return eta;
		
	}
	
}

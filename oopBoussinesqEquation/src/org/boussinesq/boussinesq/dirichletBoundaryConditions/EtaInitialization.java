package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.computationalDoman.ComputationalDomain;

public class EtaInitialization extends IsNoValue {

	public double[] etaInitialization(double[] eta){
		
		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(ComputationalDomain.etaDirichlet[i], ComputationalDomain.NOVALUE)) {

				// not Dirichlet cells
				eta[i] = eta[i];
			} else {

				// Dirichlet cells
				eta[i] = ComputationalDomain.etaDirichlet[i];
			}
		}
		
		return eta;
		
	}
	
}

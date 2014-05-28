package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;

public class EtaInitialization extends IsNoValue {

	public double[] etaInitialization(double[] eta){
		
		int endForLoop = eta.length;
		
		// initialize eta array
		for (int i = 0; i < endForLoop; i++) {
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

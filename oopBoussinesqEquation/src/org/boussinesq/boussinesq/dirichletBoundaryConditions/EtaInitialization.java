package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import org.boussinesq.boussinesq.computationalDomain.AbstractDomain;

public class EtaInitialization {

	public double[] etaInitialization(double[] eta, IsNoValue verifyDirichlet) {

		int endForLoop = eta.length;

		// initialize eta array
		for (int i = 0; i < endForLoop; i++) {
			if (verifyDirichlet.isNoValue(AbstractDomain.etaDirichlet[i],
					AbstractDomain.NOVALUE)) {

				// not Dirichlet cells
				eta[i] = eta[i];
			} else {

				// Dirichlet cells
				eta[i] = AbstractDomain.etaDirichlet[i];
			}
		}

		return eta;

	}

}

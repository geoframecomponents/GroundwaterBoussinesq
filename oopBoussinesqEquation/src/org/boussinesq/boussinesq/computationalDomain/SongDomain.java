package org.boussinesq.boussinesq.computationalDomain;

public class SongDomain extends AbstractDomain {
	
	int dim;
	
	public SongDomain() {

		AbstractDomain.Np = 1000;
		dim = Np;

		computeAdjacencyMatrixFeatures();
		computeSidesProperties();
		computePolygonsProperties();
	}

	public void computeAdjacencyMatrixFeatures() {

		AbstractDomain.Mp = new int[dim + 1];
		AbstractDomain.Mi = new int[dim * 3 - 2];
		AbstractDomain.Ml = new double[dim * 3 - 2];
		AbstractDomain.Mp[0] = 0;
		AbstractDomain.Mp[1] = 2;

		for (int i = 2; i <= (dim); i++) {
			AbstractDomain.Mp[i] = AbstractDomain.Mp[i - 1] + 3;
		}

		AbstractDomain.Mp[dim] = AbstractDomain.Mp[dim - 1] + 2;

		AbstractDomain.Mi[0] = 0;
		AbstractDomain.Mi[1] = 1;
		int index = 0;
		for (int i = 1; i < (dim); i++) {

			for (int j = AbstractDomain.Mp[i]; j < AbstractDomain.Mp[i + 1]; j++) {

				AbstractDomain.Mi[j] = index;
				index++;
			}
			index = index - 2;
		}

		AbstractDomain.Mi[AbstractDomain.Mi.length - 1] = AbstractDomain.Mi[AbstractDomain.Mi.length - 2] + 1;

		AbstractDomain.Ml[0] = -1;
		AbstractDomain.Ml[1] = 2;

		int ind = 2;

		for (int i = 1; i < dim; i++) {
			/*
			 * nested for-loop to analyze diagonal entries, which are identified
			 * by a negative number
			 */

			for (int j = AbstractDomain.Mp[i]; j < AbstractDomain.Mp[i + 1]; j++) {

				if (AbstractDomain.Mi[j] == i) {
					AbstractDomain.Ml[j] = -1;
				} else {
					Ml[j] = ind;
					ind++;
				}

			}
			ind = (int) AbstractDomain.Ml[AbstractDomain.Mp[i + 1] - 1];
		}

		AbstractDomain.Ml[AbstractDomain.Ml.length - 1] = -1;

	}

	public void computeSidesProperties() {

		AbstractDomain.lengthSides = new double[dim + 1];
		AbstractDomain.euclideanDistance = new double[dim + 1];
		AbstractDomain.hydrConductivity = new double[dim + 1];

		for (int i = 0; i < (dim + 1); i++) {
			AbstractDomain.lengthSides[i] = 1;
			AbstractDomain.euclideanDistance[i] = 1;
			AbstractDomain.hydrConductivity[i] = 0.01;

		}

	}

	public void computePolygonsProperties() {

		AbstractDomain.planArea = new double[dim];
		AbstractDomain.source = new double[dim];
		AbstractDomain.eta = new double[dim];
		AbstractDomain.etaDirichlet = new double[dim];
		AbstractDomain.bedRockElevation = new double[dim];
		AbstractDomain.porosity = new double[dim];
		AbstractDomain.c = new double[dim];
		AbstractDomain.m = new double[dim];

		for (int i = 0; i < dim; i++) {
			
			AbstractDomain.planArea[i] = 1;
			AbstractDomain.source[i] = 0;
			AbstractDomain.eta[i] = 0;
			AbstractDomain.etaDirichlet[i] = -9999;
			AbstractDomain.bedRockElevation[i] = 0;
			AbstractDomain.porosity[i] = 0.4;
			AbstractDomain.c[i] = 0;
			AbstractDomain.m[i] = 1;

		}

		AbstractDomain.c[dim - 1] = 1;
		AbstractDomain.m[dim - 1] = 1;
		AbstractDomain.etaDirichlet[0] = 1;

	}

}

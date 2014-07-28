package org.boussinesq.boussinesq.computationalDomain;

public class SongDomain {
	
	
	
	SongDomain(){

		ComputationalDomain.Np = 1200;
		ComputationalDomain.NOVALUE = -999;
		
		computeAdjacencyMatrixFeatures(ComputationalDomain.Np);
		computeSidesProperties(ComputationalDomain.Np);
		computePolygonsProperties(ComputationalDomain.Np);
	}
	
	public void computeAdjacencyMatrixFeatures(int dim){
		
		ComputationalDomain.Mp = new int[dim + 1];
		ComputationalDomain.Mi = new int[dim * 3 - 2];
		ComputationalDomain.Ml = new double[dim * 3 - 2];
		ComputationalDomain.Mp[0] = 0;
		ComputationalDomain.Mp[1] = 2;

		for (int i = 2; i < (dim); i++) {
			ComputationalDomain.Mp[i] = ComputationalDomain.Mp[i - 1] + 3;
		}

		ComputationalDomain.Mp[dim] = ComputationalDomain.Mp[dim - 1] + 2;

		ComputationalDomain.Mi[0] = 0;
		ComputationalDomain.Mi[1] = 1;
		int index = 0;
		for (int i = 1; i < (dim); i++) {

			for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {

				ComputationalDomain.Mi[j] = index;
				index++;
			}
			index = index - 2;
		}

		ComputationalDomain.Mi[ComputationalDomain.Mi.length - 1] = ComputationalDomain.Mi[ComputationalDomain.Mi.length - 2] + 1;

		ComputationalDomain.Ml[0] = -1;
		ComputationalDomain.Ml[1] = 2;

		int ind = 2;

		for (int i = 1; i < dim; i++) {
			/*
			 * nested for-loop to analyze diagonal entries, which are
			 * identified by a negative number
			 */

			for (int j = ComputationalDomain.Mp[i]; j < ComputationalDomain.Mp[i + 1]; j++) {

				if (ComputationalDomain.Mi[j] == i) {
					ComputationalDomain.Ml[j] = -1;
				} else {
					ComputationalDomain.Ml[j] = ind;
					ind++;
				}

			}
			ind = (int) ComputationalDomain.Ml[ComputationalDomain.Mp[i + 1] - 1];
		}

		ComputationalDomain.Ml[ComputationalDomain.Ml.length - 1] = -1;

		
	}
	
	public void computeSidesProperties(int dim){
		
		ComputationalDomain.lengthSides = new double[dim + 1];
		ComputationalDomain.euclideanDistance = new double[dim + 1];
		ComputationalDomain.hydrConductivity = new double[dim + 1];
		
		for (int i = 0; i < (dim +1); i++) {
			ComputationalDomain.lengthSides[i] = 1;
			ComputationalDomain.euclideanDistance[i] = 1;
			ComputationalDomain.hydrConductivity[i] = 1;

		}

	}

	
	public void computePolygonsProperties(int dim){
		

		
		ComputationalDomain.planArea = new double[dim];
		ComputationalDomain.source = new double[dim];
		ComputationalDomain.eta = new double[dim];
		ComputationalDomain.etaDirichlet = new double[dim];
		ComputationalDomain.bedRockElevation = new double[dim];
		ComputationalDomain.porosity = new double[dim];
		ComputationalDomain.c = new double[dim];
		ComputationalDomain.m = new double[dim];
		
		for (int i = 0; i < dim; i++) {
			ComputationalDomain.planArea[i] = 1;
			ComputationalDomain.source[i] = 0;
			ComputationalDomain.eta[i] = 0;
			ComputationalDomain.etaDirichlet[i] = -999;
			ComputationalDomain.bedRockElevation[i] = 0;
			ComputationalDomain.porosity[i] = 1;
			ComputationalDomain.c[i] = 0;
			ComputationalDomain.m[i] = 1;

		}
		

		ComputationalDomain.c[dim - 1] = 1;
		ComputationalDomain.m[dim - 1] = 1;
		ComputationalDomain.etaDirichlet[0] = 1;
		
	}

	public static void main(String[] args) {

	}

}

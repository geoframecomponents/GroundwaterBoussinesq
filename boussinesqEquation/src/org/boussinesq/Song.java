package org.boussinesq;

import java.util.Arrays;

public class Song {

	double[] x;
	// double[] porosity;
	double hydraulicConductivity;
	int t;
	int alpha;
	static int h1;
	int nmax;

	Song(int dim, int time, double ks) {

		x = new double[dim];
		t = time;
		alpha = 0;
		h1 = 1;
		nmax = 10;
		hydraulicConductivity = ks;
		
		for (int i = 0; i<x.length;i++){
			
			x[i] = i+1;
			
		}

	}

	public double[] computeA(double[] ax, double xi0) {

		double[] a = new double[ax.length];

		for (int i = 0; i < a.length; i++) {

			a[i] = ax[i] * Math.pow(xi0, 2);

		}
		
		return a;

	}

	public double[] computeX(Grid mesh) {

		double[] xi = new double[x.length];

		for (int i = 0; i < xi.length; i++) {

			xi[i] = x[i]
					* Math.pow(2 * mesh.porosity[i] * (alpha + 1) / (h1
							* hydraulicConductivity * Math.pow(t, alpha + 1)),
							0.5);

		}
		
		//System.out.println(Arrays.toString(xi));
		return xi;
	}

	public double[] beqSong(Grid mesh) {

		double[] ax = new double[nmax];
		double[] solutionDimensionless = new double[x.length];
		double[] solution = new double[x.length];

		ax = SongCoefficient.CoefficientSongSolution(nmax, alpha / (alpha + 1));

		//System.out.println("Ax" + Arrays.toString(ax));
		
		double xi0 = 0;
		double sum = 0;

		for (int i = 0; i < ax.length; i++) {

			sum += ax[i];
			

		}

		xi0 = Math.pow(sum, -0.5);
		
		//System.out.println(xi0);
		
		solutionDimensionless = SongDimensionless.beqSongDimensionless(
				computeX(mesh), xi0, computeA(ax, xi0));

		for (int i = 0; i < solutionDimensionless.length; i++) {

			solution[i] = h1 * Math.pow(t, alpha) * solutionDimensionless[i];

		}

		return solution;

	}

	public static void main(String[] args) {
		int time = 3600 * 24;
		Grid mesh = new Grid("Song");
		Song s = new Song(mesh.numberSidesPolygon.length, time,
				0.01);
		System.out.println(Arrays.toString(s.beqSong(mesh)));
		
		System.exit(1);

	}

}

package org.boussinesq;

//import java.util.Arrays;

public class SongDimensionless {

	public static double[] beqSongDimensionless(double[] xi, double xi0,
			double[] a) {

		// initialization of double in java is equal to zero
		double[] solution = new double[xi.length];
		/*
		 * double[] temp = new double[xi.length];
		 * 
		 * for (int i = 0; i < xi.length; i++) {
		 * 
		 * temp[i] = Math.max(0, 1 - xi[i] / xi0);
		 * 
		 * }
		 */

		for (int i = 0; i < a.length; i++) {

			for (int j = 0; j < xi.length; j++) {

				solution[j] = solution[j] + a[i]
						* Math.pow(Math.max(0, 1 - xi[j] / xi0), i);

			}

		}

		// System.out.println("Song dimensionless2\n" +
		// Arrays.toString(solution));

		return solution;

	}

	public static void main(String[] args) {

	}

}

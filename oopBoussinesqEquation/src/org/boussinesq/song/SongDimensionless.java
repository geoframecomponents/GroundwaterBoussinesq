package org.boussinesq.song;

//import java.util.Arrays;

public class SongDimensionless {

	public static double[] beqSongDimensionless(double[] xi, double xi0,
			double[] a) {

		// initialization of double in java is equal to zero
		double[] solution = new double[xi.length];

		int endForLoop1 = a.length;
		int endForLoop2 = xi.length;
		
		for (int i = 0; i < endForLoop1; i++) {

			for (int j = 0; j < endForLoop2; j++) {

				solution[j] = solution[j] + a[i]
						* Math.pow(Math.max(0, 1 - xi[j] / xi0), i);

			}

		}

		return solution;

	}

	public static void main(String[] args) {

	}

}

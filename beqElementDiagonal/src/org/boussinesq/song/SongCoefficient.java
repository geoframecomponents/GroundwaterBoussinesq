package org.boussinesq.song;

import java.util.Arrays;

public class SongCoefficient {

	public static double[] CoefficientSongSolution(int nmax, double lambda) {

		double[] a = new double[nmax];

		// System.out.println(lambda);

		a[0] = 0;
		a[1] = 0.25;
		a[2] = 0.0625 * (2 * lambda - 1);

		for (int i = 3; i < nmax; i++) {

			a[i] = (2 * lambda + 1 - i) / (Math.pow(i, 2)) * a[i - 1];

			/*
			 * if (i==4){ System.out.println(c); }
			 */

			for (int k = 2; k < i; k++) {

				a[i] += -2 * (i + 1) * Math.pow(i, -1) * a[k] * a[i + 1 - k];

				/*
				 * if (i == 4){ System.out.println(a[i]);
				 * System.out.println(a[k]); System.out.println(a[i+1-k]); }
				 */

			}

		}

		System.out.println(Arrays.toString(a));

		return a;

	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

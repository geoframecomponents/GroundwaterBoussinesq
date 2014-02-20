package org.boussinesq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * The Class TripletToRowCompressedForm.
 */
public class TripletToRowCompressedForm {

	/**
	 * Compute Mp.
	 * 
	 * @desc
	 * 
	 * @param Mj
	 *            the mj
	 * @param Np
	 *            the np
	 * @return the int[]
	 */
	public static int[] computeMp(int[] pMi, int[] pMj, double[] pMl, int Np) {

		int[] Mp = new int[Np + 1];
		int[] Mi = new int[pMi.length];
		int[] Mj = new int[pMj.length];
		double[] Ml = new double[pMl.length];

		int index = 0;
		int min = 0;
		boolean intoLoop = false;
		// Mp[index] = 0;
		//index++;

		ArrayList<Integer> alMi = new ArrayList<Integer>();
		ArrayList<Integer> alMj = new ArrayList<Integer>();
		ArrayList<Double> alMl = new ArrayList<Double>();

		for (int i = 0; i < Np; i++) {

			for (int j = 0; j < pMi.length; j++) {

				if (pMi[j] == i) {

					alMi.add(pMi[j]);
					alMl.add(pMl[j]);
					alMj.add(pMj[j]);

					intoLoop = true;

				}
			}

			if (intoLoop) {

				Mp[i] = index;
				intoLoop = false;

				//System.out.println(index);

				while (alMi.size() > 0) {

					min = alMj.indexOf(Collections.min(alMj));
					Mi[index] = alMi.get(min);
					Mj[index] = alMj.get(min);
					Ml[index] = alMl.get(min);

					index++;
					
					alMi.remove(min);
					alMj.remove(min);
					alMl.remove(min);

				}

//				alMi.clear();
//				alMj.clear();
//				alMl.clear();

			}

		}

		Mp[Np] = pMi.length;
		
		System.out.println(Arrays.toString(Mj));
		System.out.println(Arrays.toString(Ml));
		
		return Mp;

	}

	/**
	 * The main method.
	 * 
	 * @desc This main is only a test method.
	 * 
	 * @param args
	 *            the arguments
	 */
	public static void main(String[] args) {

		int[] Mi = new int[] { 0, 0, 1, 1, 1, 2, 2, 0, 4, 3, 1, 5, 4, 3, 4, 2,
				4, 5, 6, 5, 3, 5, 6, 8, 6, 4, 6, 7, 7, 9, 8, 5, 8, 9, 9, 10, 6,
				9, 10, 10, 11, 11, 7, 10, 11, 12, 9, 11, 12, 13, 12, 10, 11,
				13, 14, 15, 13, 14, 14, 15, 15, 14 };

		int[] Mj = new int[] { 0, 1, 0, 1, 2, 1, 2, 4, 0, 3, 5, 1, 3, 4, 4, 6,
				5, 4, 2, 5, 8, 6, 5, 3, 6, 9, 7, 7, 6, 4, 8, 10, 9, 9, 8, 5,
				11, 10, 9, 10, 6, 11, 12, 11, 10, 7, 13, 12, 11, 9, 12, 14, 15,
				13, 10, 11, 14, 13, 14, 15, 14, 15 };

		double[] Ml = new double[] { -1, 23, 23, -1, 24, 24, -1, 5, 5, -1, 6,
				6, 27, 27, -1, 7, 28, 28, 7, -1, 9, 29, 29, 9, -1, 10, 30, -1,
				30, 10, -1, 11, 33, -1, 33, 11, 12, 34, 34, -1, 12, -1, 13, 35,
				35, 13, 15, 36, 36, 15, -1, 16, 17, -1, 16, 17, 39, 39, -1, -1,
				40, 40 };

		// number of polygons in the mesh
		int Np = 16;

		// compute Mp from a triplet form
		int[] Mp = computeMp(Mi, Mj, Ml, Np);
		System.out.println(Arrays.toString(Mp));
	}
}

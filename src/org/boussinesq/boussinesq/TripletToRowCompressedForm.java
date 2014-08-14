package org.boussinesq.boussinesq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

/**
 * The Class TripletToRowCompressedForm.
 */
public class TripletToRowCompressedForm {

	
	static int[] Mp;
	static int[] Mj;
	static double[] Ml;
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
	public static void computeMp(int[] pMi, int[] pMj, double[] pMl, int Np) {

		Mp = new int[Np + 1];
		int[] Mi = new int[pMi.length];
		Mj = new int[pMj.length];
		Ml = new double[pMl.length];

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

				if (pMi[j] == (i+1)) {

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
					Mi[index] = alMi.get(min)-1;
					Mj[index] = alMj.get(min)-1;
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
		
		//return Mp;

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

//		int[] Mi = new int[] {2, 1, 4, 1, 3, 2, 5, 2, 5, 4, 6, 4, 7, 5, 7, 6, 8, 7, 1, 2, 3, 4, 5, 6, 7,8};
//
//		int[] Mj = new int[] {1, 2, 1, 4, 2, 3, 2, 5, 4, 5, 4, 6, 5, 7, 6, 7, 7, 8, 1, 2, 3, 4, 5, 6, 7,8};
//
//		double[] Ml = new double[] {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, -1, -1, -1, -1, -1, -1, -1,-1};

		
		int[] Mi = new int[] {1, 2, 4, 1, 2, 3, 5, 2, 3, 1, 4, 5, 6, 2, 4, 5, 7, 4, 6, 7, 5, 6, 7, 8, 7, 8};

		int[] Mj = new int[] {1, 2, 4, 1, 2, 3, 5, 2, 3, 1, 4, 5, 6, 2, 4, 5, 7, 4, 6, 7, 5, 6, 7, 8, 7, 8};

		double[] Ml = new double[] {-1, 1, 2, 1, -1, 3, 4, 3, -1, 2, -1, 5, 6, 4, 5, -1, 7, 6, -1, 8, 7, 8, -1, 9, 9,-1};
		
		// number of polygons in the mesh
		int Np = 8;

		// compute Mp from a triplet form
		computeMp(Mi, Mj, Ml, Np);
		System.out.println(Arrays.toString(Mp));
	}
}

package org.boussinesq.boussinesq;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public class DeleteRowColumnNullDiagonalEntry {

	LinkedList<Integer> alNewMj, alNewMp, alNewMp2, alNullDiag, alNullDiag2;
	ArrayList<Integer> indices;
	LinkedList<Double> alT;

	static double[] newMatT;
	static int[] newMj, newMp, newDeleteIndex, newIndices;

	public DeleteRowColumnNullDiagonalEntry() {

		alNewMj = new LinkedList<Integer>();
		alNewMp = new LinkedList<Integer>();
		alNewMp2 = new LinkedList<Integer>();
		alNullDiag = new LinkedList<Integer>();
		alNullDiag2 = new LinkedList<Integer>();
		indices = new ArrayList<Integer>();
		alT = new LinkedList<Double>();

	}

	public static double[] convertDoubles(LinkedList<Double> doubles) {
		double[] ret = new double[doubles.size()];
		Iterator<Double> iterator = doubles.iterator();
		int i = 0;
		while (iterator.hasNext()) {
			ret[i] = iterator.next().doubleValue();
			i++;
		}
		return ret;
	}

	public static int[] convertInteger(List<Integer> integers) {
		int[] ret = new int[integers.size()];
		Iterator<Integer> iterator = integers.iterator();
		int i = 0;
		while (iterator.hasNext()) {
			ret[i] = iterator.next().intValue();
			i++;
		}
		return ret;
	}

	public void computationalDomain(int size, double[] matT, int[] indexDiag,
			int[] Mp, int[] Mj) {

		for (int i = 0; i < size; i++) {

			// alNewMp.add(Mp[i]);

			if (matT[indexDiag[i]] == 0) {
				alNullDiag.add(i);
				alNullDiag2.add(i);
			}

		}

		for (int i = 0; i < matT.length; i++) {

			alNewMj.add(Mj[i]);
			alT.add(matT[i]);

		}

		int cont = 0;

		while (alNullDiag.size() != 0) {

			for (int i = 0; i < size; i++) {

				for (int j = Mp[i]; j < Mp[i + 1]; j++) {

					if (i == alNullDiag.get(cont)) {

						indices.add(j);

					} else if (Mj[j] == alNullDiag.get(cont)) {

						indices.add(j);

					}

				}

			}

			alNullDiag.remove(cont);

		}

		Collections.sort(indices);

		newIndices = new int[indices.size()];
		newIndices = convertInteger(indices);

		cont = 0;
		int contatore = indices.size() - 1;

		int prova;

		while (contatore >= 0) {

			prova = indices.get(contatore);
			alT.remove(prova);
			alNewMj.remove(prova);
			contatore = contatore - 1;

		}

		newMj = new int[alNewMj.size()];
		newMj = convertInteger(alNewMj);
		newMp = new int[size - alNullDiag2.size()];

		cont = 0;

		newMp = new int[Mp.length];
		System.arraycopy(Mp, 0, newMp, 0, Mp.length);

		cont = 0;

		contatore = 0;

		while (contatore < alNullDiag2.size()) {

			for (int i = 0; i < size; i++) {

				if (i == alNullDiag2.get(contatore)) {

					for (int j = Mp[i]; j < Mp[i + 1]; j++) {

						cont++;

					}

				}

				newMp[i] = newMp[i] - cont;

			}

			newMp[newMp.length - 1] = newMp[newMp.length - 1] - cont;

			cont = 0;
			contatore++;

		}

		cont = 0;
		contatore = 0;

		while (contatore < alNullDiag2.size()) {

			for (int i = 0; i < size; i++) {

				newMp[i] = newMp[i] - cont;

				if (i != alNullDiag2.get(contatore)) {

					for (int j = Mp[i]; j < Mp[i + 1]; j++) {

						if (Mj[j] == alNullDiag2.get(contatore)) {

							cont++;

						}

					}

				}

			}

			newMp[newMp.length - 1] = newMp[newMp.length - 1] - cont;

			cont = 0;
			contatore++;

		}

		for (int i = 0; i < Mp.length; i++) {

			if (cont < alNullDiag2.size() && i != alNullDiag2.get(cont)) {

				// System.out.println(alNullDiag2.get(cont));

				alNewMp.add(newMp[i]);

			} else if (cont < alNullDiag2.size() && i == alNullDiag2.get(cont)) {

				cont++;

			} else {

				alNewMp.add(newMp[i]);

			}

		}

		cont = 0;

		
//		contatore = ;
		
		System.out.println(Arrays.toString(newMj));
		
		contatore = 0;
		
		while (contatore < alNullDiag2.size()) {

//			System.out.println("I'm here");
			
			for (int i = 0; i < newMj.length; i++){
				
				if (newMj[i] >= alNullDiag2.get(contatore)){
					
					newMj[i] = newMj[i] - 1;
					
				}
				
			}

			contatore++;
			
		}

		newMp = convertInteger(alNewMp);

		System.out.println(Arrays.toString(newMp));
		System.out.println(Arrays.toString(newMj));

		newMatT = convertDoubles(alT);
		newMj = convertInteger(alNewMj);

	}

	public static void main(String[] args) {

		double[] matT = { 72, -36, -36, -36, 93.6, -28.8, -28.8, -28.8, 28.8,
				0, 0, 0, 0, 0, -36, 61.2, -25.2, -28.8, -25.2, 93.6, 0, -39.6,
				0, 0, 0, 0, 0, 0, 0, 36, -36, 0, -36, 79.2, -43.2, -39.6, 79.2,
				-39.6, 0, -39.6, 39.6, 0, 0, 0, 43.2, -43.2, -43.2, -43.2, 86.4 };

		int[] Mp = { 0, 3, 7, 11, 14, 17, 22, 27, 32, 35, 38, 42, 46, 49 };
		int[] Mj = { 0, 1, 4, 0, 1, 2, 5, 1, 2, 3, 6, 2, 3, 7, 0, 4, 5, 1, 4,
				5, 6, 9, 2, 5, 6, 7, 10, 3, 6, 7, 8, 11, 7, 8, 12, 5, 9, 10, 6,
				9, 10, 11, 7, 10, 11, 12, 8, 11, 12 };

		int size = 13;

		int[] indexDiag = new int[size];

		/* for-loop to analyze the matrix cell by cell */
		for (int i = 0; i < size; i++) {
			/*
			 * nested for-loop to analyze diagonal entries, which are identified
			 * by a negative number
			 */
			for (int j = Mp[i]; j < Mp[i + 1]; j++) {

				if (Mj[j] == i) {
					indexDiag[i] = j;
				}

			}
		}

		DeleteRowColumnNullDiagonalEntry prova = new DeleteRowColumnNullDiagonalEntry();

		prova.computationalDomain(size, matT, indexDiag, Mp, Mj);

		System.out.println(Arrays.toString(newMatT));

		System.exit(0);

	}

}

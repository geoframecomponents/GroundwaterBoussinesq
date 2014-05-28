package org.boussinesq.RowCompressedForm;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

public class DeleteRowColumnNullDiagonalEntry {

	public static ArrayList<Integer> alNewMj, alNullDiag;
	private ArrayList<Integer> indices;
	private ArrayList<Double> alT;

	public double[] newMatT;
	public static int[] newMj, newMp, newIndices;

	public double[] convertDoubles(List<Double> doubles) {
		double[] ret = new double[doubles.size()];
		Iterator<Double> iterator = doubles.iterator();
		int i = 0;
		while (iterator.hasNext()) {
			ret[i] = iterator.next().doubleValue();
			i++;
		}
		return ret;
	}

	public int[] convertInteger(List<Integer> integers) {
		int[] ret = new int[integers.size()];
		Iterator<Integer> iterator = integers.iterator();
		int i = 0;
		while (iterator.hasNext()) {
			ret[i] = iterator.next().intValue();
			i++;
		}
		return ret;
	}

	@SuppressWarnings("unchecked")
	public void searchRowColumnNullDiagonalEntries(int size, double[] matT,
			int[] indexDiag, int[] Mp, int[] Mj) {

		// search the null diagonal entries and store the indices in two linked
		// lists

		for (int i = 0; i < size; i++) {

			if (matT[indexDiag[i]] == 0) {
				alNullDiag.add(i);
				// alNullDiag2.add(i);
			}

		}

		int contatore = 0;

		while (contatore < alNullDiag.size()) {

			for (int i = 0; i < size; i++) {

				for (int j = Mp[i]; j < Mp[i + 1]; j++) {

					if (i == alNullDiag.get(contatore)) {

						indices.add(j);

					} else if (Mj[j] == alNullDiag.get(contatore)) {

						indices.add(j);

					}

				}

			}

			// alNullDiag.remove(cont);
			contatore++;

		}

		
		@SuppressWarnings("rawtypes")
		HashSet hs = new HashSet();
		hs.addAll(indices);
		indices.clear();
		indices.addAll(hs);

		// sort indices in ascending order

		Collections.sort(indices);

		newIndices = new int[indices.size()];
		newIndices = convertInteger(indices);
		
	}
	
	public void computeNewarrayMp(int size, int[] Mp){
		
//		newMp = new int[size - alNullDiag.size()];
		

		int contatore = 1;

		newMp = new int[Mp.length - alNullDiag.size()];

		newMp[0] = 0;

		int endForLoop = newMj.length - 1;
		
		// compute Mp from Mj

		for (int i = 0; i < endForLoop; i++) {

			if (newMj[i] > newMj[i + 1]) {

				newMp[contatore] = i + 1;
				contatore++;
				// cont = cont + 1;
			}

		}

		// newMp = convertInteger(alNewMp);

		newMp[newMp.length - 1] = newMj.length;
		
	}

	public void computeNewArrayMj(int[] Mj) {

		int endForLoop = Mj.length;
		
		for (int i = 0; i < endForLoop; i++) {

			alNewMj.add(Mj[i]);

		}
		
		int contatore = indices.size() - 1;

		int prova;

		// remove rows and columns from T and Mj linked lists

		while (contatore >= 0) {

			prova = indices.get(contatore);
			alNewMj.remove(prova);
			contatore = contatore - 1;

		}
		
		newMj = new int[alNewMj.size()];
		newMj = convertInteger(alNewMj);
		
		int cont = 0;

		contatore = 0;

		endForLoop = newMj.length;
		
		// compute new Mj

		for (int i = 0; i < endForLoop; i++) {

			while (contatore < alNullDiag.size()) {

				if (newMj[i] > alNullDiag.get(contatore))
					cont++;

				contatore++;

			}

			newMj[i] = newMj[i] - cont;
			cont = 0;
			contatore = 0;

		}

	}

	public double[] computeNewArrayT(double[] matT) {

		alT = new ArrayList<Double>();
		
		int endForLoop = matT.length;
		
		for (int i = 0; i < endForLoop; i++) {

			alT.add(matT[i]);

		}
		
		
		int contatore = indices.size() - 1;

		int prova;

		// remove rows and columns from T and Mj linked lists

		while (contatore >= 0) {

			prova = indices.get(contatore);
			alT.remove(prova);
			contatore = contatore - 1;

		}
		
		newMatT = convertDoubles(alT);

		matT = new double[newMatT.length];

		System.arraycopy(newMatT, 0, matT, 0, newMatT.length);
		
		return matT;

	}

	public void computationalDomain(int size, double[] matT,
			int[] indexDiag, int[] Mp, int[] Mj) {

		alNewMj = new ArrayList<Integer>();
		alNullDiag = new ArrayList<Integer>();
		// alNullDiag2 = new LinkedList<Integer>();
		indices = new ArrayList<Integer>();
		alT = new ArrayList<Double>();

		searchRowColumnNullDiagonalEntries(size, matT, indexDiag, Mp, Mj);

	}

	public double[] computeCellsArray(double[] eta) {

		LinkedList<Double> alEta = new LinkedList<Double>();

		int cont = 0;
		
		int endForLoop = eta.length;

		for (int i = 0; i < endForLoop; i++) {

			if (cont >= alNullDiag.size()) {

				alEta.add(eta[i]);

			} else if (i != alNullDiag.get(cont)) {

				alEta.add(eta[i]);

			} else {

				cont++;

			}

		}

		double[] newEta = new double[alEta.size()];
		newEta = convertDoubles(alEta);

		return newEta;

	}

	public double[] addRemovedElements(double[] x, double[] element) {

		double[] y = new double[x.length + alNullDiag.size()];

		LinkedList<Double> llX = new LinkedList<Double>();

		int endForLoop = x.length;
		
		for (int i = 0; i < endForLoop; i++) {

			llX.add(x[i]);

		}

		int cont = 0;

		int endWhileLoop = DeleteRowColumnNullDiagonalEntry.alNullDiag.size();
		
		while (cont < endWhileLoop) {

			llX.add(DeleteRowColumnNullDiagonalEntry.alNullDiag.get(cont),
					element[cont]);

			cont++;

		}

		y = convertDoubles(llX);

		return y;

	}

	public static void main(String[] args) {

		// prova 1

		int size = 13;

		double[] matT = { 72, -36, -36, -36, 93.6, -28.8, -28.8, -28.8, 28.8,
				0, 0, 0, 0, 0, -36, 61.2, -25.2, -28.8, -25.2, 93.6, 0, -39.6,
				0, 0, 0, 0, 0, 0, 0, 36, -36, 0, -36, 79.2, -43.2, -39.6, 79.2,
				-39.6, 0, -39.6, 39.6, 0, 0, 0, 43.2, -43.2, -43.2, -43.2, 86.4 };

		int[] Mp = { 0, 3, 7, 11, 14, 17, 22, 27, 32, 35, 38, 42, 46, 49 };
		int[] Mj = { 0, 1, 4, 0, 1, 2, 5, 1, 2, 3, 6, 2, 3, 7, 0, 4, 5, 1, 4,
				5, 6, 9, 2, 5, 6, 7, 10, 3, 6, 7, 8, 11, 7, 8, 12, 5, 9, 10, 6,
				9, 10, 11, 7, 10, 11, 12, 8, 11, 12 };

		double[] eta = { 1, 0.8, 0, 0, 0.7, 0, 0, 0, 1, 1.1, 0, 0, 1.2 };
		double[] newEta = null;

		// prova 2

		// int size = 5;
		//
		// double[] matT = { 36, -36, -36, 36, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
		//
		// int[] Mp = { 0, 2, 5, 8, 11, 13 };
		// int[] Mj = { 0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4 };
		//
		// double[] eta = { 1, 0, 0, 0, 0 };
		// double[] newEta = null;

		// prova 3

		// int size = 20;
		//
		// double[] matT = { 86.4, -43.2, -43.2, -43.2, 122.4, -39.6, -39.6,
		// -39.6, 97.2, -28.8, -28.8, -28.8, 28.8, 0, 0, 0, 0, 0, -43.2,
		// 122.4, -39.6, -39.6, -39.6, -39.6, 136.8, -28.8, -28.8,
		// -28.8, -28.8, 100.8, -21.6, -21.6, 0, -21.6, 50.4, 0, -28.8, 0, 0,
		// 36, -36, -39.6, 39.6, 0, 0, -28.8, 0, 50.4, -21.6, 0, -21.6,
		// -21.6, 93.6, -28.8, -21.6, -28.8, -28.8, 122.4, -36, -28.8,
		// -36, -36, 108, -36, 0, 0, 0, 0, 0, 0, 0, -21.6, 0, 50.4, -28.8,
		// -28.8, -28.8, 93.6, -36, -36, -36, 72 };
		//
		// int[] Mp = {0, 3, 7, 11, 15, 18, 22, 27, 32, 37, 41 ,45, 50, 55, 60,
		// 64, 67, 71, 75, 79, 82};
		//
		// int[] Mj = { 0, 1, 5, 0, 1, 2, 6, 1, 2, 3, 7, 2, 3, 4, 8, 3, 4, 9, 0,
		// 5, 6, 10, 1, 5, 6, 7, 11, 2, 6, 7, 8, 12, 3, 7, 8, 9, 13, 4, 8,
		// 9, 14, 5, 10, 11, 15, 6, 10, 11, 12, 16, 7, 11, 12, 13, 17, 8,
		// 12, 13, 14, 18, 9, 13, 14, 19, 10, 15, 16, 11, 15, 16, 17, 12,
		// 16, 17, 18, 13, 17, 18, 19, 14, 18, 19 };
		//
		// double[] eta = { 1, 0.8, 0, 0, 0.7, 0, 0, 0, 1, 1.1, 0, 0, 1.2 };
		// double[] newEta = null;

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

		// System.out.println(Arrays.toString(indexDiag));

		DeleteRowColumnNullDiagonalEntry prova = new DeleteRowColumnNullDiagonalEntry();

		prova.computationalDomain(size, matT, indexDiag, Mp, Mj);

		prova.computeNewArrayMj(Mj);
		prova.computeNewarrayMp(size, Mp);
		prova.newMatT = prova.computeNewArrayT(matT);
		
		System.out.println("Mp: " + Arrays.toString(newMp));
		System.out.println("Mj: " + Arrays.toString(newMj));
		System.out.println("T : " + Arrays.toString(prova.newMatT));

		newEta = prova.computeCellsArray(eta);

		System.out.println("eta: " + Arrays.toString(newEta));

		System.exit(0);

	}

}

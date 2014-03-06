package org.francescoS.usefulClasses;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Scanner;

public class FileRead {

	int row = 0;
	int col = 0;
	double[] array = null;
	double[][] matrix = null;

	public double[][] readDoubleMatrix(String filePath) throws IOException {

		readMatrixDimension(filePath);
		readMatrix(filePath);

		return matrix;

	}

	public double[] readDoubleArray(String filePath)
			throws FileNotFoundException {

		readMatrixDimension(filePath);

		if (col == 1) {

			readColumnArray(filePath);

		} else if (row == 1) {

			readRowArray(filePath);
		}

		return array;

	}

	public void readMatrixDimension(String filePath)
			throws FileNotFoundException {

		Scanner inputFile = null;
		inputFile = new Scanner(new File(filePath));

		// determine the number of rows/columns

		while (inputFile.hasNextLine()) {

			row++;

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inputFile.nextLine());

			while (matCol.hasNextDouble())

			{

				col++;
				matCol.nextDouble();

			}

		}

		col = col / row;

		System.out.println("NUmber of rows: " + row);
		System.out.println("Number of columns: " + col);
		inputFile.close();

	}

	public void readMatrix(String filePath) throws FileNotFoundException {

		int r = 0;
		int c = 0;
		matrix = new double[row][col];

		Scanner inFile = null;

		inFile = new Scanner(new File(filePath));

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextDouble())

			{

				matrix[r][c] = matCol.nextDouble();
				c++;

			}

			c = 0;
			r++;

		}

		inFile.close();

	}

	public void readRowArray(String filePath) throws FileNotFoundException {

		int c = 0;
		array = new double[col];

		Scanner inFile = null;

		inFile = new Scanner(new File(filePath));

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextDouble())

			{

				array[c] = matCol.nextDouble();
				c++;

			}

			c = 0;

		}

		inFile.close();

	}

	public void readColumnArray(String filePath) throws FileNotFoundException {

		int r = 0;
		array = new double[row];

		Scanner inFile = null;

		inFile = new Scanner(new File(filePath));

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextDouble())

			{

				array[r] = matCol.nextDouble();

			}

			r++;

		}

		inFile.close();

	}

	public static void main(String[] args) throws IOException {

		// double[][] matrix;
		double[] array;

		String path = "/home/francesco/b.txt";
		FileRead file = new FileRead();

		// matrix = file.readDoubleMatrix(path);
		array = file.readDoubleArray(path);

		/*
		 * for (int r = 0; r < matrix.length; r++) {
		 * 
		 * for (int c = 0; c < matrix[0].length; c++) {
		 * 
		 * System.out.print(matrix[r][c] + " ");
		 * 
		 * }
		 * 
		 * System.out.println("");
		 * 
		 * }
		 */

		System.out.println(Arrays.toString(array));

		System.exit(0);

	}

}

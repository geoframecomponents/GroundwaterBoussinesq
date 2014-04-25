package org.francescoS.usefulClasses;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
//import java.util.Arrays;
import java.util.Scanner;

// TODO: Auto-generated Javadoc
/**
 * The Class FileRead.
 */
public class FileRead {

	/** The row. */
	int row = 0;

	/** The col. */
	int col = 0;
	
	int[] intArray = null;

	/** The array. */
	double[] doubleArray = null;

	/** The matrix. */
	double[][] doubleMatrix = null;

	/**
	 * Read double matrix.
	 * 
	 * @param filePath
	 *            the file path
	 * @return the double[][]
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public double[][] readDoubleMatrix(File filePath) throws IOException {

		readDoubleMatrixDimension(filePath);
		loadDoubleMatrix(filePath);

		return doubleMatrix;

	}

	/**
	 * Read double array.
	 * 
	 * @param filePath
	 *            the file path
	 * @return the double[]
	 * @throws FileNotFoundException
	 *             the file not found exception
	 */
	public double[] readDoubleArray(File filePath)
			throws FileNotFoundException {

		readDoubleMatrixDimension(filePath);

		if (col == 1) {

			readColumnDoubleArray(filePath);

		} else if (row == 1) {

			readRowDoubleArray(filePath);
		}

		return doubleArray;

	}
	
	public int[] readIntArray(File filePath) throws FileNotFoundException{
		
		readIntMatrixDimension(filePath);
		
		if (col == 1) {

			readColumnIntArray(filePath);

		} else if (row == 1) {

			readRowIntArray(filePath);
		}
		
		return intArray;
	}

	/**
	 * Read matrix dimension.
	 * 
	 * @param filePath
	 *            the file path
	 * @throws FileNotFoundException
	 *             the file not found exception
	 */
	public void readDoubleMatrixDimension(File filePath)
			throws FileNotFoundException {

		Scanner inputFile = null;
		inputFile = new Scanner(filePath);

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

//		System.out.println("NUmber of rows: " + row);
//		System.out.println("Number of columns: " + col);
		inputFile.close();

	}
	
	public void readIntMatrixDimension(File filePath)
			throws FileNotFoundException {

		Scanner inputFile = null;
		inputFile = new Scanner(filePath);

		// determine the number of rows/columns

		while (inputFile.hasNextLine()) {

			row++;

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inputFile.nextLine());

			while (matCol.hasNextInt())

			{

				col++;
				matCol.nextInt();

			}

		}

		col = col / row;

//		System.out.println("NUmber of rows of integers: " + row);
//		System.out.println("Number of columns of integers: " + col);
		inputFile.close();

	}

	/**
	 * Read matrix.
	 * 
	 * @param filePath
	 *            the file path
	 * @throws FileNotFoundException
	 *             the file not found exception
	 */
	public void loadDoubleMatrix(File filePath) throws FileNotFoundException {

		int r = 0;
		int c = 0;
		doubleMatrix = new double[row][col];

		Scanner inFile = null;

		inFile = new Scanner(filePath);

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextDouble())

			{

				doubleMatrix[r][c] = matCol.nextDouble();
				c++;

			}

			c = 0;
			r++;

		}

		inFile.close();

	}

	/**
	 * Read row array.
	 * 
	 * @param filePath
	 *            the file path
	 * @throws FileNotFoundException
	 *             the file not found exception
	 */
	public void readRowDoubleArray(File filePath) throws FileNotFoundException {

		int c = 0;
		doubleArray = new double[col];

		Scanner inFile = null;

		inFile = new Scanner(filePath);

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextDouble())

			{

				doubleArray[c] = matCol.nextDouble();
				c++;

			}

			c = 0;

		}

		inFile.close();

	}

	/**
	 * Read column array.
	 * 
	 * @param filePath
	 *            the file path
	 * @throws FileNotFoundException
	 *             the file not found exception
	 */
	public void readColumnDoubleArray(File filePath) throws FileNotFoundException {

		int r = 0;
		doubleArray = new double[row];

		Scanner inFile = null;

		inFile = new Scanner(filePath);

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextDouble())

			{

				doubleArray[r] = matCol.nextDouble();

			}

			r++;

		}

		inFile.close();

	}

	
	
	public void readRowIntArray(File filePath) throws FileNotFoundException {

		int c = 0;
		intArray = new int[col];

		Scanner inFile = null;

		inFile = new Scanner(filePath);

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextInt())

			{

				intArray[c] = matCol.nextInt();
				c++;

			}

			c = 0;

		}

		inFile.close();

	}

	/**
	 * Read column array.
	 * 
	 * @param filePath
	 *            the file path
	 * @throws FileNotFoundException
	 *             the file not found exception
	 */
	public void readColumnIntArray(File filePath) throws FileNotFoundException {

		int r = 0;
		intArray = new int[row];

		Scanner inFile = null;

		inFile = new Scanner(filePath);

		while (inFile.hasNextLine()) {

			@SuppressWarnings("resource")
			Scanner matCol = new Scanner(inFile.nextLine());

			while (matCol.hasNextInt())

			{

				intArray[r] = matCol.nextInt();

			}

			r++;

		}

		inFile.close();

	}
	
	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args) throws IOException {

		// double[][] matrix;
//		double[] array;

//		String path = "/home/francesco/b.txt";
//		FileRead file = new FileRead();

		// matrix = file.readDoubleMatrix(path);
//		array = file.readDoubleArray(path);

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


		System.exit(0);

	}

}

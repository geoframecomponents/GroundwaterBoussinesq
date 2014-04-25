/*
 * 
 */
package org.francescoS.usefulClasses;
/**
 * A class to write your solution on text file
 * 
 * @desc	In this code is implemented a robust manner to write computational
 * 			data on text file. This code verifies if the path of text file is
 * 			correct and if overwrite the file or not (it takes in input the
 * 			boolean about overwriting)
 * 
 * @author	F. Serafin, 2014
 * Copyright GPL v. 3 (http://www.gnu.org/licenses/gpl.html)
 * */
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

/**
 * The Class FileWrite.
 */
public class FileWrite {

	/** The errestat. */
	static PrintWriter writeData;
	
	/** The Rstatfile. */
	static FileWriter outputFile;

	
	
	
	
	
	
	
	
	
	/**
	 * Open txt file.
	 *
	 * @param path the path
	 * @param overwrite the overwrite
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void openTxtFile(String name, File dirPath, boolean overwrite)
			throws IOException {


//		File file = dirPath + new File(name);
//		String path = dirPath.getPath() + name;
		
		
		try {

			File file = new File(dirPath, name);

			if (overwrite) {

				file.createNewFile();
				
//				System.out.println("Overwrite option selected");
//				System.out.println("File name is " + file.getName());

				outputFile = new FileWriter(file, false); // overwrites file
				writeData = new PrintWriter(outputFile);

			} else if (!file.exists()) {
				
				file.createNewFile();

//				System.out.println("Overwrite option not selected");
//				System.out.println("File name is " + file.getName());

				outputFile = new FileWriter(file, true); // no overwrites file
				writeData = new PrintWriter(outputFile);

			} else {

				System.err.println("Path\n" + file.getAbsolutePath()
						+ "\nexists and overwrite option not selected\n");
				System.err.println("\nEnd of code");
				System.exit(0);

			}

			
			
		} catch (IOException e) {

			File file = new File(dirPath, name);

			System.err.println("Path\n" + file.getAbsolutePath());
			// Catch exception if any
			System.err.println("Error: " + e.getMessage());
			throw new RuntimeException(e);

		}

		
		
	}

	
	
	
	
	
	
	
	
	
	/**
	 * Write double1 column.
	 *
	 * @param data the data
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void writeOneDoubleColumn(double[] data) throws IOException {

		for (int j = 0; j < data.length; j++) {
			writeData.println(data[j]);
		}
		
	}
	
	
	
	
	
	
	
	
	
	
	public static void writeFourDoubleColumn(double[] data1, double[] data2, double[] data3, double[] data4){
		
		for (int j = 0; j < data1.length; j++){
			
			writeData.print(data1[j] + "\t");
			writeData.print(data2[j] + "\t");
			writeData.print(data3[j] + "\t");
			writeData.println(data4[j]);
			
		}
		
	}
	
	
	
	
	
	
	
	
	
	
	public static void writeStringDoubleString(String description, double value, String unitMeasure){
		
		writeData.println(description + ": " + value + " " + unitMeasure);
		
	}
	
	
	
	
	
	
	
	
	
	
	public static void writeStringIntString(String description, int value, String unitMeasure){
		
		writeData.println(description + ": " + value + " " + unitMeasure);
		
	}
	
	
	
	
	
	
	
	
	
	
	public static void writeFourStringColumn(String string1, String string2, String string3, String string4){
		
		writeData.print(string1 + "\t");
		writeData.print(string2 + "\t");
		writeData.print(string3 + "\t");
		writeData.println(string4);
		
	}
	
	
	
	
	
	
	
	
	
	
	public static void closeTxtFile() throws IOException{

		writeData.println();
		System.out.println();
		outputFile.close();

	}
	
	
	
	
	
	
	
	
	
	
	public static File makeDirectory(String dirName){
		
		String defaultPath = new File(dirName).getAbsolutePath(); 
		
		File newPath = new File(defaultPath);
		newPath.mkdir();
		
		TextIO.putln("Directory of output files: " + newPath);
		
		return newPath;

		
	}
}

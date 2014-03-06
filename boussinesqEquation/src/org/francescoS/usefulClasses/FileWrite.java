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

// TODO: Auto-generated Javadoc
/**
 * The Class FileWrite.
 */
public class FileWrite {

	/** The errestat. */
	static PrintWriter errestat;
	
	/** The Rstatfile. */
	static FileWriter Rstatfile;

	/**
	 * Open txt file.
	 *
	 * @param path the path
	 * @param overwrite the overwrite
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void openTxtFile(String path, boolean overwrite)
			throws IOException {

		try {

			File file = new File(path);

			if (overwrite) {

				file.createNewFile();
				
				System.out.println("Overwrite option selected");
				System.out.println("File name is " + file.getName());

				Rstatfile = new FileWriter(file, false); // overwrites file
				errestat = new PrintWriter(Rstatfile);

			} else if (!file.exists()) {
				
				file.createNewFile();

				System.out.println("Overwrite option not selected");
				System.out.println("File name is " + file.getName());

				Rstatfile = new FileWriter(file, true); // no overwrites file
				errestat = new PrintWriter(Rstatfile);

			} else {

				System.err.println("Path\n" + file.getAbsolutePath()
						+ "\nexists and overwrite option not selected\n");
				System.err.println("\nEnd of code");
				System.exit(0);

			}

		} catch (IOException e) {

			File file = new File(path);

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
	public static void writeDouble1Column(double[] data) throws IOException {

		for (int j = 0; j < data.length; j++) {

			errestat.println(data[j]);

		}

		errestat.println();
		System.out.println();
		Rstatfile.close();

	}
}

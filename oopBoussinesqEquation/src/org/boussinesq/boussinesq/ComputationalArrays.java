package org.boussinesq.boussinesq;

import org.boussinesq.RowCompressedForm.DeleteRowColumnNullDiagonalEntry;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;

public class ComputationalArrays {

	public static int numberOfPolygon;
	public static int[] arrayMp, arrayMj, indexDiagonal;
	public static double[] bedRockElevation, porosity, planarArea, source,
			coeffC, coeffM, etaDirichlet, arrayEta;
	
	RCIndexDiagonalElement newIndexDiag;
	public static boolean matrixReduction = false;
	
	public double[] DefineDomainProperties(double[] arrayT, double[] eta, int[] indexDiag, DeleteRowColumnNullDiagonalEntry cleanMatT){
		
		newIndexDiag = new RCIndexDiagonalElement();
		
		double[] matT;
		
		if (ComputeT.unlockDeleteRowColumn) {
			
			System.out.println("Remove row column");

			matT = cleanMatT.computeNewArrayT(arrayT);
			
			cleanMatT.computeNewArrayMj(ComputationalDomain.Mi);
			arrayMj = new int[DeleteRowColumnNullDiagonalEntry.newMj.length];
			System.arraycopy(DeleteRowColumnNullDiagonalEntry.newMj, 0,
					arrayMj, 0, DeleteRowColumnNullDiagonalEntry.newMj.length);
			
			cleanMatT.computeNewarrayMp(ComputationalDomain.Np, ComputationalDomain.Mp);
			arrayMp = new int[DeleteRowColumnNullDiagonalEntry.newMp.length];
			System.arraycopy(DeleteRowColumnNullDiagonalEntry.newMp, 0,
					arrayMp, 0, DeleteRowColumnNullDiagonalEntry.newMp.length);

			numberOfPolygon = DeleteRowColumnNullDiagonalEntry.newMp.length-1;

			bedRockElevation = cleanMatT
					.computeCellsArray(ComputationalDomain.bedRockElevation);

			porosity = cleanMatT
					.computeCellsArray(ComputationalDomain.porosity);

			planarArea = cleanMatT
					.computeCellsArray(ComputationalDomain.planArea);

			source = cleanMatT.computeCellsArray(ComputationalDomain.source);

			coeffC = cleanMatT.computeCellsArray(ComputationalDomain.c);

			coeffM = cleanMatT.computeCellsArray(ComputationalDomain.m);

			indexDiagonal = newIndexDiag.computeIndexDiag(numberOfPolygon, arrayMp,
					arrayMj);

			etaDirichlet = cleanMatT
					.computeCellsArray(ComputationalDomain.etaDirichlet);
			
			arrayEta = cleanMatT
					.computeCellsArray(eta);
			
			matrixReduction = true;
			ComputeT.unlockDeleteRowColumn = false;

		} else {

			matT = new double[arrayT.length];
			System.arraycopy(arrayT, 0, matT, 0, arrayT.length);

			numberOfPolygon = ComputationalDomain.Np;

			arrayMp = new int[ComputationalDomain.Mp.length];
			System.arraycopy(ComputationalDomain.Mp, 0, arrayMp, 0,
					ComputationalDomain.Mp.length);

			arrayMj = new int[ComputationalDomain.Mi.length];
			System.arraycopy(ComputationalDomain.Mi, 0, arrayMj, 0,
					ComputationalDomain.Mi.length);

			bedRockElevation = new double[ComputationalDomain.bedRockElevation.length];
			System.arraycopy(ComputationalDomain.bedRockElevation, 0,
					bedRockElevation, 0,
					ComputationalDomain.bedRockElevation.length);

			porosity = new double[ComputationalDomain.porosity.length];
			System.arraycopy(ComputationalDomain.porosity, 0, porosity, 0,
					ComputationalDomain.porosity.length);

			planarArea = new double[ComputationalDomain.planArea.length];
			System.arraycopy(ComputationalDomain.planArea, 0, planarArea, 0,
					ComputationalDomain.planArea.length);

			source = new double[ComputationalDomain.source.length];
			System.arraycopy(ComputationalDomain.source, 0, source, 0,
					ComputationalDomain.source.length);

			coeffC = new double[ComputationalDomain.c.length];
			System.arraycopy(ComputationalDomain.c, 0, coeffC, 0,
					ComputationalDomain.c.length);

			coeffM = new double[ComputationalDomain.m.length];
			System.arraycopy(ComputationalDomain.m, 0, coeffM, 0,
					ComputationalDomain.m.length);

			indexDiagonal = new int[ComputationalDomain.Np];
			System.arraycopy(indexDiag, 0, indexDiagonal, 0,
					indexDiag.length);

			etaDirichlet = new double[ComputationalDomain.etaDirichlet.length];
			System.arraycopy(ComputationalDomain.etaDirichlet, 0, etaDirichlet,
					0, ComputationalDomain.etaDirichlet.length);

			arrayEta = new double[ComputationalDomain.eta.length];
			System.arraycopy(eta, 0, arrayEta,
					0, ComputationalDomain.eta.length);
			
		}
		
		return matT;
		
	}
	
}

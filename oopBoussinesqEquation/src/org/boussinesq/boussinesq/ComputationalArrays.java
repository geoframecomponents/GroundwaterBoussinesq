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

	public void wetDomain(double[] eta, int[] indexDiag,
			DeleteRowColumnNullDiagonalEntry removeDryCell) {

		newIndexDiag = new RCIndexDiagonalElement();

		System.out.println("Remove row column");

		removeDryCell.computeNewArrayMj(ComputationalDomain.Mi);
		arrayMj = new int[DeleteRowColumnNullDiagonalEntry.newMj.length];
		System.arraycopy(DeleteRowColumnNullDiagonalEntry.newMj, 0, arrayMj, 0,
				DeleteRowColumnNullDiagonalEntry.newMj.length);

		removeDryCell.computeNewarrayMp(ComputationalDomain.Np,
				ComputationalDomain.Mp);
		arrayMp = new int[DeleteRowColumnNullDiagonalEntry.newMp.length];
		System.arraycopy(DeleteRowColumnNullDiagonalEntry.newMp, 0, arrayMp, 0,
				DeleteRowColumnNullDiagonalEntry.newMp.length);

		numberOfPolygon = DeleteRowColumnNullDiagonalEntry.newMp.length - 1;

		bedRockElevation = removeDryCell
				.computeCellsArray(ComputationalDomain.bedRockElevation);

		porosity = removeDryCell.computeCellsArray(ComputationalDomain.porosity);

		planarArea = removeDryCell.computeCellsArray(ComputationalDomain.planArea);

		source = removeDryCell.computeCellsArray(ComputationalDomain.source);

		coeffC = removeDryCell.computeCellsArray(ComputationalDomain.c);

		coeffM = removeDryCell.computeCellsArray(ComputationalDomain.m);

		indexDiagonal = newIndexDiag.computeIndexDiag(numberOfPolygon, arrayMp,
				arrayMj);

		etaDirichlet = removeDryCell
				.computeCellsArray(ComputationalDomain.etaDirichlet);

		arrayEta = removeDryCell.computeCellsArray(eta);

	}

	public void completeDomain(double[] eta, int[] indexDiag) {

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
		System.arraycopy(indexDiag, 0, indexDiagonal, 0, indexDiag.length);

		etaDirichlet = new double[ComputationalDomain.etaDirichlet.length];
		System.arraycopy(ComputationalDomain.etaDirichlet, 0, etaDirichlet, 0,
				ComputationalDomain.etaDirichlet.length);

		arrayEta = new double[ComputationalDomain.eta.length];
		System.arraycopy(eta, 0, arrayEta, 0, ComputationalDomain.eta.length);

	}

}

package org.boussinesq.boussinesq.computationalDomain;

import org.boussinesq.rowCompressedForm.DeleteRowColumnNullDiagonalEntry;
import org.boussinesq.rowCompressedForm.RCIndexDiagonalElement;

public class ActualComputationalDomain {

	public static int numberOfPolygons;
	public static int[] array_Mp, array_Mj, indexDiagonal;
	public static double[] bedRockElevation, porosity, planarArea, source,
			ratingCurveCoeff_c, ratingCurveCoeff_m, etaDirichlet, eta;

	RCIndexDiagonalElement rcActualIndecesDiagonalTerms;

	public void wetDomain(double[] eta, int[] indexDiag,
			DeleteRowColumnNullDiagonalEntry removeDryCell) {

		rcActualIndecesDiagonalTerms = new RCIndexDiagonalElement();

		System.out.println("Remove row column");

		removeDryCell.computeNewArrayMj(AbstractDomain.Mi);
		ActualComputationalDomain.array_Mj = new int[DeleteRowColumnNullDiagonalEntry.newMj.length];
		System.arraycopy(DeleteRowColumnNullDiagonalEntry.newMj, 0, array_Mj,
				0, DeleteRowColumnNullDiagonalEntry.newMj.length);

		removeDryCell.computeNewArrayMp(AbstractDomain.Np, AbstractDomain.Mp);
		ActualComputationalDomain.array_Mp = new int[DeleteRowColumnNullDiagonalEntry.newMp.length];
		System.arraycopy(DeleteRowColumnNullDiagonalEntry.newMp, 0, array_Mp,
				0, DeleteRowColumnNullDiagonalEntry.newMp.length);

		ActualComputationalDomain.numberOfPolygons = DeleteRowColumnNullDiagonalEntry.newMp.length - 1;

		ActualComputationalDomain.bedRockElevation = removeDryCell
				.reduceToActualArray(AbstractDomain.bedRockElevation);

		ActualComputationalDomain.porosity = removeDryCell
				.reduceToActualArray(AbstractDomain.porosity);

		ActualComputationalDomain.planarArea = removeDryCell
				.reduceToActualArray(AbstractDomain.planArea);

		ActualComputationalDomain.source = removeDryCell
				.reduceToActualArray(AbstractDomain.source);

		ActualComputationalDomain.ratingCurveCoeff_c = removeDryCell
				.reduceToActualArray(AbstractDomain.c);

		ActualComputationalDomain.ratingCurveCoeff_m = removeDryCell
				.reduceToActualArray(AbstractDomain.m);

		ActualComputationalDomain.indexDiagonal = rcActualIndecesDiagonalTerms
				.computeIndexDiag(numberOfPolygons, array_Mp, array_Mj);

		ActualComputationalDomain.etaDirichlet = removeDryCell
				.reduceToActualArray(AbstractDomain.etaDirichlet);

		ActualComputationalDomain.eta = removeDryCell.reduceToActualArray(eta);

	}

	public void completeDomain(double[] eta, int[] indexDiag) {

		ActualComputationalDomain.numberOfPolygons = AbstractDomain.Np;

		ActualComputationalDomain.array_Mp = new int[AbstractDomain.Mp.length];
		System.arraycopy(AbstractDomain.Mp, 0, array_Mp, 0,
				AbstractDomain.Mp.length);

		ActualComputationalDomain.array_Mj = new int[AbstractDomain.Mi.length];
		System.arraycopy(AbstractDomain.Mi, 0, array_Mj, 0,
				AbstractDomain.Mi.length);

		ActualComputationalDomain.bedRockElevation = new double[AbstractDomain.bedRockElevation.length];
		System.arraycopy(AbstractDomain.bedRockElevation, 0, bedRockElevation,
				0, AbstractDomain.bedRockElevation.length);

		ActualComputationalDomain.porosity = new double[AbstractDomain.porosity.length];
		System.arraycopy(AbstractDomain.porosity, 0, porosity, 0,
				AbstractDomain.porosity.length);

		ActualComputationalDomain.planarArea = new double[AbstractDomain.planArea.length];
		System.arraycopy(AbstractDomain.planArea, 0, planarArea, 0,
				AbstractDomain.planArea.length);

		ActualComputationalDomain.source = new double[AbstractDomain.source.length];
		System.arraycopy(AbstractDomain.source, 0, source, 0,
				AbstractDomain.source.length);

		ActualComputationalDomain.ratingCurveCoeff_c = new double[AbstractDomain.c.length];
		System.arraycopy(AbstractDomain.c, 0, ratingCurveCoeff_c, 0,
				AbstractDomain.c.length);

		ActualComputationalDomain.ratingCurveCoeff_m = new double[AbstractDomain.m.length];
		System.arraycopy(AbstractDomain.m, 0, ratingCurveCoeff_m, 0,
				AbstractDomain.m.length);

		ActualComputationalDomain.indexDiagonal = new int[AbstractDomain.Np];
		System.arraycopy(indexDiag, 0, indexDiagonal, 0, indexDiag.length);

		ActualComputationalDomain.etaDirichlet = new double[AbstractDomain.etaDirichlet.length];
		System.arraycopy(AbstractDomain.etaDirichlet, 0, etaDirichlet, 0,
				AbstractDomain.etaDirichlet.length);

		ActualComputationalDomain.eta = new double[AbstractDomain.eta.length];
		System.arraycopy(eta, 0, eta, 0, AbstractDomain.eta.length);

	}

}

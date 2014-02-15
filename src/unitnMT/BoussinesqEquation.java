package unitnMT;

/**
 * Mass-Conservative groundwater equation integration
 * 
 * @desc	This class contains 7 methods:
 * 				1. computeIndexDiag: array of indices of diagonal entries
 * 				2. computeWetArea: compute the wet area of a cell at every
 * 					Newton iteration
 * 				3. computeWaterVolume: compute the water volume stored in
 * 					every cell at every Newton iteration
 * 				4. computeT: compute the array T at every time step
 * 				5. computeB: compute the array b at every time step
 * 				6. computeR: compute the array of residual function R at
 * 					every Newton iteration
 * 				7. computeJr: compute the array of Jacobian of the water
 * 					volume stored at every Newton iteration
 * 
 * @author	E. Cordano, G. Formetta, R. Rigon, F. Serafin, 2014
 * Copyright GPL v. 3 (http://www.gnu.org/licenses/gpl.html)
 * */

import cern.colt.Arrays;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation {

	/** The sol. */
	double[] solOld;
	double[] solNew;
	double[][] sol;

	/** The deltat. */
	int deltat = 3600;

	/** legth of the simulation */
	int simTime = 3600*24;

	double tol = 10 ^ (-1);

	double tolerance = 0;

	BoussinesqEquation(int Np, int SIZE) {

		solOld = new double[Np];
		solNew = new double[Np];

		tolerance = computeMachineEpsilonDouble();

	}

	/**
	 * Calculate machine epsilon double.
	 * 
	 * @desc this method compute the tolerance of the machine. For more info go
	 *       to
	 *       https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java
	 *       . In c/c++ section there's write that:
	 * 
	 *       In such languages as C or C++ when you do something like while( 1.0
	 *       + eps > 1.0 ) the expression is calculated not with 64 bits
	 *       (double) but with the processor precision (80 bits or more depends
	 *       on the processor and compile options). Below program calculates
	 *       exactly on 32 bits (float) and 64 bits (double)
	 * 
	 * 
	 * @return the tolerance of the machine
	 */
	private static double computeMachineEpsilonDouble() {

		// machine tolerance
		double machEps = 1.0d;

		do
			machEps /= 2.0d;
		while ((double) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	/**
	 * Compute index of diagonal.
	 * 
	 * @desc this method computes the indices of the diagonal of the adjacency
	 *       matrix; this matrix and all the sparse matrices are stored in Row
	 *       Compressed Form. More information at the web site
	 *       https://en.wikipedia.org/wiki/Sparse_matrix
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * 
	 * @return the array that holds the indices of the diagonal entries of the
	 *         sparse adjacency matrix in Row Compressed Form
	 */
	public int[] computeIndexDiag(Grid mesh) {

		// declaration of the array that holds the indices of diagonal entries
		int[] indexDiag = new int[mesh.numberSidesPolygon.length];

		/* for-loop to analyze the matrix cell by cell */
		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {
			/*
			 * nested for-loop to analyze diagonal entries, which are identified
			 * by a negative number
			 */
			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

				if (mesh.Mi[j] == i) {
					indexDiag[i] = j;
				}

			}
		}

		return indexDiag;
	}

	/**
	 * Compute wet area.
	 * 
	 * @desc this method computes the wet area of the cell; if the piezometric
	 *       head is less than the bedrock elevation, the wet area is null,
	 *       otherwise is equal to the porosity for planimetric area of the
	 *       cell.
	 * 
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the bedrock elevation
	 * @param porosity
	 *            the porosity
	 * @param planimetricArea
	 *            the planimetric area of the cell
	 * 
	 * @return the wet area of the cell
	 */
	public double computeWetArea(double eta, double zetaBedrock,
			double porosity, double planimetricArea) {

		// the wet area is the variable that the method returns
		double wetArea = 0;

		if (eta > zetaBedrock) {
			wetArea = porosity * planimetricArea;
		} else {
			wetArea = 0;
		}

		return wetArea;
	}

	/**
	 * Compute water volume.
	 * 
	 * @desc this method computes the volume of stored water in the cell like
	 *       multiplication between the wet area and the thickness of the water
	 *       table
	 * 
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the bedrock elevation
	 * @param planimetricArea
	 *            the planimetric area of the cell
	 * 
	 * @return the volume of the stored water in the cell
	 */
	public double computeWaterVolume(double eta, double zetaBedrock,
			double wetArea) {

		double volume = wetArea * (eta - zetaBedrock);

		return volume;

	}

	/**
	 * Compute T.
	 * 
	 * @desc this method computes the elements of the matrix T in Row Compressed
	 *       Form according to the equations (20) e (21) of [Cordano & Rigon,
	 *       Mass-conservative method for the integration of the two-dimensional
	 *       groundwater (Boussinesq) equation, Water Resources Research 2012]
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the piezometric head
	 * 
	 * @return the matrix T like array in Row Compressed Form
	 */
	public double[] computeT(Grid mesh, double[] eta) {

		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */
		double rowSum = 0;

		/* to identify the diagonal entry of matrix T in row-compressed form */
		int index = 0;

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

				if (mesh.Mi[j] != i) {
					// equation (21)
					arrayT[j] = -deltat
							* (1 / mesh.euclideanDistance[(int) mesh.Ml[j]])
							* mesh.hydrConductivity[(int) mesh.Ml[j]]
							* mesh.lengthSides[(int) mesh.Ml[j]]
							* Math.max(
									Math.max(0, eta[mesh.Mi[j]]
											- mesh.bottomElevation[mesh.Mi[j]]),
									Math.max(0, eta[i]
											- mesh.bottomElevation[i]));

					rowSum += -arrayT[j];

				} else {
					index = j;
				}

			}
			// equation (20)
			arrayT[index] = rowSum;
			rowSum = 0;
		}
		
		return arrayT;
	}

	/**
	 * Compute b.
	 * 
	 * @desc this method computes the elements of the array b of the known terms
	 *       of the linear system for every cell according to the equations (19)
	 *       of [Cordano & Rigon, Mass-conservative method for the integration
	 *       of the two-dimensional groundwater (Boussinesq) equation, Water
	 *       Resources Research 2012]
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the the piezometric head
	 * 
	 * @return the array b
	 */
	public double[] computeB(Grid mesh, double[] eta, double[] Tdrichelet,
			double[] etaDrichelet) {

		// declaration of the array that holds the known terms of the linear
		// system
		double[] arrB = new double[mesh.numberSidesPolygon.length];

		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {
			// compute the wet area
			double wetArea = computeWetArea(eta[i], mesh.bottomElevation[i],
					mesh.porosity[i], mesh.planArea[i]);
			// compute the water volume stored in the cell
			double volume = computeWaterVolume(eta[i], mesh.bottomElevation[i],
					wetArea);
			// equation (19)
			double sum = 0;
			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
				sum += Tdrichelet[j] * etaDrichelet[mesh.Mi[j]];
			}

			arrB[i] = volume + deltat * mesh.planArea[i] * mesh.source[i] - sum;
		}

		return arrB;
	}

	/**
	 * Compute R.
	 * 
	 * @desc this method computes the values of the residual function at every
	 *       iteration of the Newton's method for every cell according to the
	 *       equation (A3) of [Cordano & Rigon, Mass-conservative method for the
	 *       integration of the two-dimensional groundwater (Boussinesq)
	 *       equation, Water Resources Research 2012]
	 * 
	 * @param arrT
	 *            the array of T in Row Compressed Form
	 * @param arrb
	 *            the array of known terms
	 * @param zetaBedrock
	 *            the bedrock elevation
	 * @param porosity
	 *            the porosity
	 * @param Np
	 *            the number of polygons
	 * @param Mp
	 *            the array that holds the number of non-zero entries in
	 *            adjacency matrix
	 * @param Mi
	 *            the array that holds the column indices of non-zero entries
	 * @param eta
	 *            the piezometric head
	 * @param planimetricArea
	 *            the planimetric area of the cell
	 * 
	 * @return the array of the residual function
	 */
	public double[] computeR(double[] arrT, double[] arrb,
			double[] zetaBedrock, double[] porosity, int Np, int[] Mp,
			int[] Mi, double[] eta, double[] planimetricArea,double[] etaDrichelet, double NOVALUE) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[Np];

		for (int i = 0; i < Np; i++) {
			if (isNoValue(etaDrichelet[i],NOVALUE)){
			for (int j = Mp[i]; j < Mp[i + 1]; j++) {
				sum += arrT[j] * eta[Mi[j]];
			}

			double wetArea = computeWetArea(eta[i], zetaBedrock[i],
					porosity[i], planimetricArea[i]);
			double waterVolume = computeWaterVolume(eta[i], zetaBedrock[i],
					wetArea);
			// equation (A3)
			arrR[i] = waterVolume + sum - arrb[i];

			sum = 0;
			}else{
				arrR[i] = 0;
			}
		}

		return arrR;
	}

	/**
	 * Compute Jr.
	 * 
	 * @desc this method computes the Jacobian matrix of the water volume stored
	 *       into every cell. In this case the array Jr, in Row Compressed Form,
	 *       is evaluated like sum between array T and the wet area, according
	 *       the equation (A6) and (A7) of [Cordano & Rigon, Mass-conservative
	 *       method for the integration of the two-dimensional groundwater
	 *       (Boussinesq) equation, Water Resources Research 2012]. The array Jr
	 *       is a copy of T where only diagonal entries are summed to P, because
	 *       P is a diagonal matrix in Row Compressed Form too.
	 * 
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param arrT
	 *            the array of T in Row Compressed Form
	 * @param eta
	 *            the piezometric head
	 * @param zetaBed
	 *            the bedrock elevation
	 * @param porosity
	 *            the porosity
	 * @param planimetricArea
	 *            the planimetric area of a cell
	 * 
	 * @return the Jacobian array of water volume stored in Row Compressed Form
	 */
	public double[] computeJr(int[] indexDiag, double[] arrT, double[] eta,
			double[] zetaBedrock, double[] porosity, double[] planimetricArea,double[] etaDrichelet,double NOVALUE) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		double[] arrJr = new double[arrT.length];

		System.arraycopy(arrT, 0, arrJr, 0, arrT.length);

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < indexDiag.length; i++) {
			if (isNoValue(etaDrichelet[i],NOVALUE)){
				arrJr[indexDiag[i]] = arrT[indexDiag[i]]
						+ computeWetArea(eta[i], zetaBedrock[i], porosity[i],
								planimetricArea[i]);

			}else {
				arrJr[indexDiag[i]] = arrT[indexDiag[i]];
			}

			// equation (A6)
			

		}

		return arrJr;
	}

	public double[] computeTDrichelet(Grid mesh, double[] T) {

		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */

		/* to identify the diagonal entry of matrix T in row-compressed form */

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */

		/* to identify the diagonal entry of matrix T in row-compressed form */

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			if (!isNoValue(mesh.etaDrichelet[i], mesh.NOVALUE)) {// is drichelet
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
					arrayT[j] = T[j];
				}
			} else {// is not drichelet
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (!isNoValue(mesh.etaDrichelet[mesh.Mi[j]], mesh.NOVALUE)) {// is drichelet
						arrayT[j] = T[j];
					} else {
						arrayT[j] = 0.0;
					}
				}

			}

		}

		return arrayT;
	}

	public double[] computeTNoDrichelet(Grid mesh, double[] T) {

		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */

		/* to identify the diagonal entry of matrix T in row-compressed form */

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		/*
		 * variable to which sum the terms of matrix T (T is an array because is
		 * in RC-F) that are outside the diagonal; after investigation of the
		 * row of the matrix the value is stored in the diagonal of matrix T
		 */

		/* to identify the diagonal entry of matrix T in row-compressed form */

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			if (isNoValue(mesh.etaDrichelet[i], mesh.NOVALUE)) {
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (isNoValue(mesh.etaDrichelet[mesh.Mi[j]], mesh.NOVALUE)) {
						arrayT[j] = T[j];
					} else {
						arrayT[j] = 0.0;
					}
				}

			}

		}

		return arrayT;
	}

	public boolean isNoValue(double x, double noValue) {

		if (x <= noValue) {
			return true;
		} else {
			return false;
		}
	}

	public double[] newtonIteration(double[] arrb, double[] matT,
			int[] indexDiag, Grid grid1, double[] eta)
			throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;

		double maxResidual = 10;

		do {

			double[] jr = computeJr(indexDiag, matT, eta,
					grid1.bottomElevation, grid1.porosity, grid1.planArea,grid1.etaDrichelet,grid1.NOVALUE);
			matrixJr = new SparseRCDoubleMatrix2D(
					grid1.numberSidesPolygon.length,
					grid1.numberSidesPolygon.length, grid1.Mp, grid1.Mi, jr);

			double[] r = computeR(matT, arrb, grid1.bottomElevation,
					grid1.porosity, grid1.numberSidesPolygon.length, grid1.Mp,
					grid1.Mi, eta, grid1.planArea,grid1.etaDrichelet,grid1.NOVALUE);

			matrixr = new SparseDoubleMatrix1D(r);

			// System.out.println("Index prova " + indexProva);
			RCConjugateGradient cg = new RCConjugateGradient(
					grid1.numberSidesPolygon.length);

			cg.solverCG(matrixr, matrixJr);

			for (int i = 0; i < grid1.eta.length; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
				// System.out.println(cg.matSol.get(i));
				// System.out.println(solNew[i]);
			}
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));
			//System.out.println(maxResidual);

			/*
			 * for (int i = 1; i < grid1.numberSidesPolygon.length; i++) { if
			 * (Math.abs(cg.matSol.get(i)) > maxResidual) { maxResidual =
			 * Math.abs(cg.matSol.get(i)); }
			 * 
			 * }
			 */

		} while (maxResidual > tolerance*100);// arrR.getMaxLocation()[0] > tol |
		return eta;

	}

	public void runBEq(Grid grid1)
			throws IterativeSolverDoubleNotConvergedException {
		// the values of first attempt are the
		
		double[] eta = new double[grid1.numberSidesPolygon.length];
		for(int i=0; i<eta.length;i++){
			if(isNoValue(grid1.etaDrichelet[i],grid1.NOVALUE)){
				eta[i] = grid1.eta[i];
			} else {
				eta[i] = grid1.etaDrichelet[i];
			}
		}

		int[] indexDiag = computeIndexDiag(grid1);

		for (int t = 0; t < simTime; t += deltat) {

			long start=System.nanoTime();
			
			double[] matT = computeT(grid1, eta);
			double[] matTDrichelet = computeTDrichelet(grid1, matT);
			double[] matTNoDrichelet = computeTNoDrichelet(grid1, matT);

			double[] arrb = computeB(grid1, eta, matTDrichelet,grid1.etaDrichelet);
			
			long end = System.nanoTime();
			System.out.println("End time: " + (end-start));

			eta = newtonIteration(arrb, matTNoDrichelet, indexDiag, grid1, eta);
			System.out.println(Arrays.toString(eta));
			System.out.println(t/3600);
			//System.out.println(Arrays.toString(eta));
		}
		System.out.println("Exit code");
	}

	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException {

		long start=System.nanoTime();
		Grid grid1 = new Grid("Song");
		BoussinesqEquation beq = new BoussinesqEquation(
				grid1.numberSidesPolygon.length, grid1.Ml.length);
		beq.runBEq(grid1);
		long end=System.nanoTime();
		System.out.println("End time: " + (end-start));
		
		System.exit(1);

	}

}

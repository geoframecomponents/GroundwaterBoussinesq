package org.boussinesq.boussinesq;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

public class SetDirichletBC {

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
	 *       2012]
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the piezometric head
	 * 
	 * @return the matrix T like array in Row Compressed Form
	 */
	public double[] computeT(double[] eta) {

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
		double[] arrayT = new double[Mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Mesh.numberSidesPolygon.length; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
			for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {

				if (Mesh.Mi[j] != i) {
					// equation (21)
					arrayT[j] = -BoussinesqEquation.deltat
							* (1 / Mesh.euclideanDistance[(int) Mesh.Ml[j]])
							* Mesh.hydrConductivity[(int) Mesh.Ml[j]]
							* Mesh.lengthSides[(int) Mesh.Ml[j]]
							* Math.max(
									Math.max(0, eta[Mesh.Mi[j]]
											- Mesh.bottomElevation[Mesh.Mi[j]]),
									Math.max(0, eta[i]
											- Mesh.bottomElevation[i]));

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
	 * Compute B.
	 * 
	 * @desc this method computes the elements of the array b of the known terms
	 *       of the linear system for every cell according to the equations (19)
	 *       of [Cordano & Rigon, 2012]. By considering head-based boundary
	 *       condi- tions (Dirichlet), it's necessary to subtract the product
	 *       between T matrix of Dirichlet cells and eta array of Dirichlet
	 *       cells. Thus the right side of equation (32) is solved. It's
	 *       considered the outflow scale of the closure point too, thus it is
	 *       subtract the outgoing flow from of last cell at the right side of
	 *       equation (32)
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the the piezometric head
	 * @param Tdrichelet
	 *            the matrix T computes for the Dirichlet cells
	 * @param etaDrichelet
	 *            the array of known eta in the Dirichlet cells
	 * 
	 * @return the array b
	 */
	public double[] computeB(double[] eta, double[] Tdrichelet,
			double[] etaDrichelet) {

		// declaration of the array that holds the known terms of the linear
		// system
		double[] arrB = new double[Mesh.numberSidesPolygon.length];

		for (int i = 0; i < Mesh.numberSidesPolygon.length; i++) {
			// compute the wet area
			double wetArea = computeWetArea(eta[i], Mesh.bottomElevation[i],
					Mesh.porosity[i], Mesh.planArea[i]);
			// compute the water volume stored in the cell
			double volume = computeWaterVolume(eta[i], Mesh.bottomElevation[i],
					wetArea);
			// equation (19)
			double sum = 0;
			for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {
				sum += Tdrichelet[j] * etaDrichelet[Mesh.Mi[j]];
			}

			// delta t deve essere minore di 1/c
			arrB[i] = volume + BoussinesqEquation.deltat * Mesh.planArea[i] * Mesh.source[i] - sum
					- BoussinesqEquation.deltat * Mesh.planArea[i] * Mesh.c[i]
					* Math.pow(volume / Mesh.planArea[i], Mesh.m[i]);

		}

		return arrB;
	}

	/**
	 * Compute R.
	 * 
	 * @desc this method computes the values of the residual function at every
	 *       iteration of the Newton's method for every cell according to the
	 *       equation (A3) of [Cordano & Rigon, 2012]. By analyzing every cell,
	 *       the residual function is computed only for polygons that are not
	 *       Dirichlet cells, otherwise the residual is imposed equal to zero,
	 *       so the known piezometric head of Dirichlet cells remains constant
	 *       during the Newton's loop.
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
	public double[] computeR(double[] arrT, double[] arrb, double[] eta) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[Mesh.numberSidesPolygon.length];

		for (int i = 0; i < Mesh.numberSidesPolygon.length; i++) {
			if (isNoValue(Mesh.etaDrichelet[i], Mesh.NOVALUE)) {

				// non Dirichlet cells
				for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {
					sum += arrT[j] * eta[Mesh.Mi[j]];
				}

				double wetArea = computeWetArea(eta[i],
						Mesh.bottomElevation[i], Mesh.porosity[i],
						Mesh.planArea[i]);
				double waterVolume = computeWaterVolume(eta[i],
						Mesh.bottomElevation[i], wetArea);
				// equation (A3)
				arrR[i] = waterVolume + sum - arrb[i];

				sum = 0;
			} else {

				// Dirichlet cells
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
	 *       the equation (A6) and (A7) of [Cordano & Rigon, 2012]. The array Jr
	 *       is a copy of T where only diagonal entries are summed to P, because
	 *       P is a diagonal matrix in Row Compressed Form too. These operations
	 *       are made only in case the Jacobian is computed only in a non
	 *       Dirichlet cell. Otherwise the volume of water stored is constant
	 *       with eta and P is equal to zero.
	 * 
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param arrT
	 *            the array of T in Row Compressed Form
	 * @param eta
	 *            the piezometric head
	 * @param zetaBedrock
	 *            the zeta bedrock
	 * @param porosity
	 *            the porosity
	 * @param planimetricArea
	 *            the planimetric area of a cell
	 * @param etaDirichlet
	 *            the eta of Dirichlet cells
	 * @param NOVALUE
	 *            the novalue
	 * 
	 * @return the Jacobian array of water volume stored in Row Compressed Form
	 */
	public double[] computeJr(int[] indexDiag, double[] arrT, double[] eta) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		double[] arrJr = new double[arrT.length];

		System.arraycopy(arrT, 0, arrJr, 0, arrT.length);

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < indexDiag.length; i++) {

			if (isNoValue(Mesh.etaDrichelet[i], Mesh.NOVALUE)) {
				// non Dirichlet cells
				// equation (A6)
				arrJr[indexDiag[i]] = arrT[indexDiag[i]]
						+ computeWetArea(eta[i], Mesh.bottomElevation[i],
								Mesh.porosity[i], Mesh.planArea[i]);

			} else {
				// Dirichlet cells
				arrJr[indexDiag[i]] = arrT[indexDiag[i]];
			}
		}

		return arrJr;
	}

	/**
	 * Compute T for Dirichlet cells.
	 * 
	 * @desc according the head-based Boundary Conditions (Dirichlet) at the
	 *       equation (30) of [Cordano & Rigon, 2012], the matrix T in row
	 *       compressed form for a Dirichlet cell is computed only if the cell
	 *       that I'm observing or the adjacency cell are a Dirichlet cell.
	 *       Other T is imposed equal to zero.
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param T
	 *            the array of T in Row Compressed Form
	 * 
	 * @return the array of T in RC-F for Dirichlet cells
	 */
	public double[] computeTDirichlet(double[] T) {

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[Mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Mesh.numberSidesPolygon.length; i++) {

			if (!isNoValue(Mesh.etaDrichelet[i], Mesh.NOVALUE)) {

				// Dirichlet cells
				for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {
					arrayT[j] = T[j];
				}
			} else {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
				for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {

					if (!isNoValue(Mesh.etaDrichelet[Mesh.Mi[j]], Mesh.NOVALUE)) {

						// adjacent Dirichlet cell
						arrayT[j] = T[j];
					} else {

						// adjacent non Dirichlet cell
						arrayT[j] = 0.0;
					}
				}

			}

		}

		return arrayT;
	}

	/**
	 * Compute T for non Dirichlet cells.
	 * 
	 * @desc according the head-based Boundary Conditions (Dirichlet) at the
	 *       equation (29) of [Cordano & Rigon, 2012], the matrix T in row
	 *       compressed form for a non Dirichlet cell is computed only if
	 *       neither the cell that I'm observing nor the adjacency cell are a
	 *       Dirichlet cell. Other T is imposed equal to zero.
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param T
	 *            the array of T in Row Compressed Form
	 * 
	 * @return the array of T in RC-F for non Dirichlet cells
	 */
	public double[] computeTNoDirichlet(double[] T) {

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[Mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < Mesh.numberSidesPolygon.length; i++) {

			if (isNoValue(Mesh.etaDrichelet[i], Mesh.NOVALUE)) {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
				for (int j = Mesh.Mp[i]; j < Mesh.Mp[i + 1]; j++) {

					if (isNoValue(Mesh.etaDrichelet[Mesh.Mi[j]], Mesh.NOVALUE)) {

						// adjacent non Dirichlet cells
						arrayT[j] = T[j];
					} else {

						// adjacent Dirichlet cells
						arrayT[j] = 0.0;
					}
				}

			}

		}

		return arrayT;
	}

	/**
	 * Checks if is no value.
	 * 
	 * @desc this method evaluates if the entry value is a NaN or not. If is
	 *       NaN, the cell I'm observing is a non Dirichlet cell, otherwise is a
	 *       Dirichlet cell
	 * 
	 * @param x
	 *            the eta of the cell that I'm analyzing
	 * @param noValue
	 *            the no value
	 * 
	 * @return boolean true if the cell isn't a Dirichlet cell, false otherwise
	 */
	public boolean isNoValue(double x, double noValue) {

		if (x <= noValue) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Newton iteration.
	 * 
	 * @desc this method compute the Newton iteration.
	 * 
	 * @param arrb
	 *            the array of known terms
	 * @param arrT
	 *            the array of T for non Dirichlet cells
	 * @param indexDiag
	 *            the array of the indices of the diagonal entries
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @param eta
	 *            the array of eta at the previous time step
	 * @param cg
	 *            the conjugate gradient class
	 * 
	 * @return the array of eta at the following time step
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public double[] newtonIteration(double[] arrb, double[] arrT,
			int[] indexDiag, double[] eta, RCConjugateGradient cg)
			throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;

		double maxResidual = 10;

		do {

			// compute Jr
			double[] jr = computeJr(indexDiag, arrT, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixJr = new SparseRCDoubleMatrix2D(
					Mesh.numberSidesPolygon.length,
					Mesh.numberSidesPolygon.length, Mesh.Mp, Mesh.Mi, jr);

			// compute the residual function
			double[] r = computeR(arrT, arrb, eta);

			// convert array in sparse matrix for DoubleCG class
			matrixr = new SparseDoubleMatrix1D(r);

			cg.solverCG(matrixr, matrixJr);

			// compute the new eta for every cell
			for (int i = 0; i < Mesh.eta.length; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
			}

			// compute the max residual
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));

		} while (maxResidual > BoussinesqEquation.tolerance * 100);

		return eta;

	}

}

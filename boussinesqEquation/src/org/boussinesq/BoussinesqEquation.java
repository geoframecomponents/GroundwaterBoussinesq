package org.boussinesq;

/**
 * Mass-Conservative groundwater equation integration
 * 
 * @desc	In this code is implemented a conservative finite-volume
 * 			numerical solution for the two-dimensional groundwater flow
 * 			(Boussinesq) equation, which can be used for investigation
 * 			of hillslope subsurface flow processes and simulations of
 * 			catchment hydrology.
 * 
 * 			The idea is taken from:
 * 
 * 			"A mass-conservative method for the integration of the
 * 			 two dimensional groundwater (Boussinesq) equation"
 * 			E.Cordano, R.Rigon 2012 - Water Resources Research
 * 
 * @author	E. Cordano, G. Formetta, R. Rigon, F. Serafin, 2014
 * Copyright GPL v. 3 (http://www.gnu.org/licenses/gpl.html)
 * */

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation {

	/** The deltat. */
	int deltat = 360;

	/** legth of the simulation */
	int simTime = 3600 * 24 * 5;

	double tolerance = 0;

	BoussinesqEquation() {

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

			// delta t deve essere minore di 1/c
			arrB[i] = volume + deltat * mesh.planArea[i] * mesh.source[i] - sum
					- deltat * mesh.planArea[i] * mesh.c[i]
					* Math.pow(volume / mesh.planArea[i], mesh.m[i]);

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
	public double[] computeR(double[] arrT, double[] arrb,
			double[] zetaBedrock, double[] porosity, int Np, int[] Mp,
			int[] Mi, double[] eta, double[] planimetricArea,
			double[] etaDrichelet, double NOVALUE) {

		// variable where allocate the matrix-vector multiplication
		double sum = 0;
		// declaration of the array that holds the residual function for every
		// cell
		double[] arrR = new double[Np];

		for (int i = 0; i < Np; i++) {
			if (isNoValue(etaDrichelet[i], NOVALUE)) {

				// non Dirichlet cells
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
	public double[] computeJr(int[] indexDiag, double[] arrT, double[] eta,
			double[] zetaBedrock, double[] porosity, double[] planimetricArea,
			double[] etaDrichelet, double NOVALUE) {

		// declaration of the array that holds the Jacobian of water volume
		// stored
		double[] arrJr = new double[arrT.length];

		System.arraycopy(arrT, 0, arrJr, 0, arrT.length);

		// cicle only in the cells, because it's necessary to inspect only
		// diagonal entries
		for (int i = 0; i < indexDiag.length; i++) {

			if (isNoValue(etaDrichelet[i], NOVALUE)) {
				// non Dirichlet cells
				// equation (A6)
				arrJr[indexDiag[i]] = arrT[indexDiag[i]]
						+ computeWetArea(eta[i], zetaBedrock[i], porosity[i],
								planimetricArea[i]);

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
	public double[] computeTDirichlet(Grid mesh, double[] T) {

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {

			if (!isNoValue(mesh.etaDrichelet[i], mesh.NOVALUE)) {

				// Dirichlet cells
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
					arrayT[j] = T[j];
				}
			} else {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (!isNoValue(mesh.etaDrichelet[mesh.Mi[j]], mesh.NOVALUE)) {

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
	public double[] computeTNoDirichlet(Grid mesh, double[] T) {

		/*
		 * the matrix T is an array because this code uses the Row Compressed
		 * Form to stored sparse matrix
		 */
		double[] arrayT = new double[mesh.Ml.length];

		/* for-loop to analyze the mesh cell by cell */
		for (int i = 0; i < mesh.numberSidesPolygon.length; i++) {

			if (isNoValue(mesh.etaDrichelet[i], mesh.NOVALUE)) {

				// non Dirichlet cells
				/*
				 * nested for-loop to analyze shared edges between the i-th cell
				 * and the Mi[j]-th cell
				 */
				for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

					if (isNoValue(mesh.etaDrichelet[mesh.Mi[j]], mesh.NOVALUE)) {

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
			int[] indexDiag, Grid mesh, double[] eta, RCConjugateGradient cg)
			throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;

		double maxResidual = 10;

		do {

			// compute Jr
			double[] jr = computeJr(indexDiag, arrT, eta, mesh.bottomElevation,
					mesh.porosity, mesh.planArea, mesh.etaDrichelet,
					mesh.NOVALUE);

			// convert array in sparse matrix for DoubleCG class
			matrixJr = new SparseRCDoubleMatrix2D(
					mesh.numberSidesPolygon.length,
					mesh.numberSidesPolygon.length, mesh.Mp, mesh.Mi, jr);

			// compute the residual function
			double[] r = computeR(arrT, arrb, mesh.bottomElevation,
					mesh.porosity, mesh.numberSidesPolygon.length, mesh.Mp,
					mesh.Mi, eta, mesh.planArea, mesh.etaDrichelet,
					mesh.NOVALUE);

			// convert array in sparse matrix for DoubleCG class
			matrixr = new SparseDoubleMatrix1D(r);

			cg.solverCG(matrixr, matrixJr);

			// compute the new eta for every cell
			for (int i = 0; i < mesh.eta.length; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
			}

			// compute the max residual
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));

		} while (maxResidual > tolerance * 100);

		return eta;

	}

	/**
	 * Compute the Boussinesq Equation.
	 * 
	 * @desc in this method the temporal loop is implemented. Before start the
	 * 		 loop, the eta array is initialized with eta of Dirichlet if the cell
	 * 		 is a Dirichlet cells, otherwise it's inizialized with first attempt
	 * 		 value.
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void computeBEq(Grid mesh)
			throws IterativeSolverDoubleNotConvergedException, IOException {

		// allocate the memory for eta array
		double[] eta = new double[mesh.numberSidesPolygon.length];
		
		// initialize eta array
		for (int i = 0; i < eta.length; i++) {
			if (isNoValue(mesh.etaDrichelet[i], mesh.NOVALUE)) {
				
				// not Dirichlet cells
				eta[i] = mesh.eta[i];
			} else {
				
				// Dirichlet cells
				eta[i] = mesh.etaDrichelet[i];
			}
		}

		// new conjugate gradient object
		RCConjugateGradient cg = new RCConjugateGradient(
				mesh.numberSidesPolygon.length);

		int[] indexDiag = computeIndexDiag(mesh);
		FileWriter Rstatfile = new FileWriter(mesh.outputPathBeq);
		PrintWriter errestat = new PrintWriter(Rstatfile);

		for (int t = 0; t < simTime; t += deltat) {

			double[] matT = computeT(mesh, eta);
			double[] matTDrichelet = computeTDirichlet(mesh, matT);
			double[] matTNoDrichelet = computeTNoDirichlet(mesh, matT);

			double[] arrb = computeB(mesh, eta, matTDrichelet,
					mesh.etaDrichelet);

			eta = newtonIteration(arrb, matTNoDrichelet, indexDiag, mesh, eta,
					cg);

			System.out.println("Simulation time: " + (double) t / 3600);

		}
		for (int j = 0; j < eta.length; j++) {

			errestat.println(eta[j]);

		}

		errestat.println();
		System.out.println();
		Rstatfile.close();

		System.out.println("Exit code");
	}

	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException, IOException {
		String simulationType = "Song";
		// long start=System.nanoTime();
		Grid grid1 = new Grid(simulationType);
		BoussinesqEquation beq = new BoussinesqEquation();
		beq.computeBEq(grid1);
		
		if (simulationType == "Song"){
			double[] songSol = new double[grid1.numberSidesPolygon.length];
			Song s = new Song(beq.simTime, grid1.numberSidesPolygon.length,
					grid1.hydrConductivity[0]);
			
			songSol = s.beqSong(grid1.porosity);
			
			FileWriter Rstatfile = new FileWriter(grid1.outputPathSong);
			PrintWriter errestat = new PrintWriter(Rstatfile);
			
			for (int j = 0; j < songSol.length; j++) {

				errestat.println(songSol[j]);

			}

			errestat.println();
			System.out.println();
			Rstatfile.close();
			
		}
		// long end=System.nanoTime();
		// System.out.println("End time: " + (end-start));

		System.exit(1);

	}

}

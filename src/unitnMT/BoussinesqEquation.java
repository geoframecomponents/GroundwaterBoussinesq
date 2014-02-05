package unitnMT;

/**
 * Mass-Conservative groundwater equation integration
 * 
 * @desc	This class contains 4 methods:
 * 				1. estimateT: code to define T, from [Cordano Rigon 2013]
 * 				2. estimateR: code to define R, from [Cordano Rigon 2013]
 * 				3. estimateJr: code to define Jr, from [Cordano Rigon 2013]
 * 				4. newtonBEq: method that call the solver
 * 
 * 			The constructor method instantiates the object of BEq:
 * 				1. arrb: b from [Cordano Rigon 2013]
 * 				2. matT: T from [Cordano Rigon 2013]
 * 				3. matJr: Jr from [Cordano Rigon 2013]
 * 				4. arrR: R from [cordano Rigon 2013]
 * 				5. deltat: time step
 * 
 * @author	Francesco Serafin, 2014
 * Copyright GPL v. 3 (http://www.gnu.org/licenses/gpl.html)
 * */

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
//import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

// TODO: Auto-generated Javadoc
/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation {

	/** The arrb. */

	/** The mat t. */

	/** The sol. */
	double[] solOld;
	double[] solNew;
	double[][] sol;

	/** The arr r. */

	/** The deltat. */
	int deltat = 1;

	/** legth of the simulation */
	int simTime = 2;

	double tol = 10 ^ (-1);

	/**
	 * Instantiates a new boussinesq equation object.
	 * 
	 * @param Np
	 *            : the number of polygons that make up the mesh
	 * @param SIZE
	 *            : the number of non-zero in the Mij adjacency matrix
	 */
	BoussinesqEquation(int Np, int SIZE) {

		System.out.println("Number of polygons:" + Np);
		System.out.println("Number of elements of T:" + SIZE);

		solOld = new double[Np];
		solNew = new double[Np];
		sol = new double[Np][(int) simTime / deltat];

	}

	public double computePAplanar(double eta, double zetaBed, double porosity,
			double pArea) {

		double area = 0;

		if (eta - zetaBed > 0) {
			area = porosity * pArea;
		} else {
			area = 0;
		}

		return area;

	}

	public double computeVplanar(double eta, double zetaBed, double pArea) {

		double v = pArea * (eta - zetaBed);

		return v;

	}

	/**
	 * Estimate T.
	 * 
	 * @desc This method estimates the value of T, from [Cordano Rigon 2013]
	 *       equations (18), (20) e (21). Because the for-loop is going, it's
	 *       possible to fill the b array from equation (19).
	 * 
	 * @param grid1
	 *            : the object grid1 holds all the properties of the mesh. In
	 *            this class are useful these variable: - numberSidesPolygon: to
	 *            cycling 1..Np number of polygons - eta: water-table elevation
	 *            (piezometric head) - bottomElevation: bedrock elevation -
	 *            planArea: area of each cell - sourceSink: source term -
	 *            euclideanDistance: euclidean distance between each center of
	 *            every polygon (it's in array form) - hydrConductivity:
	 *            saturated hydraulic conductivity as pro- perties of every edge
	 *            - lengthSides: length of every edge of the polygon - Mp, Mi,
	 *            Ml: row compressed form of adjacency matrix
	 */
	public double[] estimateT(Grid grid1, double[] eta) {

		/*
		 * variable to add T terms outside the diagonal, that will be stored in
		 * the diagonal position
		 */
		double colSum = 0;
		/* to identify the diagonal entry in row-compressed form */
		int index = 0;

		double[] matT = new double[grid1.Ml.length];

		for (int i = 0; i < grid1.numberSidesPolygon.length; i++) {
			for (int j = grid1.Mp[i]; j < grid1.Mp[i + 1]; j++) {
				if (grid1.Mi[j] != j) {
					matT[j] = -deltat
							* (1 / grid1.euclideanDistance[(int) grid1.Ml[j]])
							* grid1.hydrConductivity[(int) grid1.Ml[j]]
							* grid1.lengthSides[(int) grid1.Ml[j]]
							* Math.max(Math.max(0, eta[grid1.Mi[j]]
									- grid1.bottomElevation[grid1.Mi[j]]), Math
									.max(0, eta[i] - grid1.bottomElevation[i]));
					colSum += -matT[j];
				} else {
					index = j;
				}
			}
			matT[index] = colSum;
			colSum = 0;
		}

		return matT;
	}

	public int[] computeIndexDiag(Grid grid1) {

		/*
		 * variable to add T terms outside the diagonal, that will be stored in
		 * the diagonal position
		 */
		double colSum = 0;
		/* to identify the diagonal entry in row-compressed form */

		int[] indexDiag = new int[grid1.numberSidesPolygon.length];

		for (int i = 0; i < grid1.numberSidesPolygon.length; i++) {
			for (int j = grid1.Mp[i]; j < grid1.Mp[i + 1]; j++) {
				if (grid1.Mi[j] == j) {
					indexDiag[i] = j;
				}
			}
		}

		return indexDiag;
	}

	public double[] estimateB(Grid grid1, double[] eta) {

		/*
		 * variable to add T terms outside the diagonal, that will be stored in
		 * the diagonal position
		 */
		double colSum = 0;
		/* to identify the diagonal entry in row-compressed form */
		// int index = 0;

		double[] arrb = new double[grid1.numberSidesPolygon.length];

		for (int i = 0; i < grid1.numberSidesPolygon.length; i++) {
			double wetArea = computePAplanar(eta[i], grid1.bottomElevation[i],
					grid1.porosity[i], grid1.planArea[i]);
			double volume = computeVplanar(eta[i], grid1.bottomElevation[i],
					wetArea);
			arrb[i] = volume + deltat * grid1.planArea[i] * grid1.source[i];
		}

		return arrb;
	}

	/**
	 * Estimate R.
	 * 
	 * @desc This method estimates the value of R, from [Cordano Rigon 2013]
	 *       equation (A3).
	 * 
	 * @param Np
	 *            : number of polygons in the mesh
	 * @param Mp
	 *            : array of row pointers in Row-Compressed Form
	 * @param Mi
	 *            : array of column indices of entries in row j
	 * @param eta
	 *            : may holds the value of the first attempt or the old values
	 *            if the method is called into the Newton-while-loop
	 * @param p
	 *            : planimetric area of each polygon
	 * @param z
	 *            : bedrock elevation
	 */
	public double[] estimateR(double[] matT, double[] arrb, double[] zetaBed,
			double[] porosity, int Np, int[] Mp, int[] Mi, double[] eta,
			double[] p) {

		double sum = 0;
		double[] arrR = new double[Np];

		for (int i = 0; i < Np; i++) {
			for (int j = Mp[i]; j < Mp[i + 1]; j++) {
				sum += matT[j] * eta[Mi[j]];
			}
			double aaa = computePAplanar(eta[i], zetaBed[i], porosity[i], p[i]);
			double vvv = computeVplanar(eta[i], zetaBed[i], aaa);
			arrR[i] = vvv + sum - arrb[i];
			// System.out.println("ArrR " + arrR.get(i));
			sum = 0;
		}

		return arrR;
	}

	/**
	 * Estimate Jr.
	 * 
	 * @desc this method needs to define the values of the matrix Jr in
	 *       Row-Compressed Form. From [Cordano Rigon 2013] equation (A6), Jr is
	 *       the sum of T (evaluated at the value of first attempt in the time
	 *       step) and the Jacobian matrix of V. The Jacobian matrix of V is a
	 *       diagonal matrix with non-null entries that are all positive, wich
	 *       correspond to the horizontal wet area of each cell. From [Casulli
	 *       2009], P is evaluted as: 1. p[i] if etaOld[i] = etaNew[i] 2.
	 *       0.5*p[i] if etaOld[i] or etaNew[i] are null 3. 0 if etaOld[i] and
	 *       etaNew[i] are null
	 * 
	 * @param Np
	 *            : number of polygons in the mesh
	 * @param Mp
	 *            : array of row pointers in Row-Compressed Form
	 * @param Mi
	 *            : array of column indices of entries in row j
	 * @param etaOld
	 *            : holds the values of the eta at the old time step
	 * @param etaNew
	 *            : holds the values of the eta at the new time step
	 * @param p
	 *            : planimetric area of each polygon
	 * @param z
	 *            : bedrock elevation
	 */
	public double[] estimateJr(int[] indexDiag, double[] matT, double[] eta,
			double[] zetaBed, double[] porosity, double[] pArea) {

		double[] matJr = new double[matT.length];

		System.arraycopy(matT, 0, matJr, 0, matT.length);

		for (int i = 0; i < indexDiag.length; i++) {

			matJr[indexDiag[i]] = matT[indexDiag[i]]
					+ computePAplanar(eta[i], zetaBed[i], porosity[i], pArea[i]);

		}

		return matJr;
	}

	public double[] newtonIteration(double[] arrb, double[] matT,
			int[] indexDiag, Grid grid1, double[] eta)  throws IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;
		int indexProva = 0;
		

		do {
			
			double[] jr = estimateJr(indexDiag, matT, eta,
					grid1.bottomElevation, grid1.porosity, grid1.planArea);
			matrixJr = new SparseRCDoubleMatrix2D(
					grid1.numberSidesPolygon.length,
					grid1.numberSidesPolygon.length, grid1.Mp, grid1.Mi, jr);

			double[] r = estimateR(matT, arrb, grid1.bottomElevation,
					grid1.porosity, grid1.numberSidesPolygon.length, grid1.Mp,
					grid1.Mi, eta, grid1.planArea);

			matrixr = new SparseDoubleMatrix1D(r);

			// System.out.println("Index prova " + indexProva);
			RCConjugateGradient cg = new RCConjugateGradient(
					grid1.numberSidesPolygon.length);
			
			cg.solverCG(matrixr, matrixJr);
			
			for (int i = 0; i < grid1.eta.length; i++) {
				eta[i] = eta[i] - cg.matSol.get(i);
				// System.out.println(cg.matSol.get(i));
				System.out.println(solNew[i]);
			}


		} while (indexProva < 10);// arrR.getMaxLocation()[0] > tol |
		return eta;

	}

	/**
	 * Newton b eq.
	 * 
	 * @param grid1
	 *            the grid1
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public void newtonBEq(Grid grid1)
			throws IterativeSolverDoubleNotConvergedException {

		// the values of first attempt are the
		System.arraycopy(grid1.eta, 0, solOld, 0, grid1.eta.length);

		for (int t = 0; t < simTime; t += deltat) {

			estimateTB(grid1, solOld);

			int indexProva = 0;

			newtonIteration();

			for (int i = 0; i < solNew.length; i++) {
				sol[i][t] = solNew[i];
			}

			System.out.println("Exit code");

		}
	}

	/**
	 * The main method.
	 * 
	 * @desc main method calls the different classes
	 * 
	 * @param args
	 *            the arguments
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException {

		Grid grid1 = new Grid();
		BoussinesqEquation beq = new BoussinesqEquation(
				grid1.numberSidesPolygon.length, grid1.Ml.length);
		beq.newtonBEq(grid1);

		System.exit(1);

	}

}

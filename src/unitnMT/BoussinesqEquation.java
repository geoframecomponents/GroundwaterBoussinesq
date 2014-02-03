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
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;

// TODO: Auto-generated Javadoc
/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation {
	
	/** The arrb. */
	double[] arrb;
	
	/** The mat t. */
	double[] matT;
	
	double[] matJr;
	
	/** The sol. */
	double[] sol;
	
	/** The arr r. */
	SparseDoubleMatrix1D arrR;
	
	/** The deltat. */
	int deltat = 1;
	
	/** legth of the simulation */
	int simTime = 2;
	
	
	
	/**
	 * Instantiates a new boussinesq equation object.
	 *
	 * @param Np: the number of polygons in the mesh
	 * @param SIZE: the number of non-zero in the Mij adjacency matrix
	 */
	BoussinesqEquation(int Np, int SIZE){
		
		System.out.println("Number of polygons:" + Np);
		System.out.println("Number of elements of T:" + SIZE);
		
		arrb = new double[Np];
		sol = new double[Np];
		matT = new double[SIZE];
		matJr = new double[SIZE];
		arrR = new SparseDoubleMatrix1D(Np);
		
	}

	
	
	/**
	 * Estimate T.
	 * 
	 * @desc This class estimates the value of T, from [Cordano Rigon 2013]
	 * 		 equations (18), (20) e (21). Because the for-loop is going,
	 * 		 it's possible to fill the b array from equation (19).
	 *
	 * @param grid1: the object grid1 holds all the properties of the mesh.
	 * 		  In this class are useful these variable:
	 * 			- numberSidesPolygon: to cycling 1..Np number of polygons
	 * 			- eta: water-table elevation (piezometric head)
	 * 			- bottomElevation: bedrock elevation
	 * 			- planArea: area of each cell
	 * 			- sourceSink: source term
	 * 			- euclideanDistance: euclidean distance between each center
	 * 				of every polygon (it's in array form)
	 * 			- hydrConductivity: saturated hydraulic conductivity as pro-
	 * 				perties of every edge
	 * 			- lengthSides: length of every edge of the polygon
	 * 			- Mp, Mi, Ml: row compressed form of adjacency matrix
	 */
	public void estimateT(Grid grid1, double[] eta){
		
		/* variable to add T terms outside the diagonal, that will be
			stored in the diagonal position*/
		double colSum = 0;
		/* to identify the diagonal entry in row-compressed form */
		int index = 0;
		
		for (int i = 0; i < grid1.numberSidesPolygon.length; i++){
			arrb[i] = (eta[i]-grid1.bottomElevation[i])*grid1.planArea[i]
					+deltat*grid1.planArea[i]*grid1.sourceSink[i];
			for (int j = grid1.Mp[i]; j < grid1.Mp[i+1]; j++){
				if (grid1.Ml[j] != -1){
					matT[j] = -deltat*(1/grid1.euclideanDistance[(int) grid1.Ml[j]])
							*grid1.hydrConductivity[(int) grid1.Ml[j]]
									*grid1.lengthSides[(int) grid1.Ml[j]]
											*Math.max(eta[grid1.Mi[j]]-grid1.bottomElevation[grid1.Mi[j]],
													eta[i]-grid1.bottomElevation[i]);
					colSum += -matT[j];
				} else {index = j;}
			}
			matT[index] = colSum;
			colSum = 0;
		}
	}
	
	
	
	/**
	 * Estimate r.
	 *
	 * @param Np the np
	 * @param Mp the mp
	 * @param eta the eta
	 * @param p the p
	 * @param z the z
	 */
	public void estimateR(int Np, int[] Mp, int[] Mi, double[] eta, double[] p, double[] z){
		double sum = 0;

		for (int i = 0; i < Np; i++){
			for (int j = Mp[i]; j < Mp[i+1]; j++){
				sum += matT[Mi[j]]*eta[Mi[j]];
			}
			arrR.set(i,p[i]*(eta[i]-z[i])+sum-arrb[i]);
			sum = 0;
		}
		
	}
	
	
	
	/**
	 * Estimate p.
	 */
	public void estimateJr(int Np, int[] Mp, int[] Mi, double[] eta, double[] p, double[] z){
	/*	for (int i = 0; i < Np; i++){
			for (int j = Mp[i]; j < Mp[i+1]; j++){
				if (j == i){
					
					matJr[j] = matT[j] + p[j]*(eta[j] - z[j]);
					
				} else {
					matJr[j] = matT[j];
				}
			}
		}*/
	}
	
	/**
	 * Newton b eq.
	 *
	 * @param grid1 the grid1
	 * @throws IterativeSolverDoubleNotConvergedException the iterative solver double not converged exception
	 */
	public void newtonBEq(Grid grid1) throws IterativeSolverDoubleNotConvergedException {
		
		System.arraycopy(grid1.eta, 0, sol, 0, grid1.eta.length);
		
		System.out.println(sol[6]);
		
		for (int t = 0; t < simTime; t += deltat){
			estimateT(grid1,sol);
			estimateR(grid1.numberSidesPolygon.length,
					grid1.Mp,
					grid1.Mi,
					grid1.topElevation,
					grid1.planArea,
					grid1.bottomElevation);
			estimateJr(grid1.numberSidesPolygon.length,
					grid1.Mp,
					grid1.Mi,
					grid1.topElevation,
					grid1.planArea,
					grid1.bottomElevation);
			while (){
				RCConjugateGradient cg = new RCConjugateGradient(grid1.numberSidesPolygon.length,
						grid1.Mp,
						grid1.Mi,matJr);
				cg.solverCG(arrR);
				for (int i = 0; i < grid1.eta.length; i++){
					sol[i] = grid1.eta[i] - cg.matSol.get(i);
					System.out.println(cg.matSol.get(i));
					System.out.println(sol[i]);
				}


				estimateR(grid1.numberSidesPolygon.length,
						grid1.Mp,
						grid1.Mi,
						sol,
						grid1.planArea,
						grid1.bottomElevation);
				estimateJr(grid1.numberSidesPolygon.length,
						grid1.Mp,
						grid1.Mi,
						sol,
						grid1.planArea,
						grid1.bottomElevation);
			}
			
		}
	}
	
	
	
	/**
	 * The main method.
	 *
	 * @param args the arguments
	 * @throws IterativeSolverDoubleNotConvergedException the iterative solver double not converged exception
	 */
	public static void main(String[] args) throws IterativeSolverDoubleNotConvergedException {
		
		Grid grid1 = new Grid();
		BoussinesqEquation beq = new BoussinesqEquation(grid1.numberSidesPolygon.length, grid1.Ml.length);
		beq.newtonBEq(grid1);
		
		
		
		System.exit(1);

	}

}

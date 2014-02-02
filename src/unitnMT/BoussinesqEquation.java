package unitnMT;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;

public class BoussinesqEquation {
	
	double[] arrb;
	double[] matT;
	SparseDoubleMatrix1D arrR;
	int deltat = 1;
	
	BoussinesqEquation(int Np, int SIZE){
		
		System.out.println("Number of polygons:" + Np);
		System.out.println("Number of elements of T:" + SIZE);
		
		arrb = new double[Np];
		matT = new double[SIZE];
		arrR = new SparseDoubleMatrix1D(Np);
		
	}

	public void estimateT(Grid grid1){
		
		double colSum = 0;
		int index = 0;
		
		for (int i = 0; i < grid1.numberSidesPolygon.length; i++){
			arrb[i] = (grid1.eta[i]-grid1.bottomElevation[i])*grid1.planArea[i]+deltat*grid1.planArea[i]*grid1.sourceSink[i];
			for (int j = grid1.Mp[i]; j < grid1.Mp[i+1]; j++){
				if (grid1.Ml[j] != -1){
					matT[j] = -deltat*(1/grid1.euclideanDistance[(int) grid1.Ml[j]])*grid1.hydrConductivity[(int) grid1.Ml[j]]*grid1.lengthSides[(int) grid1.Ml[j]]*Math.max(grid1.eta[grid1.Mi[j]+1]-grid1.bottomElevation[grid1.Mi[j]+1], grid1.eta[i+1]-grid1.bottomElevation[i+1]);
					colSum += -matT[j];
				} else {index = j;}
			}
			matT[index] = colSum;
			colSum = 0;
		}
	}
	
	public void estimateR(int Np, int[] Mp, double[] eta, double[] p, double[] z){
		double sum = 0;

		for (int i = 0; i < Np; i++){
			for (int j = Mp[i]; j < Mp[i+1]; j++){
				sum += matT[j]*eta[i];
			}
			arrR.set(i,p[i]*(eta[i]-z[i])+sum-arrb[i]);
			sum = 0;
		}
		
	}
	
	public void estimateP(){
		//calculate P
	}
	
	public void newtonBEq(Grid grid1) throws IterativeSolverDoubleNotConvergedException {
		
		estimateT(grid1);
		estimateR(grid1.numberSidesPolygon.length, grid1.Mp, grid1.topElevation, grid1.planArea, grid1.bottomElevation);
		//estimate P
		//while with R
		RCConjugateGradient cg = new RCConjugateGradient(grid1.numberSidesPolygon.length,grid1.Mp,grid1.Mi,grid1.Ml);
		cg.solverCG(arrR);
		
		estimateR(grid1.numberSidesPolygon.length, grid1.Mp, cg.matSol.toArray(), grid1.planArea, grid1.bottomElevation);
		
	}
	
	public static void main(String[] args) throws IterativeSolverDoubleNotConvergedException {
		
		Grid grid1 = new Grid();
		BoussinesqEquation beq = new BoussinesqEquation(grid1.numberSidesPolygon.length, grid1.Ml.length);
		beq.newtonBEq(grid1);
		
		

	}

}

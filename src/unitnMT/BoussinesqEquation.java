package unitnMT;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;

public class BoussinesqEquation {
	
	double[] arrb;
	double[] matT;
	int deltat = 1;
	
	BoussinesqEquation(int Np, int SIZE){
		
		arrb = new double[Np];
		matT = new double[SIZE];
		
	}

	public void estimateT(Grid grid1){
		
		double colSum = 0;
		int index = 0;
		
		for (int i = 0; i < grid1.numberSidesPolygon.length; i++){
			arrb[i] = (grid1.eta[i]-grid1.bottomElevation[i])*grid1.planArea[i]+deltat*grid1.planArea[i]*grid1.sourceSink[i];
			for (int j = grid1.Mp[i]; j < grid1.Mp[i+1]; j++){
				if (grid1.Ml[j] != -1){
					matT[j] = -deltat*(1/grid1.euclideanDistance[grid1.Ml[j]])*grid1.hydrConductivity[grid1.Ml[j]]*grid1.lengthSides[grid1.Ml[j]]*Math.max(grid1.eta[grid1.Mi[j]+1]-grid1.bottomElevation[grid1.Mi[j]+1], grid1.eta[i+1]-grid1.bottomElevation[i+1]);
					colSum += -matT[j];
				} else {index = j;}
			}
			matT[index] = colSum;
			colSum = 0;
		}
	}
	
	public static void main(String[] args) {
		
		Grid grid1 = new Grid();
		BoussinesqEquation beq = new BoussinesqEquation(grid1.lengthSides.length, grid1.Ml.length);
		beq.estimateT(grid1);
		
		

	}

}

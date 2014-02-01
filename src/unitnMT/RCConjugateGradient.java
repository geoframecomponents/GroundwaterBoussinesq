package unitnMT;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

public class RCConjugateGradient {
	
	SparseRCDoubleMatrix2D matrix_A;
	DoubleMatrix1D matrix_b;
	DoubleMatrix1D matrix_x;
	
	RCConjugateGradient(int SIZE, int[] Mp, int[] Mi, double[] Ml, DoubleMatrix1D matrix_b) throws IterativeSolverDoubleNotConvergedException{
		
		matrix_A = new SparseRCDoubleMatrix2D(SIZE,SIZE,Mp,Mi,Ml);
		
		DoubleCG conjugateGradient = new DoubleCG(matrix_x);
		DoubleMatrix1D matSol = conjugateGradient.solve(matrix_A,
				matrix_b,matrix_x);
		
	}

}

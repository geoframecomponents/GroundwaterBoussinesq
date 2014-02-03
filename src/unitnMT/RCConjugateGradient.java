package unitnMT;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.algo.solver.DoubleCG;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

// TODO: Auto-generated Javadoc
/**
 * The Class RCConjugateGradient.
 */
public class RCConjugateGradient {
	
	/** The matrix_ a. */
	SparseRCDoubleMatrix2D matrix_A;
	
	/** The matrix_x. */
	DoubleMatrix1D matrix_x;
	
	/** The mat sol. */
	DoubleMatrix1D matSol;
	
	/**
	 * Instantiates a new rC conjugate gradient.
	 *
	 * @param SIZE the size
	 * @param Mp the mp
	 * @param Mi the mi
	 * @param Ml the ml
	 */
	RCConjugateGradient (int SIZE, int[] Mp, int[] Mi, double[] Ml) {
		
		matrix_A = new SparseRCDoubleMatrix2D(SIZE,SIZE,Mp,Mi,Ml);
		matrix_x = new SparseDoubleMatrix1D(SIZE);
		
	}
	
	/**
	 * Solver cg.
	 *
	 * @param matrix_b the matrix_b
	 * @throws IterativeSolverDoubleNotConvergedException the iterative solver double not converged exception
	 */
	public void solverCG(SparseDoubleMatrix1D matrix_b) throws IterativeSolverDoubleNotConvergedException {
		
		DoubleCG conjugateGradient = new DoubleCG(matrix_x);
		matSol = conjugateGradient.solve(matrix_A,matrix_b,matrix_x);
		
	}
	
	/**
	 * The main method.
	 *
	 * @param args the arguments
	 */
	public static void main (String[] args){
		
		//Grid grid1 = new Grid();
		//RCConjugateGradient cg = new RCConjugateGradient(grid1.numberSidesPolygon.length,grid1.Mp,grid1.Mi,grid1.Ml);
		
	}

}

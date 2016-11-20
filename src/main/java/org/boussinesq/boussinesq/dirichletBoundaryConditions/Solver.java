package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AbstractRCAdjacencyMatrixBasedMesh;
import org.boussinesq.RowCompressedForm.RCConjugateGradient;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;
import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;

public class Solver {
	
	ComputeJr cJr;
	ComputeR cR;
	RCConjugateGradient cg;
	
	Solver(AbstractRCAdjacencyMatrixBasedMesh mesh){
		
		cJr = new ComputeJr();
		cR = new ComputeR();
		
		cg = new RCConjugateGradient(mesh.polygonsNumber);
		
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
	 *
	 * @return the array of eta at the following time step
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 */
	public double[] newtonIteration(double[] arrb, double[] arrT,
			int[] indexDiag, double[] eta,
			double tolerance, AbstractRCAdjacencyMatrixBasedMesh mesh) throws
			IterativeSolverDoubleNotConvergedException {

		SparseRCDoubleMatrix2D matrixJr;
		SparseDoubleMatrix1D matrixr;

		double maxResidual = 10;

		do {

			// compute Jr
			double[] jr = cJr.computeJr(indexDiag, arrT, eta,
					(CatchmentDomain) mesh);

			// convert array in sparse matrix for DoubleCG class
			matrixJr = new SparseRCDoubleMatrix2D(mesh.polygonsNumber, mesh.polygonsNumber, mesh.Mp,
					mesh.Mi, jr);

			// compute the residual function
			double[] r = cR.computeR(arrT, arrb, eta, (CatchmentDomain) mesh);

			// convert array in sparse matrix for DoubleCG class
			matrixr = new SparseDoubleMatrix1D(r);
						
			cg.solverCG(matrixr, matrixJr);
			
			maxResidual = Math.max(Math.abs(cg.matSol.getMaxLocation()[0]),
					Math.abs(cg.matSol.getMinLocation()[0]));

			// compute the new eta for every cell
			for (int i = 0; i < ((CatchmentDomain) mesh).eta.length; i++) {
				
				eta[i] = eta[i] - cg.matSol.get(i);
				
			}

			
			System.out.println("Residual: " + maxResidual);
			
			
		} while (maxResidual > tolerance * 100);

		
		return eta;

	}

}

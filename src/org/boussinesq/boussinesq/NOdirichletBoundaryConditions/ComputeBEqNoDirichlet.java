package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.PdeTermT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.PdeTermB;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq {

	double[] eta;
	double[] matT;
	double[] arrb;

	Solver newton;

	PdeTermT computeT;
	PdeTermB computeB;

	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;

	public ComputeBEqNoDirichlet(AbstractRCAdjacencyMatrixBased mesh) {

		eta = new double[mesh.polygonsNumber];
		newton = new Solver(mesh.polygonsNumber);

		computeT = new PdeTermT();
		computeB = new PdeTermB();

		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();

	}

	public void computeBEqArrays(double[] eta,
			AbstractRCAdjacencyMatrixBased mesh) {

		matT = assemblePdeTerm(eta, mesh, computeT);
		arrb = assemblePdeTerm(eta, mesh, computeB);

		double rowSum = 0;

		/* to identify the diagonal entry of matrix T in row-compressed form */
		int index = 0;

		for (int i = 0; i < mesh.polygonsNumber; i++) {

			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {

				if (mesh.Mi[j] != i) {
					// equation (21)

					rowSum += -matT[j];

				} else {
					index = j;
				}

			}
			// equation (20)
			if (rowSum == 0) {

				matT[index] = 1;

			} else {

				matT[index] = rowSum;
				rowSum = 0;

			}
		}

	}

	public double[] solutionMethod(double[] etaOld, double[] matT,
			double[] arrb, AbstractRCAdjacencyMatrixBased mesh)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, tolerance,
				mesh);
		
		for (int i = 0; i < mesh.polygonsNumber; i++){
			
			if (eta[i] < mesh.bedRockElevation[i]){
				
				eta[i] = mesh.bedRockElevation[i];
				
			}
			
		}

		return eta;

	}

	public void computeBEq(AbstractRCAdjacencyMatrixBased mesh) {
		// allocate the memory for eta array

		// firstThings(mesh);
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
				mesh.polygonsNumber, mesh.Mp, mesh.Mi);

		tolerance = cMEd.computeMachineEpsilonDouble();

		// initialize eta array
		System.arraycopy(mesh.eta, 0, eta, 0, mesh.eta.length);

		mesh.source = new double[mesh.polygonsNumber];
		mesh.outflow = new double[mesh.polygonsNumber];

		temporalLoop(mesh);

	}

	public void temporalLoop(AbstractRCAdjacencyMatrixBased mesh) {

		int contatore = 0;

		for (int t = 0; t < TimeSimulation.SIMULATIONTIME; t += TimeSimulation.TIMESTEP) {

			for (int i = 0; i < mesh.polygonsNumber; i++) {

				mesh.source[i] = mesh.rainHour[contatore];

			}

			contatore++;

			TextIO.putln("Time step " + (double) t / 3600);

			try {
				openTxtFile(t);
			} catch (IOException e1) {

				e1.printStackTrace();
			}

			computeBEqArrays(eta, mesh);

			try {
				eta = solutionMethod(eta, matT, arrb, mesh);
			} catch (IterativeSolverDoubleNotConvergedException e) {

				e.printStackTrace();
			}

			// computeOutputFeatures(eta, mesh);

			try {
				writeSolution(t, eta, mesh);
			} catch (IOException e) {

				e.printStackTrace();

			}

		}

	}

}

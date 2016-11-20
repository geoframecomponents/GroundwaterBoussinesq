package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.IOException;

import com.blogspot.geoframe.differentialEquation.partialDifferentialEquation
		.AbstractPdeTerm;
import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AbstractRCAdjacencyMatrixBasedMesh;
import com.blogspot.geoframe.mesh.unstructured.adjacencyMatrixBased
		.AdjacencyMatrixBased;
import it.blogspot.geoframe.machineEpsilon.MachineEpsilon;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeBEq;
import org.boussinesq.boussinesq.PdeTermT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.computationalDomain.CatchmentDomain;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqNoDirichlet extends ComputeBEq {

	double[] eta;
	double[] matT;
	double[] arrb;

	Solver newton;

	AbstractPdeTerm computeT;
	AbstractPdeTerm computeB;

	RCIndexDiagonalElement rcIndexDiagonalElement;

	public ComputeBEqNoDirichlet(AbstractRCAdjacencyMatrixBasedMesh mesh) {

		eta = new double[mesh.polygonsNumber];
		newton = new Solver(mesh.polygonsNumber);

		computeT = new PdeTermT();
		computeB = new PdeTermB();

		rcIndexDiagonalElement = new RCIndexDiagonalElement();

	}

	public void computeBEqArrays(double[] eta,
			AbstractRCAdjacencyMatrixBasedMesh mesh) {

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
			double[] arrb, AbstractRCAdjacencyMatrixBasedMesh mesh)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, tolerance,
				mesh);
		
		for (int i = 0; i < mesh.polygonsNumber; i++){
			
			if (eta[i] < ((CatchmentDomain) mesh).bedRockElevation[i]){
				
				eta[i] = ((CatchmentDomain) mesh).bedRockElevation[i];
				
			}
			
		}

		return eta;

	}

	public void computeBEq(AdjacencyMatrixBased mesh) {
		// allocate the memory for eta array

		// firstThings(mesh);
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
				mesh.polygonsNumber, mesh.Mp, ((CatchmentDomain) mesh).Mi);
		
		tolerance = MachineEpsilon.doublePrecision();

		// initialize eta array
		System.arraycopy(((CatchmentDomain) mesh).eta, 0, eta, 0, (
				(CatchmentDomain) mesh).eta.length);

		((CatchmentDomain) mesh).source = new double[mesh.polygonsNumber];
		((CatchmentDomain) mesh).outflow = new double[mesh.polygonsNumber];

		temporalLoop(mesh);

	}

	public void temporalLoop(AdjacencyMatrixBased mesh) {

		int contatore = 0;

		for (int t = 0; t < TimeSimulation.SIMULATIONTIME; t += TimeSimulation.TIMESTEP) {

			for (int i = 0; i < mesh.polygonsNumber; i++) {

				((CatchmentDomain) mesh).source[i] = ((CatchmentDomain) mesh).rainHour[contatore];

			}

			contatore++;

			TextIO.putln("Time step " + (double) t / 3600);

			try {
				openTxtFile(t);
			} catch (IOException e1) {

				e1.printStackTrace();
			}

			computeBEqArrays(eta, (AbstractRCAdjacencyMatrixBasedMesh) mesh);

			try {
				eta = solutionMethod(eta, matT, arrb,
						(AbstractRCAdjacencyMatrixBasedMesh) mesh);
			} catch (IterativeSolverDoubleNotConvergedException e) {

				e.printStackTrace();
			}

			// computeOutputFeatures(eta, mesh);

			try {
				writeSolution(t, eta, (AbstractRCAdjacencyMatrixBasedMesh)
						mesh);
			} catch (IOException e) {

				e.printStackTrace();

			}

		}

	}

}

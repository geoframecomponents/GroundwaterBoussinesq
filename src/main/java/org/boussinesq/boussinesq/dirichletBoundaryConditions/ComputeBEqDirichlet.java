package org.boussinesq.boussinesq.dirichletBoundaryConditions;

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

public class ComputeBEqDirichlet extends ComputeBEq {

	double[] eta;
	double[] matT;
	double[] arrb;

	double[] matTDirichlet;
	double[] matTNoDirichlet;

	double volumeDirichlet;

	Solver newton;
	AbstractPdeTerm computeT;
	AbstractPdeTerm computeB;
	ComputeTDirichlet cTDirichlet;
	ComputeTNoDirichlet cTNoDirichlet;

	int[] indexDiag;

	RCIndexDiagonalElement rcIndexDiagonalElement;
	double tolerance;

	public ComputeBEqDirichlet(AbstractRCAdjacencyMatrixBasedMesh mesh) {

		eta = new double[mesh.polygonsNumber];
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		newton = new Solver(mesh);

		computeT = new PdeTermT();
		computeB = new PdeTermB();
		cTDirichlet = new ComputeTDirichlet();
		cTNoDirichlet = new ComputeTNoDirichlet();

	}

	public void computeBEqArrays(double[] eta,
			AbstractRCAdjacencyMatrixBasedMesh mesh) {

		double rowSum = 0;

		/* to identify the diagonal entry of matrix T in row-compressed form */
		int index = 0;

		for (int i = 0; i < mesh.polygonsNumber; i++) {
			/*
			 * nested for-loop to analyze shared edges between the i-th cell and
			 * the Mi[j]-th cell
			 */
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

		matTDirichlet = cTDirichlet.computeTDirichlet(matT, (CatchmentDomain)
				mesh);
		matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT, indexDiag,
				(CatchmentDomain) mesh);

		arrb = assemblePdeTerm(eta, mesh, computeB);

		for (int i = 0; i < mesh.polygonsNumber; i++) {

			double sum = 0;
			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
				sum += matTDirichlet[j] * ((CatchmentDomain) mesh).etaDirichlet[mesh.Mi[j]];
			}

			arrb[i] = arrb[i] - sum;

		}

	}

	public double[] solutionMethod(double[] etaOld, double[] matT,
			double[] arrb, CatchmentDomain mesh)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, tolerance,
				mesh);

		for (int i = 0; i < mesh.polygonsNumber; i++) {

			if (eta[i] < mesh.bedRockElevation[i]) {

				eta[i] = mesh.bedRockElevation[i];

			}

		}

		return eta;

	}

	public void firstThings(AbstractRCAdjacencyMatrixBasedMesh mesh) {

		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
				mesh.polygonsNumber, mesh.Mp, mesh.Mi);

		tolerance = MachineEpsilon.doublePrecision();

	}

	public void computeBEq(AbstractRCAdjacencyMatrixBasedMesh mesh) {

		firstThings(mesh);

		System.arraycopy(((CatchmentDomain) mesh).eta, 0, eta, 0,
				((CatchmentDomain) mesh).eta.length);

		((CatchmentDomain) mesh).source = new double[mesh.polygonsNumber];

		temporalLoop(mesh);

	}

	public void temporalLoop(AdjacencyMatrixBased mesh) {

		EtaInitialization etaInit = new EtaInitialization();

		int contatore = 0;

		for (int t = 0; t < TimeSimulation.SIMULATIONTIME; t += TimeSimulation.TIMESTEP) {

			for (int i = 0; i < mesh.polygonsNumber; i++) {

				((CatchmentDomain) mesh).source[i] = ((CatchmentDomain) mesh).rainHour[contatore];

			}

			contatore++;

			eta = etaInit.etaInitialization(eta, (CatchmentDomain) mesh);

			TextIO.putln("Time step " + (double) t / 3600);

			try {
				openTxtFile(t);
			} catch (IOException e1) {

				e1.printStackTrace();
			}

			matT = assemblePdeTerm(eta, mesh, computeT);

			computeBEqArrays(eta, (AbstractRCAdjacencyMatrixBasedMesh) mesh);

			try {
				eta = solutionMethod(eta, matTNoDirichlet, arrb,
						(CatchmentDomain) mesh);
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

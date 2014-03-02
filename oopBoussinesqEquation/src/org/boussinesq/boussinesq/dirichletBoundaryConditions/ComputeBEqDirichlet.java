package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
//import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.Mesh;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.machineEpsilon.MachineEpsilon;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEqDirichlet extends ComputeT implements TimeSimulation {

	// allocate the memory for eta array
	double[] eta;

	double[] matT;
	double[] matTDirichlet;
	double[] matTNoDirichlet;
	double[] arrb;

	public ComputeBEqDirichlet() {
		super();

	}

	public void computeBEqArrays(ComputeTDirichlet cTDirichlet,
			ComputeTNoDirichlet cTNoDirichlet, ComputeB cB) {

		matT = computeT(eta);
		matTDirichlet = cTDirichlet.computeTDirichlet(matT);
		matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT);
		arrb = cB.computeB(eta, matTDirichlet);

	}

	public void computeBEq() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		EtaInitialization etaInit = new EtaInitialization();

		ComputeTDirichlet cTDirichlet = new ComputeTDirichlet();
		ComputeTNoDirichlet cTNoDirichlet = new ComputeTNoDirichlet();
		ComputeB cB = new ComputeB();
		Solver newton = new Solver();

		RCConjugateGradient cg = new RCConjugateGradient(Mesh.Np);
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();

		int[] indexDiag = rcIndexDiagonalElement.computeIndexDiag(Mesh.Np,
				Mesh.Mp, Mesh.Mi);

		MachineEpsilon cMEd = new MachineEpsilon();
		double tolerance = cMEd.computeMachineEpsilonDouble();

		// allocate the memory for eta array
		eta = new double[Mesh.Np];

		FileWriter Rstatfile = new FileWriter(Mesh.outputPathBeqDirichlet);
		PrintWriter errestat = new PrintWriter(Rstatfile);

		System.arraycopy(Mesh.eta, 0, eta, 0, Mesh.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {

			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(cTDirichlet, cTNoDirichlet, cB);

			eta = newton.newtonIteration(arrb, matTNoDirichlet, indexDiag, eta,
					cg, tolerance);

			System.out.println("Simulation time: " + t / 3600);

		}
		for (int j = 0; j < eta.length; j++) {

			errestat.println(eta[j]);

		}

		errestat.println();
		System.out.println();
		Rstatfile.close();

		System.out.println("Exit code");
	}
}

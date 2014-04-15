package org.boussinesq.boussinesq.NOdirichletBoundaryConditions;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.Mesh;
import org.boussinesq.boussinesq.TimeSimulation;
//import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeB;
import org.boussinesq.machineEpsilon.MachineEpsilon;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class ComputeBEq extends ComputeT implements TimeSimulation {

	// allocate the memory for eta array
	double[] eta;

	double[] matT;
	double[] arrb;

	public ComputeBEq() {
		super();
	}

	public void computeBEqArrays(ComputeB cB) {

		matT = computeT(eta);
		arrb = cB.computeB(eta);

	}

	public void computeBEq() throws IOException,
			IterativeSolverDoubleNotConvergedException {

		ComputeB cB = new ComputeB();
		Solver newton = new Solver();

		RCConjugateGradient cg = new RCConjugateGradient(Mesh.Np);
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();

		int[] indexDiag = rcIndexDiagonalElement.computeIndexDiag(Mesh.Np,
				Mesh.Mp, Mesh.Mi);

		MachineEpsilon cMEd = new MachineEpsilon();
		double tolerance = cMEd.computeMachineEpsilonDouble();

		FileWriter Rstatfile = new FileWriter(Mesh.outputPathBeqNoDirichlet);
		PrintWriter errestat = new PrintWriter(Rstatfile);

		// allocate the memory for eta array
		eta = new double[Mesh.Np];

		// initialize eta array
		System.arraycopy(Mesh.eta, 0, eta, 0, Mesh.eta.length);

		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {

			computeBEqArrays(cB);

			eta = newton.newtonIteration(arrb, matT, indexDiag, eta, cg,
					tolerance);

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

package org.boussinesq.boussinesq;

import java.io.IOException;

import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeBEq;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.ComputeBEqDirichlet;
import org.boussinesq.song.Song;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

public class BoussinesqEquation implements TimeSimulation {

	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException, IOException {
		
		String simulationType = "Song";
		String boundaryConditions = "Dirichlet";
		// long start=System.nanoTime();
		@SuppressWarnings("unused")
		Mesh mesh = new Mesh("Song");

		if (boundaryConditions.equals("Dirichlet")) {

			ComputeBEqDirichlet cBEqD = new ComputeBEqDirichlet();
			cBEqD.computeBEq();

		} else {

			ComputeBEq cBEq = new ComputeBEq();
			cBEq.computeBEq();

		}

		if (simulationType.equals("Song")) {

			Song s = new Song(SIMULATIONTIME, Mesh.Np, Mesh.hydrConductivity[0]);

			s.beqSong(Mesh.porosity);

			// long end=System.nanoTime();
			// System.out.println("End time: " + (end-start));

			System.exit(1);
		}

	}

}

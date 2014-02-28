package org.boussinesq.boussinesq;

/**
 * Mass-Conservative groundwater equation integration
 * 
 * @desc	In this code is implemented a conservative finite-volume
 * 			numerical solution for the two-dimensional groundwater flow
 * 			(Boussinesq) equation, which can be used for investigation
 * 			of hillslope subsurface flow processes and simulations of
 * 			catchment hydrology.
 * 
 * 			The idea is taken from:
 * 
 * 			"A mass-conservative method for the integration of the
 * 			 two dimensional groundwater (Boussinesq) equation"
 * 			E.Cordano, R.Rigon 2012 - Water Resources Research
 * 
 * @author	E. Cordano, G. Formetta, R. Rigon, F. Serafin, 2014
 * Copyright GPL v. 3 (http://www.gnu.org/licenses/gpl.html)
 * */

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.boussinesq.RowCompressedForm.*;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.boussinesq.song.Song;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;
//import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix1D;
//import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix2D;

/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation {

	/** The deltat. */
	static int deltat = 3600;

	/** legth of the simulation */
	static int simTime = 3600 * 24;

	static double tolerance;
	
	/**
	 * Compute the Boussinesq Equation.
	 * 
	 * @desc in this method the temporal loop is implemented. Before start the
	 *       loop, the eta array is initialized with eta of Dirichlet if the
	 *       cell is a Dirichlet cells, otherwise it's inizialized with first
	 *       attempt value.
	 * 
	 * @param mesh
	 *            the object mesh is passed so every field of the mesh class is
	 *            available
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void computeBEq(String boundaryConditions)
			throws IterativeSolverDoubleNotConvergedException, IOException {

		// allocate the memory for eta array
		double[] eta = new double[Mesh.Np];
		// new conjugate gradient object
		RCConjugateGradient cg = new RCConjugateGradient(Mesh.Np);
		RCIndexDiagonalElement rcIndexDiagonalElement = new RCIndexDiagonalElement();
		
		int[] indexDiag = rcIndexDiagonalElement.computeIndexDiag(Mesh.Np, Mesh.Mp, Mesh.Mi);
		
		tolerance = MachineEpsilon.computeMachineEpsilonDouble();

		if (boundaryConditions.equals("Dirichlet")) {

			SetDirichletBC dirichletBC = new SetDirichletBC();

			FileWriter Rstatfile = new FileWriter(Mesh.outputPathBeqDirichlet);
			PrintWriter errestat = new PrintWriter(Rstatfile);

			System.arraycopy(Mesh.eta, 0, eta, 0, Mesh.eta.length);
			
			for (int t = 0; t < simTime; t += deltat) {

				// initialize eta array
				for (int i = 0; i < eta.length; i++) {
					if (dirichletBC.isNoValue(Mesh.etaDirichlet[i],
							Mesh.NOVALUE)) {

						// not Dirichlet cells
						eta[i] = eta[i];
					} else {

						// Dirichlet cells
						eta[i] = Mesh.etaDirichlet[i];
					}
				}

				double[] matT = ComputeT.computeT(eta);
				double[] matTDrichelet = dirichletBC.computeTDirichlet(matT);
				double[] matTNoDrichelet = dirichletBC
						.computeTNoDirichlet(matT);

				double[] arrb = dirichletBC.computeB(eta, matTDrichelet,
						Mesh.etaDirichlet);

				eta = dirichletBC.newtonIteration(arrb, matTNoDrichelet,
						indexDiag, eta, cg);

				System.out.println("Simulation time: " + t / 3600);

			}
			for (int j = 0; j < eta.length; j++) {

				errestat.println(eta[j]);

			}

			errestat.println();
			System.out.println();
			Rstatfile.close();

			System.out.println("Exit code");

		} else {

			FileWriter Rstatfile = new FileWriter(Mesh.outputPathBeqNoDirichlet);
			PrintWriter errestat = new PrintWriter(Rstatfile);

			SetNoDirichletBC noDirichletBC = new SetNoDirichletBC();

			// initialize eta array
			System.arraycopy(Mesh.eta, 0, eta, 0, Mesh.eta.length);

			for (int t = 0; t < simTime; t += deltat) {

				double[] matT = ComputeT.computeT(eta);

				double[] arrb = noDirichletBC.computeB(eta);

				eta = noDirichletBC.newtonIteration(arrb, matT, indexDiag, eta,
						cg);

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

	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException, IOException {
		String simulationType = "Song";
		String boundaryConditions = "Dirichlet";
		// long start=System.nanoTime();
		@SuppressWarnings("unused")
		Mesh mesh = new Mesh(simulationType);
		BoussinesqEquation beq = new BoussinesqEquation();
		beq.computeBEq(boundaryConditions);

		if (simulationType == "Song") {
			double[] songSol = new double[org.boussinesq.boussinesq.Mesh.Np];
			Song s = new Song(simTime, Mesh.Np, Mesh.hydrConductivity[0]);

			songSol = s.beqSong(Mesh.porosity);

			FileWriter Rstatfile = new FileWriter(
					org.boussinesq.boussinesq.Mesh.outputPathSong);
			PrintWriter errestat = new PrintWriter(Rstatfile);

			for (int j = 0; j < songSol.length; j++) {

				errestat.println(songSol[j]);

			}

			errestat.println();
			System.out.println();
			Rstatfile.close();

		}
		// long end=System.nanoTime();
		// System.out.println("End time: " + (end-start));

		System.exit(1);

	}

}

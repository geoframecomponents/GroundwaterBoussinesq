package org.boussinesq.boussinesq;

import java.io.File;
import java.io.IOException;
import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeBEqNoDirichlet;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.ComputeBEqDirichlet;
import org.boussinesq.song.Song;
import org.francescoS.GUI.SelectOptions;
import org.francescoS.usefulClasses.GUIpathFileRead;
import org.francescoS.usefulClasses.TextIO;
import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation implements TimeSimulation {

	/** The boundary conditions. */
	String boundaryConditions;

	/** Directory of the solution. */
	public static File solutionDir;

	/**
	 * Define simulation type.
	 * 
	 * @desc defineSimulationType this method instantiate a new object of type
	 *       SelectOptions. That object create a new Window where select the
	 *       type of simulation. Once the simulation is defined, the
	 *       corresponding computational domain is loaded
	 * 
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void defineSimulationType() throws IOException {

		// prepare a new GUI
		SelectOptions viewOption = new SelectOptions();
		// add the different options
		viewOption.addOptions("Song simulation");
		viewOption.addOptions("Catchment basin simulation");
		// view the GUI and select between the option
		viewOption.showComboboxDemo();

		// load the corresponding computational domain
		if (SelectOptions.name.equals("Song simulation")) {

			ComputationalDomain.callSongDomain();

			// run the Song analytical solution
			Song s = new Song(SIMULATIONTIME, ComputationalDomain.Np,
					ComputationalDomain.hydrConductivity[0]);
			s.beqSong(ComputationalDomain.porosity);

		} else {

			ComputationalDomain.callCatchmentDomain();
		}

	}

	/**
	 * Define solution print location.
	 * 
	 * @return the file
	 */
	public File defineSolutionPrintLocation() {

		GUIpathFileRead guiDir = new GUIpathFileRead();
		File path = guiDir.saveDialog("Input path of BEq solution");

		return path;
	}

	/**
	 * Define boundary conditions type.
	 * 
	 * @param beq
	 *            the beq
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void defineBoundaryConditionsType(BoussinesqEquation beq)
			throws IterativeSolverDoubleNotConvergedException, IOException {

		beq.boundaryConditions = "NoDirichlet";

		for (int i = 0; i < ComputationalDomain.etaDirichlet.length; i++) {

			if (ComputationalDomain.etaDirichlet[i] != ComputationalDomain.NOVALUE) {

				System.out.println(ComputationalDomain.etaDirichlet[i]);

				beq.boundaryConditions = "Dirichlet";
				break;
			}

		}

		TextIO.putln("Simulation boundary conditions: "
				+ beq.boundaryConditions);

		if (beq.boundaryConditions.equals("Dirichlet")) {

			ComputeBEqDirichlet cBEqD = new ComputeBEqDirichlet();
			cBEqD.computeBEq();

		} else {

			ComputeBEqNoDirichlet cBEq = new ComputeBEqNoDirichlet();
			cBEq.computeBEq();

		}

	}

	/**
	 * The main method.
	 * 
	 * @param args
	 *            the arguments
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException, IOException {

		BoussinesqEquation beq = new BoussinesqEquation();
		beq.defineSimulationType();
		solutionDir = beq.defineSolutionPrintLocation();
		beq.defineBoundaryConditionsType(beq);

		System.exit(1);

	}

}

package org.boussinesq.boussinesq;

import java.io.File;
import java.io.IOException;

import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeBEqNoDirichlet;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.ComputeBEqDirichlet;
import org.boussinesq.song.Song;
import org.wordpress.growworkinghard.GUI.SelectOptions;
import org.wordpress.growworkinghard.usefulClasses.GUIpathFileRead;
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

/**
 * Boussinesq Equation class
 * 
 * @desc	This class is the MAIN CLASS to compute the piezometric head in 
 * 			every cell of a domain. This program implements the boussinesq
 * 			equation to obtain the result
 * 
 * @author	F. Serafin, 2014
 * Copyright GPL v. 3 (http://www.gnu.org/licenses/gpl.html)
 * */
public class BoussinesqEquation implements TimeSimulation {

	/** The boundary conditions. */
	public static String boundaryConditions;
	
	String simulationType;

	/** Directory of the solution. */
	public static File solutionDir;

	
	
	
	
	
	
	
	
	
	/**
	 * Define simulation type.
	 * 
	 * @desc this method instantiate a new object of type SelectOptions. That
	 *       object create a new Window where select the type of simulation.
	 *       Once the simulation is defined, the corresponding computational
	 *       domain is loaded
	 * 
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void defineSimulationType() throws IOException {

		SelectOptions viewOption = new SelectOptions();// prepare new GUI

		// add the different options
		viewOption.addOptions("Song simulation");
		viewOption.addOptions("Catchment basin simulation");

		// view the GUI and select between the option
		viewOption.showComboboxDemo();

		simulationType = SelectOptions.name;
		
		// load the corresponding computational domain
		if (SelectOptions.name.equals("Song simulation")) {

			ComputationalDomain.callSongDomain();

			// run the Song analytical solution
			Song s = new Song(SIMULATIONTIME, ComputationalDomain.Np,
					ComputationalDomain.hydrConductivity[0]);
			s.beqSong(ComputationalDomain.porosity);

		} else {

			// load data of catchment basin
			ComputationalDomain.callCatchmentDomain();
		}

	}

	
	
	
	
	
	
	
	
	
	/**
	 * Define solution print location.
	 * 
	 * @desc this method instantiate new object of type GUIpathFileRead. That
	 *       object open a window manager where you can select the path to save
	 *       data. Like a complete windows manager you can create new folder,
	 *       rename the folder and move into your hard disk.
	 * 
	 * @return the path where to save files with value at every time step
	 */
	public File defineSolutionPrintLocation() {

		GUIpathFileRead guiDir = new GUIpathFileRead();
		File path = guiDir.saveDialog("Input path of BEq solution");

		return path;
	}

	
	
	
	
	
	
	
	
	
	/**
	 * Define boundary conditions type.
	 * 
	 * @desc this method define boundary conditions type, observing the array of
	 *       the eta of Dirichlet cells. In fact if this array is not null the
	 *       compute of the solution follows BEqDirichlet code, slower than
	 *       BEqNoDirichlet code.
	 * 
	 * @param beq
	 *            the object boussinesq equation is composed by two variables:
	 *            the type of boundary condition and the path where save the
	 *            solution
	 * @throws IterativeSolverDoubleNotConvergedException
	 *             the iterative solver double not converged exception
	 * @throws IOException
	 *             Signals that an I/O exception has occurred.
	 */
	public void defineBoundaryConditionsType(BoussinesqEquation beq)
			throws IterativeSolverDoubleNotConvergedException, IOException {

		// initialize BC like case without Dirichlet cells
		BoussinesqEquation.boundaryConditions = "NoDirichlet";

		// search a Dirichlet cell into the array of eta of Dirichlet
		for (int i = 0; i < ComputationalDomain.etaDirichlet.length; i++) {

			if (ComputationalDomain.etaDirichlet[i] != ComputationalDomain.NOVALUE) {

				System.out.println(ComputationalDomain.etaDirichlet[i]);
				BoussinesqEquation.boundaryConditions = "Dirichlet";
				break;// go out the loop at the first Dirichlet cell
			}

		}

		// the simulation type is shown by the video output
		TextIO.putln("Simulation boundary conditions: "
				+ BoussinesqEquation.boundaryConditions);

		// choose the type of simulation at run time
		if (BoussinesqEquation.boundaryConditions.equals("Dirichlet")) {

			ComputeBEqDirichlet cBEqD = new ComputeBEqDirichlet();
			cBEqD.computeBEq(BoussinesqEquation.boundaryConditions,beq.simulationType);

		} else {

			ComputeBEqNoDirichlet cBEq = new ComputeBEqNoDirichlet();
			cBEq.computeBEq(BoussinesqEquation.boundaryConditions,beq.simulationType);

		}

	}

	
	
	
	
	
	
	
	
	
	/**
	 * The main method.
	 * 
	 * @desc the main calls only the methods useless to start the simulation
	 * 			1- define the simulation type (song or catchment)
	 * 			2- define the location where save the results
	 * 			3- define the type of simulation that must be run
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

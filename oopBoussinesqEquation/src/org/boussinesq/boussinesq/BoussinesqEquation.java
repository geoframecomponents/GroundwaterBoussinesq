package org.boussinesq.boussinesq;

//import java.io.File;
import java.io.File;
import java.io.IOException;
//import java.util.Arrays;







import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeBEq;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.ComputeBEqDirichlet;
import org.boussinesq.song.Song;
import org.francescoS.usefulClasses.GUIpathFileRead;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation implements TimeSimulation {

	String boundaryConditions;
	public static File solutionPath;
	public static File solutionDir;
	
	public void defineBoundaryConditionsType(BoussinesqEquation beq){
		
		beq.boundaryConditions = "NoDirichlet";
		
		for (int i = 0; i < ComputationalDomain.etaDirichlet.length; i++){
			
			if (ComputationalDomain.etaDirichlet[i] != ComputationalDomain.NOVALUE){
				
				System.out.println(ComputationalDomain.etaDirichlet[i]);
				
				beq.boundaryConditions = "Dirichlet";
				break;
			}
			
		}
		
		TextIO.putln("Simulation boundary conditions: " + beq.boundaryConditions);
		
	}
	
	/**
	 * The main method.
	 *
	 * @param args the arguments
	 * @throws IterativeSolverDoubleNotConvergedException the iterative solver double not converged exception
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException, IOException {
		
		String sep = System.getProperty("file.separator");
				

		String simulationType = "Song";

		// long start=System.nanoTime();
		
//		solutionDir = FileWrite.makeDirectory(ReadFromScreen.readText("Write the name of the solution folder\n(it's better without space)"));
		
		GUIpathFileRead guiDir = new GUIpathFileRead();
		
		solutionPath = guiDir.saveDialog();
		
		solutionDir = new File(solutionPath, sep);
		
		BoussinesqEquation beq = new BoussinesqEquation();
		
		if (simulationType.equals("Song")){
			
			ComputationalDomain.callSongDomain();
			
		} else {
			
			ComputationalDomain.callCatchmentDomain();
		}

		beq.defineBoundaryConditionsType(beq);
		
		if (beq.boundaryConditions.equals("Dirichlet")) {

			ComputeBEqDirichlet cBEqD = new ComputeBEqDirichlet();
			cBEqD.computeBEq();

		} else {

			ComputeBEq cBEq = new ComputeBEq();
			cBEq.computeBEq();

		}

		if (simulationType.equals("Song")) {

			Song s = new Song(SIMULATIONTIME, ComputationalDomain.Np, ComputationalDomain.hydrConductivity[0]);

			s.beqSong(ComputationalDomain.porosity);

			// long end=System.nanoTime();
			// System.out.println("End time: " + (end-start));

			
		}

		System.exit(1);
		
	}

}

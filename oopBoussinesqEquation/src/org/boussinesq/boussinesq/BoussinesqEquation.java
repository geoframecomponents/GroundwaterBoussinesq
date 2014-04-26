package org.boussinesq.boussinesq;

import java.io.File;
import java.io.IOException;

import org.boussinesq.boussinesq.NOdirichletBoundaryConditions.ComputeBEqNoDirichlet;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.ComputeBEqDirichlet;
import org.boussinesq.song.Song;
import org.francescoS.GUI.selectOptions;
import org.francescoS.usefulClasses.GUIpathFileRead;
import org.francescoS.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

/**
 * The Class BoussinesqEquation.
 */
public class BoussinesqEquation implements TimeSimulation {

	String boundaryConditions;
	public static File solutionDir;
	
	
	public void defineSimulationType() throws IOException{
		
		selectOptions viewOption = new selectOptions();
		viewOption.addOptions("Song simulation");
		viewOption.addOptions("Catchment basin simulation");
		viewOption.showComboboxDemo();
		
		if (selectOptions.name.equals("Song simulation")){
			
			ComputationalDomain.callSongDomain();
			
			Song s = new Song(SIMULATIONTIME, ComputationalDomain.Np, ComputationalDomain.hydrConductivity[0]);
			s.beqSong(ComputationalDomain.porosity);
			
		} else {
			
			ComputationalDomain.callCatchmentDomain();
		}
		
	}
	
	public File defineSolutionPrintLocation(){
		
		GUIpathFileRead guiDir = new GUIpathFileRead();
		File path = guiDir.saveDialog("Input path of BEq solution");
				
		return path;
	}
	
	public void defineBoundaryConditionsType(BoussinesqEquation beq) throws IterativeSolverDoubleNotConvergedException, IOException{
		
		beq.boundaryConditions = "NoDirichlet";
		
		for (int i = 0; i < ComputationalDomain.etaDirichlet.length; i++){
			
			if (ComputationalDomain.etaDirichlet[i] != ComputationalDomain.NOVALUE){
				
				System.out.println(ComputationalDomain.etaDirichlet[i]);
				
				beq.boundaryConditions = "Dirichlet";
				break;
			}
			
		}
		
		TextIO.putln("Simulation boundary conditions: " + beq.boundaryConditions);
		
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
	 * @param args the arguments
	 * @throws IterativeSolverDoubleNotConvergedException the iterative solver double not converged exception
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public static void main(String[] args)
			throws IterativeSolverDoubleNotConvergedException, IOException {
		
		BoussinesqEquation beq = new BoussinesqEquation();
		
		beq.defineSimulationType();
		
		solutionDir = beq.defineSolutionPrintLocation();

		// long start=System.nanoTime();
		
		beq.defineBoundaryConditionsType(beq);

		System.exit(1);
		
	}

}

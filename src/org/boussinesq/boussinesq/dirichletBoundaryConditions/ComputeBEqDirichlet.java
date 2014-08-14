package org.boussinesq.boussinesq.dirichletBoundaryConditions;

import java.io.IOException;
<<<<<<< HEAD
import java.text.DecimalFormat;
=======
>>>>>>> thesis_structure

import org.boussinesq.RowCompressedForm.RCConjugateGradient;
import org.boussinesq.RowCompressedForm.RCIndexDiagonalElement;
import org.boussinesq.boussinesq.ComputeBEq;
<<<<<<< HEAD
import org.boussinesq.boussinesq.ComputeT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.Solver;
import org.boussinesq.boussinesq.computationalDomain.ComputationalDomain;
import org.boussinesq.machineEpsilon.MachineEpsilon;
=======
import org.boussinesq.boussinesq.PdeTermT;
import org.boussinesq.boussinesq.TimeSimulation;
import org.boussinesq.boussinesq.dirichletBoundaryConditions.Solver;
import org.boussinesq.machineEpsilon.MachineEpsilon;
import org.meshNumericalMethods.unstructuredMesh.adjacencyMatrixBased.AbstractRCAdjacencyMatrixBased;
>>>>>>> thesis_structure
import org.wordpress.growworkinghard.usefulClasses.TextIO;

import cern.colt.matrix.tdouble.algo.solver.IterativeSolverDoubleNotConvergedException;

<<<<<<< HEAD
public class ComputeBEqDirichlet extends ComputeBEq implements TimeSimulation {
=======
public class ComputeBEqDirichlet extends ComputeBEq {
>>>>>>> thesis_structure

	double[] eta;
	double[] matT;
	double[] arrb;
<<<<<<< HEAD
	
=======

>>>>>>> thesis_structure
	double[] matTDirichlet;
	double[] matTNoDirichlet;

	double volumeDirichlet;

	Solver newton;
	RCConjugateGradient cg;
<<<<<<< HEAD
	ComputeT computeT;
	ComputeTDirichlet cTDirichlet;
	ComputeTNoDirichlet cTNoDirichlet;
	
	
=======
	PdeTermT computeT;
	PdeTermB computeB;
	ComputeTDirichlet cTDirichlet;
	ComputeTNoDirichlet cTNoDirichlet;

>>>>>>> thesis_structure
	int[] indexDiag;

	RCIndexDiagonalElement rcIndexDiagonalElement;
	MachineEpsilon cMEd;
	double tolerance;
<<<<<<< HEAD
	
	public ComputeBEqDirichlet() {
		
		eta = new double[ComputationalDomain.Np];
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
		newton = new Solver();
		cg = new RCConjugateGradient(ComputationalDomain.Np);
		
		computeT = new ComputeT();
		cTDirichlet = new ComputeTDirichlet();
		cTNoDirichlet = new ComputeTNoDirichlet();
	

	}
	

	
	
	
	
	public void computeBEqArrays(double[] eta) {

		
		matT = computeT.computeT(eta);
		matTDirichlet = cTDirichlet.computeTDirichlet(matT);
		matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT, indexDiag);
		ComputeB cB = new ComputeB();
		arrb = cB.computeB(eta, matTDirichlet);

	}
	
	
	
	
	
	
	
	
	public double[] solutionMethod(double[] etaOld, double[] matT, double[] arrb) throws IterativeSolverDoubleNotConvergedException{
		
		double[] eta = new double[etaOld.length];
		
		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, cg,
				tolerance);
		
		return eta;
		
	}

	
	
	
	public void firstThings(){
		
		indexDiag = rcIndexDiagonalElement.computeIndexDiag(ComputationalDomain.Np,
				ComputationalDomain.Mp, ComputationalDomain.Mi);

		tolerance = cMEd.computeMachineEpsilonDouble();
		
	}
	
	
	
	
	
	
	public void computeBEq(String boundaryCondition, String simulationType) throws IOException,
			IterativeSolverDoubleNotConvergedException {

		firstThings();
		
		EtaInitialization etaInit = new EtaInitialization();

		computeInitialVolume();

		System.arraycopy(ComputationalDomain.eta, 0, eta, 0, ComputationalDomain.eta.length);

//		DecimalFormat df = new DecimalFormat("#.##");
		
		int contatore = 0;
		
		for (int t = 0; t < SIMULATIONTIME; t += TIMESTEP) {
			
			for (int i = 0; i<ComputationalDomain.Np; i++ ){
				
				ComputationalDomain.source[i] = ComputationalDomain.rainHour[contatore];
				
			}
			
//			t = Double.parseDouble(df.format(t));
			
			eta = etaInit.etaInitialization(eta);

			computeBEqArrays(eta);

			eta = solutionMethod(eta, matTNoDirichlet, arrb);
					
			computeOutputFeatures(eta);
			
//			if (((int)Math.round(t*100)) % 100 == 0){
				
				TextIO.putln("Time step " + (double) t/3600);
				
				openTxtFile(t);
				writeSolution(t, eta, boundaryCondition, simulationType);
				
//			}
			
			computeVolumeConservation();
			
			contatore++;
			
=======

	public ComputeBEqDirichlet(AbstractRCAdjacencyMatrixBased mesh) {

		eta = new double[mesh.polygonsNumber];
		rcIndexDiagonalElement = new RCIndexDiagonalElement();
		cMEd = new MachineEpsilon();
		newton = new Solver();
		cg = new RCConjugateGradient(mesh.polygonsNumber);

		computeT = new PdeTermT();
		computeB = new PdeTermB();
		cTDirichlet = new ComputeTDirichlet();
		cTNoDirichlet = new ComputeTNoDirichlet();

	}

	public void computeBEqArrays(double[] eta,
			AbstractRCAdjacencyMatrixBased mesh) {
		
		matT = assemblePdeTerm(eta, mesh, computeT);

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


		matTDirichlet = cTDirichlet.computeTDirichlet(matT, mesh);
		matTNoDirichlet = cTNoDirichlet.computeTNoDirichlet(matT, indexDiag,
				mesh);

		arrb = assemblePdeTerm(eta, mesh, computeB);

		for (int i = 0; i < mesh.polygonsNumber; i++) {

			double sum = 0;
			for (int j = mesh.Mp[i]; j < mesh.Mp[i + 1]; j++) {
				sum += matTDirichlet[j] * mesh.etaDirichlet[mesh.Mi[j]];
			}

			arrb[i] = arrb[i] - sum;

		}

	}

	public double[] solutionMethod(double[] etaOld, double[] matT,
			double[] arrb, AbstractRCAdjacencyMatrixBased mesh)
			throws IterativeSolverDoubleNotConvergedException {

		double[] eta = new double[etaOld.length];

		eta = newton.newtonIteration(arrb, matT, indexDiag, etaOld, cg,
				tolerance, mesh);

		return eta;

	}

	public void firstThings(AbstractRCAdjacencyMatrixBased mesh) {

		indexDiag = rcIndexDiagonalElement.computeIndexDiag(
				mesh.polygonsNumber, mesh.Mp, mesh.Mi);

		tolerance = cMEd.computeMachineEpsilonDouble();

	}

	public void computeBEq(AbstractRCAdjacencyMatrixBased mesh) {

		firstThings(mesh);

		System.arraycopy(mesh.eta, 0, eta, 0, mesh.eta.length);

		temporalLoop(mesh);

	}

	public void temporalLoop(AbstractRCAdjacencyMatrixBased mesh) {

		int contatore = 0;

		EtaInitialization etaInit = new EtaInitialization();

		for (int t = 0; t < TimeSimulation.SIMULATIONTIME; t += TimeSimulation.TIMESTEP) {

			// for (int i = 0; i < mesh.polygonsNumber; i++) {
			//
			// mesh.source[i] = mesh.rainHour[contatore];
			//
			// }

			// t = Double.parseDouble(df.format(t));

			eta = etaInit.etaInitialization(eta, mesh);

			TextIO.putln("Time step " + (double) t / 3600);

			try {
				openTxtFile(t);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			computeBEqArrays(eta, mesh);

			try {
				eta = solutionMethod(eta, matTNoDirichlet, arrb, mesh);
			} catch (IterativeSolverDoubleNotConvergedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			computeOutputFeatures(eta, mesh);

			try {
				writeSolution(t, eta, mesh);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			contatore++;

>>>>>>> thesis_structure
		}

	}
}

package org.boussinesq.machineEpsilon;

public class MachineEpsilon {

	/**
	 * Calculate machine epsilon double.
	 * 
	 * @desc this method compute the tolerance of the machine. For more info go
	 *       to
	 *       https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java
	 *       . In c/c++ section there's write that:
	 * 
	 *       In such languages as C or C++ when you do something like while( 1.0
	 *       + eps > 1.0 ) the expression is calculated not with 64 bits
	 *       (double) but with the processor precision (80 bits or more depends
	 *       on the processor and compile options). Below program calculates
	 *       exactly on 32 bits (float) and 64 bits (double)
	 * 
	 * 
	 * @return the tolerance of the machine for double
	 */
	public static double computeMachineEpsilonDouble() {

		// machine tolerance
		double machEps = 1.0d;

		do
			machEps /= 2.0d;
		while ((double) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	/**
	 * Calculate machine epsilon float.
	 * 
	 * @desc this method compute the tolerance of the machine. For more info go
	 *       to
	 *       https://en.wikipedia.org/wiki/Machine_epsilon#Approximation_using_Java
	 *       . In c/c++ section there's write that:
	 * 
	 *       In such languages as C or C++ when you do something like while( 1.0
	 *       + eps > 1.0 ) the expression is calculated not with 64 bits
	 *       (double) but with the processor precision (80 bits or more depends
	 *       on the processor and compile options). Below program calculates
	 *       exactly on 32 bits (float) and 64 bits (double)
	 * 
	 * 
	 * @return the tolerance of the machine for float
	 */
	private static float computeMachineEpsilonFloat() {
		float machEps = 1.0f;

		do
			machEps /= 2.0f;
		while ((float) (1.0 + (machEps / 2.0)) != 1.0);

		return machEps;
	}

	public static void main(String[] args) {

		double doubleTolerance;
		float floatTolerance;

		doubleTolerance = computeMachineEpsilonDouble();
		floatTolerance = computeMachineEpsilonFloat();

		System.out.println("The machine precision for double is : "
				+ doubleTolerance);
		System.out.println("The machine precision for float is : "
				+ floatTolerance);

	}

}

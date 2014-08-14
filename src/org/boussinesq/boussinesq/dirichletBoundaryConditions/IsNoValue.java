package org.boussinesq.boussinesq.dirichletBoundaryConditions;

public class IsNoValue {

	/**
	 * Checks if is no value.
	 * 
	 * @desc this method evaluates if the entry value is a NaN or not. If is
	 *       NaN, the cell I'm observing is a non Dirichlet cell, otherwise is a
	 *       Dirichlet cell
	 * 
	 * @param x
	 *            the eta of the cell that I'm analyzing
	 * @param noValue
	 *            the no value
	 * 
	 * @return boolean true if the cell isn't a Dirichlet cell, false otherwise
	 */
	public boolean isNoValue(double x, double noValue) {

		if (x <= noValue) {
			return true;
		} else {
			return false;
		}
	}
	
}

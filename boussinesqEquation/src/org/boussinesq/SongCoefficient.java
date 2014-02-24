package org.boussinesq;

import java.util.Arrays;

public class SongCoefficient {

	public static double[] CoefficientSongSolution(int nmax,double lambda){
		
		double[] a = new double[nmax];
		
		a[0] = 0.25;
		a[1] = (2*lambda-1)/16;
		
		for (int i = 2; i < nmax; i++){
			
			a[i] = (2*lambda+1-(i+1))/(Math.pow(i+1, 2))*a[i-1];
			
			for (int k = 1; k < (i-1); k++){
				
				a[i] = a[i]-2*(i+2)/(i+1)*a[k]*a[i-k+1];
				
			}
			//System.out.println(a[i]);
		}
		
		//System.out.println(Arrays.toString(a));
		
		return a;
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}

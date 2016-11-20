package song;

import org.junit.Test;

public class Song {

	double[] x;
	// double[] porosity;
	double hydraulicConductivity;
	int t;
	int alpha;
	static int h1;
	int nmax;

	@Test
	public void test() {

		int time = 3600 * 24 * 5;
		int numberOfCells = 1000;

		double[] porosity = new double[numberOfCells];

		for (int i = 0; i < numberOfCells; i++) {

			porosity[i] = 0.4;

		}

		songInit(time, numberOfCells, 0.1);
		beqSong(porosity);

	}

	public void songInit(int time, int Np, double ks) {

		x = new double[Np];
		t = time;
		alpha = 0;
		h1 = 1;
		nmax = 10;
		hydraulicConductivity = ks;

		for (int i = 0; i < x.length; i++) {

			x[i] = i;

		}

	}

	private double[] computeA(double[] ax, double xi0) {

		double[] a = new double[ax.length];

		for (int i = 0; i < a.length; i++) {

			a[i] = ax[i] * Math.pow(xi0, 2);

		}

		double sum = 0;

		for (int i = 0; i < a.length; i++) {

			sum += a[i];

		}

		System.out.println("Somma a: " + sum);

		return a;

	}

	private double[] computeXI(double[] porosity) {

		double[] xi = new double[x.length];

		for (int i = 0; i < xi.length; i++) {

			xi[i] = x[i]
					* Math.pow(
							2
									* porosity[i]
									* (alpha + 1)
									/ (h1 * hydraulicConductivity * Math.pow(t,
											alpha + 1)), 0.5);

		}

		// System.out.println(Arrays.toString(xi));
		return xi;
	}
	
	private void beqSong(double[] porosity) {

		String song = "songks";
		song = song.concat(Double.toString(hydraulicConductivity));
		song = song.concat("days").concat(Integer.toString(t/(3600*24)));
		
		double[] ax = new double[nmax];
		double[] solutionDimensionless = new double[x.length];
		double[] solution = new double[x.length];

		ax = SongCoefficient.CoefficientSongSolution(nmax,
				(alpha / (alpha + 1)));

		// System.out.println("Ax" + Arrays.toString(ax));

		double xi0 = 0;
		double sum = 0;

		for (int i = 1; i < ax.length; i++) {

			sum += ax[i];

		}

		xi0 = Math.pow(sum, -0.5);

		// System.out.println(xi0);

		solutionDimensionless = SongDimensionless.beqSongDimensionless(
				computeXI(porosity), xi0, computeA(ax, xi0));

		for (int i = 0; i < solutionDimensionless.length; i++) {

			// System.out.println(solutionDimensionless[i]);

			solution[i] = h1 * Math.pow(t, alpha) * solutionDimensionless[i];

		}

	}

}

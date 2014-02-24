package org.boussinesq;

// TODO: Auto-generated Javadoc
/**
 * The Class Grid.
 */
public class Grid {
	
	String outputPathBeq;
	String outputPathSong;
	
	//POLYGONS PROPERTIES
	
	/** The number sides polygon. */
	int[] numberSidesPolygon = {4,4,4,4,4,4,4,4,
			4,4,4,4,4,4,4,4};

	/** The plan area. */
	double[] planArea = {6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
	
	/** The source.  per unit area of the polygon*/
	double[] source = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	/** The eta. */
	double[] eta = {12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12};
	
	/** The eta. */
	double[] etaDrichelet = {12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12};
	
	double NOVALUE = -999;
	
	double[] h;
	
	/** The top elevation. */
	double[] topElevation = {16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
	
	/** The bottom elevation. */
	double[] bottomElevation = {10,10,10,10,10,10,10,10,10,10,
			10,10,10,10,10,10};
	
	double[] porosity = {0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,
			0.10,0.10,0.10,0.10,0.10,0.10};
	
	double[] c = {0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,
			0.10,0.10,0.10,0.10,0.10,0.10};
	
	double[] m = {0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,
			0.10,0.10,0.10,0.10,0.10,0.10};
	
	
	//SIDES PROPERTIES
	
	
	/** The length sides. */
	double[] lengthSides = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
	
	/** The euclidean distance. */
	double[] euclideanDistance = {1.5,1.5,1.5,
			1.5,3,3,3,1.5,
			3,3,3,3,3,
			1.5,3,3,3,1.5,
			1.5,1.5,1.5,
			1,2,2,1,
			1,2,2,2,2,1,
			1,2,2,2,2,1,
			1,2,2,1};
	
	/** The hydr conductivity. */
	double[] hydrConductivity = {10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
				10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
				10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
				10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
				10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
				10^(-5),10^(-5),10^(-5),10^(-5),10^(-5)};


	//ADJACENCY MATRIX PROPERTIES
	
	
	/** The Mp. */
	int[] Mp = {0,3,7,10,13,18,23,28,31,34,39,44,49,52,55,59,62};
	
	/** The Mi. */
	int[] Mi = {0,1,4, 0,1,2,5, 1,2,6, 3,4,8, 0,3,4,5,9, 1,4,5,6,10,
			2,5,6,7,11, 6,7,12, 3,8,9, 4,8,9,10,13, 5,9,10,11,14,
			6,10,11,12,15, 7,11,12, 9,13,14, 10,13,14,15, 11,14,15};
	
	/** The Ml. */
	double[] Ml = {-1,23,5, 23,-1,24,6, 24,-1,7, -1,27,9, 5,27,-1,28,10,
			6,28,-1,29,11, 7,29,-1,30,12, 30,-1,13, 9,-1,33,
			10,33,-1,34,15, 11,34,-1,35,16, 12,35,-1,36,17,
			13,36,-1, 15,-1,39, 16,39,-1,40, 17,40,-1};
	
	

	
	Grid(String configuration){
		
		if (configuration.equals("test0")){
			//POLYGONS PROPERTIES
			
			/** The number sides polygon. */
			this.numberSidesPolygon = new int[]{4,4,4,4,4,4,4,4,
					4,4,4,4,4,4,4,4};

			/** The plan area. */
			this.planArea = new double[]{6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
			
			/** The source.  per unit area of the polygon*/
			this.source = new double[]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
			
			/** The eta. */
			this.eta = new double[]{12,12,12,12,12,14,12,12,23,12,12,12,12,12,12,12};
			
			/** The eta. */
			this.etaDrichelet = new double[]{13,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,
					-999,-999,-999,-999};
			
			this.NOVALUE = -999;
			
			/** The top elevation. */
			this.topElevation = new double[]{16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
			
			/** The bottom elevation. */
			this.bottomElevation = new double[]{10,10,10,10,10,10,10,10,10,10,
					10,10,10,10,10,10};
			
			this.porosity = new double[]{0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,0.10,
					0.10,0.10,0.10,0.10,0.10,0.10};
			
			
			//SIDES PROPERTIES
			
			
			/** The length sides. */
			this.lengthSides = new double[]{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
					3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
			
			/** The euclidean distance. */
			this.euclideanDistance = new double[]{1.5,1.5,1.5,
					1.5,3,3,3,1.5,
					3,3,3,3,3,
					1.5,3,3,3,1.5,
					1.5,1.5,1.5,
					1,2,2,1,
					1,2,2,2,2,1,
					1,2,2,2,2,1,
					1,2,2,1};
			
			/** The hydr conductivity. */
			this.hydrConductivity = new double[]{10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
						10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
						10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
						10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
						10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),10^(-5),
						10^(-5),10^(-5),10^(-5),10^(-5),10^(-5)};


			//ADJACENCY MATRIX PROPERTIES
			
			
			/** The Mp. */
			this.Mp = new int[]{0,3,7,10,13,18,23,28,31,34,39,44,49,52,55,59,62};
			
			/** The Mi. */
			this.Mi = new int[]{0,1,4, 0,1,2,5, 1,2,6, 3,4,8, 0,3,4,5,9, 1,4,5,6,10,
					2,5,6,7,11, 6,7,12, 3,8,9, 4,8,9,10,13, 5,9,10,11,14,
					6,10,11,12,15, 7,11,12, 9,13,14, 10,13,14,15, 11,14,15};
			
			/** The Ml. */
			this.Ml =new double[]{-1,23,5, 23,-1,24,6, 24,-1,7, -1,27,9, 5,27,-1,28,10,
					6,28,-1,29,11, 7,29,-1,30,12, 30,-1,13, 9,-1,33,
					10,33,-1,34,15, 11,34,-1,35,16, 12,35,-1,36,17,
					13,36,-1, 15,-1,39, 16,39,-1,40, 17,40,-1};
		}else if (configuration.equals("Song")){
			
			
			this.outputPathBeq = "/home/francesco/desktop/beq.txt";
			this.outputPathSong = "/home/francesco/desktop/song.txt";
			int dim = 1000;
			//POLYGONS PROPERTIES
			this.numberSidesPolygon= new int[dim];
			this.planArea= new double[dim];
			this.source = new double[dim];
			this.eta = new double[dim];
			this.etaDrichelet = new double[dim];
			this.topElevation = new double[dim];
			this.bottomElevation = new double[dim];
			this.porosity = new double[dim];
			this.c = new double[dim];
			this.m = new double[dim];
			this.lengthSides = new double[dim + 1];
			this.euclideanDistance = new double[dim + 1];
			this.hydrConductivity = new double[dim+1];

			this.NOVALUE = -999;
			
			for (int i=0; i<dim;i++){
				this.numberSidesPolygon[i] = 2;
				this.planArea[i] = 1;
				this.source[i] = 0;
				this.eta[i] = 0;
				this.etaDrichelet[i] = -999;
				this.bottomElevation[i] = 0;
				this.porosity[i] = 0.4;
				this.lengthSides[i]= 1;
				this.euclideanDistance[i] = 1;
				this.hydrConductivity[i] = 0.01;
				this.c[i] = 0;
				this.m[i] = 1;
				
			}
			
			
			this.c[dim-1] = 1;
			this.m[dim-1] = 1;
			this.etaDrichelet[0] = 1;
			this.lengthSides[dim]= 1;
			this.euclideanDistance[dim] = 1;
			this.hydrConductivity[dim] = 0.01;
			
			//ADJACENCY MATRIX PROPERTIES
			
			this.Mp = new int[dim+1];
			this.Mi = new int[dim*3-2];
			this.Ml = new double[dim*3-2];
			this.Mp[0] = 0;
			this.Mp[1] = 2;
			
;			//System.out.println("Son qui");
			
			for (int i= 2; i<(dim);i++){
				this.Mp[i] = this.Mp[i-1]+3;
			}
			
			this.Mp[dim] = this.Mp[dim-1]+1;
			
			this.Mi[0] = 0;
			this.Mi[1] = 1;
			int index = 0;
			for (int i=1;i<(dim);i++){
				
				for(int j=Mp[i];j<Mp[i+1];j++){
					
					this.Mi[j]= index;
					index ++;
				}
				index = index-2;
			}
			
			//System.out.println(this.Mi.length);
			this.Mi[this.Mi.length-1]= this.Mi[this.Mi.length-2]+1;
			
			this.Ml[0] = -1;
			this.Ml[1] = 2;
			
			int ind =2;
			
			for (int i = 1; i < dim; i++) {
				/*
				 * nested for-loop to analyze diagonal entries, which are identified
				 * by a negative number
				 */
				
				for (int j = this.Mp[i]; j < this.Mp[i + 1]; j++) {

					if (this.Mi[j] == i) {
						this.Ml[j] = -1;
					}else{
						this.Ml[j] = ind;
						ind++;
					}

				}
				ind = (int) this.Ml[Mp[i+1]-1];
			}
			
			this.Ml[this.Ml.length-1] = -1;
			
			//System.out.println(Arrays.toString(this.Mp));
			//System.out.println(Arrays.toString(this.Mi));
			//System.out.println(Arrays.toString(this.Ml));
			/** The Mp. */
			//this.Mp = new int[]{0,2,5,8,11,14,17,20,23,26,27};
			
			/** The Mi. */
			//this.Mi = new int[]{0,1,0,1,2,1,2,3,2,3,4,3,4,5,4,5,6,5,6,7,6,7,8,
				//	7,8,9,8,9};
			
			/** The Ml. */
			//this.Ml =new double[]{-1,2,2,-1,3,3,-1,4,4,-1,5,5,-1,6,6,-1,7,7,
					//-1,8,8,-1,9,9,-1,10,10,-1};
		}
			
		}
		
	
	/**
	 * The main method.
	 *
	 * @param args the arguments
	 */
	public static void main (String[] args){
		
		
		//Grid gridTest = new Grid("Song");
		/*Grid gridTest = new Grid();
		
		System.out.println("");
		System.out.println("Array POLYGONS properties length");
		System.out.println("");
		
		System.out.println("numberSidesPolygon array length: "+ gridTest.numberSidesPolygon.length);
		System.out.println("planArea array length: "+ gridTest.planArea.length);
		System.out.println("sourceSink array length: "+ gridTest.source.length);
		System.out.println("eta array length: "+ gridTest.eta.length);
		System.out.println("topElevation array length: "+ gridTest.topElevation.length);
		System.out.println("bottomElevation array length: "+ gridTest.bottomElevation.length);
		
		System.out.println("");
		System.out.println("Array SIDES properties length");
		System.out.println("");
		
		System.out.println("lengthSides array length: "+ gridTest.lengthSides.length);
		System.out.println("euclideanDistance array length: "+ gridTest.euclideanDistance.length);
		System.out.println("hydrConductivity array length: "+ gridTest.hydrConductivity.length);
		
		System.out.println("");
		System.out.println("ADJACENCY MATRICES properties length");
		System.out.println("");
		
		System.out.println("Mp array length: "+ gridTest.Mp.length);
		System.out.println("Mi array length: "+ gridTest.Mi.length);
		System.out.println("Ml array length: "+ gridTest.Ml.length);*/
		
	}
	
}
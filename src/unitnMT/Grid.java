package unitnMT;

public class Grid {
	
	//POLYGONS PROPERTIES
	
	int[] numberSidesPolygon = {4,4,4,4,4,4,4,4,
			4,4,4,4,4,4,4,4};

	double[] planArea = {6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6};
	double[] sourceSink = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	
	double[] eta = {12,12,12,12,12,12,12,12,12,12,12,12,12,12,12,12};
	double[] topElevation = {16,16,16,16,16,16,16,16,16,16,16,16,16,16,16,16};
	double[] bottomElevation = {10,10,10,10,10,10,10,10,10,10,
			10,10,10,10,10,10};
	
	
	//SIDES PROPERTIES
	
	
	double[] lengthSides = {2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
			3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
	double[] euclideanDistance = {1.5,1.5,1.5,
			1.5,3,3,3,1.5,
			3,3,3,3,3,
			1.5,3,3,3,1.5,
			1.5,1.5,1.5,
			1,2,2,1,
			1,2,2,2,2,1,
			1,2,2,2,2,1,
			1,2,2,1};
	double[] hydrConductivity = {10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),
				10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),
				10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),
				10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),
				10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),10^(-4),
				10^(-4),10^(-4),10^(-4),10^(-4),10^(-4)};


	//ADJACENCY MATRIX PROPERTIES
	
	
	int[] Mp = {0,3,7,10,13,18,23,28,31,34,39,44,49,52,55,59,62};
	int[] Mi = {0,1,4, 0,1,2,5, 1,2,6, 3,4,8, 0,3,4,5,9, 1,4,5,6,10,
			2,5,6,7,11, 6,7,12, 3,8,9, 4,8,9,10,13, 5,9,10,11,14,
			6,10,11,12,15, 7,11,12, 9,13,14, 10,13,14,15, 11,14,15};
	double[] Ml = {-1,23,5, 23,-1,24,6, 24,-1,7, -1,27,9, 5,27,-1,28,10,
			6,28,-1,29,11, 7,29,-1,30,12, 30,-1,13, 9,-1,33,
			10,33,-1,34,15, 11,34,-1,35,16, 12,35,-1,36,17,
			13,36,-1, 15,-1,39, 16,39,-1,40, 17,40,-1};
	
	

	
	//Grid(){}
	
	public static void main (String[] args){
		
		Grid gridTest = new Grid();
		
		System.out.println("");
		System.out.println("Array POLYGONS properties length");
		System.out.println("");
		
		System.out.println("numberSidesPolygon array length: "+ gridTest.numberSidesPolygon.length);
		System.out.println("planArea array length: "+ gridTest.planArea.length);
		System.out.println("sourceSink array length: "+ gridTest.sourceSink.length);
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
		System.out.println("Ml array length: "+ gridTest.Ml.length);
		
	}
	
}
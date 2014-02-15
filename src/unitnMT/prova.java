package unitnMT;

import java.util.Collection;
import java.util.LinkedList;

public class prova {
	

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		final Grid grid1 = new Grid("Song");
		final int[] indexDiag = new int[grid1.numberSidesPolygon.length];
		
		Collection<Integer> elems = new LinkedList<Integer>();
		for (int i = 0; i < grid1.numberSidesPolygon.length; ++i) {
		    elems.add(i);
		}
		
		// to explains the JVM type 64 bit or 32 bit. 
		System.getProperties().list(System.out);
		String bits = System.getProperty("sun.arch.data.model", "?");
		System.out.println("\n\n\nJVM type: " + bits + "-bit");
		
		System.out.println("JVM processors available: "
							+ Parallel.NUM_CORES);
		
		System.out.println("\nShow multi-core input order in the indexDiag array\n");

		
		Parallel.For(elems, 
		 // The operation to perform with each item
		 new Parallel.Operation<Integer>() {
			
		    public void perform(Integer params) {
		        //System.out.println(param);
		        
		    	
				//System.out.println("Number of cores: "Parallel.NUM_CORES);

					for (int j = grid1.Mp[params]; j < grid1.Mp[params + 1]; j++) {
						if (grid1.Mi[j] == params) {
							indexDiag[params] = j;
							
							System.out.println(indexDiag[params]);
						}
					}
		    };
		});
		
		System.out.println("\nNow prints the ordered indexDiag array\n");
		
		for (int i=0; i<grid1.numberSidesPolygon.length; i++){
			System.out.println(indexDiag[i]);
		}
		
		
		
		double machEps = 1.0;
		 
        do
           machEps /= 2.0;
        while ((double) (1.0 + (machEps / 2.0)) != 1.0);
 
        System.out.println("\nSystem tolerance: " + machEps);
		
		
		
		System.exit(0);
	}

}

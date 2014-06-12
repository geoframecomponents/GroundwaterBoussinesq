package org.wordpress.growworkinghard.usefulClasses;

public class CodeTitle {

	
	public void firstLines(){
		
		TextIO.putln("#---------------------------------------------#");
		TextIO.putln("#                                             #");
		TextIO.putln("#                                             #");
		
	}
	
	public void lastLines(){
		
		TextIO.putln("#                                             #");
		TextIO.putln("#                                             #");
		TextIO.putln("#---------------------------------------------#");
		
	}
	
	public void continueProgram(){
		
		String string;
		
		TextIO.putln("\n --------------------------------------------- ");
		TextIO.putln(" type yes to continue: ");
		string = TextIO.getlnString();
		
		if (string.equals("yes")){}
		else {
			
			TextIO.putln("Thank you for open this program\nI hope you will use it in a few time");
			System.exit(0);
			
		}
		
		
	}
	
	public static void main(String[] args) {
		
		CodeTitle prova = new CodeTitle();
		prova.firstLines();
		prova.lastLines();
		prova.continueProgram();
		
		System.exit(0);

	}

}

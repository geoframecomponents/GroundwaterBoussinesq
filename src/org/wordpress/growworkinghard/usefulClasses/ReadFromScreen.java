package org.wordpress.growworkinghard.usefulClasses;

public class ReadFromScreen {

	
	
	
	
	
	
	
	
	
	
	public static String readText(String outputText){
		
		String inputText = null;
		
		TextIO.putln(outputText);
		inputText = TextIO.getln();
		
		return inputText;
		
	}
	
	
	
	
	
	
	
	
	
	
	public static void main(String[] args) {
		
		String outputTxt = "Digit the name of solution folder";
		String inputTxt = readText(outputTxt);

		TextIO.putln(inputTxt);
		
	}

}

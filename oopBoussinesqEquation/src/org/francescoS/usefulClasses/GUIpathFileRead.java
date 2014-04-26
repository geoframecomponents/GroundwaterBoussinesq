package org.francescoS.usefulClasses;

import java.io.File;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JPanel;

public class GUIpathFileRead extends JPanel {
	String sep = System.getProperty("file.separator");

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	JFileChooser fileChooser = new JFileChooser();
	
	public File openDialog(String title) {
		// TODO Auto-generated constructor stub
		fileChooser.setDialogTitle(title);
		int n = fileChooser.showOpenDialog(GUIpathFileRead.this);
		
		File f = null;
		
		if (n==JFileChooser.APPROVE_OPTION){
			
			 f = fileChooser.getSelectedFile();
			
		
		}else if (n==JFileChooser.CANCEL_OPTION){
			
			TextIO.putln("WARNING!!! You have select CANCEL option\n"
					+ "You need to restart the program");
			System.exit(0);
			
			
		} else {
			
			TextIO.putln("ERROR ENCOUNTERED!!!");	
			System.exit(0);
			
		}
		
		return f;
		
	}
	
	
	public File saveDialog(){
		
		File f, solutionDir;
		
		// parent component of the dialog
		JFrame parentFrame = new JFrame();		
		
		fileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		int returnVal = fileChooser.showOpenDialog(parentFrame);
		
		f = new File(fileChooser.getSelectedFile().getAbsolutePath());
//		f = fileChooser.getCurrentDirectory();
		
		if(returnVal == JFileChooser.APPROVE_OPTION) {
		   System.out.println("You chose to open this directory: " +
		        fileChooser.getSelectedFile().getAbsolutePath());
		}
		
		solutionDir = new File(f, sep);
		
		return solutionDir;
		
	}
		
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		GUIpathFileRead gui = new GUIpathFileRead();
		gui.openDialog("Firt test");
		
		

	}

}

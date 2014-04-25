package org.francescoS.GUI;


import javax.swing.*;
import javax.swing.JPanel;
import java.awt.*;
public class inputWindow {
   JFrame Frame1 = new JFrame("Test Frame");
   JPanel j2 = new JPanel();
   JButton b1 = new JButton ("Click me");
   JTextField t1,t2;
   JPanel j1 = new JPanel (new FlowLayout());
   JLabel l1, l2;

   
   inputWindow()

    {
	   t1 = new JTextField(20);
	   t2 = new JTextField(20);
	   
	   l1 = new JLabel("Initial Condition: ");
	   l2 = new JLabel("Source");
	   
	   
	   
	   
    j1.add(l1);

    j1.add(t1);
    
    j1.add(b1);
    
    j1.add(l2);
    
    j1.add(t2);

//    j2.add(b1);

    Frame1.add(j1);

    Frame1.add(j2);

    Frame1.setLayout(new FlowLayout());

    Frame1.setSize(400,150);

    Frame1.setVisible(true);

   Frame1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

  }


public static void main (String[] arg0){
	
	inputWindow prova = new inputWindow();
	
}
}
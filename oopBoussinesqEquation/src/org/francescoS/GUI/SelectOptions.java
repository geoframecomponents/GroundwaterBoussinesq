package org.francescoS.GUI;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;
 
public class SelectOptions {
    
   private JFrame mainFrame;
   private JLabel headerLabel;
   private JLabel statusLabel;
   private static JPanel controlPanel;
   DefaultComboBoxModel<String> options;
   public static String item;
   public static String name;
   static int click = 0;
   static boolean done = false;
   
   public SelectOptions(){
	   
	   options = new DefaultComboBoxModel<String>();
	      options.addElement("--Select--");
	   prepareGUI();

   }
   
   public void addOptions(String input){
	   
	   options.addElement(input);
	   
   }
   
   private void prepareGUI(){
      mainFrame = new JFrame("Java Swing Examples");
      mainFrame.setSize(400,400);
      mainFrame.setLayout(new GridLayout(3, 1));
      mainFrame.addWindowListener(new WindowAdapter() {
         public void windowClosing(WindowEvent windowEvent){
            System.exit(0);
         }        
      });    
      headerLabel = new JLabel("", JLabel.CENTER);        
      statusLabel = new JLabel("",JLabel.CENTER);    

      statusLabel.setSize(350,100);

      controlPanel = new JPanel();
      controlPanel.setLayout(new FlowLayout());

      mainFrame.add(headerLabel);
      mainFrame.add(controlPanel);
      mainFrame.add(statusLabel);
      mainFrame.setVisible(true);  
   }

   public void showComboboxDemo(){                                    
      headerLabel.setText("Control in action: JComboBox"); 
      

      final JComboBox<String> fruitCombo = new JComboBox<String>(options);    
      fruitCombo.setSelectedIndex(0);

      JScrollPane fruitListScrollPane = new JScrollPane(fruitCombo);    

      JButton showButton = new JButton("OK");
      
      controlPanel.add(fruitListScrollPane);          
      controlPanel.add(showButton);

      mainFrame.setVisible(true); 
      
      
      
      do {
              // do stuff

     
      
      ActionListener actionListener = new ActionListener() {
      public void actionPerformed(ActionEvent actionEvent) {
    	  
    	  
      ItemSelectable is = (ItemSelectable)actionEvent.getSource();
      
      String data = "";
      if (fruitCombo.getSelectedIndex() != -1) {                     
         data = "SIMULATION SELECTED: " 
            + fruitCombo.getItemAt
              (fruitCombo.getSelectedIndex());             
      }              
      statusLabel.setText(data);
      
      name=selectedString(is);
      
//      JOptionPane.showMessageDialog(null,"You have Selected: " + name);
     
      
      
      }

      };
      
      fruitCombo.addActionListener(actionListener);


      

//    else buttonPress = false;  // insures buttonPress is false, not needed


      showButton.addActionListener(new ActionListener() {
         public void actionPerformed(ActionEvent e) { 

        	 
        	 mainFrame.dispose();
        	 click = 1;
        	 
        	 
         }
      }); 


      if (click == 1) done = true;
      // ends loop
    
    } while(!done);

  


      
   } 
   
   static private String selectedString(ItemSelectable is) {
   Object selected[] = is.getSelectedObjects();
   return ((selected.length == 0) ? "null" : (String)selected[0]);
 } 
   
   public static void main(String[] args){
	      SelectOptions swingControlDemo = new SelectOptions();   
	      swingControlDemo.addOptions("Song");
	      swingControlDemo.addOptions("Catchment");
	      swingControlDemo.showComboboxDemo();
	      
	      System.out.println(name);
	      
	      System.exit(0);
	   }
}
package IMRT_Base;

/* The function solves the problem using AMPL. To do that, the function
 * first create a script to be run from the AMPL terminal. It also creates a
 * .dat file with some specific parameter values (extraXXXXX.dat). 
 * Once the problem is solved, the function return solution  * vector "x". 
 */

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Random;
/**
 *
 * @author guille
 */
public class AMPL_Solver {
    public int organs;
    public int beams;
    public int[] bmlts;
    public int totalBmlts;
    public int[] angles;
    public int[] voxels;
    public int[] aPar;
    public int[] vPar;
    public double[] wPar; //weights
    public double[] EUD0Par;
    public double[] LB;
    public double[] UB;
    public Boolean[] isTarget;
    public double epsilon;
    //public double t;
    public double[] x;
    public String jobThreadID;
    public boolean output;
    //public double PTV_UB;
    
    /**
     * @param args the command line arguments
     */
    //public AMPL_Solver (Organs_bkp[] o,TreatmentPlan sol, double e) throws IOException {
    public AMPL_Solver (Organs[] o,TreatmentPlan sol, double e, String jobThreadID, boolean outputl) 
            throws IOException {
        
        
        this.organs = o.length;
        this.beams = sol.beams;
        this.bmlts=new int[this.beams];
        this.totalBmlts = sol.beamlets;
        this.angles= new int[this.beams];
        for(int i=0; i< this.beams; i++){
            this.angles[i] = sol.selAngles[i].index;
            this.bmlts[i] = sol.selAngles[i].beamlets;
        }
        this.voxels= new int[this.organs];
        this.aPar =new int[this.organs];
        this.vPar =new int[this.organs];
        this.EUD0Par =new double[this.organs];
        this.UB =new double[this.organs];
        this.LB =new double[this.organs];
        this.wPar =new double[this.organs];
        this.isTarget = new Boolean[this.organs];
        for(int i=0; i< o.length; i++){
            this.voxels[i]= o[i].totalVoxels;
            this.aPar[i] =  o[i].a;
            this.vPar[i] =  o[i].v;
            this.EUD0Par[i] =  o[i].desiredDose;
            this.UB[i] = o[i].doseUB;
            this.LB[i] = o[i].doseLB;
            this.isTarget[i]=o[i].isTarget;
            this.wPar[i] = o[i].weight;
            //if (o[i].isTarget){
            //    this.PTV_UB =  o[i].doseUB;
            //}
        }
        this.epsilon = e;
        //this.t= o[0].doseLB;
        this.x = new double[this.totalBmlts]; 
        this.jobThreadID = jobThreadID;
        this.output=outputl;
        
    }
    
    public void generateParametersFile() throws IOException{
        String parameterFile = this.jobThreadID + "extra.dat";
        //Deleting parameter file extra.txt
              
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(1)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
        }
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + r.nextDouble()*50 + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        
        
        
        bwParametersFile.close();
        
    }
    public void generateParametersFile(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        String parameterFile = this.jobThreadID + "extra.dat";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(2)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
    }
    
    public void generateParametersFile_logFunction() throws IOException{
        String parameterFile = this.jobThreadID + "extraLogFunction.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_logFunction(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraLogFunction.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_convexLogFunction() throws IOException{
        String parameterFile = this.jobThreadID + "extraConvexLogFunction.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_convexLogFunction(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraConvexLogFunction.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
     public void generateParametersFile_variableWeights() throws IOException{
        String parameterFile = this.jobThreadID + "extraVariableWeights.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_variableWeights(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraVariableWeights.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_adaptiveWeightedSum() throws IOException{
        String parameterFile = this.jobThreadID + "extraAdaptiveWeightedSum.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_adaptiveWeightedSum(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraAdaptiveWeightedSum.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
     public void generateParametersFile_weightedSum() throws IOException{
        String parameterFile = this.jobThreadID + "extraWeightedSum.dat";
        //Deleting parameter file extra.txt
        
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(3)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        Random r = new Random();
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        //writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) +" 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_weightedSum(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        
        String parameterFile = "./" + this.jobThreadID + "extraWeightedSum.dat";
        System.out.println(parameterFile);  
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(4)");
    		}
            }
    	}catch(Exception e){
    	}
        Random r = new Random();
        //creating the new file
        
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("var x := ", bwParametersFile);
        for (int i=0;i<this.totalBmlts;i++){
            int j = i+1;
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param w := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.wPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine((this.organs+1) + " 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            int j=i+1;
            writeLine("param R" + j + " := " + this.voxels[i] + ";\n", bwParametersFile);
            
            if (this.isTarget[i]){
                writeLine("param t := " + EUD0Par[i] + ";\n", bwParametersFile);
                writeLine("param OAR_targetUB := " + this.UB[i] + ";\n", bwParametersFile);
            }else{
                writeLine("param UB" + j + " := " + this.UB[i] + ";\n", bwParametersFile);
                //writeLine("param LB" + j + " := " + this.LB[i] + ";\n", bwParametersFile);
            }
        }
        
        //writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        //writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        //writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        //writeLine("param t := " + this.t + ";\n", bwParametersFile);
        
        bwParametersFile.close();
        
    }
    
    public void getSolution() throws FileNotFoundException, IOException{
        String dir = this.jobThreadID + "currentSol.txt";
        String[] auxReader=null;
        File f = new File(dir);
        if (f.exists()) {
            BufferedReader fileIn = new BufferedReader(new FileReader(f));
            String line = "";
            line=fileIn.readLine(); //avoid first line;
            line=fileIn.readLine();
            auxReader = line.split(" ");
            int j=0;
            double[] auxX = new double[this.x.length * 2];
            while (!";".equals(auxReader[0])){

                for (String auxReader1 : auxReader) {
                    if (!"".equals(auxReader1)) {
                        auxX[j] = Double.parseDouble(auxReader1);
                        j++;
                    }
                }
                line=fileIn.readLine();
                //System.out.println(line);
                auxReader = line.split(" ");
            }
            j=0;
            for (int i=0; i<auxX.length;i++){
                if (i%2 == 0){
                    j=(int)auxX[i];
                }else{
                    if (j>0){
                        x[j-1] = auxX[i];
                    }else{
                        System.err.print("ERROR al leer " + jobThreadID + "currentSol.txt");
                        System.err.print("auxX = ");
                        for (int k=0; k<auxX.length;k++){
                            System.err.print(auxX[k]+ " - ");
                        }
                        System.err.println();
                        System.err.print("x = ");
                        for (int k=0; k<this.x.length;k++){
                            System.err.print(this.x[k]+ " - ");
                        }
                        x[j-1] = auxX[i];
                    }
                }       
            }
        }else{
            System.out.println("ERROR: CurrentSol file wasn't generated");
        }
    }

     private BufferedWriter createBufferedWriter(String route) throws IOException {
        BufferedWriter bw;
        File auxFile = new File(route);
                if (!auxFile.exists()) {
                    bw = new BufferedWriter(new FileWriter(route));
                }else{
                    bw = new BufferedWriter(new FileWriter(route,true));
                }
        return bw;     
    }
     
    
     public void runLogFunction_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptLogFunction.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptLogFunction.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"logisticModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraLogFunction.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "logisticModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "logisticModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
            
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param v{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("- log((1+(((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])/EUD0["+(o1.index +1)+"])^v["+(o1.index +1)+"])^-1)", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptLogFunction.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving Log Function (Wu et al., 2002): ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
     public void runConvexLogFunction_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptConvexLogFunction.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptConvexLogFunction.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"convexLogModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraConvexLogFunction.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "convexLogModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "convexLogModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param v{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("- log(((((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])/EUD0["+(o1.index +1)+"])^v["+(o1.index +1)+"])^-1)", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptConvexLogFunction.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving Convex Log Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
     public void runWeightedSum_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptWeightedSum.sh");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptWeightedSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"weightedSumModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraWeightedSum.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "weightedSumModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "weightedSumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptWeightedSum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving WeightedSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void runVariableWeights_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptVariableWeights.sh with weights ( " +this.wPar[1] +" - "+this.wPar[2] + " )");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptVariableWeights.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"variableWeightsModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraVariableWeights.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
        
        /*for (Organs o1 : o) {
            writeLine("display gEUD_" + o1.name + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile);
            writeLine("display d_" + o1.name + " > " + this.jobThreadID + "dvh_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display gEUD_" + o1.name + "_UB > "+"gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile);
            }
        }*/
        
         bwParametersFile.close();
        
        /*CREATING THE LOGISTIC MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "variableWeightsModel.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "variableWeightsModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        scriptFile = this.jobThreadID + "scriptVariableWeights.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving VariableWeights Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void runAdaptiveWeightedSum_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(5)");
    		}
            }
    	}catch(Exception e){
    	}
        System.out.println("Creating " + this.jobThreadID + "scriptAdaptiveWeightedSum.sh with weights ( " +this.wPar[1] +" - "+this.wPar[2] + " )");  
        //Creating the script file
        String scriptFile = this.jobThreadID + "scriptAdaptiveWeightedSum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model " + this.jobThreadID +"adaptiveWeightedSum.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraAdaptiveWeightedSum.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        
        
        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        
         bwParametersFile.close();
        
        /*CREATING THE WEIGHTED SUM MODEL FILE FOR AMPL*/
        
        System.out.println("Creating " + this.jobThreadID + "adaptiveWeightedSum.mod");  
        //Creating the script file
        scriptFile = this.jobThreadID + "adaptiveWeightedSum.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(6)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 =0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        writeLine("param w{1 .. "+ (o.length) + "}; \n", bwParametersFile);
        writeLine("param EUD0{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        
        
        /*for (Organs o1 : o) {
            writeLine("par d_" + o1.name + " {i in 1 .. R"+(o1.index +1)+"} = (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]); \n", bwParametersFile); 
        }*/
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        for (Organs o1 : o) {
            if (!o1.isTarget){
                writeLine("+ w["+(o1.index+1)+"] * ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])", bwParametersFile);
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                writeLine("#OAR_UB"+(o1.index +1)+": 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= UB"+(o1.index +1)+"; \n", bwParametersFile);
            }
        }
        bwParametersFile.close();
        
        
        
        
        
        
        System.out.println("PASO 2");
        //Running Process
        
        scriptFile = this.jobThreadID + "scriptAdaptiveWeightedSum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving Adaptive WeightedSum Function: ");
        //System.out.print("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // EUD_0: ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.print(" // Weights: ");
        for(int i=0; i< this.wPar.length; i++){
            System.out.print(this.wPar[i]+ " - ");
        }
        System.out.println();
        
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void runLexRectum_Solver(Organs o[]) throws IOException{
        
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(7)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "scriptLexicoRectum.sh");  
        String scriptFile = this.jobThreadID + "scriptLexicoRectum.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(9.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"lexicoRectumModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create LexicoRectum.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "lexicoRectumModel.mod");  
        scriptFile = this.jobThreadID + "lexicoRectumModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(9.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("RECTUM".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "scriptLexicoRectum.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving LexicoRectum problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void run_gEUD_Rectum_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(7)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "script_gEUD.sh");  
        String scriptFile = this.jobThreadID + "script_gEUD.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"gEUDConstrainedModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create gEUDConstrainedModel.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "gEUDConstrainedModel.mod");  
        scriptFile = this.jobThreadID + "gEUDConstrainedModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("RECTUM".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                if ("BLADDER".equals(o1.name)){
                    //writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                    //    + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= epsilon;\n", bwParametersFile); 
                    writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) = epsilon;\n", bwParametersFile); 
                }
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "script_gEUD.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving gEUD_Constrained problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void run_gEUD_Bladder_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(7)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "script_gEUD.sh");  
        String scriptFile = this.jobThreadID + "script_gEUD.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"gEUDConstrainedModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create gEUDConstrainedModel.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "gEUDConstrainedModel.mod");  
        scriptFile = this.jobThreadID + "gEUDConstrainedModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(11.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("BLADDER".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }else{
                if ("RECTUM".equals(o1.name)){
                    //writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                    //    + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) <= epsilon;\n", bwParametersFile); 
                    writeLine("constraintOAR_OAR:((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) = epsilon;\n", bwParametersFile); 
                }
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "script_gEUD.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving gEUD_Constrained problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void runLexBladder_Solver(Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(10)");
    		}
            }
    	}catch(Exception e){
    	}
        
        //Create scriptLexicoRectum.sh file
        System.out.println("Creating " + this.jobThreadID + "scriptLexicoBladder.sh");  
        String scriptFile = this.jobThreadID + "scriptLexicoBladder.sh";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(10.1)");
    		}
            }
    	}catch(Exception e){
    	}
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        writeLine("model "+this.jobThreadID+"lexicoBladderModel.mod;\n", bwParametersFile);
        for (Organs o1 : o) {
            writeLine("data "+ this.jobThreadID + "DDM_" + o1.name + ".dat;\n", bwParametersFile);
        }
        
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);

        for (Organs o1 : o) {
            writeLine("display ((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"])" + " > " + this.jobThreadID + "gEUD_" + o1.name + ".txt;\n", bwParametersFile); 
            writeLine(" display {i in 1 .. R"+(o1.index +1)+"} (sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j]) > " + this.jobThreadID + "DVH_" + o1.name + ".txt;\n", bwParametersFile); 
            if (o1.isTarget){
                writeLine("display " + "((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) > " + this.jobThreadID + "gEUD_" + o1.name + "_UB.txt;\n", bwParametersFile); 
            }
        }
        writeLine("display _con > "+ this.jobThreadID +"lagrangeMultiplier.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Create LexicoBladder.mod for AMPL*/
        System.out.println("Creating " + this.jobThreadID + "lexicoBladderModel.mod");  
        scriptFile = this.jobThreadID + "lexicoBladderModel.mod";
        try{
            File file = new File(scriptFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(10.2)");
    		}
            }
    	}catch(Exception e){
    	}
        //creating the new file
        bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("option solver ipopt; \n", bwParametersFile);
        //writeLine("options ipopt_options \"linear_solver=ma57 linear_system_scaling=mc19 wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile);
        //writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile);
        if(output){
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001 output_file="+System.currentTimeMillis()+"_"+this.jobThreadID+"_"+"salida file_print_level=3\"; \n", bwParametersFile); //Cambiar luego
        }else{
            writeLine("options ipopt_options \"wantsol=8 print_level=4 tol=0.0001\"; \n", bwParametersFile); //Cambiar luego
        }
        
        for (Organs o1 : o) {
            writeLine("param R" + (o1.index +1) + "; #number of voxels of " + o1.name + "\n", bwParametersFile);
        }
        writeLine("param bmlt; 	#number of beamlets \n", bwParametersFile);
        
        for (Organs o1 : o) {
            writeLine("param ddm" + o1.name + "{1 .. R" + (o1.index +1) + ", 1 .. bmlt};\n", bwParametersFile);
        }
        for (Organs o1 : o) {
            writeLine("param UB" + (o1.index +1) + ";\n", bwParametersFile);
        }
        writeLine("param t;\n" + "param epsilon;\n" + "param OAR_targetUB;\n", bwParametersFile);
        writeLine("var x {1 .. bmlt} >= 0, <=400, default 1; \n", bwParametersFile);
        
        writeLine("param a{1 .. "+ (o.length + 1) + "}; \n", bwParametersFile);
        
        writeLine("minimize Total_Cost: ", bwParametersFile);
        
        for (Organs o1 : o) {
            if ("BLADDER".equals(o1.name)){
                writeLine("((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} "
                        + "(sum {j in 1..bmlt} x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]);\n", bwParametersFile); 
            }
        }
        writeLine(";\n s.t. \n", bwParametersFile);
        
        for (Organs o1 : o) {
            if (o1.isTarget){
                writeLine("equalityTarget: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o1.index +1)+"]))^(1/a["+(o1.index +1)+"]) >= t; \n", bwParametersFile);
                writeLine("constraintOAR_Target: 	((1/R"+(o1.index +1)+")*(sum {i in 1..R"+(o1.index +1)+"} (sum {j in 1..bmlt} "
                        + "x[j]*ddm"+o1.name+"[i,j])^a["+(o.length +1)+"]))^(1/a["+(o.length +1)+"]) <=OAR_targetUB;\n", bwParametersFile);
            }
        }
        bwParametersFile.close();
            
        
        //Running Process
        scriptFile = this.jobThreadID + "scriptLexicoBladder.sh";
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving LexicoBladder problem: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.print(" // ");
        for(int i=0; i< this.EUD0Par.length; i++){
            System.out.print(this.EUD0Par[i]+ " - ");
        }
        System.out.println();
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void writeLine(String l, BufferedWriter bw) throws IOException{
	
	String row;
	row = l; 
	bw.write(row);
	

    }
   
}

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
public class JuMP {
    public int organs;
    public int beams;
    public int[] bmlts;
    public int totalBmlts;
    public int[] angles;
    public int[] voxels;
    public int[] aPar;
    public int[] vPar;
    public double[] EUD0Par;   
    public double epsilon;
    public double t;
    public double[] x;
    public String jobThreadID;
    
    /**
     * @param args the command line arguments
     */
    //public AMPL_Solver (Organs_bkp[] o,TreatmentPlan sol, double e) throws IOException {
    public JuMP (Organs[] o,TreatmentPlan sol, double e, String jobThreadID) 
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
        for(int i=0; i< o.length; i++){
            this.voxels[i]= o[i].totalVoxels;
            this.aPar[i] =  o[i].a;
            this.vPar[i] =  o[i].v;
            if (i==0){
                this.EUD0Par[i] =  o[i].doseLB;
            }else{
                this.EUD0Par[i] =  o[i].doseUB;
            }
            
        }
        this.epsilon = e;
        this.t= o[0].doseLB;
        this.x = new double[this.totalBmlts]; 
        this.jobThreadID = jobThreadID;
    }
    
    public void generateParametersFile() throws IOException{
        //Deleting parameter file extra.txt
        String parameterFile = this.jobThreadID + "Extra.dat";
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
        writeLine("4 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        writeLine("param t := " + this.t + ";\n", bwParametersFile);
        bwParametersFile.close();
        
    }
    public void generateParametersFile(double[] x) throws IOException{
        //Deleting parameter file extra.txt
        String parameterFile = this.jobThreadID + "Extra.dat";
        try{
            File file = new File(parameterFile);
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(2)");
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
            writeLine(j + " " + x[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param a := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.aPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine("4 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        writeLine("param t := " + this.t + ";\n", bwParametersFile);
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_logFunction(Organs o[]) throws IOException{
        //Deleting parameter file extra.txt
        String parameterFile = this.jobThreadID + "JuMPScript_logFunction.txt";
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
        //creating the new file
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("using JuMP\n", bwParametersFile);
        writeLine("Package.add(\"Gurobi\")\n", bwParametersFile);
        writeLine("Package.add(\"Ipopt\")\n", bwParametersFile);
        writeLine("m = Model()\n", bwParametersFile);
        writeLine("t = "+this.t +"\n", bwParametersFile);
        writeLine("bmlt = "+this.bmlts +"\n", bwParametersFile);
        
        
        String auxString = "a = [";
        for (int i=0;i<this.organs-1;i++){
            auxString = auxString + this.aPar[i] + "\t";
        }
        writeLine(auxString + this.aPar[this.organs-1] + "\t10]\n", bwParametersFile);
        
        auxString = "v = [";
        for (int i=0;i<this.organs-1;i++){
            auxString = auxString + this.vPar[i] + "\t";
        }
        writeLine(auxString + this.vPar[this.organs-1] + "\t8]\n", bwParametersFile);
        
        auxString = "EUD0 = [";
        for (int i=0;i<this.organs-1;i++){
            auxString = auxString + this.EUD0Par[i] + "\t";
        }
        writeLine(auxString + this.EUD0Par[this.organs-1] + "t80]\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            auxString = "R" + i + " = " + this.voxels[i] + "\n";
        }
        
        writeLine("@defVar(m, x[1:bmlt] >= 0 )\n", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            writeLine("ddm" + o[i].name + " = readdlm(./" + this.jobThreadID + "DDM_" + o[i].name + ".txt, '\t')", bwParametersFile);
        }
        
        auxString = "@setNLObjective(m, Min, ";
        for (int i=0;i<this.organs;i++){
            auxString = auxString + "- log( (1+( ((1/R" + i + ") * (sum{(x[j]*ddm"+ o[i].name +"[i,j])^a[" + i + "], i = R" + i +
                    ", j = bmlt})) ^ (1/a["+ i +"])/EUD0[" + i + "])^v[" + i + "])^-1) ";
        }
        writeLine(auxString + ")\n", bwParametersFile);

        auxString = "@addNLConstraint(m, ";
        for (int i=0;i<this.organs;i++){
            if (o[i].isTarget){
                auxString = auxString + "((1/R" + i + ") * (sum{(x[j]*ddm"+ o[i].name +"[i,j])^a[" + i + "], i = R" + i +
                    ", j = bmlt})) ^ (1/a["+ i +"]) == t ";
            }
            writeLine(auxString + ")\n", bwParametersFile);
        }
        
        auxString = "@addNLConstraint(m, ( (1/R1) * (sum{(x[j]*ddmPTVHD[i,j])^a[4], i = R1, j = bmlt}))^(1/a[4]) <=OAR_targetUB)";
        writeLine(auxString + "\n", bwParametersFile);
        bwParametersFile.close();
        
    }
    
    public void generateParametersFile_logFunction(double[] x, Organs[] o) throws IOException{
        //Deleting parameter file extra.txt
        String parameterFile = this.jobThreadID + "JuMPScript_logFunction.txt";
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
        //creating the new file
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(parameterFile);
        
        writeLine("using JuMP\n", bwParametersFile);
        writeLine("Package.add(\"Gurobi\")\n", bwParametersFile);
        writeLine("Package.add(\"Ipopt\")\n", bwParametersFile);
        writeLine("m = Model()\n", bwParametersFile);
        writeLine("t = "+this.t +"\n", bwParametersFile);
        writeLine("bmlt = "+this.bmlts +"\n", bwParametersFile);
        
        
        String auxString = "a = [";
        for (int i=0;i<this.organs-1;i++){
            auxString = auxString + this.aPar[i] + "\t";
        }
        writeLine(auxString + this.aPar[this.organs-1] + "\t10]\n", bwParametersFile);
        
        auxString = "v = [";
        for (int i=0;i<this.organs-1;i++){
            auxString = auxString + this.vPar[i] + "\t";
        }
        writeLine(auxString + this.vPar[this.organs-1] + "\t8]\n", bwParametersFile);
        
        auxString = "EUD0 = [";
        for (int i=0;i<this.organs-1;i++){
            auxString = auxString + this.EUD0Par[i] + "\t";
        }
        writeLine(auxString + this.EUD0Par[this.organs-1] + "t80]\n", bwParametersFile);
        
        for (int i=0;i<this.organs;i++){
            auxString = "R" + i + " = " + this.voxels[i] + "\n";
        }
        
        writeLine("@defVar(m, x[1:bmlt] >= 0 )\n", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            writeLine("ddm" + o[i].name + " = readdlm(./" + this.jobThreadID + "DDM_" + o[i].name + ".txt, '\t')", bwParametersFile);
        }
        
        auxString = "@setNLObjective(m, Min, ";
        for (int i=0;i<this.organs;i++){
            auxString = auxString + "- log( (1+( ((1/R" + i + ") * (sum{(x[j]*ddm"+ o[i].name +"[i,j])^a[" + i + "], i = R" + i +
                    ", j = bmlt})) ^ (1/a["+ i +"])/EUD0[" + i + "])^v[" + i + "])^-1) ";
        }
        writeLine(auxString + ")\n", bwParametersFile);

        auxString = "@addNLConstraint(m, ";
        for (int i=0;i<this.organs;i++){
            if (o[i].isTarget){
                auxString = auxString + "((1/R" + i + ") * (sum{(x[j]*ddm"+ o[i].name +"[i,j])^a[" + i + "], i = R" + i +
                    ", j = bmlt})) ^ (1/a["+ i +"]) == t ";
            }
            writeLine(auxString + ")\n", bwParametersFile);
        }
        
        auxString = "@addNLConstraint(m, ( (1/R1) * (sum{(x[j]*ddmPTVHD[i,j])^a[4], i = R1, j = bmlt}))^(1/a[4]) <=OAR_targetUB)";
        writeLine(auxString + "\n", bwParametersFile);
        bwParametersFile.close();
        
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
        writeLine("4 10\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param v := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.vPar[i] + "\t", bwParametersFile);
        }
        //OAR-Target
        writeLine("4 8\t", bwParametersFile);
        writeLine(";\n", bwParametersFile);
        
        writeLine("param EUD0 := ", bwParametersFile);
        for (int i=0;i<this.organs;i++){
            int j = i+1;
            writeLine(j + " " + this.EUD0Par[i] + "\t", bwParametersFile);
        }
        writeLine(";\n", bwParametersFile);
        
        writeLine("param R1 := " + this.voxels[0] + ";\n", bwParametersFile);
        writeLine("param R2 := " + this.voxels[1] + ";\n", bwParametersFile);
        writeLine("param R3 := " + this.voxels[2] + ";\n", bwParametersFile);
        writeLine("param bmlt := " + this.totalBmlts + ";\n", bwParametersFile);
        writeLine("param epsilon := " + this.epsilon + ";\n", bwParametersFile);
        writeLine("param t := " + this.t + ";\n", bwParametersFile);
        bwParametersFile.close();
        
    }
    
    public void getSolution() throws FileNotFoundException, IOException{
        String dir = this.jobThreadID + "currentSol.txt";
        String[] auxReader=null;
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine(); //avoid first line;
        line=fileIn.readLine();
        auxReader = line.split(" ");
        int j=0;
        double[] auxX = new double[this.x.length * 2];
        while (!";".equals(auxReader[0])){
            
            for (int i=0; i<auxReader.length;i++){
                if (!"".equals(auxReader[i])){
                    auxX[j]=Double.parseDouble(auxReader[i]);
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
            e.printStackTrace();
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
            e.printStackTrace();
    	}
        
        //creating the new file
        
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model WuFormula.mod;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("data "+ this.jobThreadID + "DDM_" + o[y].name+".dat;\n", bwParametersFile);
        }
        System.out.println("PASO 1");
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extraLogFunction.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("display gEUD_" +o[y].name+ " > "  + this.jobThreadID + "gEUD_"+o[y].name+".txt;\n",bwParametersFile); 
        }
        //writeLine("display gEUD_PTV > "+ Thread.currentThread().getName()+"gEUD_PTVHD.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Rectum > "+ Thread.currentThread().getName()+"gEUD_RECTUM.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Bladder > "+ Thread.currentThread().getName()+"gEUD_BLADDER.txt;\n", bwParametersFile);
        bwParametersFile.close();
        System.out.println("PASO 2");
        //Running Process
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
        System.out.println("t = " + this.t );
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
            e.printStackTrace();
    	}
        

        //creating the new file
        String scriptFile = this.jobThreadID + "scriptLexicoRectum.sh";
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model LexicoRectum.mod;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("data "+ this.jobThreadID + "DDM_" + o[y].name+".dat;\n", bwParametersFile);
        }
        
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("display gEUD_" +o[y].name+ " > "  + this.jobThreadID + "gEUD_"+o[y].name+".txt;\n",bwParametersFile); 
        }
        //writeLine("display gEUD_PTV > "+ Thread.currentThread().getName()+"gEUD_PTVHD.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Rectum > "+ Thread.currentThread().getName()+"gEUD_RECTUM.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Bladder > "+ Thread.currentThread().getName()+"gEUD_BLADDER.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        
        //Running Process
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
        System.out.println("t = " + this.t );
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void run_gEUD_Solver (Organs[] o) throws IOException{
        //Deleting Solution File currentSol.txt
        try{
            File file = new File(this.jobThreadID + "currentSol.txt");
            if (file.exists()) {
    		if(!file.delete()){
                    System.out.println("Delete operation failed.(8)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        
        //creating the new file
        String scriptFile = this.jobThreadID + "script_gEUD.sh";
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model gEUDConstrained.mod;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("data "+ this.jobThreadID + "DDM_" + o[y].name+".dat;\n", bwParametersFile);
        }
        
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("display gEUD_" +o[y].name+ " > "  + this.jobThreadID + "gEUD_"+o[y].name+".txt;\n",bwParametersFile); 
        }
        //writeLine("display gEUD_PTV > "+ Thread.currentThread().getName()+"gEUD_PTVHD.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Rectum > "+ Thread.currentThread().getName()+"gEUD_RECTUM.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Bladder > "+ Thread.currentThread().getName()+"gEUD_BLADDER.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        
        //Running Process
        Process p = new ProcessBuilder("ampl", scriptFile).start();
        InputStream is = p.getInputStream();
        InputStreamReader isr = new InputStreamReader(is);
        BufferedReader br = new BufferedReader(isr);
        String line;
        System.out.print("Solving gEUD Constrained Model: ");
        //System.out.println("Angles: ");
        for(int i=0; i< this.beams; i++){
            System.out.print(" " +this.angles[i]+ " -- ");
        }
        System.out.println("t = " + this.t );
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
                    System.out.println("Delete operation failed.(9)");
    		}
            }
    	}catch(Exception e){
            e.printStackTrace();
    	}
        
        //creating the new file
        String scriptFile = this.jobThreadID + "scriptLexicoBladder.sh";
        BufferedWriter bwParametersFile=null;
        bwParametersFile =createBufferedWriter(scriptFile);
        
        writeLine("model LexicoBladder.mod;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("data "+ this.jobThreadID + o[y].name+".dat;\n", bwParametersFile);
        }
        
        //writeLine("data "+ this.jobID +"DDM_RECTUM.dat;\n", bwParametersFile);
        //writeLine("data "+ this.jobID +"DDM_BLADDER.dat;\n", bwParametersFile);
        writeLine("data "+ this.jobThreadID +"extra.dat;\n", bwParametersFile);
        writeLine("solve;\n", bwParametersFile);
        writeLine("display x > "+ this.jobThreadID +"currentSol.txt;\n", bwParametersFile);
        for (int y=0; y<o.length;y++){
            writeLine("display gEUD_" +o[y].name+ " > "  + this.jobThreadID + "gEUD_"+o[y].name+".txt;\n",bwParametersFile); 
        }
        //writeLine("display gEUD_PTV > "+ Thread.currentThread().getName()+"gEUD_PTVHD.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Rectum > "+ Thread.currentThread().getName()+"gEUD_RECTUM.txt;\n", bwParametersFile);
        //writeLine("display gEUD_Bladder > "+ Thread.currentThread().getName()+"gEUD_BLADDER.txt;\n", bwParametersFile);
        bwParametersFile.close();
        
        //Running Process
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
        System.out.println("t = " + this.t );
        while ((line = br.readLine()) != null) {
            System.out.println(line);
        }
    }
    
    public void writeLine(String l, BufferedWriter bw) throws IOException{
	
	String row = "";
	row = l; 
	bw.write(row);
	

    }
   
}

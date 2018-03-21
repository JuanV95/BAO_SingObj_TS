package IMRT_Base;
/* LogFunction
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Comparator;
/**
 *
 * @author gcab623
 */
public class TreatmentPlan{
    public int beams;
    public int beamlets;
    public double slope;
    public double[] intensity;
    public Beam[] selAngles;
    public double[] weights;
    public double[] gEUD;
    public double singleObjectiveValue;
    public boolean visited;
    public double Vx[][];
    public double Dx[][];
    public boolean output;
    
    public TreatmentPlan(int b, int bl, int a, int o, boolean outputl){
        this.beams=b;
        this.beamlets=bl;
        this.intensity = new double[bl];
        this.selAngles = new Beam[a];
        this.gEUD=new double[o+1]; // # of Organs + PTV UB
        this.weights = new double[o]; // # of Organs + PTV UB
        this.Vx=new double[o][]; // # of Organs
        this.Dx=new double[o][]; // # of Organs
        this.singleObjectiveValue=0;
        this.visited=false;
        this.slope=0;
        this.output=outputl;
    }
    
    public void updateSol(TreatmentPlan s){
        this.beamlets=s.beamlets;
        this.beams=s.beams;
        this.singleObjectiveValue = s.singleObjectiveValue;
        this.intensity = new double[this.beamlets];
        this.selAngles = new Beam[this.beams];
        this.visited = s.visited;
        this.gEUD=new double[s.gEUD.length];
        System.arraycopy(s.intensity, 0, this.intensity, 0, 
                Math.min(s.intensity.length, this.intensity.length));
        System.arraycopy(s.selAngles, 0, this.selAngles, 0, this.beams);
        System.arraycopy(s.gEUD, 0, this.gEUD, 0, s.gEUD.length);
        System.arraycopy(s.weights, 0, this.weights, 0, s.weights.length);
        for (int i=0;i<s.Dx.length;i++){
            if (null != s.Dx[i]){
                this.Dx[i] = new double[s.Dx[i].length];
                System.arraycopy(s.Dx[i], 0, this.Dx[i], 0, 
                Math.min(s.Dx[i].length, this.Dx[i].length));
            }
        }
        for (int i=0;i<s.Vx.length;i++){
            if (null != s.Vx[i]){
                this.Vx[i] = new double[s.Vx[i].length];
                System.arraycopy(s.Vx[i], 0, this.Vx[i], 0, 
                Math.min(s.Vx[i].length, this.Vx[i].length));
            }
        }
        this.slope=s.slope;
    }    
    
    public void generateReferencePoint(Organs[] o, int opt, String jobThreadID ) throws IOException{
        
        switch (opt){
            case 1: 
                this.solveLogFunction(o, jobThreadID); 
                break;
            case 2: 
                this.solveLexiRectum(o, jobThreadID); 
                break;
            case 4: 
                this.solveVariableWeights(o, jobThreadID); 
                break;
            case 5: 
                this.solveConvexLogFunction(o, jobThreadID); 
                break;
            case 6: 
                this.solveWeightedSum(o, jobThreadID); 
                break;  
            case 7: 
                this.solveAdaptiveWeightedSum(o, jobThreadID); 
                break;  
        }
                
            
    }
    
    //public void solveLogFunction(Organs_bkp[] o, DDM_bkp M) throws IOException{
    public void solveLogFunction(Organs[] o, String jobThreadID) throws IOException{
      double[] EUD0 = new double[o.length];
      
      for (int i=0; i<o.length; i++ ){
          if (o[i].isTarget){
              EUD0[i]=o[i].doseLB; //Target
          }else{
              EUD0[i]=o[i].doseUB; //OAR
          }
      }
      
      AMPL_Solver logFunction = new AMPL_Solver(o, this, 0, jobThreadID, output);
       /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
      double auxSumIntensity=0;
      for(int i=0;i<this.intensity.length;i++){
          auxSumIntensity =+ this.intensity[i];
      }
      if (auxSumIntensity > 0){
          System.out.println("Generating parameters File (0)");
          logFunction.generateParametersFile_logFunction(this.intensity);
      }else{
          System.out.println("Generating parameters File (1)");
          logFunction.generateParametersFile_logFunction();
      }
      logFunction.runLogFunction_Solver(o);
      logFunction.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(logFunction.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o,jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.getLogFunctionValue(o);
  }
  
   public void solveVariableWeights(Organs[] o, String jobThreadID) throws IOException{
      double[] EUD0 = new double[o.length];
      for (int i=0; i<o.length; i++ ){
          if (o[i].isTarget){
              EUD0[i]=o[i].doseLB; //Target
          }else{
              EUD0[i]=o[i].doseUB; //OAR
          }
      }
      /*******************/
      /*Determine Weights*/
      /*******************/
      
      Random r = new Random();
      this.weights[1] = 0;
      while (this.weights[1] == 0 || this.weights[1] == 1){
        this.weights[1] = Math.floor(100.0 * r.nextDouble()) /100.0;
      }
      this.weights[2] = 1 - this.weights[1];
      
      AMPL_Solver variableWeights = new AMPL_Solver(o, this, 0, jobThreadID, output);
      variableWeights.wPar[1] = this.weights[1];
      variableWeights.wPar[2] = this.weights[2];
       /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
      double auxSumIntensity=0;
      for(int i=0;i<this.intensity.length;i++){
          auxSumIntensity =+ this.intensity[i];
      }
      if (auxSumIntensity > 0){
          System.out.println("Generating parameters File (0)");
          variableWeights.generateParametersFile_variableWeights(this.intensity);
      }else{
          System.out.println("Generating parameters File (1)");
          variableWeights.generateParametersFile_variableWeights();
      }
      variableWeights.runVariableWeights_Solver(o);
      variableWeights.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(variableWeights.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o,jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.getVariableWeightsValue(o);
  }

   public void solveAdaptiveWeightedSum(Organs[] o, String jobThreadID) throws IOException{
      double[] EUD0 = new double[o.length];
      for (int i=0; i<o.length; i++ ){
          if (o[i].isTarget){
              EUD0[i]=o[i].doseLB; //Target
          }else{
              EUD0[i]=o[i].doseUB; //OAR
          }
      }
      
      AMPL_Solver adaptiveWeightedSum = new AMPL_Solver(o, this, 0, jobThreadID, output);
      adaptiveWeightedSum.wPar[1] = this.weights[1];
      adaptiveWeightedSum.wPar[2] = this.weights[2];
       /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
      double auxSumIntensity=0;
      for(int i=0;i<this.intensity.length;i++){
          auxSumIntensity =+ this.intensity[i];
      }
      if (auxSumIntensity > 0){
          System.out.println("Generating parameters File (0)");
          adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum(this.intensity);
      }else{
          System.out.println("Generating parameters File (1)");
          adaptiveWeightedSum.generateParametersFile_adaptiveWeightedSum();
      }
      adaptiveWeightedSum.runAdaptiveWeightedSum_Solver(o);
      adaptiveWeightedSum.getSolution();
      System.arraycopy(adaptiveWeightedSum.x, 0, this.intensity, 0, this.beamlets);
      //Random r = new Random();
      //this.intensity=new double[this.beamlets];
      //for(int i=0;i<this.intensity.length;i++){
      //    this.intensity[i] = r.nextDouble();
      //}
      
      //calculate gEUDs
      System.arraycopy(getEUD(o,jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      //this.gEUD[0]=70;
      //this.gEUD[1] = 51.0 + r.nextDouble() * 8.0;
      //this.gEUD[2] = 31.0 + r.nextDouble() * 5.0;
      
      
      
      this.getAdaptiveWeightedSumValue(o);
  }

  public void solveConvexLogFunction(Organs[] o, String jobThreadID) throws IOException{
      double[] EUD0 = new double[o.length];
      
      for (int i=0; i<o.length; i++ ){
          if (o[i].isTarget){
              EUD0[i]=o[i].doseLB; //Target
          }else{
              EUD0[i]=o[i].doseUB; //OAR
          }
      }
      
      AMPL_Solver convexLogFunction = new AMPL_Solver(o, this, 0, jobThreadID, output);
       /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
      double auxSumIntensity=0;
      for(int i=0;i<this.intensity.length;i++){
          auxSumIntensity =+ this.intensity[i];
      }
      if (auxSumIntensity > 0){
          System.out.println("Generating parameters File (0)");
          convexLogFunction.generateParametersFile_convexLogFunction(this.intensity);
      }else{
          System.out.println("Generating parameters File (1)");
          convexLogFunction.generateParametersFile_convexLogFunction();
      }
      convexLogFunction.runConvexLogFunction_Solver(o);
      convexLogFunction.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(convexLogFunction.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o,jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.getConvexLogFunctionValue(o);
  } 
  
  public void solveWeightedSum(Organs[] o, String jobThreadID) throws IOException{
      double[] EUD0 = new double[o.length];
      
      for (int i=0; i<o.length; i++ ){
          if (o[i].isTarget){
              EUD0[i]=o[i].doseLB; //Target
          }else{
              EUD0[i]=o[i].doseUB; //OAR
          }
      }
      
      AMPL_Solver weightedSum = new AMPL_Solver(o, this, 0, jobThreadID,output);
       /*##################################################
       # "generateParametersFile" method  per each option #
       # "generateParametersFile" method  per each option #
       #################################################*/
      double auxSumIntensity=0;
      for(int i=0;i<this.intensity.length;i++){
          auxSumIntensity =+ this.intensity[i];
      }
      if (auxSumIntensity > 0){
          System.out.println("Generating parameters File (0)");
          weightedSum.generateParametersFile_weightedSum(this.intensity);
      }else{
          System.out.println("Generating parameters File (1)");
          weightedSum.generateParametersFile_weightedSum();
      }
      weightedSum.runWeightedSum_Solver(o);
      weightedSum.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(weightedSum.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o,jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.getWeightedSumValue(o);
  }
  
  public void solveTestFunction(Organs[] o, String jobThreadID) throws IOException{
      Random r = new Random();
      double aux;
      aux=0;
      this.intensity=new double[this.beamlets];
      for(int i=0;i<this.beamlets;i++){
          this.intensity[i] = r.nextDouble();
      }
      this.gEUD[0]= o[0].desiredDose;
      aux=(1+ r.nextDouble());
      this.gEUD[1]= o[1].desiredDose * aux;
      aux=(1+ r.nextDouble());
      this.gEUD[2]= o[2].desiredDose * aux;
      aux=r.nextDouble();
      this.gEUD[3]= o[0].doseUB * aux;
      aux=1+ r.nextDouble();
      this.singleObjectiveValue = aux * 5;
  }
    
  //public void getLogFunctionValue(Organs_bkp[] o) throws IOException{
  public void getLogFunctionValue(Organs[] o) throws IOException{
      double objectiveValue=0;
      double EUD0, aux_gEUD, v;
      for (int i =1;i<o.length; i++){
          EUD0 = o[i].desiredDose;
          aux_gEUD = this.gEUD[i];
          v=o[i].v;
          objectiveValue = objectiveValue - Math.log(Math.pow(1+Math.pow(aux_gEUD/EUD0, v),-1));
      }
      this.singleObjectiveValue = objectiveValue;
      System.out.println("Log function Value = " + this.singleObjectiveValue );
  }
  
  public void getVariableWeightsValue(Organs[] o) throws IOException{
      double objectiveValue=0;
      double aux_gEUD, aux_weights;
      for (int i =1;i<o.length; i++){
          aux_gEUD = this.gEUD[i];
          aux_weights = this.weights[i];
          objectiveValue = objectiveValue + (aux_weights * aux_gEUD);
      }
      this.singleObjectiveValue = objectiveValue;
      System.out.println("Weighted Sum with variable weights = " + this.singleObjectiveValue );
  }
  
  
  public void getAdaptiveWeightedSumValue(Organs[] o) throws IOException{
      double objectiveValue=0;
      double aux_gEUD, aux_weights;
      for (int i =1;i<o.length; i++){
          aux_gEUD = this.gEUD[i];
          aux_weights = this.weights[i];
          objectiveValue = objectiveValue + (aux_weights * aux_gEUD);
      }
      this.singleObjectiveValue = objectiveValue;
      System.out.println("Adaptive Weighted Sum with variable weights = " + this.singleObjectiveValue );
  }
  
  public void getConvexLogFunctionValue(Organs[] o) throws IOException{
      double objectiveValue=0;
      double EUD0, aux_gEUD, v;
      for (int i =1;i<o.length; i++){
          EUD0 = o[i].desiredDose;
          aux_gEUD = this.gEUD[i];
          v=o[i].v;
          objectiveValue = objectiveValue - Math.log(Math.pow(Math.pow(aux_gEUD/EUD0, v),-1));
      }
      this.singleObjectiveValue = objectiveValue;
      System.out.println("Inverse function Value = " + this.singleObjectiveValue );
  }
  
  public void getWeightedSumValue(Organs[] o) throws IOException{
      double objectiveValue=0;
      double weight, aux_gEUD, v;
      for (int i =1;i<o.length; i++){
          if(!o[i].isTarget){
            weight = o[i].weight;
            aux_gEUD = this.gEUD[i];
            objectiveValue = objectiveValue + aux_gEUD * weight;
          }
      }
      this.singleObjectiveValue = objectiveValue;
      System.out.println("Weighted Sum Value = " + this.singleObjectiveValue );
  }
  
  public void getVxDx(String jobThreadID, Organs[] o, int[][] valueVx, int[][] valueDx) throws FileNotFoundException, IOException{
        // load dose vector (dose deposited at each voxel
        for (Organs o1 : o) {
            String dir = jobThreadID + "DVH_" + o1.name + ".txt";
            double[] dvh = new double[o1.totalVoxels];
            System.out.println(dir + "\t" + o1.totalVoxels);
            String[] auxReader=null;
            File f = new File(dir);
            if (f.exists()) {
                BufferedReader fileIn = new BufferedReader(new FileReader(f));
                String line = "";
                line=fileIn.readLine(); //avoid first line;
                line=fileIn.readLine();
                auxReader = line.split(" ");
                while (!";".equals(auxReader[0])){
                    for (int j=0;j<auxReader.length;j++) {
                        if (!"".equals(auxReader[j])){
                            int auxIndex = (int)Double.parseDouble(auxReader[j]);
                            while("".equals(auxReader[j+1])){
                                j++;
                            }
                            dvh[auxIndex-1] = Double.parseDouble(auxReader[j+1]);
                            j++;
                        }
                    }
                    line=fileIn.readLine();
                    auxReader = line.split(" ");
                }
            }else{
                System.out.println("ERROR: CurrentSol file wasn't generated ("+dir+")");
            }
            Arrays.sort(dvh);
            //calculate Vx indexes
            if (valueVx[o1.index] != null) {
                this.Vx[o1.index] = new double [valueVx[o1.index].length];
                for (int j = 0; j < valueVx[o1.index].length; j++) {
                    double limit = (valueVx[o1.index][j] * o1.desiredDose) / 100;
                    double count = 0;
                    for (int k = o1.totalVoxels - 1; k>=0; k--) {
                        if (dvh[k]>limit){
                            count++;
                        }else{
                            break;
                        }
                    }
                    this.Vx[o1.index][j] = (count / o1.totalVoxels) * 100;
                }
            }
            //calculate Dx indexes
            if (valueDx[o1.index] != null) {
                this.Dx[o1.index] = new double [valueDx[o1.index].length];
                for (int j = 0; j < valueDx[o1.index].length; j++) {
                    int limit = o1.totalVoxels - (int) (valueDx[o1.index][j] * o1.totalVoxels) / 100;
                    this.Dx[o1.index][j] = dvh[limit];
                }
            }
        }
    }
  
  public void getSlope(String jobID, int lineNum) throws FileNotFoundException, IOException{
      String dir = jobID + "lagrangeMultiplier.txt";
      File f = new File(dir);
      BufferedReader brLM = new BufferedReader(new FileReader(f));
      String line = "";
      String[] auxReader;  
      for (int i=0;i<lineNum;i++){
        line=brLM.readLine();
      }
      auxReader = line.split("  ");
      this.slope = 1/Double.parseDouble(auxReader[auxReader.length-1]);
  }
  
  //public void solveLexiRectum(Organs_bkp[] o, DDM_bkp M) throws IOException{
  public void solveLexiRectum(Organs[] o, String jobThreadID) throws IOException{
      AMPL_Solver lexRect = new AMPL_Solver(o, this, 0, jobThreadID,output);
      
      if (this.intensity.length > 0){
          lexRect.generateParametersFile(this.intensity);
      }else{
          lexRect.generateParametersFile();
      }
      lexRect.runLexRectum_Solver(o);
      lexRect.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(lexRect.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
      this.getSlope(jobThreadID,3);
  }
  //public void solveLexiBladder(Organs_bkp[] o, DDM_bkp M) throws IOException{
  public void solveLexiBladder(Organs[] o, String jobThreadID) throws IOException{
      AMPL_Solver lexBladder = new AMPL_Solver(o, this, 0, jobThreadID, output);
      if (this.intensity.length > 0){
          lexBladder.generateParametersFile(this.intensity);
      }else{
          lexBladder.generateParametersFile();
      }
      lexBladder.runLexBladder_Solver(o);
      lexBladder.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(lexBladder.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
      this.getSlope(jobThreadID,3);
  }
  
  public void solve_gEUD_Rectum(Organs[] o, double epsilon, String jobThreadID) throws IOException{
      AMPL_Solver gEUD_Model = new AMPL_Solver(o, this, epsilon, jobThreadID, output);
      if (this.intensity.length > 0){
          gEUD_Model.generateParametersFile(this.intensity);
      }else{
          gEUD_Model.generateParametersFile();
      }
      gEUD_Model.run_gEUD_Rectum_Solver(o);
      gEUD_Model.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(gEUD_Model.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
      this.getSlope(jobThreadID,4);
  }
  
  public void solve_gEUD_Bladder(Organs[] o, double epsilon, String jobThreadID) throws IOException{
      AMPL_Solver gEUD_Model = new AMPL_Solver(o, this, epsilon, jobThreadID, output);
      if (this.intensity.length > 0){
          gEUD_Model.generateParametersFile(this.intensity);
      }else{
          gEUD_Model.generateParametersFile();
      }
      gEUD_Model.run_gEUD_Bladder_Solver(o);
      gEUD_Model.getSolution();
      this.intensity=new double[this.beamlets];
      System.arraycopy(gEUD_Model.x, 0, this.intensity, 0, this.beamlets);
      //calculate gEUDs
      System.arraycopy(getEUD(o, jobThreadID), 0, this.gEUD, 0, this.gEUD.length);
      this.singleObjectiveValue = this.gEUD[1]; //gEUD Rectum
      this.getSlope(jobThreadID,4);
  }
  

  public double[] getEUD(Organs[] o, String jobThreadID) throws FileNotFoundException, IOException{
      File[] f= new File[o.length];
      BufferedReader[] fileIn= new BufferedReader[o.length];
      double[] aux_gEUD = new double[this.gEUD.length];
      for (int y=0;y<o.length;y++){
          f[y]= new File("./"+jobThreadID+"gEUD_" + o[y].name + ".txt");
          fileIn[y]= new BufferedReader(new FileReader(f[y]));
          String line = fileIn[y].readLine();
          line = fileIn[y].readLine();          
          String[] auxReader = null;
          if (line!=null){
              auxReader = line.split(" = ");
               aux_gEUD[y] = Double.parseDouble(auxReader[1]);
              //System.out.println(o[y].name +" : " +  aux_gEUD[y] + " [Gy] ");
          }else{
              System.out.println("error: No gEUD value for: " + o[y].name + "(./"+jobThreadID+"gEUD_" + o[y].name + ".txt)");
          }
          fileIn[y].close();
          if (o[y].isTarget){
            // get gEUD for OAR-PTV
            File g = new File("./"+jobThreadID+"gEUD_" + o[y].name + "_UB.txt");
            BufferedReader gIn= new BufferedReader(new FileReader(g));
            line = gIn.readLine();
            line = gIn.readLine();
            auxReader = null;
            if (line!=null){
              auxReader = line.split(" = ");
              aux_gEUD[this.gEUD.length-1] = Double.parseDouble(auxReader[1]);
              //System.out.println(o[y].name +" : " +  aux_gEUD[y] + " [Gy] ");
            }else{
              System.out.println("error: No gEUD value for OAR-PTV (./"+jobThreadID+"gEUD_" + o[y].name + "_UB.txt)");
            }
          }
      }
      return aux_gEUD;
  }
  
  public boolean sameBAC(Beam [] bac1){
        boolean flag=false;
            for (int j=0;j<bac1.length;j++){
                flag=false; //Assuming bac1 is not equal to bac2
                for (int k=0;k<this.selAngles.length;k++){
                    if (bac1[j].index==this.selAngles[k].index){
                        flag=true;
                        break;
                    }
                }
                if (flag==false){
                    break;
                }
            }
        return flag;
    }
  
  public void printSol(String dirFile) throws IOException{
      BufferedWriter bwFile = null ;
      File solFile = new File(dirFile);
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(solFile));
        }
        for (int i=0; i< this.beams;i++){  
            writeLine(this.selAngles[i].index + "\t", bwFile);
        }
        for (int j=0; j< this.gEUD.length;j++){
            writeLine(this.gEUD[j] + "\t", bwFile);
        }
        writeLine(this.singleObjectiveValue + "\t", bwFile);
        
        for (int j=0; j< this.weights.length;j++){
            writeLine(this.weights[j] + "\t", bwFile);
        }
        
        for (double[] Vx1 : this.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    writeLine(Vx1[j] + "\t", bwFile);
                }
            }
        }
        for (double[] Dx1 : this.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    writeLine(Dx1[j] + "\t", bwFile);
                }
            }
        }
        
        writeLine("x:\t", bwFile);
        
        for (int j=0; j< this.intensity.length;j++){
            writeLine(this.intensity[j] + "\t", bwFile);
        }
        writeLine("\n", bwFile);
        bwFile.close();
  }
  
  public void printSol(String dirFile, int option) throws IOException{
      BufferedWriter bwFile = null ;
      File solFile = new File(dirFile);
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(solFile));
        }
        for (int i=0; i< this.beams;i++){  
            writeLine(this.selAngles[i].index + "\t", bwFile);
        }
        for (int j=0; j< this.gEUD.length;j++){
            writeLine(this.gEUD[j] + "\t", bwFile);
        }
        writeLine(this.singleObjectiveValue + "\t", bwFile);
        
        for (int j=0; j< this.weights.length;j++){
            writeLine(this.weights[j] + "\t", bwFile);
        }
        
        for (double[] Vx1 : this.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    writeLine(Vx1[j] + "\t", bwFile);
                }
            }
        }
        for (double[] Dx1 : this.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    writeLine(Dx1[j] + "\t", bwFile);
                }
            }
        }
        writeLine("\n", bwFile);
        bwFile.close();
  }
  
  public void printSol(String dirFile, long time) throws IOException{
      BufferedWriter bwFile = null ;
      File solFile = new File(dirFile);
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(solFile));
        }
        for (int i=0; i< this.beams;i++){  
            writeLine(this.selAngles[i].index + "\t", bwFile);
        }
        for (int j=0; j< this.gEUD.length;j++){
            writeLine(this.gEUD[j] + "\t", bwFile);
        }
        writeLine(this.singleObjectiveValue + "\t", bwFile);
        writeLine(time + "\t", bwFile);
        //writeLine("x[]:\t", bwFile);
        
        for (int j=0; j< this.weights.length;j++){
            writeLine(this.weights[j] + "\t", bwFile);
        }
        
        for (double[] Vx1 : this.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    writeLine(Vx1[j] + "\t", bwFile);
                }
            }
        }
        for (double[] Dx1 : this.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    writeLine(Dx1[j] + "\t", bwFile);
                }
            }
        }
        
        for (int j=0; j< this.intensity.length;j++){
            writeLine(this.intensity[j] + "\t", bwFile);
        }
        writeLine("\n", bwFile);
        bwFile.close();
  }
  public boolean isDominated(TreatmentPlan s) {
      boolean isDominated=true; //means, s dominates 'this'
      //We assume that first index in gEDU  corresponds to the target
      //which is equal for all the solutions
      for (int i=1;i < this.gEUD.length-1; i++){
            if (this.gEUD[i]<s.gEUD[i]){
                isDominated=false;
                break;
            }
      }
      return (isDominated);
  }
  
  public void computeNewWeights(TreatmentPlan sol, ArrayList<TreatmentPlan> convexHull, double[] bounds, Organs[] o ){
      double[] w = new double[sol.weights.length];
      double[] vertex1 = new double[3];double[] vertex2 = new double[3];
      double[] currentPoint = new double[3];
      double[] originalWeights1 = new double[3];double[] originalWeights2 = new double[3];
      double distance = 10000, auxDist; //big number
      boolean weightsDone = false;
      Random r = new Random();
      
        switch (convexHull.size()) {
            case 1://We have only one non-dominated point --> We return same weights
                w[0] = convexHull.get(0).weights[0]; //we don't use this anyway  
                w[1] = convexHull.get(0).weights[1];
                w[2] = convexHull.get(0).weights[2];
                weightsDone = true;
                System.out.println("There is only one solution in the convex hull.");
                System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            +      "New Weights: ["      + w[1]           + ", " + w[2]           + "] ");
                break;
            case 2://We have two points in the convex hull (we might have more non-dominated points though). 
                System.arraycopy(convexHull.get(0).gEUD, 0, vertex1, 0, vertex1.length);
                System.arraycopy(convexHull.get(1).gEUD, 0, vertex2, 0, vertex2.length);
                System.arraycopy(convexHull.get(0).weights, 0, originalWeights1, 0, originalWeights1.length);
                System.arraycopy(convexHull.get(1).weights, 0, originalWeights2, 0, originalWeights2.length);
                
                //determine whether this solution is supported or not (i.e. it belongs to the convex hull)
                if (sol.sameBAC(convexHull.get(0).selAngles)){
                    double[] wLeg1 = new double[3];double[] wLeg2 = new double[3];        
                    //In this case, wLeg1 is a vertical line, i.e. [1,0]
                    wLeg1[0]= sol.weights[0]; //we don't use this anyway  
                    wLeg1[1] = bounds[1];
                    wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]); 
                    wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                    wLeg1[2] = (double) 1-wLeg1[1]; 
                    wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000) / 1000);
                    if(wLeg1[1] + wLeg1[2] != 1 ){ 
                        double dif = 1 - (wLeg1[1] + wLeg1[2]); 
                        wLeg1[2] = wLeg1[2] + dif; 
                    }
                    //wLeg1[1]= 1; wLeg1[2]= 0;  
                    if(wLeg1[1] + wLeg1[2] != 1 ){
                        double dif = 1 - (wLeg1[1] + wLeg1[2]);
                        wLeg1[2] = wLeg1[2] + dif;
                    }
                    
                    wLeg2[0]= sol.weights[0]; //we don't use this anyway
                    //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                    wLeg2[1] = Math.min(bounds[1], ((vertex2[2]-vertex1[2])/((vertex1[1]-vertex2[1])-(vertex1[2]-vertex2[2]))));
                    wLeg2[1] = Math.max(1 - bounds[2], wLeg2[1]);
                    wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                    wLeg2[2] = (double) 1-wLeg2[2];
                    wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                    if(wLeg2[1] + wLeg2[2] != 1 ){
                        double dif = 1 - (wLeg2[1] + wLeg2[2]);
                        wLeg2[2] = wLeg2[2] + dif;
                    }
                    
                    //Since this solution is a supported one (i.e. it belongs to the convex hull), we calculate the
                    //average of the weights of legs 1 and 2
                    
                    w[0]= sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], (wLeg1[1]+wLeg2[1])/2);
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                    w[2] = (double) 1-w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                    if(w[1] + w[2] != 1 ){
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    System.out.println("This is a Supported Solution (index = 0)");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "wLegs (" + wLeg1[1] + ", " + wLeg1[2] + ") and "
                                        + "(" + wLeg2[1] + ", " + wLeg2[2] + ") "
                             + "New Weights: [" + w[1] +                ", " + w[2]                + "] ");
                    weightsDone = true;
                }else if(sol.sameBAC(convexHull.get(1).selAngles)){
                double[] wLeg1 = new double[3];double[] wLeg2 = new double[3];        
                    wLeg2[0]= sol.weights[0]; //we don't use this anyway
                    //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                    wLeg2[2] = bounds[2];
                    wLeg2[2] = Math.max(1 - bounds[1], wLeg2[2]);
                    wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                    wLeg2[1] = (double) 1-wLeg2[2];
                    wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                    if(wLeg2[1] + wLeg2[2] != 1 ){
                        double dif = 1 - (wLeg2[1] + wLeg2[2]);
                        wLeg2[2] = wLeg2[2] + dif;
                    }
                    
                    wLeg1[0]= sol.weights[0]; //we don't use this anyway
                    //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                    wLeg1[1] = Math.min(bounds[1], ((vertex2[2]-vertex1[2])/((vertex1[1]-vertex2[1])-(vertex1[2]-vertex2[2]))));
                    wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                    wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000) / 1000);
                    wLeg1[2] = (double) 1-wLeg1[2];
                    wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000) / 1000);
                    if(wLeg1[1] + wLeg1[2] != 1 ){
                        double dif = 1 - (wLeg1[1] + wLeg1[2]);
                        wLeg1[2] = wLeg1[2] + dif;
                    }
                    
                    //Since this solution is a supported one (i.e. it belongs to the convex hull), we calculate the
                    //average of the weights of legs 1 and 2
                    
                    w[0]= sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], (wLeg1[1]+wLeg2[1])/2);
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000) / 1000);
                    w[2] = (double) 1-w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000) / 1000);
                    if(w[1] + w[2] != 1 ){
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    weightsDone = true;
                    System.out.println("This is a Supported Solution (index = 1)");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "wLegs (" + wLeg1[1] + ", " + wLeg1[2] + ") and "
                                        + "(" + wLeg2[1] + ", " + wLeg2[2] + ") "
                             + "New Weights: [" + w[1] +                ", " + w[2]                + "] ");
                }else{ //case for unsupported solutions
                    //Since there is only one leg ('cause there are only two supported solutions) all the
                    //unsupported solutions will use the same weights
                    w[0]= sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], ((vertex2[2]-vertex1[2])/((vertex1[1]-vertex2[1])-(vertex1[2]-vertex2[2]))));
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                    w[2] = (double) 1-w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                    if(w[1] + w[2] != 1 ){
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    System.out.println("This is an unsupported Solution, and there is only one segment in the convex hull.");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "Vertex Weights (" + originalWeights1[1] + ", " + originalWeights1[2] + ") and "
                                        + "(" + originalWeights1[1] + ", " + originalWeights1[2] + ") "
                             + "New Weights (segment weights): [" + w[1] +                ", " + w[2]                + "] ");
                    weightsDone = true;
                    
                }
                break;
            default:
                //First, check whether sol belongs to the convex hull
                for(int i=0;i<convexHull.size();i++){
                    if (sol.sameBAC(convexHull.get(i).selAngles)){
                        double[] wLeg1 = new double[3];double[] wLeg2 = new double[3];
                        if (i==0){ //if sol is the first point in the convex hull
                            System.arraycopy(convexHull.get(0).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(1).gEUD, 0, vertex2, 0, vertex2.length);
                            //System.arraycopy(convexHull.get(0).weights, 0, originalWeights1, 0, originalWeights1.length);
                            //System.arraycopy(convexHull.get(1).weights, 0, originalWeights2, 0, originalWeights2.length);
                
                            wLeg1[0]= sol.weights[0]; //we don't use this anyway
                            wLeg1[1] = bounds[1];
                            wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                            wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                            wLeg1[2] = (double) 1-wLeg1[1];
                            wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000000) / 1000000);
                            if(wLeg1[1] + wLeg1[2] != 1 ){
                                double dif = 1 - (wLeg1[1] + wLeg1[2]);
                                wLeg1[2] = wLeg1[2] + dif;
                            }
                                   
                            wLeg2[0]= sol.weights[0]; //we don't use this anyway
                            wLeg2[1] = Math.min(bounds[1], ((vertex2[2]-vertex1[2])/((vertex1[1]-vertex2[1])-(vertex1[2]-vertex2[2]))));
                            wLeg2[1] = Math.max(1 - bounds[2], wLeg2[1]);
                            wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                            wLeg2[2] = (double) 1-wLeg2[2];
                            wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                            if(wLeg2[1] + wLeg2[2] != 1 ){
                                double dif = 1 - (wLeg2[1] + wLeg2[2]);
                                wLeg2[2] = wLeg2[2] + dif;
                            }
                            System.out.println("This is a Supported Solution (index = 0)");
                            
                        }else if(i == convexHull.size()-1){//if sol is the last point in the convex hull    
                            System.arraycopy(convexHull.get(convexHull.size()-2).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(convexHull.size()-1).gEUD, 0, vertex2, 0, vertex2.length);
                            //System.arraycopy(convexHull.get(0).weights, 0, originalWeights1, 0, originalWeights1.length);
                            
                
                            wLeg2[0]= sol.weights[0]; //we don't use this anyway
                            //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                            wLeg2[2] = bounds[2];
                            wLeg2[2] = Math.max(1 - bounds[1], wLeg2[2]);
                            wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                            wLeg2[1] = (double) 1-wLeg2[2];
                            wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                            if(wLeg2[1] + wLeg2[2] != 1 ){
                                double dif = 1 - (wLeg2[1] + wLeg2[2]);
                                wLeg2[2] = wLeg2[2] + dif;
                            }
                    
                            wLeg1[0]= sol.weights[0]; //we don't use this anyway
                            //wLeg1[1] = Math.min(bounds[1], ((vertex1[2]-extreme1[2])/((extreme1[1]-vertex1[1])-(extreme1[2]-vertex1[2]))));
                            wLeg1[1] = Math.min(bounds[1], ((vertex2[2]-vertex1[2])/((vertex1[1]-vertex2[1])-(vertex1[2]-vertex2[2]))));
                            wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                            wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                            wLeg1[2] = (double) 1-wLeg1[2];
                            wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000000) / 1000000);
                            if(wLeg1[1] + wLeg1[2] != 1 ){
                                double dif = 1 - (wLeg1[1] + wLeg1[2]);
                                wLeg1[2] = wLeg1[2] + dif;
                            }
                            System.out.println("This is a Supported Solution (index = " + (convexHull.size()-1) +")");
                        }else{ //sol belongs to the convex hull but it is neither the first nor the last point in the convex hull
                               
                            System.arraycopy(convexHull.get(i-1).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(i).gEUD, 0, currentPoint, 0, currentPoint.length);
                            System.arraycopy(convexHull.get(i+1).gEUD, 0, vertex2, 0, vertex2.length);
                            //System.arraycopy(convexHull.get(i-1).weights, 0, originalWeights1, 0, originalWeights1.length);
                            //System.arraycopy(convexHull.get(i+1).weights, 0, originalWeights2, 0, originalWeights2.length);
                            
                            wLeg1[0]= sol.weights[0]; //we don't use this anyway
                            wLeg1[1] = Math.min(bounds[1], ((currentPoint[2]-vertex1[2])/((vertex1[1]-currentPoint[1])-(vertex1[2]-currentPoint[2]))));
                            wLeg1[1] = Math.max(1 - bounds[2], wLeg1[1]);
                            wLeg1[1] = (double) (Math.floor(wLeg1[1] * 1000000) / 1000000);
                            wLeg1[2] = (double) 1-wLeg1[1];
                            wLeg1[2] = (double) (Math.floor(wLeg1[2] * 1000000) / 1000000);
                            if(wLeg1[1] + wLeg1[2] != 0 ){
                                double dif = 1 - (wLeg1[1] + wLeg1[2]);
                                wLeg1[2] = wLeg1[2] + dif;
                            }
                            
                            wLeg2[0]= sol.weights[0]; //we don't use this anyway
                            wLeg2[1] = Math.min(bounds[1], ((vertex2[2]-currentPoint[2])/((currentPoint[1]-vertex2[1])-(currentPoint[2]-vertex2[2]))));
                            wLeg2[1] = Math.max(1 - bounds[2], wLeg2[1]);
                            wLeg2[1] = (double) (Math.floor(wLeg2[1] * 1000000) / 1000000);
                            wLeg2[2] = (double) 1-wLeg2[1];
                            wLeg2[2] = (double) (Math.floor(wLeg2[2] * 1000000) / 1000000);
                            if(wLeg2[1] + wLeg2[2] != 0 ){
                                double dif = 1 - (wLeg2[1] + wLeg2[2]);
                                wLeg2[2] = wLeg2[2] + dif;
                            }
                            System.out.println("This is a Supported Solution (index = " + i +")");
                        }
                        //Since this solution is a supported one (i.e. it belongs to the convex hull), we calculate the
                        //average of the weights of legs 1 and 2

                        w[0]= sol.weights[0]; //we don't use this anyway
                        w[1] = Math.min(bounds[1], (wLeg1[1]+wLeg2[1])/2);
                        w[1] = Math.max(1 - bounds[2], w[1]);
                        w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                        w[2] = (double) 1-w[1];
                        w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                        if(w[1] + w[2] != 1 ){
                            double dif = 1 - (w[1] + w[2]);
                            w[2] = w[2] + dif;
                        }
                        System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "wLegs (" + wLeg1[1] + ", " + wLeg1[2] + ") and "
                                        + "(" + wLeg2[1] + ", " + wLeg2[2] + ") "
                             + "New Weights: [" + w[1] +                ", " + w[2]                + "] ");
                        weightsDone = true;
                    }
                }        
                if (!weightsDone){ //means this solution (sol) in unsuported (i.e. it does not belong to the convex hull)
                    //locate the closest segment of the convex hull to sol
                    for(int i=0;i<convexHull.size()-1;i++){
                        double[] p0 = new double[3];double[] p1 = new double[3];double[] p2 = new double[3];
                        System.arraycopy(convexHull.get(i).gEUD, 0, p0, 0, p0.length);
                        System.arraycopy(convexHull.get(i+1).gEUD, 0, p1, 0, p1.length);
                        System.arraycopy(sol.gEUD, 0, p2, 0, p2.length);
                        auxDist = Math.abs((p2[2]-p1[2])*p0[1] -(p2[1]-p1[1])*p0[2] + p2[1]*p1[2]-p2[2]*p1[1]);
                        auxDist = auxDist / Math.sqrt(Math.pow(p2[2]-p1[2],2) + Math.pow(p2[1]-p1[1],2));
                        if (auxDist < distance){
                            distance = auxDist;
                            System.arraycopy(convexHull.get(i).gEUD, 0, vertex1, 0, vertex1.length);
                            System.arraycopy(convexHull.get(i+1).gEUD, 0, vertex2, 0, vertex1.length);
                            System.arraycopy(convexHull.get(i).weights, 0, originalWeights1, 0, originalWeights1.length);
                            System.arraycopy(convexHull.get(i+1).weights, 0, originalWeights2, 0, originalWeights2.length);
                        }
                    }
                    w[0]= sol.weights[0]; //we don't use this anyway
                    w[1] = Math.min(bounds[1], ((vertex2[2]-vertex1[2])/((vertex1[1]-vertex2[1])-(vertex1[2]-vertex2[2]))));
                    w[1] = Math.max(1 - bounds[2], w[1]);
                    w[1] = (double) (Math.floor(w[1] * 1000000) / 1000000);
                    w[2] = (double) 1-w[1];
                    w[2] = (double) (Math.floor(w[2] * 1000000) / 1000000);
                    if(w[1] + w[2] != 0 ){
                        double dif = 1 - (w[1] + w[2]);
                        w[2] = w[2] + dif;
                    }
                    System.out.println("This is an unsupported Solution.");
                    System.out.println("Original Weights: [" + sol.weights[1] + ", " + sol.weights[2] + "] "
                            + "Vertex Weights (" + originalWeights1[1] + ", " + originalWeights1[2] + ") and "
                                        + "(" + originalWeights1[1] + ", " + originalWeights1[2] + ") "
                             + "New Weights (segment weights): [" + w[1] +                ", " + w[2]                + "] ");
                }
                break;
        }
        for(int i=0;i<this.weights.length;i++){
          this.weights[i]=w[i];
        }
  }
  
  public static void writeLine(String l, BufferedWriter bw) throws IOException{
	
	bw.write(l);
    }

    public static Comparator<TreatmentPlan> gEUDRectum_Comparator = new Comparator<TreatmentPlan>() {
        public int compare(TreatmentPlan tp1,TreatmentPlan tp2){
            if(tp1.gEUD[1] < tp2.gEUD[1]){
                return -1;
            }else{ 
                return 1;
            }
        }
    };    
    
    public static Comparator<TreatmentPlan> judgementFunct_Comparator = new Comparator<TreatmentPlan>() {
        public int compare(TreatmentPlan tp1,TreatmentPlan tp2){
            if(tp1.singleObjectiveValue < tp2.singleObjectiveValue){
                return -1;
            }else{ 
                return 1;
            }
        }
    };
}


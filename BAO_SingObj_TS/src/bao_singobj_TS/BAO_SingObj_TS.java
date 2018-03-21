/*
 * This is a Local Search implementation to solve multi-objective BAO problem.
 * Solution structure consists of a vector of beams. 
 * Neighborhood move implies an increment (decrement) of one degree over
 * one of the beams in the current solution.
 * To decide among a set of neighbours we can calculate either values of sample
 * points (e.g. Lexicographic values) or gEUD-based functions (e.g. LogFunction).
 * To add any new decision criterion some few changes must be done in this file
 * as well as in the TreatmentPlan class.
 * For each iteration we obtain a local optima in term of the decision criterion.
 * At the end, the algorithm generates a list of files organized as follows:
 * 1.- Best_Solution_"CriterionName"_XX.txt: Shows the best neighbour obtained at each 
 * iteration. XX corresponds to the Local Search iteration.
 * 2.- generatedSolutions_"CriterionName"_XX: Shows all the beam angle cofigurations
 * visited.
 * The file fun_selectBAC.m is an Octave file that generates a file with all
 * the local optima. This file must be the input data for the FMO algorithm 
 * which is the one that will generate the non-dominated points for each
 * local optima obtained here.
 */

package bao_singobj_TS;

import IMRT_Base.*;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.text.DecimalFormat;
import java.util.Random;
import java.util.Scanner;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.*;
import java.util.concurrent.Callable;

/**
 *
 * @author guille
 */
public class BAO_SingObj_TS {
    public static int numOrgans;
    public static int numAngles;
    public static ArrayList<TreatmentPlan> generatedBACs; 
    public static int[][] Vx; //Vx index. % of the volume receiving more than x% of the prescribed dose
                              //One row per organ (usually only for target regions)
    public static int[][] Dx; //Dx index. Minimum dose of the hottest x% of the volume.
                              //One row per organ (usually only for OAR regions)
    
    //public static int poolSize = 3;
    public static Organs[] o;
    public static ArrayList<TreatmentPlan> visitedBACs; 
    public static double[][] ndp; 
    public static ArrayList<TreatmentPlan> localOptima;
    //public static ArrayList<Long[]> voxelIndex[]; 
    //public static DDM M;
    public static int option;
    public static boolean randomInitial, parallelNeighbourhoodGeneration;
    public static int[][] initialBACs;
    public static String pathFile = "";
    public static int[][] beamletSet;
    public static int step, totalObjFuncEval=0;
    public static int iterations, selectionCriterion;
    public static String jobID;
    public static long globalTimer;
    public static String solver="ipopt", name="";
    public static int neighbourhoodSize, globalIter, globalObjFuncEval;
    public static int[] tabuList = new int[360];
    public static int penalizacion; //  Numero de iteraciones que un angulo estara penalisado
    public static int iterationsLimit;// Numero de iteraciones de la tabu search
    public static int noImprovementLimit;// Numero maximo de iteraciones donde no se encuentra un best
    public static int tipoStructura;// Tipo de la estructura de memoria a utilizar
    public static int iteracionesParaAprender;
    public static int tipoDeTabuList;
    public static double[] performanceList = new double[360];
    public static double[][] datosLocalSearch;
    public static double[][] performanceList2 = new double[360][360];
    public static int[][] tabuList2;
    public static int[][] frequencyList;
    public static boolean busquedaRestringida = false;
    public static int anchoCono;
    public static boolean output;
    public static boolean comparacionOpcion;
    
    
        
    public static void main(String[] args) throws IOException, InterruptedException, ExecutionException{
    //agrs[] = jobID, objFunction option (1 logFunction 2 lexRectum 3 lexBladder,    
        
    /************************************************************
    * Ask user which objective function he/she wants to         * 
    * consider. Options are:                                    *
    * 1- LogFunction                                            *
    * 2- Lexico Rectum
    * 3- Lexico Bladder
    *************************************************************/
    //option = menuInit();
    readInputFile(args[1]);
    beamletSet = new int[numAngles][4];
    
    /************************************************************
    * Initialise some of the variables of the algorithm         * 
    *************************************************************/
    //Thread.currentThread().setName("0"); //0 means main program... kind of currentSol
    jobID = args[0]; // Used by NESI. Corresponds to the JOBID when submitted using file *.sl
    
    DDM M; M = new DDM(o.length);
    String beamInfoDir = pathFile + "beamsInfo.txt";
    beamletSet = M.loadNumBixels(beamInfoDir);
    
    switch (selectionCriterion){
        case 1: //random selection
            name = "TS";
            break;
        case 2: //Selection based on singObjFunc value 
            name = "TSPL";
            break;
        case 3: //Next Descent strategy
            break;
    }
  
     /******************Start Iterations***************************
    * Each Iteration leads to a best Optima                    * 
    *************************************************************/        
    
    localOptima = new ArrayList<>();
    String finalBACsFile="";
    switch (option){
        case 1: 
            finalBACsFile = "./Results/logFunction/"+jobID+"singObj_"+name+"_finalBACs_LogFunction_"+numAngles+"Beams_step"+step+".txt";
            break;
        case 2: 
            finalBACsFile = "./Results/lexicoRectum/"+jobID+"singObj_"+name+"_finalBACs_LexicoR_"+numAngles+"Beams_step"+step+".txt";
            break;
        case 3: 
            finalBACsFile = "./Results/lexicoBladder/"+jobID+"singObj_"+name+"_finalBACs_LexicoB_"+numAngles+"Beams_step"+step+".txt";
            break;
        case 4: 
            finalBACsFile = "./Results/inverseFunction/"+jobID+"singObj_"+name+"_finalBACs_inverseFunction_"+numAngles+"Beams_step"+step+".txt";
            break;
        case 6: 
            finalBACsFile = "./Results/weightedSum/"+jobID+"singObj_"+name+"_finalBACs_weightedSum_"+numAngles+"Beams_step"+step+".txt";
            break;
    }
    //neighbourhoodSize = numAngles * 2; //Two neighbours per beam angle; 
    generatedBACs = new ArrayList<>(); //Those BACs that have been generated    
    visitedBACs = new ArrayList<>(); //Those BACs whose neighbourhood has been generated
    /*************************************************************  
    * Create CheckFile. Allows us to shut down the algorithm     *
    * ***********************************************************/
    BufferedWriter bwCheckFile=null;
    File checkFile = new File("./checkFiles/"+jobID+"checkFile.txt");
    if (checkFile.exists()){
        bwCheckFile = new BufferedWriter(new FileWriter(checkFile, true));
        writeLine("Remove this file if you want the process to be stopped\n" , bwCheckFile);
        bwCheckFile.close();
    }else{
        bwCheckFile = new BufferedWriter(new FileWriter(checkFile));
        writeLine("Remove this file if you want the process to be stopped\n" , bwCheckFile);
        bwCheckFile.close();
    }
    globalObjFuncEval = 0;
    globalTimer = System.currentTimeMillis();
    for(globalIter = 0; globalIter<iterations; globalIter++){
        tabuSearch();
        localOptima.get(globalIter).printSol(finalBACsFile);
    }
    System.out.println("Total Objective Function Evaluations: "+globalObjFuncEval+"\n");
    System.out.println("Total Time (secs): "+(System.currentTimeMillis()-globalTimer)/1000+"\n");
}
    
    public static void generateDDM(TreatmentPlan sol, String processName) throws FileNotFoundException, IOException{
        DDM M = new DDM(o.length);
        System.out.println("Generating DDM...");
        for (Organs o1 : o) {
            //Read Intesities of ALL Beam Angles 'j' for organ 'y'
            for (int j = 0; j<sol.beams; j++) {
                sol.selAngles[j].readIntensities(pathFile, o1);
            }
            //Writing DDM for organ 'y'
            o1.totalVoxels=M.writeDDM(sol.beamlets, sol.selAngles, o1, processName + "DDM_", processName);
        }
        System.out.println("");    
    }
    
    public static TreatmentPlan generateInitialBAC(TreatmentPlan currentSol) throws IOException{
        Random r = new Random();
        System.out.println("First BAC (" + jobID + "):");
        ArrayList<Integer> auxAng = new ArrayList<>();
        auxAng.add(r.nextInt(360));
        while (auxAng.get(0)%step != 0){
            auxAng.set(0,auxAng.get(0)-1);
        }
        int auxCounter = 1;
        
        while(auxAng.size() < currentSol.beams){
            boolean auxFlag;
            do{
                
                auxFlag = true;
                auxAng.add(r.nextInt(360));
                while (auxAng.get(auxCounter)%step != 0){
                    auxAng.set(auxCounter,auxAng.get(auxCounter)-1);
                }
                for (int j=0;j<auxCounter;j++){
                    int max = Math.max(auxAng.get(auxCounter), auxAng.get(j));
                    int min = Math.min(auxAng.get(auxCounter), auxAng.get(j));
                    if (max-min < 20){
                        auxFlag = false;
                    }
                    if (max-min >= 175 && max-min <= 185){
                        auxFlag = false;
                    }
                }
                if (!auxFlag){
                    auxAng.remove(auxCounter);
                }
            }while (!auxFlag);
            auxCounter++;
        }
        Collections.sort(auxAng);
        
        for (int i=0;i<currentSol.beams;i++){
            currentSol.selAngles[i]=new Beam(beamletSet[auxAng.get(i)][1],auxAng.get(i), pathFile);
            currentSol.selAngles[i].intensity = (ArrayList<double[]>[])new ArrayList[o.length];
            currentSol.selAngles[i].init=currentSol.beamlets;
            currentSol.beamlets=currentSol.beamlets + currentSol.selAngles[i].beamlets;
            System.out.println(currentSol.selAngles[i].index +"("+currentSol.selAngles[i].init+")"+ " - ");
            
        }
        return (currentSol);
        //System.out.println("(TotalBeamlets: "+currentSol.beamlets + " )");
    }
    
    public static TreatmentPlan  generateInitialBAC(int[] b, TreatmentPlan currentSol) throws IOException{
        System.out.println("First BAC (" + jobID + "):");
        ArrayList<Integer> auxAng = new ArrayList<>();
        for (int i= 0; i< b.length; i++){
            auxAng.add(b[i]);
        }
        //Collections.sort(auxAng);
        for (int i=0;i<currentSol.beams;i++){
            currentSol.selAngles[i]=new Beam(beamletSet[auxAng.get(i)][1],auxAng.get(i),pathFile);
            currentSol.selAngles[i].intensity = (ArrayList<double[]>[])new ArrayList[o.length];
            currentSol.selAngles[i].init=currentSol.beamlets;
            currentSol.beamlets=currentSol.beamlets + currentSol.selAngles[i].beamlets;
            System.out.println(currentSol.selAngles[i].index +"("+currentSol.selAngles[i].init+")"+ " - ");   
        }
        return (currentSol);
        //System.out.println("(TotalBeamlets: "+currentSol.beamlets + " )");
    }
            
  public static TreatmentPlan[] generateNeighbourhood(TreatmentPlan sol, 
            int ns, int[][] bs, String cFile, String performanceListFile, String performanceList2File, int iteracion) throws IOException, InterruptedException, ExecutionException{
        
        int[][] auxAngles = new int[ns][numAngles], angulosReemplazo = new int [ns][2];
        TreatmentPlan[] nhood = new TreatmentPlan [ns];
        int neighbour = 0, angle, anguloQueSacare=0, nuevoAngulo, auxAngle = 0;
        Random  rnd = new Random();
        ArrayList<Integer> validBeams;
        
        System.out.print("GENERATING NEIGHBOURHOOD FOR BAC: " );
      
        for (int j = 0; j<numAngles;j++){
            System.out.print(sol.selAngles[j].index + " - " );
        }
        System.out.println();
        /**********Calculate BAC to be Visited**********
        * We generate (2 * #beams) neighbours of our   *
        * current solution. Each beam angle is modified*
        * in +/- "step" degrees.                       *
        ***********************************************/
        validBeams = new ArrayList<>();
        
        for(int k = 0; k < 360; k=k+step){
            if(isValid(k, sol)){
                validBeams.add(k);
            }
        }
            
        while(neighbour<neighbourhoodSize){

            for (int l=0; l<numAngles;l++){
                auxAngles[neighbour][l]=sol.selAngles[l].index;
            }

            anguloQueSacare = rnd.nextInt(numAngles);

            if(iteracion >= iteracionesParaAprender && tipoStructura == 1){

                nuevoAngulo = seleccionPorRuleta(validBeams);

                for(int i = 0; i < validBeams.size(); i++){

                    if(validBeams.get(i) == nuevoAngulo){

                        validBeams.remove(i);
                    }
                }

                angulosReemplazo[neighbour][0] = auxAngles[neighbour][anguloQueSacare];
                angulosReemplazo[neighbour][1] = nuevoAngulo;        
                auxAngles[neighbour][anguloQueSacare] = nuevoAngulo;    

                neighbour++;

            }else if(iteracion > iteracionesParaAprender && tipoStructura == 2){

                nuevoAngulo = seleccionPorRuletaMatriz(anguloQueSacare,validBeams);

                for(int i = 0; i < validBeams.size(); i++){

                    if(validBeams.get(i) == nuevoAngulo){
                        System.out.println("Se elimino el angulo "+validBeams.get(i)+" de los posibles candidatos");
                        validBeams.remove(i);
                    }
                }

                angulosReemplazo[neighbour][0] = auxAngles[neighbour][anguloQueSacare];
                angulosReemplazo[neighbour][1] = nuevoAngulo;        
                auxAngles[neighbour][anguloQueSacare] = nuevoAngulo;    

                neighbour++;

            }else{

                auxAngle = validBeams.get(rnd.nextInt(validBeams.size()));

                for(int i = 0; i < validBeams.size(); i++){

                    if(validBeams.get(i) == auxAngle){

                        validBeams.remove(i);
                    }
                }

                angulosReemplazo[neighbour][0] = auxAngles[neighbour][anguloQueSacare];
                angulosReemplazo[neighbour][1] = auxAngle;        
                auxAngles[neighbour][anguloQueSacare] = auxAngle;    

                neighbour++;

            }    
       }
        
        updateFrequencyList(angulosReemplazo, ns);
        
        /***********Sort New BACs***********/
        ArrayList<Integer> angSorter;
        for (int i=0; i<auxAngles.length;i++){
            angSorter = new ArrayList<>();
            for (int j=0; j<auxAngles[i].length;j++){
                angSorter.add(auxAngles[i][j]);
            }
            Collections.sort(angSorter);
            for (int j=0; j<auxAngles[i].length;j++){
                auxAngles[i][j]=angSorter.get(j);
            }
        }
        
        /**********Neighbours Evaluation************
        * Evaluate neighbours that must be visited *
        *******************************************/
        if (parallelNeighbourhoodGeneration){
            //getNeighbour[] gNeighbour = new getNeighbour[ns]; //Callables
            Future<TreatmentPlan>[] ft = new Future[ns]; //Callables
            ExecutorService pool = Executors.newFixedThreadPool(ns);
            //create Callables and FutureTasks
            for (int i = 0; i<ns;i++){
                Callable<TreatmentPlan> gNeighbour = new getNeighbour(auxAngles[i], sol, cFile, Integer.toString(i+1), o, angulosReemplazo, i);
                ft[i] = new FutureTask<>(gNeighbour);
                ft[i] = pool.submit(gNeighbour);
            }
            for (int i = 0; i<ft.length;i++){
                nhood[i] = new TreatmentPlan(0,0,numAngles, numOrgans, output); 
                nhood[i].updateSol(ft[i].get());
            }
            pool.shutdown();
        }else{
            for (int i = 0; i<ns;i++){
                /**********Setting Up neighbours************
                *******************************************/
                nhood[i] = new TreatmentPlan(0,0,numAngles, numOrgans, output);
                nhood[i].beams=numAngles;
                for (int j = 0; j<numAngles;j++){
                    angle = (int) auxAngles[i][j];
                    nhood[i].selAngles[j]=new Beam(bs[angle][1],angle,pathFile );
                    nhood[i].selAngles[j].intensity = (ArrayList<double[]>[])new ArrayList[o.length];
                    nhood[i].selAngles[j].init=nhood[i].beamlets;
                    //System.out.println(nhood[i].selAngles[j].index + "(" + nhood[i].selAngles[j].init + ") - " + nhood[i].selAngles[j].beamlets);
                    nhood[i].beamlets= nhood[i].beamlets+nhood[i].selAngles[j].beamlets;
                }
                int auxIndex=generatedBefore(nhood[i].selAngles); //auxIndex<0 means no generated before
                if (auxIndex<0){
                    generateDDM(nhood[i], jobID);
                    //To liberate memory space
                    for (int j=0; j<nhood[i].beams;j++){
                        nhood[i].selAngles[j].intensity=null;
                    }
                    nhood[i].intensity=new double[nhood[i].beamlets];

                    /*******************************************
                    * Update Intensity values. To speed up opt *
                    * process. Only beamlets from common beams *
                    * are updated. Beamlets belonging to the   *
                    * new beam will be set at zero.            *
                    *******************************************/
                    if (false){ // if false, all the optimization will start with the intensity vector set to 1 (very inefficient)
                        System.arraycopy(getNewIntensityVector(sol, nhood[i].selAngles, nhood[i].beamlets), 0, nhood[i].intensity, 0,nhood[i].beamlets);
                    }
                    /*******************************************
                    * Make the optimization process.           *
                    *******************************************/
                    nhood[i].generateReferencePoint(o, option, jobID);
                    nhood[i].getVxDx(jobID, o, Vx, Dx);
                    //Random r=new Random();
                    //nhood[i].singleObjectiveValue=r.nextDouble();
                    generatedBACs.add(nhood[i]);
                    
                    if(tipoStructura==1){
                        updatePerformanceList(sol.singleObjectiveValue, nhood[i].singleObjectiveValue, angulosReemplazo, i);
                    }else if(tipoStructura==2){
                        updatePerformanceList2(sol.singleObjectiveValue, nhood[i].singleObjectiveValue, angulosReemplazo, i);
                    }
                        
                    totalObjFuncEval++;
                }else{
                    nhood[i].updateSol(generatedBACs.get(auxIndex));
                }
                nhood[i].printSol(cFile);
                
                /*Only for NEXT DESCENT Local Search*/
                if (name=="nextDescent" && nhood[i].singleObjectiveValue<sol.singleObjectiveValue){
                    break; //return a neighbourhood with only one solution better than the current solution
                }
                
            }
        }
        System.out.println("*****************************************" );
        System.out.println("*    NEIGHBOURHOOD HAS BEEN GENERATED   *" );
        System.out.println("*****************************************" );
        
        if(tipoStructura==1){
            
            printPerformanceList(performanceListFile, iteracion);
        }else if(tipoStructura==2){
            
            printPerformanceList2(performanceList2File, iteracion);
        }
        
        return nhood;
    }
  
  public static boolean isValid(int angle, TreatmentPlan cSol){
    boolean flag=true; int margenInferior=0; int margenSuperior=360; int imagen;
     
    for (int i = 0; i<numAngles; i++){
        if(busquedaRestringida==true){
            imagen = cSol.selAngles[i].index+180;
            if(imagen>=360){
                imagen=imagen-360;
            }
            if((angle <= imagen + step*anchoCono) && (angle >= imagen - step*anchoCono)){
                flag = false;
                break;
            }
        }
      if(((cSol.selAngles[i].index) != margenInferior) && 
          ((cSol.selAngles[i].index) != margenSuperior-step)){
        if(busquedaRestringida==true){   
            if((angle <= (cSol.selAngles[i].index) + step*anchoCono) && 
                (angle >= (cSol.selAngles[i].index) - step*anchoCono)){
                flag = false;
                break;
            }
        }else{
            if((angle <= (cSol.selAngles[i].index) + step) && 
                (angle >= (cSol.selAngles[i].index) - step)){
                flag = false;
                break;
            }
        }
      }else if(cSol.selAngles[i].index == margenInferior){
          if(busquedaRestringida==true){ 
            if((angle >= margenSuperior-step*anchoCono || angle <= margenInferior+step*anchoCono)){
                flag = false;
                break;
            }
          }else{
              if((angle >= margenSuperior-step || angle <= margenInferior+step)){
                flag = false;
                break;
              }
          }
      }else if(cSol.selAngles[i].index == margenSuperior-step){
          if(busquedaRestringida==true){  
            if((angle >= cSol.selAngles[i].index-step*anchoCono || angle <= margenInferior)){
                flag = false;
                break;
            }
          }else{
            if((angle >= cSol.selAngles[i].index-step || angle <= margenInferior)){
                flag = false;
                break;
            } 
          }
      } 
    }   
    if (angle%step != 0){
      flag =false;
    }
    return flag;
  }
  
  public static TreatmentPlan getBestNeighbour(TreatmentPlan[] nhood, int ns, double bestFitness, TreatmentPlan currentSol){
        double bestValue= 10000;
        int bestIndex= 0;
        
        for (int i = 0; i<ns;i++){
            if (nhood[i] != null){
                
                for(int j=0; j < numAngles; j++){
                    
                    if(currentSol.selAngles[j].index!=nhood[i].selAngles[j].index){
                        
                        if(isTabu(tabuList[nhood[i].selAngles[j].index],currentSol.selAngles[j].index) && tipoDeTabuList!=0){
                            
                            if (nhood[i].singleObjectiveValue<bestFitness && nhood[i].singleObjectiveValue>0){
                                bestValue = nhood[i].singleObjectiveValue;
                                bestIndex = i;
                            }        
                        }else if(nhood[i].singleObjectiveValue<bestValue && nhood[i].singleObjectiveValue>0){
                            bestValue = nhood[i].singleObjectiveValue;
                            bestIndex = i;
                        }    
                    }  
                }
            }
        }
        return nhood[bestIndex];
    }
    
  public static void addNewVisitedBAC(TreatmentPlan sol){
      visitedBACs.add(sol);
    }
  
  
  public static int generatedBefore(Beam[] angles){
      boolean flag=true;
      int index = -1;
      for (int i=0;i<visitedBACs.size();i++){
          flag=true;
          for (int j=0;j<numAngles;j++){
            if (visitedBACs.get(i).selAngles[j].index != angles[j].index){
                flag=false;
                break;
            }
          }
          if (flag){
              index=i;
              break;
          }
      }
      return index;
    }

  public static void savePoint(Beam[] angles, double[] gEUD, int count){
             for (int j = 0; j<numAngles; j++){
                ndp[count][j]=angles[j].index;
             }
             for (int j = numAngles; j<numAngles+numOrgans; j++){
                ndp[count][j]=gEUD[j-numAngles];
             }
  }
  
  
  public static double[] getNewIntensityVector(TreatmentPlan s, Beam[] n, int nBlts) throws IOException{
  
  //Identify what angles from s are not in n. The algorithm should be modified in case 
  //the number of angles of n different from s is more than 1
  
  int indexInS=0, indexInN=0;
  boolean flag;
  
  
  for(int i=0; i<s.selAngles.length;i++ ){
      flag=false;
      for(int j=0; j<n.length;j++ ){
          if (n[j].index == s.selAngles[i].index){
              flag=true;
              break;
          }
      }
      if (!flag){
          indexInS = i;
      }
  }
  
  double[] x = new double[nBlts];
  for(int j=0; j<n.length;j++ ){
      flag=false;
      for(int a=0; a<s.selAngles.length;a++){ //busca si el beam de n esta en s
          if (n[j].index == s.selAngles[a].index){
                if (n[j].beamlets == s.selAngles[a].beamlets){ //just checking
                    System.arraycopy(s.intensity, s.selAngles[a].init, x, n[j].init, s.selAngles[a].beamlets);
                    flag = true;
                    break;
                }else{
                    System.out.println("!!!!!!!!!!!!!!!Error number of beamlets!!!!!!!!!!!!!!!!!");
                }
          }
    }
    if (!flag){
        indexInN = j;
        if (n[indexInN].beamlets == s.selAngles[indexInS].beamlets){
            for(int b=0; b < s.beams;b++){
                int sFirstBeamlet=0, sLastBeamlet=0, nFirstBeamlet=0, nLastBeamlet=0;
                sFirstBeamlet=s.selAngles[b].init;
                sLastBeamlet = s.selAngles[b].init + s.selAngles[b].beamlets;
                if (b != indexInS){
                    for(int a=0; a < n.length;a++){
                        if (n[a].index == s.selAngles[b].index){
                            nFirstBeamlet=n[a].init;
                            nLastBeamlet=n[a].init+n[a].beamlets;
                            break;
                        }
                    }
                }else{
                    nFirstBeamlet=n[indexInN].init;
                    nLastBeamlet=n[indexInN].init+n[indexInN].beamlets;
                }
                System.arraycopy(s.intensity, sFirstBeamlet, x, nFirstBeamlet, s.selAngles[b].beamlets);
            }     
        }else{
            for(int b=0; b < s.beams;b++){
                int sFirstBeamlet=0, sLastBeamlet=0, nFirstBeamlet=0, nLastBeamlet=0;
                sFirstBeamlet=s.selAngles[b].init;
                sLastBeamlet = s.selAngles[b].init + s.selAngles[b].beamlets;
                if (b != indexInS){
                    for(int a=0; a < n.length;a++){
                        if (n[a].index == s.selAngles[b].index){
                            nFirstBeamlet=n[a].init;
                            nLastBeamlet=n[a].init+n[a].beamlets;
                            break;
                        }
                    }
                    System.arraycopy(s.intensity, sFirstBeamlet, x, nFirstBeamlet, s.selAngles[b].beamlets);
                }else{
                    nFirstBeamlet=n[indexInN].init;
                    nLastBeamlet=n[indexInN].init+n[indexInN].beamlets;
                    double[] zeros= new double[Math.min(s.selAngles[b].beamlets, n[indexInN].beamlets)];
                    System.arraycopy(zeros, 0, x, nFirstBeamlet, zeros.length);
                }
                
            }
        }
    }
  }
  
  /*###############################################################*/
  BufferedWriter bwFile = null ;
  File intensityFile = new File("./intensitiesRecord");
  if (intensityFile.exists()) {
    bwFile = new BufferedWriter(new FileWriter(intensityFile, true));
  }else{
    bwFile = new BufferedWriter(new FileWriter(intensityFile));
  }
  String cBAC="", nBAC="";
  for (Beam selAngle : s.selAngles) {
    cBAC = cBAC + selAngle.index + "-";
  }
  for (Beam n1 : n) {
    nBAC = nBAC + n1.index + "-";
  }
  writeLine(cBAC + " /\t/ " + nBAC +"\n",bwFile);
  for(int j=0;j< Math.max(s.beamlets, nBlts);j++){
      if(j<s.beamlets && j< nBlts){
          writeLine(s.intensity[j] + " \t " + x[j],bwFile);
          if (s.intensity[j] != x[j]){
              writeLine("*\n",bwFile);
          }else{
              writeLine("\n",bwFile);
          }
      }else{
          if(j<s.beamlets){
              writeLine(s.intensity[j] + " \t        " + "\n",bwFile);
          }else{
              writeLine("          \t " + x[j],bwFile);
          }
      }
  }
  writeLine(" ###################### \n",bwFile);
  bwFile.close();
  
  /*###############################################################*/
  
  return x;
  }
  
  public static void writeLine(String l, BufferedWriter bw) throws IOException{	
	bw.write(l);
  }
  
  
  public static void tabuSearch() throws InterruptedException, ExecutionException, IOException {
        totalObjFuncEval=0;
        TreatmentPlan bestSol, bestNeighbour, currentSol;
        String bestFile="", currentFile = "", summaryDir="";
        String tabuListFile= "./Results/logFunction/"+jobID+"TabuList_"+globalIter+".txt";
        String tabuList2File= "./Results/logFunction/"+jobID+"TabuList2_"+globalIter+".txt";
        String frequencyListFile= "./Results/logFunction/"+jobID+"FrequencyList_"+globalIter+".txt";
        String performanceListFile= "./Results/logFunction/"+jobID+"PerformanceListFile_"+globalIter+".txt";
        String performanceList2File= "./Results/logFunction/"+jobID+"PerformanceListFile2_"+globalIter+".txt";
        int iter = 0, noImprovement = 0, contReset=0, iteracionEnQueSeSuperoLS=0;
        double minPerformance=1000, maxPerformance=0, bestFitnessTSTiempoLS=0, bestFitnessTSEvaluacionesLS= 0, primerBestMejorQueLS= 0;
        long localTimer = System.currentTimeMillis(), tsTimeAlSuperarLS= 0;
        File checkFile = new File("./checkFiles/"+jobID+"checkFile.txt");
        boolean flagTiempo = false, flagEvaluaciones = false, flagBest = false;
        
        switch (option){
            case 1: 
                bestFile = "./Results/logFunction/"+jobID+"singObj_"+name+"bestSolution_LogFunction_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                currentFile = "./Results/logFunction/"+jobID+"singObj_"+name+"generatedSolution_LogFunction"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                summaryDir = "./Results/logFunction/"+jobID+"singObj_"+name+"summary_LogFunction_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                break;
            case 2: 
                bestFile = "./Results/lexicoRectum/"+jobID+"singObj_"+name+"bestSolution_LexicoR_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                currentFile = "./Results/lexicoRectum/"+jobID+"singObj_"+name+"generatedSolution_LexicoR_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                summaryDir = "./Results/lexicoRectum/"+jobID+"singObj_"+name+"summary_LexicoR_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                break;
            case 3: 
                bestFile = "./Results/lexicoBladder/"+jobID+"singObj_"+name+"bestSolution_LexicoB_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                currentFile = "./Results/lexicoBladder/"+jobID+"singObj_"+name+"generatedSolution_LexicoB_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                summaryDir = "./Results/lexicoBladder/"+jobID+"singObj_"+name+"summary_LexicoB_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                break;
            case 4: 
                bestFile = "./Results/inverseFunction/"+jobID+"singObj_"+name+"bestSolution_inverseFunction_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                currentFile = "./Results/inverseFunction/"+jobID+"singObj_"+name+"generatedSolution_inverseFunction"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                summaryDir = "./Results/inverseFunction/"+jobID+"singObj_"+name+"summary_inverseFunction_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                break;
            case 5: 
                bestFile = "./Results/ConvexLogFunction/"+jobID+"singObj_"+name+"bestSolution_ConvexLogFunction_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                currentFile = "./Results/ConvexLogFunction/"+jobID+"singObj_"+name+"generatedSolution_ConvexLogFunction"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                summaryDir = "./Results/ConvexLogFunction/"+jobID+"singObj_"+name+"summary_ConvexLogFunction_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                break;
            case 6: 
                bestFile = "./Results/weightedSum/"+jobID+"singObj_"+name+"bestSolution_weightedSum_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                currentFile = "./Results/weightedSum/"+jobID+"singObj_"+name+"generatedSolution_weightedSum_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                summaryDir = "./Results/weightedSum/"+jobID+"singObj_"+name+"summary_weightedSum_"+numAngles+"Beams_step"+step+"_"+globalIter+".txt";
                break;
        }
        currentSol = new TreatmentPlan(0,0,numAngles, numOrgans, output);
        currentSol.beams=numAngles;
        for(int i=0;i<numAngles;i++){
                currentSol.selAngles[i] = new Beam(0,0);
                currentSol.selAngles[i].intensity = (ArrayList<double[]>[])new ArrayList[o.length];
        }
        
        try {
            /*************First Beam Angle Configuration*************
            * Set (randomly) the first beam angle configuration     *
            * We select the beam angles and calculate the total     *
            * number of beamlets                                    *
            ********************************************************/
            TreatmentPlan auxCurrentSol = new TreatmentPlan(0,0,numAngles, numOrgans, output);
            if (randomInitial){
                auxCurrentSol.updateSol(generateInitialBAC(currentSol));
            }else{
                auxCurrentSol.updateSol(generateInitialBAC(initialBACs[globalIter], currentSol));
            }
            currentSol.updateSol(auxCurrentSol);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        
        try {
        /**************Local Search**************
        * Create DDM by reading the intensities *
        * and writing the DDM files. Also we    *
        * create a voxelIndex for each Organ    *
        ****************************************/ 
            generateDDM(currentSol, jobID);
        } catch (FileNotFoundException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        //To liberate memory space
        for (int j=0; j<currentSol.beams;j++){
            currentSol.selAngles[j].intensity=null;
        }
        
        /**********Local Search**********/
        
        
        //Initialize LS parameters
   

        /******Setting up currentSol******
        * Create a current solution based*
        * on previously selected beams   *
        *********************************/ 
        currentSol.intensity=new double[currentSol.beamlets];
        for (int i=0; i<currentSol.intensity.length; i++){
            currentSol.intensity[i]=1;
        }

        /*##########################################################
        # Generate the reference point. If we want to add a new    # 
        # objective function as an option, we need to modify       #
        # the method "generateReferencePoint" in the treatmentPlan #
        # class.                                                   # 
        ###########################################################*/ 
        
        System.out.print("Generating Reference point....");
        

        currentSol.generateReferencePoint(o, option, jobID);
        currentSol.getVxDx(jobID, o, Vx, Dx);
        //Random r=new Random();
        //currentSol.singleObjectiveValue=r.nextDouble();
        
        
        
        totalObjFuncEval++;
        generatedBACs.add(currentSol);
        System.out.println("done!");
        /******Update List of Visited BACs******
        * The visited BAC is included into the *
        * visited BAC list and is printed out  *
        * to the corresponding files           *
        ***************************************/ 
        
        addNewVisitedBAC(currentSol); //Add Angles and Reference Point
        bestSol = new TreatmentPlan(0,0,numAngles, numOrgans, output);
        bestSol.updateSol(currentSol);
        currentSol.printSol(currentFile);
        bestSol.printSol(bestFile);
        tabuList = new int[360];
        tabuList2 = new int[360][360]; 
        performanceList = new double[360];
        performanceList2 = new double[360][360];
        frequencyList = new int[360][360];      
        
        if(tipoDeTabuList==1){  
            for(int i=0; i< numAngles && tipoDeTabuList!=0; i++){
                tabuList[currentSol.selAngles[i].index] = penalizacion;
            } 
            printTabuList(tabuListFile);
        }
        if(iteracionesParaAprender==0){
            readPerformances();
        }
   
        while (iter<iterationsLimit){
            
            
            /****Generating neighbourhood list****
            * For each generated neighbour, the *
            * reference point is calculated and *
            * the corresponding BAC is included *
            * in the visited BAC list           *
            *************************************/
            TreatmentPlan[]neighbourhood = new TreatmentPlan[neighbourhoodSize];
            neighbourhood = generateNeighbourhood(currentSol, neighbourhoodSize, beamletSet, currentFile, performanceListFile, performanceList2File, iter);
            
            /****Select Next Current Solution****
            ************************************/
            bestNeighbour = new TreatmentPlan(0,0,numAngles, numOrgans, output);
            bestNeighbour.updateSol(getBestNeighbour(neighbourhood,neighbourhoodSize,bestSol.singleObjectiveValue,currentSol));
            if(tipoDeTabuList==1){
                updateTabuList(bestNeighbour, currentSol,tabuListFile);  
            }else if(tipoDeTabuList==2){
                updateTabuList2(bestNeighbour, currentSol,tabuList2File);
            }  
            currentSol.updateSol(bestNeighbour);
            currentSol.printSol(currentFile);
            if(currentSol.singleObjectiveValue<bestSol.singleObjectiveValue && currentSol.singleObjectiveValue>0){
                noImprovement=0;
                bestSol.updateSol(currentSol);
                System.out.println("New Best Sol Found : " + bestSol.singleObjectiveValue);
                System.out.print("Beam Angle Configuration : ");
                for (int i = 0; i<numAngles;i++){
                    System.out.print(" " + bestSol.selAngles[i].index + " - "  );
                }
                System.out.println();
                System.out.print("gEUD Values : ");
                for (int i = 0; i<numOrgans;i++){
                    System.out.print(" " + bestSol.gEUD[i] + " - "  );
                }
                System.out.println();
                try {
                    bestSol.printSol(bestFile);
                } catch (IOException ex) {
                    ex.printStackTrace();
                }    
            }else{
                noImprovement++;
            }
            if (noImprovement > noImprovementLimit){
                System.out.print("Limit for movements without improvement has been reached. Restarting Current Solution");
                currentSol = new TreatmentPlan(0,0,numAngles, numOrgans, output);
                currentSol.beams=numAngles;
                TreatmentPlan auxCurrentSol = new TreatmentPlan(0,0,numAngles, numOrgans, output);
                auxCurrentSol.updateSol(generateInitialBAC(currentSol));
                currentSol.updateSol(auxCurrentSol);
                generateDDM(currentSol, jobID);
                //To liberate memory space
                for (int j=0; j<currentSol.beams;j++){
                    currentSol.selAngles[j].intensity=null;
                }
                currentSol.intensity=new double[currentSol.beamlets];
                for (int i=0; i<currentSol.intensity.length; i++){
                    currentSol.intensity[i]=1;
                }
                
                
                System.out.print("Generating Reference point....");
                currentSol.generateReferencePoint(o, option, jobID);
                currentSol.getVxDx(jobID, o, Vx, Dx);
                //currentSol.singleObjectiveValue=r.nextDouble();
                
                
                totalObjFuncEval++;
                generatedBACs.add(currentSol);
                System.out.println("done!");
                addNewVisitedBAC(currentSol); //Add Angles and Reference Point
                printResetSol(currentFile, currentSol);
                noImprovement=0;
                contReset++;
            }
            
            //Se compara el tiempo que lleva tabuSearch vs el tiempo total de localSearch
            if(((System.currentTimeMillis()-localTimer)/1000) >= datosLocalSearch[globalIter][0] && flagTiempo == false){
                
                bestFitnessTSTiempoLS = bestSol.singleObjectiveValue;
                flagTiempo = true;
            }
            
            if(bestSol.singleObjectiveValue <= datosLocalSearch[globalIter][2] && flagBest == false){
            
                primerBestMejorQueLS = bestSol.singleObjectiveValue; 
                tsTimeAlSuperarLS = ((System.currentTimeMillis()-localTimer)/1000);
                iteracionEnQueSeSuperoLS = iter;
                flagBest = true;
            }
            
            if (!checkFile.exists()) {
                System.out.println("ALGORITHM HAS BEEN SHUT DOWN. CHECKFILE NOT FOUND");
                break;
            }
            
            iter++;
        }
        System.out.println ("Local Optima Found");
        System.out.println("Final Time: " + (System.currentTimeMillis()-localTimer)/1000);
        System.out.println("Best Obj Function Value: " + bestSol.singleObjectiveValue);
        System.out.print("Beam Angle Configuration : ");
        for (int i = 0; i<numAngles;i++){
            System.out.print(" " + bestSol.selAngles[i].index + " - "  );
        }
        System.out.println();
        System.out.print("gEUD Values : ");
        for (int i = 0; i<numOrgans;i++){
            System.out.print(" " + bestSol.gEUD[i] + " - "  );
        }
        globalObjFuncEval = globalObjFuncEval + totalObjFuncEval;
        
        //Se compara la cantidad de evaluaciones de la funcion objetivo con LocalSearch
        if(totalObjFuncEval >= datosLocalSearch[globalIter][1] && flagEvaluaciones == false){
            
            bestFitnessTSEvaluacionesLS = bestSol.singleObjectiveValue;
            flagEvaluaciones = true;
        }
        
        for(int i = 0; i < performanceList.length; i++){
            
            if(performanceList[i] < minPerformance){
                
                minPerformance = performanceList[i];
            }
            
            if(performanceList[i] > maxPerformance){
                
                maxPerformance = performanceList[i];
            }
        }
        
        printFrequencyList(frequencyListFile);
        
        printSummary(summaryDir, localTimer, Integer.toString(globalIter) , bestSol, contReset, 
            minPerformance, maxPerformance, bestFitnessTSTiempoLS, bestFitnessTSEvaluacionesLS, 
            primerBestMejorQueLS, tsTimeAlSuperarLS, iteracionEnQueSeSuperoLS);
        
        localOptima.add(bestSol);
}

  public static void readInputFile(String input) throws IOException{
        
        String sp="\t";
        String dir = "./"+input;
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine();
        //First read numbero of Organs and angles
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                numOrgans = Integer.parseInt(auxReader[0]);
                numAngles = Integer.parseInt(auxReader[1]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //Go to the next input line
        while(line != null){
            if (!line.contains("%")){
                break;
            }
            line=fileIn.readLine();
        }
        //Info Organs
        o = new Organs[numOrgans];
        
        while(line != null){
            if (!line.contains("%")){
                for (int y=0;y<numOrgans;y++){
                    String[] auxReader = line.split(sp);
                    o[y]=new Organs(
                        auxReader[0], 
                        Integer.parseInt(auxReader[1]), 
                        Double.parseDouble(auxReader[2]),
                        Integer.parseInt(auxReader[3]),
                        Integer.parseInt(auxReader[4]),
                        Integer.parseInt(auxReader[5]),
                        Integer.parseInt(auxReader[6]),
                        Integer.parseInt(auxReader[7]),
                        Integer.parseInt(auxReader[8]),
                        Integer.parseInt(auxReader[9]), 
                        Integer.parseInt(auxReader[10]),
                        Integer.parseInt(auxReader[11]), 
                        Integer.parseInt(auxReader[12]), 
                        Boolean.parseBoolean(auxReader[13]));
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        //get filepath
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                pathFile=auxReader[0];
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get option
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                option=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get selectionCriterion (1 = TS, 2 = nextDescent, 3= gradientBas)
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                selectionCriterion=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get stepsize
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                step=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get iterations
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                iterations=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get numero de iteraciones de tabusearch
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                iterationsLimit=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get penalizacion
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                penalizacion=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get noImprovementLimit
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                noImprovementLimit=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get largo del vecindario
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                neighbourhoodSize=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get tipo de structura de memoria
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                tipoStructura=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get interaciones para aprender
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                iteracionesParaAprender=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get tipo de TabuList
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                tipoDeTabuList=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get busqueda restringida
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                busquedaRestringida=Boolean.parseBoolean(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get ancho del cono
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                anchoCono=Integer.parseInt(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get output de ipopt
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                output=Boolean.parseBoolean(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get solver
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                solver=auxReader[0];
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        //get Parallel neighbour generation flag
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                parallelNeighbourhoodGeneration=Boolean.parseBoolean(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //get Vx Index. % of the volume receiving more than x% of the prescribed dose
        Vx = new int[numOrgans][];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                    String[] auxReader = line.split(sp);
                    Vx[Integer.parseInt(auxReader[0])]=new int[auxReader.length-1];
                    for(int i=1;i<auxReader.length;i++){
                        Vx[Integer.parseInt(auxReader[0])][i-1]=Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        
        //get Dx Index. Minimum dose of the hottest x% of the volume
        Dx = new int[numOrgans][];
        while(line != null){
            if (!line.contains("%")){
                while(!line.contains("%")){
                    String[] auxReader = line.split(sp);
                    Dx[Integer.parseInt(auxReader[0])]=new int[auxReader.length-1];
                    for(int i=1;i<auxReader.length;i++){
                        Dx[Integer.parseInt(auxReader[0])][i-1]=Integer.parseInt(auxReader[i]);
                    }
                    line=fileIn.readLine();
                }
                break;
            }
            line=fileIn.readLine();
        }
        
        //get initial BACs
        initialBACs = new int[iterations][numAngles];
        
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                if (!auxReader[0].equals("none")){
                    randomInitial=false;
                    for (int i=0; i<iterations;i++){
                        for (int j=0; j<numAngles;j++){
                            initialBACs[i][j]= Integer.parseInt(auxReader[j]);
                        }
                        line=fileIn.readLine();
                        auxReader = line.split(sp);
                    }
                }else{
                    randomInitial = true;
                }
                break;
            }
            line=fileIn.readLine();
        }
        
         //get Opcion de comparacion
        while(line != null){
            if (!line.contains("%")){
                String[] auxReader = line.split(sp);
                comparacionOpcion=Boolean.parseBoolean(auxReader[0]);
                line=fileIn.readLine();
                break;
            }
            line=fileIn.readLine();
        }
        
        //get Datos LocalSearch para comparar
        datosLocalSearch = new double[iterations][3];
        
        if(comparacionOpcion==true){
            
            while(line != null){
                if (!line.contains("%")){
                    String[] auxReader = line.split(sp);

                    for (int i=0; i<iterations;i++){
                        for (int j=0; j<3;j++){

                            if(j==2){

                                datosLocalSearch[i][j]= Double.parseDouble(auxReader[j]);
                            }else{

                                datosLocalSearch[i][j]= Integer.parseInt(auxReader[j]);
                            }

                        }
                        line=fileIn.readLine();
                        auxReader = line.split(sp);
                    }

                    break;
                }
                line=fileIn.readLine();
            } 
        }
                       
        fileIn.close();
}
  

 public static void printSummary(String dirFile, long localTimer, String iter, TreatmentPlan bestSol, int totalReset
         , double minPerformance, double maxPerformance, double bestFitnessTSTiempoLS
         , double bestFitnessTSEvaluacionesLS, double primerBestMejorQueLS, long tsTimeAlSuperarLS, int iteracionEnQueSeSuperoLS) throws IOException{
       
        BufferedWriter bwFile = null ;
        File summaryFile = new File(dirFile);
        if (summaryFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(summaryFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(summaryFile));
        }
        
        if(comparacionOpcion==true){
            
            bwFile.write("Iter\tbestSol\tTotalTime\tEvals\tresetCalls\tMinPL\tMaxPL\tLS_Time\tLS_Evals\timprovementF(x)\timprovementTime\timprovementIter\n");
            bwFile.write(iter+"\t"+String.format("%.8f",bestSol.singleObjectiveValue)+"\t"+(System.currentTimeMillis()- localTimer)/1000+
                "\t"+totalObjFuncEval+"\t"+totalReset+"\t"+String.format("%.5f",minPerformance)+"\t"+String.format("%.5f",maxPerformance)+"\t"
                +String.format("%.8f",bestFitnessTSTiempoLS)+"\t"+String.format("%.8f",bestFitnessTSEvaluacionesLS)+"\t"+String.format("%.8f",primerBestMejorQueLS)+"\t"+tsTimeAlSuperarLS+"\t"+iteracionEnQueSeSuperoLS);
        }else{
            
            bwFile.write("Iter\tbestSol\tTotalTime\tEvals\tresetCalls\n");
            bwFile.write(iter+"\t"+String.format("%.8f",bestSol.singleObjectiveValue)+"\t"+(System.currentTimeMillis()- localTimer)/1000+
                "\t"+totalObjFuncEval+"\t"+totalReset);
        }
        //writeLine("############## ITER " +iter+ "##############\n",bwFile);
        //writeLine("Total Obj. Func. Eval.      :\t"+totalObjFuncEval+"\n",bwFile);
        //writeLine("Total Time (secs)           :\t"+(System.currentTimeMillis()- localTimer)/1000+"\n",bwFile);
        //writeLine("Total Objective Function Evaluations (Cummulative):\t"+globalObjFuncEval+"\n",bwFile);
        //writeLine("Total Time (cummulative):\t"+(System.currentTimeMillis()- globalTimer)+"\n",bwFile);
        //writeLine("Counter:\t"+counter[0]+ "\t" +counter[1]+ "\t" +counter[2]+ "\t" 
        //     +counter[3]+ "\t" +counter[4]+ "\t" +counter[5]+ "\n",bwFile);
        bwFile.close();
  } 

 public static class getNeighbour implements Callable<TreatmentPlan> {
    private int[] auxAngles;
    private TreatmentPlan sol;
    private final String cFile;
    private final String processName;
    private Organs[] org;
    private int[][] angulosReemplazo;
    private int indice;
    
    public getNeighbour(int[] auxAng, TreatmentPlan currentSol, String cF, 
            String neighbourNum, Organs[] auxOrg, int[][] angReempl, int j){
        indice=Integer.parseInt(neighbourNum)-1;
        Thread.currentThread().setName(neighbourNum);
        this.processName = jobID + "_" + Thread.currentThread().getName();
        System.out.println("Generating Neighbour Num: " + processName);
        auxAngles = new int[numAngles];
        System.arraycopy(auxAng, 0, auxAngles, 0, auxAngles.length);
        sol = new TreatmentPlan(0,0,numAngles, numOrgans, output);
        sol.updateSol(currentSol);
        cFile = cF;
        this.org = new Organs[auxOrg.length];
        for (int y=0;y<auxOrg.length;y++){
            this.org[y] = new Organs(auxOrg[y].name,auxOrg[y].index, auxOrg[y].weight, auxOrg[y].voxelEnd, auxOrg[y].totalDose, auxOrg[y].actualMinDose,
            auxOrg[y].actualMaxDose, auxOrg[y].doseUB, auxOrg[y].doseLB, auxOrg[y].desiredDose, auxOrg[y].a, auxOrg[y].v, auxOrg[y].totalVoxels, auxOrg[y].isTarget);
        }
        angulosReemplazo=new int[neighbourhoodSize][2];
        for (int i=0;i<neighbourhoodSize;i++){
            angulosReemplazo[i][0]=angReempl[i][0];
            angulosReemplazo[i][1]=angReempl[i][1];
        }
        indice=j;
    }
    
    @Override
    public TreatmentPlan call() throws Exception {
        int angle;    
        TreatmentPlan neighbour = new TreatmentPlan(0,0,numAngles, numOrgans, output);
        neighbour.beams=numAngles;
                
        for (int j = 0; j<numAngles;j++){
            angle = (int) auxAngles[j];
            neighbour.selAngles[j]=new Beam(beamletSet[angle][1],angle,pathFile );
            neighbour.selAngles[j].intensity = (ArrayList<double[]>[])new ArrayList[o.length];
            neighbour.selAngles[j].init=neighbour.beamlets;
            neighbour.beamlets= neighbour.beamlets+neighbour.selAngles[j].beamlets;
        }
        int auxIndex=generatedBefore(neighbour.selAngles); //auxIndex<0 means no generated before
        if (auxIndex<0){
     
            synchronized(o){
                generateDDM(neighbour, this.processName);
                for (int y=0;y<o.length;y++){
                    this.org[y].updateOrgans(o[y]);
                }
            }
            //To liberate memory space
            for (int j=0; j<neighbour.beams;j++){
                neighbour.selAngles[j].intensity=null;
            }
            neighbour.intensity=new double[neighbour.beamlets];

            /*******************************************
            * Update Intensity values. To speed up opt *
            * process. Only beamlets from common beams *
            * are updated. Beamlets belonging to the   *
            * new beam will be set at zero.            *
            *******************************************/
            System.arraycopy(getNewIntensityVector(sol, neighbour.selAngles, neighbour.beamlets), 0, neighbour.intensity, 0,neighbour.beamlets);
            /*******************************************
            * Run the optimization process.           *
            *******************************************/
            neighbour.generateReferencePoint(this.org, option, this.processName);
            neighbour.getVxDx(this.processName, this.org, Vx, Dx);
            //Random r=new Random();
            //neighbour.singleObjectiveValue=r.nextDouble();
            synchronized(performanceList){
             if(tipoStructura==1){
                updatePerformanceList(sol.singleObjectiveValue, 
                        neighbour.singleObjectiveValue, angulosReemplazo, indice);
             }
            }
            
            synchronized(visitedBACs){
                //addNewVisitedBAC(neighbour); //Add Angles and Reference Point
                generatedBACs.add(neighbour);
            }
            totalObjFuncEval++;
        }else{
            neighbour.updateSol(generatedBACs.get(auxIndex));
        }
        neighbour.printSol(cFile);
            return neighbour;
    }

}
 
public static void updateTabuList(TreatmentPlan bSol, TreatmentPlan cSol, String tabuListFile) throws IOException{
    
    for(int i=0; i<360; i=i+step){
        
        if(tabuList[i]>0){
            
            tabuList[i] = tabuList[i]-1;
        }
    }
  
    for(int i=0; i<numAngles; i++){
        boolean flag=true;
        for(int j=0; j<numAngles; j++){
            if(bSol.selAngles[i].index==cSol.selAngles[j].index){
                flag=false;
                break;
            }
        }
        if (flag){
            tabuList[bSol.selAngles[i].index]=penalizacion;
        }
    }

    printTabuList(tabuListFile);

} 

public static void printTabuList(String tabuListFile) throws IOException{
      BufferedWriter bwFile = null ;
      File listFile = new File(tabuListFile);
      
      if (listFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(tabuListFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(tabuListFile));
        }
      
      for(int i=0;i<360;i=i+step){
          
            if(tabuList[i]>0){
                
                 bwFile.write("Angulo " +i+ " = " +tabuList[i]+"\t");
            }        
      }
      
      bwFile.write("\n");
      
      bwFile.close();
}

public static void printResetSol(String dirFile , TreatmentPlan sol) throws IOException{
    
        BufferedWriter bwFile = null ;
        File solFile = new File(dirFile);
        
        if (solFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(solFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(solFile));
        }

        bwFile.write("Se reseteo currentSol \n");
        
        for (int i=0; i< sol.beams;i++){  
            bwFile.write(sol.selAngles[i].index +"\t");
        }
        for (int j=0; j< sol.gEUD.length;j++){
            bwFile.write(sol.gEUD[j] + "\t");
        }
        bwFile.write(sol.singleObjectiveValue + "\t");
        
        for (int j=0; j< sol.weights.length;j++){
            bwFile.write(sol.weights[j] + "\t");
        }
        
        for (double[] Vx1 : sol.Vx) {
            if (Vx1 != null) {
                for (int j = 0; j < Vx1.length; j++) {
                    bwFile.write(Vx1[j] + "\t");
                }
            }
        }
        for (double[] Dx1 : sol.Dx) {
            if (Dx1 != null) {
                for (int j = 0; j < Dx1.length; j++) {
                    bwFile.write(Dx1[j] + "\t");
                }
            }
        }
        
        bwFile.write("x:\t");
        
        for (int j=0; j< sol.intensity.length;j++){
            bwFile.write(sol.intensity[j] + "\t");
        }
        bwFile.write("\n");
        bwFile.close();   
}

public static void updatePerformanceList(double fitnessS, double fitnessE, int[][] angulosReemplazo, int i) throws IOException{
    
    performanceList[angulosReemplazo[i][0]] = performanceList[angulosReemplazo[i][0]] + ((fitnessE - fitnessS)/fitnessS);
    performanceList[angulosReemplazo[i][1]] = performanceList[angulosReemplazo[i][1]] + ((fitnessS - fitnessE)/fitnessS);     
}

public static void updatePerformanceList2(double fitnessS, double fitnessE, int[][] angulosReemplazo, int i) throws IOException{
    
    performanceList2[angulosReemplazo[i][1]][angulosReemplazo[i][0]] = performanceList2[angulosReemplazo[i][1]][angulosReemplazo[i][0]] + ((fitnessS - fitnessE)/fitnessS);  
}

public static void printPerformanceList(String performanceListFile, int iteracion) throws IOException{
    
    BufferedWriter bwFile = null ;
    File listFile = new File(performanceListFile);
      
        if (listFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(performanceListFile, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(performanceListFile));
        }
        
        bwFile.write("\nIteracion n° "+(iteracion+1)+"\n\n");
            
        for(int i = 0; i < performanceList.length; i++){
            
            if(performanceList[i]!=0){
                
                bwFile.write(i+"\t"+performanceList[i]+"\n");
            }
        }
        
      bwFile.close();
    
}

public static void printPerformanceList2(String performanceList2File, int iteracion) throws IOException{
    
    BufferedWriter bwFile = null ;
    File listFile = new File(performanceList2File);
      
        if (listFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(performanceList2File, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(performanceList2File));
        }
        
        bwFile.write("\nIteracion n° "+(iteracion+1)+"\n\n");
            
        for(int i = 0; i < 360; i=i+step){
            for(int j = 0; j < 360; j=j+step){
                bwFile.write(performanceList2[i][j]+"\t");
            }
            bwFile.write("\n");
        }
        
      bwFile.close();
    
}

public static int seleccionPorRuleta(ArrayList<Integer> angulos){ // Funcion para escoger el nuevo angulo usando la ruleta
    double[] probabilidad = new double [angulos.size()];
    int nuevoAngulo=0;
    double total=0, random, acum = 0;  
    Random  rnd = new Random();
    
    random = rnd.nextDouble();
    
    double menorPerformance = 1000;
   
    System.out.println("\nComenzo Busqueda guiada\n");
    
    for(int i=0;i<angulos.size();i++){
        if(performanceList[angulos.get(i)] < menorPerformance){
            menorPerformance = performanceList[angulos.get(i)];  
        }
    }
    total=0;
    for(int i=0;i<angulos.size();i++){
   
        probabilidad[i]=performanceList[angulos.get(i)]+Math.abs(menorPerformance); 
           
        total = total + probabilidad[i];
    }
    
    acum=0;
    for(int i=0;i<angulos.size();i++){
        probabilidad[i] = acum + (probabilidad[i]/total);
        acum = probabilidad[i];
    }
    
    /*ESTO ES INEFICIENTE. SE PUEDE HACER ANTES USANDO EL ACUMUADOR*/
    for(int i=0;i<angulos.size();i++){
        if(random<probabilidad[i]){
            nuevoAngulo=angulos.get(i);
            break;
        }
    }
 
    return nuevoAngulo;
}

public static int seleccionPorRuletaMatriz(int anguloQueSacare, ArrayList<Integer> angulos){
   
    double[] probabilidad = new double [angulos.size()];
    int nuevoAngulo=0, randomInt;
    double total=0, randomDouble, acum = 0, menorPerformance = 1000;  
    Random  rnd = new Random();
    randomDouble = rnd.nextDouble(); 
   
    System.out.println("\nComenzo Busqueda guiada con matriz\n");
    
    for(int i=0;i<angulos.size();i++){
        if(performanceList2[angulos.get(i)][anguloQueSacare] < menorPerformance){
            menorPerformance = performanceList2[angulos.get(i)][anguloQueSacare];  
        }
    }
    total=0;
    
    for(int i=0;i<angulos.size();i++){
        probabilidad[i]=performanceList2[angulos.get(i)][anguloQueSacare]+Math.abs(menorPerformance); 
        total = total + probabilidad[i];
    }
    
    if(total==0){
        randomInt = rnd.nextInt(angulos.size());
        return angulos.get(randomInt);
    }
    
    acum=0;
    for(int i=0;i<angulos.size();i++){
        probabilidad[i] = acum + (probabilidad[i]/total);
        acum = probabilidad[i];
    }
    
    for(int i=0;i<angulos.size();i++){
        System.out.println("Angulo "+angulos.get(i)+" probabilidad "+probabilidad[i]);
        if(randomDouble<probabilidad[i]){
            nuevoAngulo=angulos.get(i);
            break;
        }
    }
 
    return nuevoAngulo;
}

public static void updateTabuList2(TreatmentPlan bSol, TreatmentPlan cSol, String tabuList2File) throws IOException{
    
    int anguloEntra = 0, anguloSale = 0;
    boolean[][] flag = new boolean[2][neighbourhoodSize];
    
    for(int i=0; i<360; i=i+step){
        for(int j=0; j<360; j=j+step){
            if(tabuList2[i][j]>0){
                tabuList2[i][j] = tabuList2[i][j]-1;
            }
        }
    }
  
    for(int i=0; i<numAngles; i++){
        for(int j=0; j<numAngles; j++){
            if(bSol.selAngles[i].index==cSol.selAngles[j].index){
                flag[0][j]=true;
                flag[1][i]=true;
                break;
            }
        }
    }
    
    for(int i=0; i<numAngles; i++){
            
        if(flag[0][i]==false){
            
            anguloSale=cSol.selAngles[i].index;
        }
        
        if(flag[1][i]==false){
            
            anguloEntra=bSol.selAngles[i].index;
        }
    }
    
    tabuList2[anguloEntra][anguloSale] = penalizacion; 

    printTabuList2(tabuList2File);

} 

public static void printTabuList2(String tabuList2File) throws IOException{
        BufferedWriter bwFile = null ;
        File listFile = new File(tabuList2File);
      
        if (listFile.exists()) {
            bwFile = new BufferedWriter(new FileWriter(tabuList2File, true));
        }else{
            bwFile = new BufferedWriter(new FileWriter(tabuList2File));
        }
      
        for(int i=0;i<360;i=i+step){
          
            for(int j=0;j<360;j=j+step){
              
                if(tabuList2[i][j]>0){
                
                    bwFile.write("Angulos " +i+ " - "+j+" = "+tabuList2[i][j]+"\t");
                }     
            }             
        }
      
      bwFile.write("\n");
      
      bwFile.close();
}

public static boolean isTabu(int anguloEntra, int anguloSale){
   
    boolean flag=false;
    
    if(tipoDeTabuList==1){
        
        if(tabuList[anguloEntra]>0){
            flag=true;
        }     
    }else if(tipoDeTabuList==2){
        
        if(tabuList2[anguloEntra][anguloSale]>0){
            flag=true;
        }
    }
    
    return flag;
}

public static void printFrequencyList(String printFrequencyListFile) throws IOException{
    BufferedWriter bwFile = null ;
    File listFile = new File(printFrequencyListFile);

    if (listFile.exists()) {
        bwFile = new BufferedWriter(new FileWriter(printFrequencyListFile, true));
    }else{
        bwFile = new BufferedWriter(new FileWriter(printFrequencyListFile));
    }

    for(int i=0;i<360;i=i+step){

        for(int j=0;j<360;j=j+step){

            bwFile.write(frequencyList[i][j]+"\t");     
        } 

        bwFile.write("\n");
    }
      
    bwFile.close();
}

public static void updateFrequencyList(int[][] angulosReemplazo, int ns){
    
    for(int i=0; i<ns; i++){
        frequencyList[angulosReemplazo[i][1]][angulosReemplazo[i][0]]++;
    }  
}

public static void readPerformances() throws FileNotFoundException, IOException{
    
    String dir = "./Performances.txt";
    File f = new File(dir);
    BufferedReader fileIn = new BufferedReader(new FileReader(f));
    String line = "";
    line=fileIn.readLine();
    String sp="\t";
    
    while(line != null){
            
        String[] auxReader = line.split(sp);

        for (int i=0; i<360;i=i+step){
                performanceList[i]= Double.parseDouble(auxReader[0]); 

                line=fileIn.readLine();
                auxReader = line.split(sp);
        }
        break;
    }  
    
    fileIn.close();
}

}




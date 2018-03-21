package IMRT_Base;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.*;
import java.util.*;
        
/**
 *
 * @author gcab623
 */

public class DDM {
    public ArrayList<Long[]> voxelIndex[];
    
    public  double[][][] ddm0;
    public  boolean[] tumorVoxelFlag;
    public int[] OrganVoxels;
    public int  totalVoxels;
    public int  totalAngles;
    public int  totalFakeAngles;
    public int  totalBeamlets;
    public int  totalFakeBeamlets;
    public int  nzElements;
    
    public DDM(int bl, int o) {
        this.totalAngles = o;
        this.totalBeamlets=bl;
        this.totalFakeBeamlets=0;
        this.totalFakeAngles = 36;
        //this.ddm0 = new double[10000][o][bl]; //Voxels-Organ-Beamlets
        //this.ddm0 = new double[3][bl]; //Voxels-Organ-Beamlets
        this.nzElements=0;
        
    }
    public DDM(int numOrgans) {
        voxelIndex = (ArrayList<Long[]>[])new ArrayList[numOrgans];
        
    }
    
    public int writeDDM(int numBlts, Beam[] beams , Organs org, String fileName, String processName) throws IOException{
     //numbero of Voxels of organ 'y'
     ArrayList<double[]> DDM = new ArrayList <>();
     voxelIndex[org.index] = new ArrayList<>();
     int vIndex, col;
     double voxelNum, intensity;
     for (int j=0;j<beams.length;j++){
         for (int i=0;i<beams[j].intensity[org.index].size();i++){
             voxelNum = beams[j].intensity[org.index].get(i)[0];
             col = (int)beams[j].intensity[org.index].get(i)[1]-1+ beams[j].init; //beamlet Index;
             intensity = beams[j].intensity[org.index].get(i)[2];
             //vIndex=voxelIndex[org.index].indexOf(voxelNum);
             vIndex= (int)getIndexOf(voxelNum, org.index);
             if (vIndex==-1){ //voxel is not in the voxelIndex (new voxel)
                 Long[] auxVoxInd = new Long [2];
                 vIndex = voxelIndex[org.index].size();
                 auxVoxInd[0] = (long)voxelNum;
                 auxVoxInd[1] = (long)vIndex;
                 voxelIndex[org.index].add(auxVoxInd);
                 Collections.sort(voxelIndex[org.index], new Comparator<Long[]>() {
                    public int compare(Long[] array1, Long[] array2) {
                        return Long.compare(array1[0], array2[0]);
                    }
                });
                 double[] auxDDM = new double[numBlts+1];
                 DDM.add(auxDDM);
             }
             if (DDM.get(vIndex)[col]==0){
                DDM.get(vIndex)[col] = intensity;
                DDM.get(vIndex)[numBlts] += intensity; 
             }else{
                System.err.println("Error 4. Se sobre-escribio DDM, beamlet" + col + ", voxel " + vIndex);
             }
        }
     }
     //o[org.index].totalVoxels=voxelIndex[org.index].size();
     //int numVoxels = org.totalVoxels;
     int numVoxels = voxelIndex[org.index].size();
     
     BufferedWriter bwDDM = null ;
     fileName = fileName + org.name + ".dat";
     BufferedWriter bwVoxelIndex = null ;
     String vIndexName = "./"+processName+"voxelIndex" + org.name + ".txt";
     
     
     File fileAMPL = new File(fileName);
     if (!fileAMPL.exists()) {
        bwDDM = new BufferedWriter(new FileWriter(fileName));
     }else{
        fileAMPL.delete();
        bwDDM = new BufferedWriter(new FileWriter(fileName));
     }

     File fileVoxelIndex = new File(vIndexName);
     if (!fileVoxelIndex.exists()) {
        bwVoxelIndex = new BufferedWriter(new FileWriter(vIndexName));
     }else{
        fileVoxelIndex.delete();
        bwVoxelIndex = new BufferedWriter(new FileWriter(vIndexName));
     }
     
     writeLine("param ddm" + org.name + ":\t", bwDDM);
     for (int auxIndex=1;auxIndex<numBlts;auxIndex++){
        writeLine(auxIndex + "\t", bwDDM);
     }
     writeLine(numBlts + " :=\n", bwDDM);
     
     ArrayList<Long> deletedVoxels = new ArrayList <>();
     int deletedCounter = 0;
     int row = 1;
     for (int v=0;v<numVoxels;v++){
         if (DDM.get(v)[numBlts] > 0 || org.isTarget){
             writeLine(row + "\t" + voxelIndex[org.index].get(v)[0] + "\t" + voxelIndex[org.index].get(v)[1] + "\n", bwVoxelIndex);
             writeLine(row + "\t", bwDDM);
             for (int j=0;j<numBlts;j++){
                writeLine(DDM.get(v)[j] + "\t", bwDDM);
             }
             writeLine("\n", bwDDM);
             row++;
         }else{
             deletedVoxels.add((long)v);
             deletedCounter++;
         }
    }
    writeLine(";\n", bwDDM);
    bwDDM.close();
    bwVoxelIndex.close();
    System.out.println(deletedCounter + " voxels have been removed from " + org.name);
    for (int v=0;v<deletedVoxels.size();v++){
        for (int i=0;i<voxelIndex[org.index].size();i++){
            if (voxelIndex[org.index].get(i)[1].equals(deletedVoxels.get(v))){
                voxelIndex[org.index].remove(i);
                break;
            }
        }
    }
    //o[org.index].totalVoxels=voxelIndex[org.index].size();
    //System.out.println("Voxels : " + voxelIndex[org.index].size() + " -- Removed: " + deletedVoxels.size());    
   return voxelIndex[org.index].size();  
 }
    
 public long getIndexOf(double origVoxel, int orgIndex) {
        int posLeft, posRight, pivot=0;
        posLeft=0;posRight = voxelIndex[orgIndex].size();
                
        if (posRight != 0){
            if (voxelIndex[orgIndex].get(0)[0] == origVoxel){
                return voxelIndex[orgIndex].get(0)[1];
            }else{
                while (posLeft < posRight) {
                    pivot = posLeft + (int) (posRight- posLeft)/2 ;
                    if (voxelIndex[orgIndex].get(pivot)[0] < origVoxel){ //move posLeft to pivot
                        posLeft = pivot;
                    }else if (voxelIndex[orgIndex].get(pivot)[0] > origVoxel){ //move posLeft to pivot
                        posRight = pivot;
                    }else if (voxelIndex[orgIndex].get(pivot)[0] == origVoxel){
                        break;
                    }
                    if (posRight-posLeft == 1){
                        posRight=-1;
                        posLeft = 0;
                    }
                }
                if (posLeft > posRight){
                    return -1;
                }
                return voxelIndex[orgIndex].get(pivot)[1];
            }
        }else{
            return -1;
        }
  }
  
 public int[][] loadNumBixels(String beamsInfoDir) throws IOException{
        int[][] beamlets = new int[360][4];
        //int[][] beamlets = new int[numBeams][4];
        int auxBeamletIndex=0;
        
        File beamInfo = new File(beamsInfoDir);
        if (! beamInfo.exists()) {
            System.err.println("Couldn't find 'beamsInfo.txt' file");
        }else{
            String line ="";
            String[] auxReader=null;   
            File f= new File(beamsInfoDir);
            BufferedReader fileIn= new BufferedReader(new FileReader(f));
            for (int i=0;i<360;i++){
                line=fileIn.readLine();
                auxReader = line.split("\t");
                beamlets[i][0]=(int) Double.parseDouble(auxReader[0]); //beamIndex
                beamlets[i][1]=(int) Double.parseDouble(auxReader[1]); //numBeamlets
                beamlets[i][2]= auxBeamletIndex;                       //firstBeamletIndex
                beamlets[i][3]=(int) Double.parseDouble(auxReader[2]) - 1; //lastBeamletIndex
                auxBeamletIndex = auxBeamletIndex + (int) Double.parseDouble(auxReader[1]);
            }
            fileIn.close();
        }
        return (beamlets);
 }
    
    /*   
    
    public void generateDDM_big(Beam angles[], Organs org[], String pathFile[]) throws FileNotFoundException, IOException{
        //This algorithm does the same than generateDDM but uses 
        //some different (smaller) data structures. It takes more time
        //to produce the same files but avoid to run out of memory
        
        
        //this.ddm0 = new double[10000][1][this.totalBeamlets]; 
        String sp = "\t";
        String[] auxReader;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
        String dir;
        int index=0;
        int i = 0;                                                                                                                                                                                                                           
        //initialize Number of Voxels
        this.OrganVoxels = new int[org.length];
        
        //Set init beams per angle
        this.totalBeamlets = 0;
        for (int z=0;z<angles.length;z++){
            angles[z].init=i;
            this.totalBeamlets=this.totalBeamlets+angles[z].beamlets;
            //System.out.println(angles[z].index+"---"+angles[z].beamlets);
            i=i+angles[z].beamlets;
        }
        
        //determine number of fake beamlets
        
        int[][]fakeAngles = new int[this.totalFakeAngles][2]; // fakeAngle - num of fakeBeamlets per fakeAngle
        int step = 360/this.totalFakeAngles;
        
        for (i=0;i<step;i=i++){
            
            
            
            
            
        }
       
        i=0;
        this.totalAngles = angles.length;
        auxReader = new String[this.totalBeamlets];
        
    	for (int y=0;y<org.length;y++){
            BufferedWriter bwAMPL = null ;
            String dirAMPL = "DDM_AMPL_"+y+".dat";
            File fileAMPL = new File(dirAMPL);
            if (!fileAMPL.exists()) {
                bwAMPL = new BufferedWriter(new FileWriter(dirAMPL));
            }else{
                fileAMPL.delete();
                bwAMPL = new BufferedWriter(new FileWriter(dirAMPL));
            }
            this.ddm0=new double[10000][1][this.totalBeamlets];
            double voxelIndex[]=new double[30000];
            int[][] searchIndex= new int[30000][2]; // Helps in getting a (more) efficient search
            i=0;
            for (int z=0;z<angles.length;z++){
                dir = pathFile[1]+org[y].name+"_"+angles[z].index+".txt";
                File f = new File(dir);
                BufferedReader fileIn = new BufferedReader(new FileReader(f));
                String line = "";
                line=fileIn.readLine();
                auxReader = line.split(sp);
                int ang=0;
                line=fileIn.readLine(); // All values are separated by \t
                line=fileIn.readLine(); // All values are separated by \t
                auxReader = line.split(sp);
                //Create the first register
                //calculate beamlet
                ang=angles[z].init+(int)Double.parseDouble(auxReader[1])-1;
                line=fileIn.readLine(); // All values are separated by \t
                auxReader = line.split(sp);
                voxelIndex[i]=Double.parseDouble(auxReader[0]);
                this.ddm0[i][1][ang]=Double.parseDouble(auxReader[2]);
                this.nzElements++;
                i++;
                line=fileIn.readLine(); // All values are separated by \t
                searchIndex[ang][1]=i-1;
                while (line!=null){ // end of the structure
                    auxReader = line.split(sp);
                    boolean isNew=true;
                    //UpdateCurrentBeamlet
                    if(ang!=angles[z].init+(int)Double.parseDouble(auxReader[1])-1){
                        searchIndex[ang][1]=i-1;
                        if(searchIndex[ang][1]<searchIndex[ang][0]){
                            searchIndex[ang][1]=searchIndex[ang][0];
                        }
                        dir = pathFile[1]+org[y].name+"_"+angles[z].index+".txt";
                        ang=angles[z].init+(int)Double.parseDouble(auxReader[1])-1;
                        searchIndex[ang][0]=i;
                        searchIndex[ang][1]=i;
                    }
                    int j =0;
                    while (j<=ang){
                        index=Arrays.binarySearch(voxelIndex, searchIndex[j][0], searchIndex[j][1], Double.parseDouble(auxReader[0]));
                        if (index>=0){
                            isNew=false;
                            break;
                        }
                        j++;
                    }
                    if (isNew){//New voxel
                        voxelIndex[i]=Double.parseDouble(auxReader[0]);
                        this.ddm0[i][1][ang]=Double.parseDouble(auxReader[2]);
                        this.nzElements++;
                        i++;
                        searchIndex[ang][1]=i;
                    }else{
                        this.ddm0[index][0][ang]=Double.parseDouble(auxReader[2]);
                        this.nzElements++;
                    }
                    line=fileIn.readLine(); // All values are separated by \t
                    if (line!=null){
                        if (line.split(sp).length <3){
                            line = null;
                        }
                    }
                    //Arrays.sort(voxelIndex, 0, i-1);
                }
                fileIn.close();
                
                
            }
            this.OrganVoxels[y]=i;
            this.totalVoxels = this.totalVoxels+ i;
            writeDDM_AMPL(org, y);
            writeDDM_Plain(org, y);
        }
        
   }
   
   public void generateDDM(Beam angles[], Organs org[], String pathFiles[]) throws FileNotFoundException, IOException{
        double voxelIndex[][]=new double[org.length][10000];
        String sp = "\t";
        String[] auxReader;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
        String dir;
        int index=0;
        int i = 0;                                                                                                                                                                                                                           
        //initialize Number of Voxels
        this.OrganVoxels = new int[org.length];
        //Set init beams per angle
        this.totalBeamlets = 0;
        for (int z=0;z<angles.length;z++){
            angles[z].init=i;
            this.totalBeamlets=this.totalBeamlets+angles[z].beamlets;
            //System.out.println(angles[z].index+"---"+angles[z].beamlets);
            i=i+angles[z].beamlets;
        }
        
        i=0;
        this.totalAngles = angles.length;
        this.ddm0=new double[10000][org.length][this.totalBeamlets];
        
    	for (int y=0;y<org.length;y++){
            int[][] searchIndex= new int[1200][2]; // Helps in getting a (more) efficient search
            i=0;
            for (int z=0;z<angles.length;z++){
                //dir = "/media/ubuntuFiles/Data/Input/CERRprostate/"+org[y].name+"_"+angles[z].index+".txt";
                dir = "../CERRprostate/"+org[y].name+"_"+angles[z].index+".txt";
                //dir = "/home/guille/Data/Input/CERRprostate/"+org[y].name+"_"+angles[z].index+".txt";
                File f = new File(dir);
                BufferedReader fileIn = new BufferedReader(new FileReader(f));
                String line = "";
                line=fileIn.readLine();
                auxReader = line.split(sp);
                int ang=0;
                line=fileIn.readLine(); // All values are separated by \t
                line=fileIn.readLine(); // All values are separated by \t
                auxReader = line.split(sp);
                //Create the first register
                //calculate beamlet
                ang=angles[z].init+(int)Double.parseDouble(auxReader[1])-1;
                line=fileIn.readLine(); // All values are separated by \t
                auxReader = line.split(sp);
                voxelIndex[y][i]=Double.parseDouble(auxReader[0]);
                this.ddm0[i][y][ang]=Double.parseDouble(auxReader[2]);
                this.nzElements++;
                i++;
                line=fileIn.readLine(); // All values are separated by \t
                searchIndex[ang][1]=i-1;
                while (line!=null){ // end of the structure
                    auxReader = line.split(sp);
                    boolean isNew=true;
                    //UpdateCurrentBeamlet
                    if(ang!=angles[z].init+(int)Double.parseDouble(auxReader[1])-1){
                        searchIndex[ang][1]=i-1;
                        if(searchIndex[ang][1]<searchIndex[ang][0]){
                            searchIndex[ang][1]=searchIndex[ang][0];
                        }dir = "../CERRprostate/"+org[y].name+"_"+angles[z].index+".txt";
                        ang=angles[z].init+(int)Double.parseDouble(auxReader[1])-1;
                        searchIndex[ang][0]=i;
                        searchIndex[ang][1]=i;
                    }
                    int j =0;
                    while (j<=ang){
                        index=Arrays.binarySearch(voxelIndex[y], searchIndex[j][0], searchIndex[j][1], Double.parseDouble(auxReader[0]));
                        if (index>=0){
                            isNew=false;
                            break;
                        }
                        j++;
                    }
                    if (isNew){//New voxel
                        voxelIndex[y][i]=Double.parseDouble(auxReader[0]);
                        this.ddm0[i][y][ang]=Double.parseDouble(auxReader[2]);
                        this.nzElements++;
                        i++;
                        searchIndex[ang][1]=i;
                    }else{
                        this.ddm0[index][y][ang]=Double.parseDouble(auxReader[2]);
                        this.nzElements++;
                    }
                    line=fileIn.readLine(); // All values are separated by \t
                    if (line!=null){
                        if (line.split(sp).length <3){
                            line = null;
                        }
                    }
                    //Arrays.sort(voxelIndex, 0, i-1);
                }
                fileIn.close();
                
                
            }
            this.OrganVoxels[y]=i;
            
            this.totalVoxels = this.totalVoxels+ i;
        }
        this.tumorVoxelFlag= new boolean[this.OrganVoxels[0]];
        BufferedWriter bwAMPL = null ;
        String dirAMPL = "DDM_AMPL.dat";
        File fileAMPL = new File(dirAMPL);
        if (!fileAMPL.exists()) {
            bwAMPL = new BufferedWriter(new FileWriter(dirAMPL));
        }else{
            fileAMPL.delete();
            bwAMPL = new BufferedWriter(new FileWriter(dirAMPL));
        }
        
        for (int y=0;y<org.length;y++){
            writeLine("param ddm" + org[y].name+ ":\t", bwAMPL);
            for (int auxIndex=1;auxIndex<this.totalBeamlets;auxIndex++){
                writeLine(auxIndex + "\t", bwAMPL);
            }
            writeLine(this.totalBeamlets + " :=\n", bwAMPL);
            
            int auxVoxels=1;
            for (i=0; i<this.OrganVoxels[y];i++){
                double auxSum=0;
                if (y==0){ //we only evaluate voxels from tumor
                    for (int z=0;z<this.totalBeamlets;z++){
                        auxSum = auxSum + this.ddm0[i][y][z];
                    }
                }
                //only "relevant" voxels from the tumor and all voxels from OARs
                if (auxSum>0.01 || y>0){
                    if (y==0)
                        this.tumorVoxelFlag[i]=true;
                    writeLine(auxVoxels + "\t", bwAMPL);
                    for (int z=0;z<this.totalBeamlets;z++){
                        writeLine(this.ddm0[i][y][z]+ "\t", bwAMPL);                
                    }
                    auxVoxels++;
                    writeLine("\n", bwAMPL);
                }else{
                    this.tumorVoxelFlag[i]=false;
                }
            }
            this.OrganVoxels[y]=auxVoxels-1;
            writeLine(";\n", bwAMPL);
        }
        //bwBest.close();
        bwAMPL.close();
   } 
    
    
    
    
    
   public void writeDDM_AMPL(Organs[] organ, int orgIndex) throws IOException{
        this.tumorVoxelFlag= new boolean[this.OrganVoxels[0]];
        BufferedWriter bwAMPL = null ;
        String dirAMPL = "DDM_AMPL_"+orgIndex+".dat";
        File fileAMPL = new File(dirAMPL);
        if (!fileAMPL.exists()) {
            bwAMPL = new BufferedWriter(new FileWriter(dirAMPL));
        }else{
            fileAMPL.delete();
            bwAMPL = new BufferedWriter(new FileWriter(dirAMPL));
        }
        writeLine("param ddm" + organ[orgIndex].name+ ":\t", bwAMPL);
        for (int auxIndex=1;auxIndex<this.totalBeamlets;auxIndex++){
            writeLine(auxIndex + "\t", bwAMPL);
        }
        writeLine(this.totalBeamlets + " :=\n", bwAMPL);
            
        int auxVoxels=1;
        for (int i=0; i<this.OrganVoxels[orgIndex];i++){
            double auxSum=0;
            if (orgIndex==0){ //we only evaluate voxels from tumor
                for (int z=0;z<this.totalBeamlets;z++){
                    auxSum = auxSum + this.ddm0[i][1][z];
                }
            }
            //only "relevant" voxels from the tumor and all voxels from OARs
            if (auxSum>0.01 || orgIndex>0){
                if (orgIndex==0){
                    this.tumorVoxelFlag[i]=true;
                    writeLine(auxVoxels + "\t", bwAMPL);
                    for (int z=0;z<this.totalBeamlets;z++){
                        writeLine(this.ddm0[i][1][z]+ "\t", bwAMPL);                
                    }
                    auxVoxels++;
                    writeLine("\n", bwAMPL);
                }else{
                    this.tumorVoxelFlag[i]=false;
                }
            }
        }    
        this.OrganVoxels[orgIndex]=auxVoxels-1;
        writeLine(";\n", bwAMPL);
        bwAMPL.close();
    }
    
    public void writeDDM_Plain(Organs[] organ, int orgIndex) throws IOException{
        this.tumorVoxelFlag= new boolean[this.OrganVoxels[0]];
        BufferedWriter bw = null ;
        String dir = "DDM_"+orgIndex+".txt";
        File file = new File(dir);
        if (!file.exists()) {
            bw = new BufferedWriter(new FileWriter(dir));
        }else{
            file.delete();
            bw = new BufferedWriter(new FileWriter(dir));
        }
        int auxVoxels=1;
        for (int i=0; i<this.OrganVoxels[orgIndex];i++){
            double auxSum=0;
            if (orgIndex==0){ //we only evaluate voxels from tumor
                for (int z=0;z<this.totalBeamlets;z++){
                    auxSum = auxSum + this.ddm0[i][1][z];
                }
            }
            //only "relevant" voxels from the tumor and all voxels from OARs
            if (auxSum>0.01 || orgIndex>0){
                if (orgIndex==0){
                    for (int z=0;z<this.totalBeamlets;z++){
                        writeLine(this.ddm0[i][1][z]+ "\t", bw);                
                    }
                    auxVoxels++;
                    writeLine("\n", bw);
                }
            }
        }    
        this.OrganVoxels[orgIndex]=auxVoxels-1;
        bw.close();
    }*/
    
    public static void writeLine(String l, BufferedWriter bw) throws IOException{
	
	String row = "";
	row = l ; 
	bw.write(row );
	

    }
    
    
    
   
    
}

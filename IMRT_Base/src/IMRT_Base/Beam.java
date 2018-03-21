package IMRT_Base;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

/**
 *
 * @author gcab623
 */
public class Beam {
    public int index;
    public int beamlets;
    public double[][] beamletsCoord;
    public ArrayList<double[]>[] intensity; 
    //public ArrayList<ArrayList<double[]>> intensity; //1stCol= VoxelIndex -- 2ndCol = beamletIndex -- 3rdCol= intensity
    public int init;
    public int minAbsX;
    public int maxAbsX;
    public int minAbsY;
    public int maxAbsY;
    
    
    
    public Beam(int bl, int ang, String path) throws IOException{
        this.beamlets=bl;
        this.index=ang;
        this.init=0;//depends on the set of angles
        this.beamletsCoord = new double[bl][5];
        readCoordinates(path);
    }
    
    public Beam(int bl, int ang) throws IOException{
        this.beamlets=bl;
        this.index=ang;
        this.init=0;//depends on the set of angles
    }
    
    public void readCoordinates(String s) throws FileNotFoundException, IOException{
        if (this.beamletsCoord==null){
            this.beamletsCoord= new double[this.beamlets][5];
        }
        String[] auxReader = new String[3];
        String sp="\t";
        String dir = s + "beamletCoord" +"_"+this.index+".txt";
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine();
        int i=0;
        do{
            auxReader = line.split(sp);
            this.beamletsCoord[i][0] = i; //beamletIndex
            this.beamletsCoord[i][1] = Double.parseDouble(auxReader[6]); //Relative xCoord
            this.beamletsCoord[i][2] = Double.parseDouble(auxReader[7]); //Relative yCoord
            this.beamletsCoord[i][3] = Double.parseDouble(auxReader[8]); //xCoord
            this.beamletsCoord[i][4] = Double.parseDouble(auxReader[9]); //yCoord
            i++;
            line=fileIn.readLine();
        }while(line != null);  
    }
    public void readIntensities(String s, Organs org) throws FileNotFoundException, IOException{
            
        this.intensity[org.index] = new ArrayList<>();
        String[] auxReader = new String[3];
        String sp="\t";
        String dir = s + org.name +"_"+this.index+".txt";
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        String line = "";
        line=fileIn.readLine();
        while(line != null){
            double[] auxIntensity = new double[3];
            auxReader = line.split(sp);
            auxIntensity[0] = Double.parseDouble(auxReader[0]); //voxelIndex
            auxIntensity[1] = Double.parseDouble(auxReader[1]); //beamletIndex
            auxIntensity[2] = Double.parseDouble(auxReader[2]); //intensity
            this.intensity[org.index].add(auxIntensity); 
            line=fileIn.readLine();
        }
        fileIn.close();
    }
    
    public void readIntensities(String[] s, Organs org, ArrayList<Long> voxelIndex) throws FileNotFoundException, IOException{
        
        
        String[] auxReader = new String[3]; //VoxelIndex - BeamletIndex - Intensity
        String sp="\t";
        String dir = new String();
        dir = s[1] + org.name +"_"+this.index+".txt";
        File f = new File(dir);
        BufferedReader fileIn = new BufferedReader(new FileReader(f));
        this.intensity[org.index] = new ArrayList<>();
        dir = s[1] + org.name +"_"+this.index+".txt";
        String line = "";
        line=fileIn.readLine();
        //this.intensity[y] = new ArrayList<>();
        while(line != null){
            double[] auxIntensity = new double[3];
            auxReader = line.split(sp);
            auxIntensity[0] = Double.parseDouble(auxReader[0]); //voxelIndex
            auxIntensity[1] = Double.parseDouble(auxReader[1]); //beamletIndex
            auxIntensity[2] = Double.parseDouble(auxReader[2]); //intensity
            this.intensity[org.index].add(auxIntensity);
            //Generate Voxel Index of organ 'y' using beam angle 'j'
            /*if (!voxelIndex.contains((long)auxIntensity[0])){
                    voxelIndex.add((long)auxIntensity[0]);
            }*/
            line=fileIn.readLine();
        }
        Collections.sort(voxelIndex);
        fileIn.close();
    }
    
    public int getIndexOf(long origVoxel, int orgIndex) {
        int posLeft, posRight, pivot=0;
        posLeft=0;posRight = this.intensity[orgIndex].size();
        while (posLeft < posRight) {
            pivot = posLeft + (int) (posRight- posLeft)/2 ;
            if (this.intensity[orgIndex].get(pivot)[0] < origVoxel){ //move posLeft to pivot
                posLeft = pivot;
            }else if (this.intensity[orgIndex].get(pivot)[0] > origVoxel){ //move posLeft to pivot
                posRight = pivot;
            }else if (this.intensity[orgIndex].get(pivot)[0] == origVoxel){
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
        //We need to find the first index corresponding to the 'origVoxel'
        while (this.intensity[orgIndex].get(pivot)[0] == origVoxel){
            pivot--;
             if (pivot <0){
                 break;
             }
        }
        
        pivot++;
        return pivot; 
    }

}

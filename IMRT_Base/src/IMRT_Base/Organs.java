package IMRT_Base;
/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author gcab623
 */
public class Organs {
    public String name;
    public int index;
    public double weight; //In case weighted sum method is considered 0<weight<1
    public double voxelEnd;
    public double totalDose;
    public double actualMinDose;
    public double actualMaxDose;
    public double doseUB;
    public double doseLB;
    public double desiredDose; //Maximum dose for OARs and Minimum Dose for Target
    public int a;
    public int v;
    public int totalVoxels;
    public boolean isTarget;

    public  Organs(String n, int i, double vI, double vE, double TD, double mD, double MD, double UB, double LB, double dD, int aPar, int vPar, int tV, boolean isTar){
        this.name=n; this.index= i; this.weight = vI; this.voxelEnd=vE;
        this.totalDose=TD; this.actualMinDose=mD; this.actualMaxDose=MD;
        this.doseUB = UB; this.doseLB=LB; this.desiredDose=dD; this.a=aPar; this.v = vPar;
        this.totalVoxels = tV;
        this.isTarget = isTar;  
    }
    public void updateOrgans(Organs o){
        this.name = o.name;
        this.index = o.index;
        this.weight = o.weight;
        this.voxelEnd = o.voxelEnd;
        this.totalDose = o.totalDose;
        this.actualMinDose = o.actualMinDose;
        this.actualMaxDose = o.actualMaxDose;
        this.doseUB = o.doseUB;
        this.doseLB = o.doseLB;
        this.desiredDose = o.desiredDose;
        this.a = o.a;
        this.v = o.v;
        this.totalVoxels = o.totalVoxels;
        this.isTarget = o.isTarget;
    }
}


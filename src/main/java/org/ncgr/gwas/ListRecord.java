package org.ncgr.gwas;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;

/**
 * Encapsulates a single PLINK LIST record, with the various case and control counts, plus stats.
 */
public class ListRecord {

    public String contig;
    public String id;
    public String hetGenotype;
    public int caseREF, controlREF;
    public int caseHET, controlHET;
    public int caseHOM, controlHOM;
    public int caseNC, controlNC;
    public double stdStat;
    public double pValue;
    public double oddsRatio;
    
    /**
     * Construct from individual fields, with genotypes, cases, and controls in single strings with | separator.
     */
    // ListRecord(
    public ListRecord(String contig, String id, String hetGenotype,
		      int caseREF, int controlREF,
		      int caseHET, int controlHET,
		      int caseHOM, int controlHOM,
		      int caseNC, int controlNC,
		      double stdStat, double pValue) {
	this.contig = contig;
	this.id = id;
	this.hetGenotype = hetGenotype;
	this.caseREF = caseREF;
	this.controlREF = controlREF;
	this.caseHET = caseHET;
	this.controlHET = controlHET;
	this.caseHOM = caseHOM;
	this.controlHOM = controlHOM;
	this.caseNC = caseNC;
	this.controlNC = controlNC;
	this.stdStat = stdStat;
	this.pValue = pValue;
	this.oddsRatio = getOddsRatio();
    }

    /**
     * Construct from an output line.
     */
    public ListRecord(String line) {
	String[] fields = line.split("\t");
	this.contig = fields[0];
	this.id = fields[1];
	this.hetGenotype = fields[2];
	this.caseREF = Integer.parseInt(fields[3]);
	this.controlREF = Integer.parseInt(fields[4]);
	this.caseHET = Integer.parseInt(fields[5]);
	this.controlHET = Integer.parseInt(fields[6]);
	this.caseHOM = Integer.parseInt(fields[7]);
	this.controlHOM = Integer.parseInt(fields[8]);
	this.caseNC = Integer.parseInt(fields[9]);
	this.controlNC = Integer.parseInt(fields[10]);
	this.stdStat = Double.parseDouble(fields[11]);
	this.pValue = Double.parseDouble(fields[12]);
	this.oddsRatio = Double.parseDouble(fields[13]);
    }

    /**
     * Return an output line.
     */
    public String toString() {
	return contig+"\t"+id+"\t"+hetGenotype+
	    "\t"+caseREF+"\t"+controlREF+
	    "\t"+caseHET+"\t"+controlHET+
	    "\t"+caseHOM+"\t"+controlHOM+
	    "\t"+caseNC+"\t"+controlNC+
	    "\t"+stdStat+"\t"+pValue+"\t"+oddsRatio;
    }

    /**
     * Return the header line.
     */
    public static String getHeader() {
	return "Contig\tID\tHET\tCaseREF\tControlREF\tCaseHET\tControlHET\tCaseHOM\tControlHOM\tCaseNC\tControlNC\tStdStat\tp\tOR";
    }

    /**
     * Return the odds ratio
     */
    public double getOddsRatio() {
	int caseAlleles = caseHET + 2*caseHOM;
	int controlAlleles = controlHET + 2*controlHOM;
	return ((double) caseAlleles / (double) controlAlleles) /  ((double) caseREF / (double) controlREF);
    }
}

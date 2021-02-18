package org.ncgr.gwas;

import java.util.Map;
import java.util.HashMap;
import java.util.List;
import java.util.LinkedList;

/**
 * Encapsulates a single segregation record, with all genotypes at a locus and case and control counts, plus stats.
 */
public class SegRecord implements Comparable {

    public String contig;
    public int start;
    public String id;
    public List<String> genotypes = new LinkedList<>();
    public double maf;
    public int noCallCount;
    public Map<String,Integer> cases = new HashMap<>();    // keyed by genotype
    public Map<String,Integer> controls = new HashMap<>(); // keyed by genotype
    public double stdStat;
    public double pValue;

    /**
     * Construct from individual fields, with genotypes, cases, and controls in single strings with | separator.
     */
    public SegRecord(String contig, int start, String id, String genotypeString, double maf, int noCallCount, String caseString, String controlString, double stdStat, double pValue) {
	this.contig = contig;
	this.start = start;
	this.id = id;
	this.maf = maf;
	this.noCallCount = noCallCount;
	this.stdStat = stdStat;
	this.pValue = pValue;
	setGenotypes(genotypeString, caseString, controlString);
    }

    /**
     * Construct from an output line.
     * 0      1     2  3         4   5       6     7        8    9
     * contig start id genotypes maf nocalls cases controls stat p
     */
    public SegRecord(String line) {
	String[] fields = line.split("\t");
	this.contig = fields[0];
	this.start = Integer.parseInt(fields[1]);
	this.id = fields[2];
	String genotypeString = fields[3];
	this.maf = Double.parseDouble(fields[4]);
	this.noCallCount = Integer.parseInt(fields[5]);
	String caseString = fields[6];
	String controlString = fields[7];
	this.stdStat = Double.parseDouble(fields[8]);
	this.pValue = Double.parseDouble(fields[9]);
	setGenotypes(genotypeString, caseString, controlString);
    }

    /**
     * Compare based on contig then start.
     */
    public int compareTo(Object o) {
	SegRecord that = (SegRecord) o;
	if (!this.contig.equals(that.contig)) {
	    return this.contig.compareTo(that.contig);
	} else {
	    return Integer.compare(this.start, that.start);
	}
    }

    /**
     * Set the genotype values from strings separated by ":".
     */
    public void setGenotypes(String genotypeString, String caseString, String controlString) {
	String[] genotypeArray = genotypeString.split(":");
	String[] caseArray = caseString.split(":");
	String[] controlArray = controlString.split(":");
	for (int i=0; i<genotypeArray.length; i++) {
	    String genotype = genotypeArray[i];
	    genotypes.add(genotype);
	    try {
		cases.put(genotype, Integer.parseInt(caseArray[i]));
		controls.put(genotype, Integer.parseInt(controlArray[i]));
	    } catch (NumberFormatException ex) {
		System.err.println(ex);
		System.err.println("genotypeString="+genotypeString);
		System.err.println("caseString="+caseString);
		System.err.println("controlString="+controlString);
		System.exit(1);
	    }
	}
    }

    /**
     * Return an output line.
     * 0      1     2  3         4   5       6     7        8    9
     * contig start id genotypes maf nocalls cases controls stat p
     */
    public String toString() {
	String genotypeString = "";
	String caseString = "";
	String controlString = "";
	for (String genotype : genotypes) {
	    if (genotypeString.length()>0) {
		genotypeString += "|";
		caseString += "|";
		controlString += "|";
	    }
	    genotypeString += genotype;
	    caseString += cases.get(genotype);
	    controlString += controls.get(genotype);
	}
	//     0           1          2       3                   4        5                6               7                  8            9
	return contig+"\t"+start+"\t"+id+"\t"+genotypeString+"\t"+maf+"\t"+noCallCount+"\t"+caseString+"\t"+controlString+"\t"+stdStat+"\t"+pValue;
    }

    /**
     * Return a map of genotype to odds ratio
     */
    public Map<String,Double> getOddsRatios() {
	// calculate the odds ratio for each genotype
	Map<String,Double> oddsRatios = new HashMap<>();
	for (String genotype : genotypes) {
	    int otherCases = 0;
	    int otherControls = 0;
	    for (String otherGenotype : genotypes) {
		if (otherGenotype.equals(genotype)) continue;
		otherCases += cases.get(otherGenotype);
		otherControls += controls.get(otherGenotype);
	    }
	    double oddsRatio = ((double) cases.get(genotype) / (double) controls.get(genotype)) /  ((double) otherCases / (double) otherControls);
	    oddsRatios.put(genotype, oddsRatio);
	}
	return oddsRatios;
    }
}

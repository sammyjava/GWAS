package org.ncgr.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;
import java.util.TreeSet;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.ConcurrentSkipListSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.samtools.util.CloseableIterator;

/**
 * Compute probability risk scores for the subjects in a study with VCF segregation data already computed.
 */
public class SegPRS {
    static int DEFAULT_MAX_NOCALLS = 1000;
    static double DEFAULT_MIN_MAF = 0.01;

    /**
     * Main class outputs a tab-delimited list of subjects and risk scores.
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

	Option segFileOption = new Option("sf", "segfile", true, "VCFSegregation output file");
	segFileOption.setRequired(true);
	options.addOption(segFileOption);
	//					 
        Option vcfFileOption = new Option("vf", "vcffile", true, "VCF file");
        vcfFileOption.setRequired(true);
        options.addOption(vcfFileOption);
	//
	Option labelFileOption = new Option("lf", "labelfile", true, "case/control labels file");
	labelFileOption.setRequired(true);
	options.addOption(labelFileOption);
	//
	Option minMAFOption = new Option("maf", "minmaf", true, "minimum MAF for a locus to be output ("+DEFAULT_MIN_MAF+")");
	minMAFOption.setRequired(false);
	options.addOption(minMAFOption);
        //
        Option maxNoCallsOption = new Option("mnc", "maxnocalls", true, "maximum number of no-calls for a locus to be output ("+DEFAULT_MAX_NOCALLS+")");
        maxNoCallsOption.setRequired(false);
        options.addOption(maxNoCallsOption);
	//
	Option contigOption = new Option("c", "contig", true, "contig/chromosome for analysis region (all)");
	contigOption.setRequired(false);
	options.addOption(contigOption);
	//
	Option startOption = new Option("s", "start", true, "starting coordinate of analysis region (full contigs)");
	startOption.setRequired(false);
	options.addOption(startOption);
	//
	Option endOption = new Option("e", "end", true, "ending coordinate of analysis region (full contigs)");
	endOption.setRequired(false);
	options.addOption(endOption);
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("SegPRS", options);
            System.exit(1);
            return;
        }

        // spit out help if nothing supplied
        if (cmd.getOptions().length==0) {
            formatter.printHelp("SegPRS", options);
            System.exit(1);
            return;
        }

        // filtering parameters
        double minMAF = DEFAULT_MIN_MAF;
        if (cmd.hasOption("minmaf")) {
            minMAF = Double.parseDouble(cmd.getOptionValue("minmaf"));
        }
        int maxNoCalls = DEFAULT_MAX_NOCALLS;
        if (cmd.hasOption("maxnocalls")) {
            maxNoCalls = Integer.parseInt(cmd.getOptionValue("maxnocalls"));
        }

	// region parameters
	String contig = null;
	if (cmd.hasOption("contig")) {
	    contig = cmd.getOptionValue("contig");
	}
	int start = 0;
	int end = 0;
	if (cmd.hasOption("start") && cmd.hasOption("end")) {
	    start = Integer.parseInt(cmd.getOptionValue("start"));
	    end = Integer.parseInt(cmd.getOptionValue("end"));
	}

	// get the labels for the subjects
	Map<String,String> sampleLabels = new HashMap<>();
	String labelFilename = cmd.getOptionValue("labelfile");
	BufferedReader labelReader = new BufferedReader(new FileReader(labelFilename));
	String labelLine = null;
	while ((labelLine=labelReader.readLine())!=null) {
	    String[] fields = labelLine.split("\t");
	    if (fields.length==2) {
		sampleLabels.put(fields[0], fields[1]);
	    }
	}

	// spin through the segregation file and build PRS for each subject at each locus
	ConcurrentSkipListMap<String,Double> samplePRS = new ConcurrentSkipListMap<>(); // keyed by sample name
	ConcurrentSkipListMap<String,Integer> sampleN = new ConcurrentSkipListMap<>();  // keyed by sample name
	VCFFileReader vcfReader = new VCFFileReader(new File(cmd.getOptionValue("vcffile")));
	String segFilename = cmd.getOptionValue("segfile");
	BufferedReader segReader = new BufferedReader(new FileReader(segFilename));
	String segLine = null;
	while ((segLine=segReader.readLine())!=null) {
	    if (segLine.startsWith("#")) continue;
	    SegRecord rec = new SegRecord(segLine);
	    if (contig!=null && !rec.contig.equals(contig)) continue;
	    if (start!=0 && rec.start<start) continue;
	    if (end!=0 && rec.start>end) continue;
	    if (rec.noCallCount>maxNoCalls) continue;
	    // exclude zero and +-infinity odds ratio genotypes
	    Map<String,Double> oddsRatiosIncludingInfinity = rec.getOddsRatios();
	    Map<String,Double> logOddsRatios = new HashMap<>();
	    for (String genotype : oddsRatiosIncludingInfinity.keySet()) {
		double or = oddsRatiosIncludingInfinity.get(genotype);
		if (or!=Double.POSITIVE_INFINITY && or!=Double.NEGATIVE_INFINITY && Math.abs(or)>0.0) {
		    logOddsRatios.put(genotype, Math.log10(or));
		}
	    }
	    // DEBUG
	    System.err.println(segLine+logOddsRatios);
	    //
	    // we're supposed to close the CloseableIterater when done
	    CloseableIterator ci = vcfReader.query(rec.contig, rec.start, rec.start);
	    try {
		for (Object o : ci.toList()) {
		    VariantContext vc = (VariantContext) o;
		    // require at least two genotypes
		    if (vc.getGenotypes().size()<2) continue;
		    // minimum MAF criterion
		    int calledCount = vc.getCalledChrCount();
		    int numAboveMAF = 0;
		    for (Allele a : vc.getAlleles()) {
			int calledAlleleCount = vc.getCalledChrCount(a);
			double maf = (double)calledAlleleCount / (double)calledCount;
			if (maf>minMAF) numAboveMAF++; // includes max allele, usually REF
		    }
		    if (numAboveMAF<2) continue;
		    // update the PRS for each genotype/sample
		    ConcurrentSkipListSet concurrentGenotypes = new ConcurrentSkipListSet<>(vc.getGenotypes());
		    concurrentGenotypes.parallelStream().forEach(obj -> {
			    Genotype g = (Genotype) obj;
			    if (!g.isNoCall()) {
				String gString = g.getGenotypeString();
				String sampleName = g.getSampleName();
				if (sampleLabels.containsKey(sampleName) && logOddsRatios.containsKey(gString)) {
				    double logOR = logOddsRatios.get(gString);
				    if (samplePRS.containsKey(sampleName)) {
					double prs = samplePRS.get(sampleName);
					prs += logOR;
					samplePRS.put(sampleName, prs);
					int n = sampleN.get(sampleName);
					n++;
					sampleN.put(sampleName, n);
				    } else {
					samplePRS.put(sampleName, logOR);
					sampleN.put(sampleName, 1);
				    }
				}
			    }
			});
		}
	    } catch (Exception e) {
		System.err.println(e);
	    } finally {
		ci.close();
	    }
	}
	vcfReader.close();
	segReader.close();

	// output results
	for (String sampleName : samplePRS.keySet()) {
	    System.out.println(sampleName+"\t"+sampleLabels.get(sampleName)+"\t"+samplePRS.get(sampleName)/sampleN.get(sampleName));
	}
    }
}

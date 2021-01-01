package org.ncgr.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.Map;
import java.util.HashMap;
import java.util.Set;
import java.util.HashSet;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * Loads a PLINK list file and computes the Cochran-Armitage Test p-value for segregation between genotypes in a case/control experiment.
 *
 * Cases and controls are given by a labels file.
 *
 * @author Sam Hokin
 */
public class ListSegregation {
    static int DEFAULT_MAX_NOCALLS = 1000;
    static double DEFAULT_MIN_MAF = 0.01;

    /**
     * Main class outputs a tab-delimited list of the contingency matrix for each locus, plus the Cochran-Armitage trend test p value.
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

        Option listFileOption = new Option("list", "listfile", true, "PLINK list output file");
        listFileOption.setRequired(true);
        options.addOption(listFileOption);
	//
	Option minMAFOption = new Option("maf", "minmaf", true, "minimum MAF for a locus to be output ("+DEFAULT_MIN_MAF+")");
	minMAFOption.setRequired(false);
	options.addOption(minMAFOption);
        //
        Option maxNoCallsOption = new Option("mnc", "maxnocalls", true, "maximum number of no-calls for a locus to be output ("+DEFAULT_MAX_NOCALLS+")");
        maxNoCallsOption.setRequired(false);
        options.addOption(maxNoCallsOption);
        //
        Option labelFileOption = new Option("lf", "labelfile", true, "label file containing case/control labels for each subject");
        labelFileOption.setRequired(false);
        options.addOption(labelFileOption);
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("ListSegregation", options);
            System.exit(1);
            return;
        }

        // spit out help if nothing supplied
        if (cmd.getOptions().length==0) {
            formatter.printHelp("ListSegregation", options);
            System.exit(1);
            return;
        }

        // some general parameters
        double minMAF = DEFAULT_MIN_MAF;
        if (cmd.hasOption("minmaf")) {
            minMAF = Double.parseDouble(cmd.getOptionValue("minmaf"));
        }
        int maxNoCalls = DEFAULT_MAX_NOCALLS;
        if (cmd.hasOption("maxnocalls")) {
            maxNoCalls = Integer.parseInt(cmd.getOptionValue("maxnocalls"));
        }

        // true if case, false if control, keyed by sample ID used in LIST file
        Map<String,Boolean> subjectStatus = new HashMap<>();
        if (cmd.hasOption("labelfile")) {
            // read sample labels from a tab-delimited file. Comment lines start with #.
            // 28304	case
            // 60372	ctrl
            String labelFilename = cmd.getOptionValue("labelfile");
            BufferedReader labelReader = new BufferedReader(new FileReader(labelFilename));
            String labelLine;
            while ((labelLine=labelReader.readLine())!=null) {
                if (labelLine.startsWith("#")) continue;
		String[] fields = labelLine.split("\t");
		if (fields.length==2) {
		    String sampleId = fields[0];
		    boolean isCase = fields[1].equals("case");
		    subjectStatus.put(sampleId, isCase);
		}
            }
            labelReader.close();
        }

        // find the desired sample names as they appear in the LIST file
	// 6	AA_A_9_30018537_FS	AP	42	42	58	58	76	76	107	107	114
	Set<String> sampleNames = new HashSet<>();
	BufferedReader listReader = new BufferedReader(new FileReader(cmd.getOptionValue("listfile")));
        String listLine;
	while ((listLine=listReader.readLine())!=null) {
	    String[] fields = listLine.split("\\t");
	    String contig = fields[0];
	    String id = fields[1];
	    String genotype = fields[2];
	    for (int i=3; i<fields.length; i++) sampleNames.add(fields[i]);
	}
	listReader.close();
	Set<String> caseSampleNames = new HashSet<>();     // case subjects in the LIST file
        Set<String> controlSampleNames = new HashSet<>();  // control subjects in the LIST file
	int nCases = 0;
	int nControls = 0;
        for (String sampleName : subjectStatus.keySet()) {
            boolean found = false;
            if (sampleNames.contains(sampleName)) {
                found = true;
                if (subjectStatus.get(sampleName)) {
                    caseSampleNames.add(sampleName);
		    nCases++;
                } else {
                    controlSampleNames.add(sampleName);
		    nControls++;
                }
            }
            if (!found) System.err.println("Subject "+sampleName+" NOT FOUND in LIST file.");
	}
	System.err.println("Found "+nCases+" cases and "+nControls+" controls in LIST file.");
	
	// spin through the LIST file and get stats on qualified loci
	// sample names appear as family individual
	// there are four lines per id: minority HOM genotype, HET, majority (REF), no-calls
	// 6	AA_A_9_30018537_FS	AA	174	174	509	509	1099	1099	1360
	// 6	AA_A_9_30018537_FS	AP	42	42	58	58	76	76	107
	// 6	AA_A_9_30018537_FS	PP	45	45	55	55	57	57	59
	// 6	AA_A_9_30018537_FS	00
	// 6	AA_A_9_30018537_S	PP	800	800	1760	1760	2468	2468	3996
	// 6	AA_A_9_30018537_S	PA	55	55	58	58	70	70	76
	// 6	AA_A_9_30018537_S	AA	42	42	45	45	57	57	59
	// 6	AA_A_9_30018537_S	00
	System.out.println(ListRecord.getHeader());
	int ncRejects = 0;
	int mafRejects = 0;
	listReader = new BufferedReader(new FileReader(cmd.getOptionValue("listfile")));
	listLine = null;                                                                                                                                                                                    
        while ((listLine=listReader.readLine())!=null) {
	    // count HOM cases and controls
	    String[] fields = listLine.split("\\t");
	    String contig = fields[0];
	    String id = fields[1];
	    if (!fields[1].equals(id)) {
		System.err.println("Error: expected id="+id+" but found "+fields[1]+" instead. Aborting.");
		System.exit(1);
	    }
	    String homGenotype = fields[2];
	    int caseHOM = 0;
	    int controlHOM = 0;
	    for (int i=3; i<fields.length; i++) {
		String sampleName = fields[i];
		if (caseSampleNames.contains(sampleName)) {
		    caseHOM++;
		} else if (controlSampleNames.contains(sampleName)) {
		    controlHOM++;
		}
	    }
	    // count HET cases and controls
	    listLine = listReader.readLine();
	    fields = listLine.split("\\t");
	    if (!fields[1].equals(id)) {
		System.err.println("Error: expected id="+id+" but found "+fields[1]+" instead. Aborting.");
		System.exit(1);
	    }
	    String hetGenotype = fields[2];
	    int caseHET = 0;
	    int controlHET = 0;
	    for (int i=3; i<fields.length; i++) {
		String sampleName = fields[i];
		if (caseSampleNames.contains(sampleName)) {
		    caseHET++;
		} else if (controlSampleNames.contains(sampleName)) {
		    controlHET++;
		}
	    }
	    // count REF cases and controls
	    listLine = listReader.readLine();
	    fields = listLine.split("\\t");
	    if (!fields[1].equals(id)) {
		System.err.println("Error: expected id="+id+" but found "+fields[1]+" instead. Aborting.");
		System.exit(1);
	    }
	    String refGenotype = fields[2];
	    int caseREF = 0;
	    int controlREF = 0;
	    for (int i=3; i<fields.length; i++) {
		String sampleName = fields[i];
		if (caseSampleNames.contains(sampleName)) {
		    caseREF++;
		} else if (controlSampleNames.contains(sampleName)) {
		    controlREF++;
		}
	    }
	    // count NC cases and controls
	    listLine = listReader.readLine();
	    fields = listLine.split("\\t");
	    if (!fields[1].equals(id)) {
		System.err.println("Error: expected id="+id+" but found "+fields[1]+" instead. Aborting.");
		System.exit(1);
	    }
	    if (!fields[2].equals("00")) {
		System.err.println("Error: expected a 00 line but found a "+fields[2]+" line. Aborting.");
		System.exit(1);
	    }
	    int caseNC = 0;
	    int controlNC = 0;
	    for (int i=3; i<fields.length; i++) {
		String sampleName = fields[i];
		if (caseSampleNames.contains(sampleName)) {
		    caseNC++;
		} else if (controlSampleNames.contains(sampleName)) {
		    controlNC++;
		}
	    }
	    // divide by two for the double entry of sample names
	    caseHOM = caseHOM/2;
	    controlHOM = controlHOM/2;
	    caseHET = caseHET/2;
	    controlHET = controlHET/2;
	    caseREF = caseREF/2;
	    controlREF = controlREF/2;
	    caseNC = caseNC/2;
	    controlNC = controlNC/2;
	    // no-calls filter
	    if ((caseNC+controlNC)>maxNoCalls) {
		ncRejects++;
		continue;
	    }
	    // MAF filter
	    int totalAlleles = 2*(caseHOM + controlHOM + caseHET + controlHET + caseREF + controlREF);
	    int caseAlleles = caseHET + 2*caseHOM;
	    int controlAlleles = controlHET + 2*controlHOM;
	    double maf = (double) (caseAlleles+controlAlleles) / (double) (totalAlleles);
	    if (maf<minMAF) {
		mafRejects++;
		continue;
	    }
	    // Cochran-Armitage test
	    int[][] countTable = new int[2][3]; // 2 labels x 3 genotypes
	    int[] weights = new int[3];
	    // straight allelic association (not additive)
	    weights[0] = 0;
	    weights[1] = 1;
	    weights[2] = 1;
	    CochranArmitage ca = new CochranArmitage(weights);
	    countTable[0][0] = controlREF;
	    countTable[1][0] = caseREF;
	    countTable[0][1] = controlHET;
	    countTable[1][1] = caseHET;
	    countTable[0][2] = controlHOM;
	    countTable[1][2] = caseHOM;
	    double pValue = ca.test(countTable);
	    // we can still get a few cases with 0 alternative counts
	    if (Double.isNaN(pValue)) pValue = 1.0;
	    // output the line
	    ListRecord rec = new ListRecord(contig, id, hetGenotype, caseREF, controlREF, caseHET, controlHET, caseHOM, controlHOM, caseNC, controlNC, ca.standardStatistic, pValue);
	    System.out.println(rec);
	}
	listReader.close();
	System.err.println("maxNC rejects="+ncRejects+" minMAF rejects="+mafRejects);
    }
}

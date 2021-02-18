package org.ncgr.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.Arrays;
import java.util.List;
import java.util.LinkedList;
import java.util.Map;
import java.util.HashMap;
import java.util.Optional;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeType;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;

import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;

import org.mskcc.cbio.portal.stats.FisherExact;

/**
 * Spits out the subject IDs with case or control status per disease value, from a dbGaP study.
 * Disease status are given by a phenotype file in dbGaP format.
 * Sex and Race may also be specified as filters.
 * Subjects with neither case nor control status are marked "unkn".
 * If maxcases is given, then only maxcases cases and controls will be stored without any "unkn".
 *
 * @author Sam Hokin
 */
public class PhenoSubjects {

    String phenofilename;
    String labelsfilename;
    String samplefilename;
    String sampleVar;

    String caseValue, controlValue;
    String diseaseVar, diseaseName;
    String desiredSexValue, desiredRaceValue;

    int maxCases;

    String[] ccVars;
    List<String> ignoredSamples = new LinkedList<>();
    Map<String,String> sampleIdMap = new HashMap<>(); // keyed by dbGaPSubjectId=dbGaP_Subject_ID
    Map<String,String[]> subjectStatus = new HashMap<>(); // status array (case, ctrl, unkn) keyed by study ID

    /**
     * Main class outputs a tab-delimited list of the subject IDs and disease or case/control status.
     */
    public static void main(String[] args) throws FileNotFoundException, IOException {
        Options options = new Options();
        CommandLineParser parser = new DefaultParser();
        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd;

	Option sampleFileOption = new Option("sf", "samplefile", true, "dbGaP samples file (needed if contains mapping from dbGaP_Subject_ID to sample ID used in VCF file)");
	sampleFileOption.setRequired(false);
	options.addOption(sampleFileOption);
        //
	Option sampleVarOption = new Option("sv", "samplevar", true, "study sample ID variable in dbGaP samples file (e.g. SAMPID; required if -sf)");
	sampleVarOption.setRequired(false);
	options.addOption(sampleVarOption);
	//
        Option phenoFileOption = new Option("pf", "phenofile", true, "dbGaP phenotype file");
        phenoFileOption.setRequired(true);
        options.addOption(phenoFileOption);
        //
        Option ccVarOption = new Option("ccv", "casecontrolvar", true, "(comma-separated) case/control variables in dbGaP phenotype file (e.g. ANALYSIS_CAT or pheno_241.2)");
	ccVarOption.setRequired(true);
        options.addOption(ccVarOption);
	// 
        // NOTE: this only allows a single value of case or control in the segregating variable! (Some files have control=1, say, and several case values.)
        Option caseValueOption = new Option("caseval", true, "case value in dbGaP phenotype file (e.g. Case)");
        caseValueOption.setRequired(true);
        options.addOption(caseValueOption);
        //
        Option controlValueOption = new Option("controlval", true, "control value in dbGaP phenotype file (e.g. Control)");
        controlValueOption.setRequired(true);
        options.addOption(controlValueOption);
	//
	Option diseaseVarOption = new Option("dv", "diseasevar", true, "disease variable in dbGaP phenotype file (e.g. PRIMARY_DISEASE; required if -dn)");
        diseaseVarOption.setRequired(false);
        options.addOption(diseaseVarOption);
	//
	Option diseaseNameOption = new Option("dn", "diseasename", true, "desired case disease name in dbGaP phenotype file (e.g. Schizophrenia; required if -dv)");
        diseaseNameOption.setRequired(false);
        options.addOption(diseaseNameOption);
	//
	Option desiredSexOption = new Option("ds", "desiredsex", true, "phenotype value of desired sex (e.g. M or 1)");
	desiredSexOption.setRequired(false);
	options.addOption(desiredSexOption);
	//
	Option desiredRaceOption = new Option("dr", "desiredrace", true, "phenotype value of desired race (e.g. W)");
	desiredRaceOption.setRequired(false);
	options.addOption(desiredRaceOption);
	//
	Option maxCasesOption = new Option("mc", "maxcases", true, "maxiumum number of cases for cases=controls output (and no unknowns)");
	maxCasesOption.setRequired(false);
	options.addOption(maxCasesOption);
	//
	Option otherSubjectsOption = new Option("os", "othersubjects", true, "output case/control subjects other than those in the given labels file");
	otherSubjectsOption.setRequired(false);
	options.addOption(otherSubjectsOption);
	
        try {
            cmd = parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println(e.getMessage());
            formatter.printHelp("PhenoSubjects", options);
            System.exit(1);
            return;
        }

        // spit out help if nothing supplied
        if (cmd.getOptions().length==0) {
            formatter.printHelp("PhenoSubjects", options);
            System.exit(1);
            return;
        }

	PhenoSubjects ps = new PhenoSubjects();

        // comma-separated case-control variables
	ps.ccVars = cmd.getOptionValue("casecontrolvar").split(",");
        
	ps.caseValue = cmd.getOptionValue("caseval");
	ps.controlValue = cmd.getOptionValue("controlval");
	ps.sampleVar = cmd.getOptionValue("samplevar");

	if (cmd.hasOption("diseasevar")) {
	    ps.diseaseVar = cmd.getOptionValue("diseasevar");
	    ps.diseaseName = cmd.getOptionValue("diseasename");
	}

	if (cmd.hasOption("desiredsex")) {
	    ps.desiredSexValue = cmd.getOptionValue("desiredsex");
	}

	if (cmd.hasOption("desiredrace")) {
	    ps.desiredRaceValue = cmd.getOptionValue("desiredrace");
	}

	// max case/controls to output
	if (cmd.hasOption("maxcases")) {
	    ps.maxCases = Integer.parseInt(cmd.getOptionValue("maxcases"));
	}

	// can only do the maxCases thing with a single disease
	if (ps.ccVars.length>1 && ps.maxCases>0) {
	    System.err.println("ERROR: you may only specify a single disease with maxcases>0");
	    System.exit(1);
	}

	// get the List of ignored samples
	if (cmd.hasOption("othersubjects")) {
	    ps.labelsfilename = cmd.getOptionValue("othersubjects");
	    ps.loadIgnoredSamples();
	}

	// get the sample ID map
	if (cmd.hasOption("samplefile")) {
	    ps.samplefilename = cmd.getOptionValue("samplefile");
	    ps.loadSampleIdMap();
	}

	// get the subject status map from the pheno file
	ps.phenofilename = cmd.getOptionValue("phenofile");
	ps.loadSubjectStatus();

        // two-column labels file output
	ps.printLabelsFile();
    }

    /**
     * The optional sample file relates dbGaP_Subject_ID in the phenotypes file to the sample ID.
     *
     * dbGaP_Subject_ID dbGaP_Sample_ID	BioSample Accession  SUBJID	SAMPID	SAMP_SOURCE SOURCE_SAMPID SAMPLE_USE
     * 1284423	    1836728	        SAMN03897975	     PT-1S8D	28278	KAROLINSKA  28278	  Seq_DNA_WholeExome; Seq_DNA_SNP_CNV
     *
     * phs001949.v1.pht009705.v1.p1.P3DT_Sample.MULTI.txt
     * dbGaP_Subject_ID  dbGaP_Sample_ID  BioSample Accession  SUBJECT_ID  SAMPLE_ID
     * 3110589           3768892          SAMN13583397         1           42
     *
     * NOTE: there may be MORE THAN ONE LINE for the same dbGaP_Subject_ID! We'll assume that SAMPLE IDs are unique.
     */
    void loadSampleIdMap() throws FileNotFoundException, IOException {
	BufferedReader sampleReader = new BufferedReader(new FileReader(samplefilename));
	int sampleVarOffset = -1;
	boolean headerLine = true;
	String line = null;
	while ((line=sampleReader.readLine())!=null) {
	    if (line.startsWith("#")) {
		continue; // comment
	    } else if (line.trim().length()==0) {
		continue; // blank
	    } else if (headerLine) {
		// variable header
		String[] vars = line.split("\t");
		for (int i=0; i<vars.length; i++) {
		    if (vars[i].equals(sampleVar)) sampleVarOffset = i;
		}
		headerLine = false;
	    } else {
		String[] data = line.split("\t");
		String dbGaPSubjectId = data[0]; // assume first column is dbGaP_Subject_ID, which I hope is always true -- not necessarily unique!!
		String sampleId = data[sampleVarOffset]; // presume this is unique
		sampleIdMap.put(sampleId, dbGaPSubjectId);
	    }
	}
	sampleReader.close();
    }

    /**
     * The required phenotypes file provides case/control information per sample (per disease).
     * 
     * with DISEASE column:
     * dbGaP_Subject_ID SUBJID  SEX PRIMARY_DISEASE   ANALYSIS_CAT SITE  Coverage_Pass
     * 1287483          PT-FJ7E M   Bipolar_Disorder  Case         BROAD N
     *
     * without DISEASE column:
     * dbGaP_SubjID HASGENSP AMDSTAT CATARACT COR ID2  LPSCBASE LCOR  LCORBASE LCORSCORE LNUC  LNUCBASE LNUCSCORE LPSC  LPSCSCORE NUC PSC ...
     * 53181        Y        7       5        1   1890 0        COR-C 0.6      8.83      NUC-B 1.95     4.34      PSC-C 0         1   2   ...
     *
     * ##                            phv00421304.v1.p1.c1  phv00421305.v1.p1.c1  phv00421306.v1.p1.c1  phv00421307.v1.p1.c1  phv00421308.v1.p1.c1  phv00421309.v1.p1.c1  phv00421310.v1.p1.c1  ...
     * dbGaP_Subject_ID  SUBJECT_ID  age_alt               SEX                   RACE                  PC1                   PC2                   pheno_008             pheno_008.5 ...
     * 3110589           1           74                    M                     W                     -0.01                 0.00                  1                     1           ...
     */
    void loadSubjectStatus() throws FileNotFoundException, IOException {
        String line = "";
        boolean headerLine = true;
	int sexVarOffset = -1;
	int raceVarOffset = -1;
	int diseaseVarOffset = -1;
        int[] ccVarOffsets = new int[ccVars.length];
	BufferedReader phenoReader = new BufferedReader(new FileReader(phenofilename));
        while ((line=phenoReader.readLine())!=null) {
            if (line.startsWith("#")) {
                continue; // comment
            } else if (line.trim().length()==0) {
                continue; // blank
            } else if (headerLine) {
                // variable header
                String[] vars = line.split("\t");
                for (int i=0; i<vars.length; i++) {
                    for (int j=0; j<ccVars.length; j++) {
                        if (vars[i].equals(ccVars[j])) ccVarOffsets[j] = i;
                    }
		    if (vars[i].equals("SEX")) sexVarOffset = i;
		    if (vars[i].equals("RACE")) raceVarOffset = i;
		    if (diseaseVar!=null && vars[i].equals(diseaseVar)) diseaseVarOffset = i;
                }
		// DEBUG
		System.err.println("diseaseVar="+diseaseVar+" offset="+diseaseVarOffset);
		for (int i=0; i<ccVars.length; i++) System.err.println("ccVar="+ccVars[i]+"\toffset="+ccVarOffsets[i]);
		System.err.println("SEX offset="+sexVarOffset);
		System.err.println("RACE offset="+raceVarOffset);
		//
                headerLine = false;
            } else {
                // subject/sample line
                String[] data = line.split("\t");
		String dbGaPSubjectId = data[0]; // assume first column is dbGaP_Subject_ID, which I hope is always true
		List<String> sampleIds = new LinkedList<>(); // we may have more than one sample ID per dbGaP_Subject_ID!
		if (sampleIdMap.size()==0) {
		    // no samples file, so assume ID in the second column is study ID
		    sampleIds.add(data[1]);
		} else {
		    // spin through the records to get all the sample IDs for this dbGaPSubjectId
		    for (String sampleId : sampleIdMap.keySet()) {
			String dgsId = sampleIdMap.get(sampleId);
			if (dgsId.equals(dbGaPSubjectId)) sampleIds.add(sampleId);
		    }
		}
                // disease, sex, race
		boolean hasDisease = diseaseVarOffset>0;
		boolean hasSex = sexVarOffset>0;
		boolean hasRace = raceVarOffset>0;
		String diseaseValue = null;
		String sexValue = null;
		String raceValue = null;
		if (hasDisease) diseaseValue = data[diseaseVarOffset];
		if (hasSex) sexValue = data[sexVarOffset];
		if (hasRace) raceValue = data[raceVarOffset];
		boolean blankDisease = hasDisease && diseaseValue.trim().length()==0;
                boolean isDisease = hasDisease && diseaseValue.contains(diseaseName);
		boolean isDesiredSex = hasSex && sexValue.equals(desiredSexValue);
		boolean isDesiredRace = hasRace && raceValue.equals(desiredRaceValue);
		// subject status for each of the ccVars
		String status[] = new String[ccVars.length];
		for (int j=0; j<ccVars.length; j++) {
		    String ccValue = data[ccVarOffsets[j]];
		    boolean isCase = isDisease && ccValue.equals(caseValue);  // case for this ccVar
		    boolean isControl = (isDisease || blankDisease) && ccValue.equals(controlValue); // control for this ccVar
		    if (isCase) {
			status[j] = "case";
		    } else if (isControl) {
			status[j] = "ctrl";
		    } else {
			status[j] = "unkn";
		    }
		}
		// potentially multiple sampleIds for same subject
		for (String sampleId : sampleIds) {
		    subjectStatus.put(sampleId, status);
		}
            }
        }
	phenoReader.close();
    }

    /**
     * Return a list of samples to ignore.
     */
    void loadIgnoredSamples() throws FileNotFoundException, IOException {
	BufferedReader labelsReader = new BufferedReader(new FileReader(labelsfilename));
	String line = null;
	while ((line=labelsReader.readLine())!=null) {
	    if (line.startsWith("#") || line.startsWith("sample") || line.trim().length()==0) {
		continue;
	    }
	    String[] fields = line.split("\t");
	    String sampleId = fields[0];
	    ignoredSamples.add(sampleId);
	}
	labelsReader.close();
    }

    /**
     * Print a two-column labels file to STDOUT.
     */
    void printLabelsFile() {
	String header = "sample";
        for (int j=0; j<ccVars.length; j++)  {
            header += "\t"+ccVars[j];
        }
        System.out.println(header);
	// output samples to reach # cases = maxCases = # controls
	if (maxCases>0) {
	    // usually far less cases than controls, so let's pick them out from each batch individually
	    List<String> cases = new LinkedList<>();
	    List<String> controls = new LinkedList<>();
	    for (String sampleId : subjectStatus.keySet()) {
		if (!ignoredSamples.contains(sampleId)) {
		    String[] status = subjectStatus.get(sampleId);
		    if (status[0].equals("case")) {
			cases.add(sampleId);
		    } else if (status[0].equals("ctrl")) {
			controls.add(sampleId);
		    }
		}
	    }
	    // check that we've got >= maxCases
	    if (cases.size()<maxCases) {
		System.err.println("ERROR: cases.size()="+cases.size()+" < maxCases="+maxCases);
		System.exit(1);
	    }
	    if (controls.size()<maxCases) {
		System.err.println("ERROR: controls.size()="+controls.size()+" < maxCases="+maxCases);
		System.exit(1);
	    }
	    // randomly select and output cases
	    int caseCount = 0;
	    while (caseCount<maxCases) {
		Optional<String> optional = cases.stream().skip((int)(cases.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    caseCount++;
		    String sampleId = optional.get();
		    System.out.println(sampleId+"\t"+"case");
		    boolean removed = cases.remove(sampleId);
		}
	    }
	    // randomly select and output cases
	    int ctrlCount = 0;
	    while (ctrlCount<maxCases) {
		Optional<String> optional = controls.stream().skip((int)(controls.size()*Math.random())).findFirst();
		if (optional.isPresent()) {
		    ctrlCount++;
		    String sampleId = optional.get();
		    System.out.println(sampleId+"\t"+"ctrl");
		    boolean removed = controls.remove(sampleId);
		}
	    }
	} else {
	    // output all samples except those to be ignored
	    for (String sampleId : subjectStatus.keySet()) {
		if (!ignoredSamples.contains(sampleId)) {
		    String[] status = subjectStatus.get(sampleId);
		    String output = sampleId;
		    for (int j=0; j<status.length; j++) {
			output += "\t"+status[j];
		    }
		    System.out.println(output);
		}
	    }
	}
    }
}

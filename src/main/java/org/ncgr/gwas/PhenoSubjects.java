package org.ncgr.gwas;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;
import java.io.IOException;

import java.util.List;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;

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
 * Spits out the subject IDs with case or control status from a dbGaP study.
 * Cases and control status for the given disease are given by a phenotype file in dbGaP format.
 * Sex may also be specified.
 *
 * @author Sam Hokin
 */
public class PhenoSubjects {

    /**
     * Main class outputs a tab-delimited list of the subject IDs and case/control status.
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
        Option ccVarOption = new Option("ccv", "casecontrolvar", true, "case/control variable in dbGaP phenotype file (e.g. ANALYSIS_CAT)");
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
	Option debugOption = new Option("d", "debug", false, "enable debug mode");
	debugOption.setRequired(false);
	options.addOption(debugOption);
        //
	Option desiredSexOption = new Option("ds", "desiredsex", true, "phenotype value of desired sex (e.g. M or 1)");
	desiredSexOption.setRequired(false);
	options.addOption(desiredSexOption);
	//
	Option desiredRaceOption = new Option("dr", "desiredrace", true, "phenotype value of desired race (e.g. W)");
	desiredRaceOption.setRequired(false);
	options.addOption(desiredRaceOption);
	
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

	String ccVar = cmd.getOptionValue("casecontrolvar");
	String caseValue = cmd.getOptionValue("caseval");
	String controlValue = cmd.getOptionValue("controlval");
	String sampleVar = cmd.getOptionValue("samplevar");

	boolean debug = cmd.hasOption("debug");

	String diseaseVar = null;
	String diseaseName = null;
	if (cmd.hasOption("diseasevar")) {
	    diseaseVar = cmd.getOptionValue("diseasevar");
	    diseaseName = cmd.getOptionValue("diseasename");
	}

        String desiredSexValue = null;
	if (cmd.hasOption("desiredsex")) {
	    desiredSexValue = cmd.getOptionValue("desiredsex");
	}

	String desiredRaceValue = null;
	if (cmd.hasOption("desiredrace")) {
	    desiredRaceValue = cmd.getOptionValue("desiredrace");
	}

	// the optional sample file relates dbGaP_Subject_ID in the phenotypes file to the sample ID.
	//
	// dbGaP_Subject_ID dbGaP_Sample_ID	BioSample Accession  SUBJID	SAMPID	SAMP_SOURCE SOURCE_SAMPID SAMPLE_USE
	// 1284423	    1836728	        SAMN03897975	     PT-1S8D	28278	KAROLINSKA  28278	  Seq_DNA_WholeExome; Seq_DNA_SNP_CNV
	//
	// phs001949.v1.pht009705.v1.p1.P3DT_Sample.MULTI.txt
	// dbGaP_Subject_ID  dbGaP_Sample_ID  BioSample Accession  SUBJECT_ID  SAMPLE_ID
	// 3110589           3768892          SAMN13583397         1           42
	//
	// NOTE: there may be MORE THAN ONE LINE for the same dbGaP_Subject_ID! We'll assume that SAMPLE IDs are unique.
	Map<String,String> sampleIdMap = new HashMap<>(); // keyed by dbGaPSubjectId=dbGaP_Subject_ID
	if (cmd.hasOption("samplefile")) {
	    BufferedReader sampleReader = new BufferedReader(new FileReader(cmd.getOptionValue("samplefile")));
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
	}

        // the required phenotypes file provides case/control information per sample.
	// 
	// with DISEASE column:
	// dbGaP_Subject_ID SUBJID  SEX PRIMARY_DISEASE   ANALYSIS_CAT SITE  Coverage_Pass
	// 1287483          PT-FJ7E M   Bipolar_Disorder  Case         BROAD N
	//
	// without DISEASE column:
	// dbGaP_SubjID HASGENSP AMDSTAT CATARACT COR ID2  LPSCBASE LCOR  LCORBASE LCORSCORE LNUC  LNUCBASE LNUCSCORE LPSC  LPSCSCORE NUC PSC ...
	// 53181        Y        7       5        1   1890 0        COR-C 0.6      8.83      NUC-B 1.95     4.34      PSC-C 0         1   2   ...
	//
	// ##                            phv00421304.v1.p1.c1  phv00421305.v1.p1.c1  phv00421306.v1.p1.c1  phv00421307.v1.p1.c1  phv00421308.v1.p1.c1  phv00421309.v1.p1.c1  phv00421310.v1.p1.c1  ...
	// dbGaP_Subject_ID  SUBJECT_ID  age_alt               SEX                   RACE                  PC1                   PC2                   pheno_008             pheno_008.5 ...
	// 3110589           1           74                    M                     W                     -0.01                 0.00                  1                     1           ...
        Map<String,Boolean> subjectStatus = new HashMap<>(); // true if case, false if control, keyed by study ID
        String line = "";
        boolean headerLine = true;
        int ccVarOffset = -1;
	int sexVarOffset = -1;
	int raceVarOffset = -1;
	int diseaseVarOffset = -1;
        int nCases = 0;
        int nControls = 0;
	BufferedReader phenoReader = new BufferedReader(new FileReader(cmd.getOptionValue("phenofile")));
        while ((line=phenoReader.readLine())!=null) {
            if (line.startsWith("#")) {
                continue; // comment
            } else if (line.trim().length()==0) {
                continue; // blank
            } else if (headerLine) {
                // variable header
                String[] vars = line.split("\t");
                for (int i=0; i<vars.length; i++) {
                    if (vars[i].equals(ccVar)) ccVarOffset = i;
		    if (vars[i].equals("SEX")) sexVarOffset = i;
		    if (vars[i].equals("RACE")) raceVarOffset = i;
		    if (diseaseVar!=null && vars[i].equals(diseaseVar)) diseaseVarOffset = i;
                }
		if (debug) {
		    System.err.println("diseaseVar="+diseaseVar+" offset="+diseaseVarOffset);
		    System.err.println("ccVar="+ccVar+" offset="+ccVarOffset);
		    System.err.println("SEX offset="+sexVarOffset);
		    System.err.println("RACE offset="+raceVarOffset);
		}
                headerLine = false;
            } else {
                String[] data = line.split("\t");
		String dbGaPSubjectId = data[0]; // assume first column is dbGaP_Subject_ID, which I hope is always true
		List<String> sampleIds = new ArrayList<>(); // we may have more than one sample ID per dbGaP_Subject_ID!
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
                String ccValue = data[ccVarOffset];
		String sexValue = null;
		String raceValue = null;
		String diseaseValue = null;
		if (sexVarOffset>0) sexValue = data[sexVarOffset];
		if (raceVarOffset>0) raceValue = data[raceVarOffset];
		if (diseaseVarOffset>0) diseaseValue = data[diseaseVarOffset];
                boolean isCase = ccValue.equals(caseValue);
                boolean isControl = ccValue.equals(controlValue);
		boolean isDisease = diseaseValue==null || diseaseValue.contains(diseaseName);
		boolean isDesiredSex = desiredSexValue==null || sexValue==null || sexValue.equals(desiredSexValue);
		boolean isDesiredRace = desiredRaceValue==null || raceValue==null || raceValue.equals(desiredRaceValue);
                if (((isDisease && isCase) || isControl) && isDesiredSex && isDesiredRace) {
		    for (String sampleId : sampleIds) {
			subjectStatus.put(sampleId, isCase); // true = case
		    }
		    if (isCase) {
			nCases++;
		    } else {
			nControls++;
		    }
		}
            }
        }

        for (String sampleId : subjectStatus.keySet()) {
            boolean isCase = subjectStatus.get(sampleId);
            if (isCase) {
                System.out.println(sampleId+"\tcase");
            } else {
                System.out.println(sampleId+"\tctrl");
            }
        }
    }
}

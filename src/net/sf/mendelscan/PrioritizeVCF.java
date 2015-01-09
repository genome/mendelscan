/**
 * @(#)PrioritizeVCF.java
 *
 * Copyright (c) 2013 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.mendelscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.BitSet;
import java.util.HashMap;
import java.util.TreeMap;
import java.lang.Math;

/**
 * A class for prioritizing variants in a VCF file
 *
 * @version	1.1
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */
public class PrioritizeVCF {

	public PrioritizeVCF(String[] args, HashMap<String, String> params)
	{
		String usage = "USAGE: java -jar MendelScan.jar score [VCF] OPTIONS\n" +
		"\tOPTIONS:\n" +
		"\t--vep-file\tVariant annotation in VEP format\n" +
		"\t--ped-file\tPedigree file in 6-column tab-delimited format\n" +
		"\t--gene-file\tA list of gene expression values for tissue of interest\n" +
		"\t--output-file\tOutput file to contain human-friendly scored variants\n" +
		"\t--output-vcf\tOutput file to contain scored variants in VCF format\n" +
		"\t--inheritance\tPresumed model of inheritance: dominant, recessive, x-linked [dominant]\n\n" +
		"\tSegregation Scoring: Segregation score multiplied by these values for dominant/recessive\n" +
		"\t--seg-score-case-ref\tA case sample was called reference/wild-type (0.50/0.10)\n" +
		"\t--seg-score-case-het\tA case sample was called heterozygous (NA/0.50)\n" +
		"\t--seg-score-case-hom\tA case sample was called homozygous variant (0.80/NA)\n" +
		"\t--seg-score-control-het\tA case sample was called heterozygous (0.10/NA)\n" +
		"\t--seg-score-control-hom\tA case sample was called homozygous variant (0.01/0.10)\n\n" +
		"\tPopulation Scoring: Population score for these classes defined by dbSNP information\n" +
		"\t--pop-score-novel\tVariant is not present in dbSNP according to VCF (1.00)\n" +
		"\t--pop-score-mutation\tVariant from mutation (OMIM) or locus-specific databases (0.95)\n" +
		"\t--pop-score-known\tVariant known to dbSNP but no mutation or GMAF info (0.60)\n" +
		"\t--pop-score-rare\tVariant in dbSNP with GMAF < 0.01 (0.20)\n" +
		"\t--pop-score-uncommon\tVariant in dbSNP with GMAF 0.01-0.05 (0.02)\n" +
		"\t--pop-score-common\tVariant in dbSNP with GMAF >= 0.05 (0.01)\n\n" +
		"\tAnnotation Scoring: Annotation score based on canonical or most-damaging VEP consequence\n" +
		"\t--anno-score-1\tScore for intergenic mutations [0.01]\n" +
		"\t--anno-score-2\tScore for intronic mutations [0.01]\n" +
		"\t--anno-score-3\tScore for downstream mutations [0.01]\n" +
		"\t--anno-score-4\tScore for upstream mutations [0.01]\n" +
		"\t--anno-score-5\tScore for UTR mutations [0.01]\n" +
		"\t--anno-score-6\tScore for ncRNA mutations [0.01]\n" +
		"\t--anno-score-7\tScore for miRNA mutations [0.01]\n" +
		"\t--anno-score-8\tScore for synonymous coding mutations [0.05]\n" +
		"\t--anno-score-9\tScore for splice region mutations [0.20]\n" +
		"\t--anno-score-10\tScore for nonstop mutations [1.00]\n" +
		"\t--anno-score-11\tScore for missense mutations not predicted damaging [0.80]\n" +
		"\t--anno-score-12\tScore for missense mutations damaging by 1/3 algorithms [0.95]\n" +
		"\t--anno-score-13\tScore for missense mutations damaging by 2/3 algorithms [0.95]\n" +
		"\t--anno-score-14\tScore for missense mutations damaging by 3/3 algorithms [0.95]\n" +
		"\t--anno-score-15\tScore for essential splice site mutations [1.00]\n" +
		"\t--anno-score-16\tScore for frameshift mutations [1.00]\n" +
		"\t--anno-score-17\tScore for nonsense mutations [1.00]\n";

		String vepFile = null;
		String pedFile = null;
		String geneFile = null;
		String outFile = null;
		String outVCF = null;
		String inheritanceModel = "dominant";
		Integer minDepth = 20;

		// Print usage if -h or --help invoked //
		if(params.containsKey("help") || params.containsKey("h"))
		{
			System.err.println(usage);
			return;
		}
		else
		{
			if(params.containsKey("vep-file"))
				vepFile = params.get("vep-file");

			if(params.containsKey("ped-file"))
				pedFile = params.get("ped-file");

			if(params.containsKey("gene-file"))
				geneFile = params.get("gene-file");

			if(params.containsKey("output-file"))
				outFile = params.get("output-file");

			if(params.containsKey("output-vcf"))
				outVCF = params.get("output-vcf");

			if(params.containsKey("inheritance"))
				inheritanceModel = params.get("inheritance");
		}


	    try
	    {
	    	// First get input file //
	    	BufferedReader in = MendelScan.getInfile(args);
	    	// If no input, print usage //

	    	if(in == null)
	    	{
	    		System.err.println(usage);
				return;
	    	}


	    	HashMap<String, String> pedInfo = null;

	    	// Parse out sample information //
	    	HashMap<String, Boolean> controlSamples = new HashMap();
	    	HashMap<String, Boolean> maleSamples = new HashMap();

	    	if(params.containsKey("ped-file"))
	    	{
		    	System.err.println("Loading sample information from " + pedFile + "...");
		    	try {
		    		int numCases = 0;
		    		int numControls = 0;
		    		int numMales = 0;
			    	pedInfo = MendelScan.loadPED(pedFile);

		    		for (String sample : pedInfo.keySet())
		    		{
		    			String pedLine = pedInfo.get(sample);
		    			String[] pedLineContents = pedLine.split("\t");
		    			// Family, PaternalID, MaternalID, Sex, Status //
		    			String familyID = pedLineContents[0];
		    			String sex = pedLineContents[3];
		    			String status = pedLineContents[4];

		    			if(sex.equals("1"))
		    			{
		    				maleSamples.put(sample, true);
		    				numMales++;
		    			}

		    			if(status.equals("1"))
		    			{
		    				numControls++;
		    				controlSamples.put(sample, true);
		    			}
		    			else if(status.equals("2"))
		    			{
		    				numCases++;
		    			}
		    		}

		    		System.err.println(numMales + " males, " + numCases + " cases, " + numControls + " controls");
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Warning: Exception thrown while loading pedigree information: " + e.getMessage());
		    		System.err.println(e.getLocalizedMessage());
		    		e.printStackTrace(System.err);
		    	}
	    	}
	    	else
	    	{
	    		System.err.println("WARNING: No pedigree file provided, so all samples assumed to be female affected!");
	    	}



	    	HashMap<String, Double> geneRank = new HashMap();

	    	boolean geneRankingProvided = false;
	    	if(params.containsKey("gene-file"))
	    	{
	    		geneRankingProvided = true;
		    	System.err.println("Loading gene expression information from " + geneFile + "...");
		    	try {
			    	geneRank = MendelScan.loadGeneExpression(geneFile);
			    	System.err.println("Expression rank loaded for " + geneRank.size() + " genes");
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Warning: Exception thrown while loading gene expression information: " + e.getMessage());
		    		System.err.println(e.getLocalizedMessage());
		    		e.printStackTrace(System.err);
		    	}
	    	}
	    	else
	    	{
	    		System.err.println("No gene expression file provided, so all expression scores will be 1.00");
	    	}


	    	HashMap<String, String> vepAnnot = null;
	    	if(params.containsKey("vep-file"))
	    	{
		    	System.err.println("Loading VEP from " + vepFile + "...");
		    	try {

		    	}
		    	catch (Exception e)
		    	{
		    		System.err.println("Exception thrown while loading VEP: " + e.getMessage() + " : " + e.getLocalizedMessage());
		    		e.printStackTrace(System.err);
		    	}

	    		vepAnnot = MendelScan.loadVEP(vepFile);
	    		System.err.println(vepAnnot.size() + " variants had VEP annotation");
	    	}
	    	else
	    	{
	    		System.err.println("ERROR: You must provide a VEP annotation file!");
	    		return;
	    	}

	    	System.err.println("Scoring variants under " + inheritanceModel + " disease model");



	    	// Declare file-parsing variables //

	    	String line;

	    	long numSamples = 0;
	    	long numSamplesMale = 0;
	    	long numSamplesCase = 0;
	    	long numSamplesControl = 0;
	    	long numVariants = 0;
	    	long numVariantsAnnot = 0;
	    	TreeMap<String, Integer> stats = new TreeMap();

			// FOrmat scores for printing //
			DecimalFormat dfScore = new DecimalFormat("0.000");
			DecimalFormat dfScoreLong = new DecimalFormat("0.00000");

	    	// Save the sample names by column number //

	    	HashMap<Integer, String> samplesByColumn = new HashMap();


	    	// Try to proceed with the input file //

	    	if(in != null && in.ready())
	    	{
	    		// Open the output files //
	    		PrintStream outFileHandle = null;
	    		PrintStream errFileHandle = null;
	    		PrintStream outVCFHandle = null;
	    		if(params.containsKey("output-file"))
	    		{
	    			outFileHandle = new PrintStream( new FileOutputStream(outFile) );
	    			errFileHandle = new PrintStream( new FileOutputStream(outFile + ".excluded") );
	    		}

	    		if(params.containsKey("output-vcf"))
	    		{
	    			outVCFHandle = new PrintStream( new FileOutputStream(outVCF) );
	    		}

	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{
	    			String[] lineContents = line.split("\t");

	    			if(line.startsWith("#"))
	    			{
	    				// Print it to output file if necessary //
	    				if(params.containsKey("output-vcf"))
	    				{
	    					// Just prior to chromosome/sample line, append MendelScan information lines
	    					if(line.startsWith("#CHROM"))
	    					{
	    						outVCFHandle.println("##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"Overall score calculated from product of the five MendelScan scores\">");
	    						outVCFHandle.println("##INFO=<ID=SEGSCORE,Number=1,Type=Float,Description=\"Segregation score based on model of inheritance\">");
	    						outVCFHandle.println("##INFO=<ID=POPSCORE,Number=1,Type=Float,Description=\"Population score based on rareness of variant\">");
	    						outVCFHandle.println("##INFO=<ID=ANNOSCORE,Number=1,Type=Float,Description=\"Annotation score based on predicted variant consequence\">");
	    						outVCFHandle.println("##INFO=<ID=EXPRSCORE,Number=1,Type=Float,Description=\"Expression score reflecting the relative expression of gene\">");
	    						outVCFHandle.println("##INFO=<ID=VARGENE,Number=1,Type=String,Description=\"The HUGO gene symbol of affected gene\">");
	    						outVCFHandle.println("##INFO=<ID=VARCLASS,Number=1,Type=String,Description=\"The predicted consequence of the variant\">");

	    					}
	    					outVCFHandle.println(line);
	    				}
	    				// Parse the VCF header //

	    				if(line.startsWith("#CHROM"))
	    				{
	    					// Check for incorrect VCF header format
	    					if(lineContents.length < 9)
	    					{
	    						System.err.println("Error: Insufficient number of columns in VCF header: " + line);
	    						in.close();
	    						return;
	    					}
	    					else
	    					{
	    						String fileHeader = "CHROM\tPOSITION\tREF\tALT\tOVERALL_SCORE\tSEGREGATION_SCORE\tPOPULATION_SCORE\tANNOTATION_SCORE\tEXPRESSION_SCORE\tDBSNP_STATUS\tDBSNP_ID\tDBSNP_INFO\tCASES\tCASES_REF\tCASES_HET\tCASES_HOM\tCASES_MISSING\tCONTROLS\tCONTROLS_REF\tCONTROLS_HET\tCONTROLS_HOM\tCONTROLS_MISSING\tHUGO_GENE\tENSEMBL_GENE\tCANONICAL\tCLASS\tTX_POS\tAA_POS\tAA_CHANGE\tPOLYPHEN\tSIFT\tCONDEL";

	    						for(int colCounter = 9; colCounter < lineContents.length; colCounter++)
	    						{
	    							numSamples++;
	    							String sample = lineContents[colCounter];

	    							String sampleLogLine = sample;

	    							String sampleGender = "";
	    							String sampleStatus = "";

	    							// Save sample name //
	    							samplesByColumn.put(colCounter, sample);

	    							if(controlSamples.containsKey(sample))
	    							{
	    								sampleLogLine += "\tControl";
	    								numSamplesControl++;
	    								sampleStatus = "0";
	    							}
	    							else
	    							{
	    								sampleLogLine += "\tCase";
	    								numSamplesCase++;
	    								sampleStatus = "1";
	    							}

	    							if(maleSamples.containsKey(sample))
	    							{
	    								sampleLogLine += "\tMale";
	    								numSamplesMale++;
	    								sampleGender = "M";
	    							}
	    							else
	    							{
	    								sampleLogLine += "\tFemale";
	    								sampleGender = "F";
	    							}

	    							fileHeader += "\t" + sample + ":" + sampleGender + ":" + sampleStatus;

	    							System.err.println("SAMPLE:\t" + sampleLogLine);
	    						}

	    						// Print file header //
	    						if(!outFileHandle.equals(null))
	    						{
	    							outFileHandle.println(fileHeader);
	    						}
	    					}

	    				}
	    			}

	    			// PARSE VARIANTS FROM BODY OF VCF //
	    			else
	    			{
    					// Check for incorrect VCF header format
    					if(lineContents.length < 9)
    					{
    						System.err.println("Error: Insufficient number of columns in VCF header: " + line);
    						in.close();
    						return;
    					}
    					else
    					{

    						try {
    	    					// Parse the VCF main file //
    	    					String chrom = lineContents[0];
    	    					String position = lineContents[1];
    	    					String id = lineContents[2];
    	    					String ref = lineContents[3];
    	    					String altAlleles = lineContents[4];
    	    					String qual = lineContents[5];
    	    					String filter = lineContents[6];
    	    					String info = lineContents[7];
    	    					String format = lineContents[8];


    	    					if(filter.equals("PASS") || filter.equals("."))
    	    					{
    	    						numVariants++;

    		    					String[] alts = altAlleles.split(",");
    		    					String vepKey = chrom + "_" + position + "_" + ref + "/" + altAlleles.replace(",", "/");
      		    					// Go through VCF one alternative allele at a time, finding a VEP key that works //
    		    					int altNumber = 0;

    		    					for(int altCounter = 0; altCounter < alts.length; altCounter++)
    		    					{
    		    						String alt = alts[altCounter];
    		    						altNumber = altCounter + 1;

    		    						if(alt.length() < ref.length())
    		    						{
    		    							// Deletion //
    		    							String thisRef = ref.replaceFirst(alt, "");
    		    							String thisAlt = "-";
    		    							int thisPos = Integer.parseInt(position) + 1;
    		    							vepKey = chrom + "_" + thisPos + "_" + thisRef + "/" + thisAlt;
    		    						}
    		    						else if(alt.length() > ref.length())
    		    						{
    		    							// Insertion //
    		    							String thisAlt = alt.replaceFirst(ref, "");
    		    							String thisRef = "-";
    		    							int thisPos = Integer.parseInt(position) + 1;
    		    							vepKey = chrom + "_" + thisPos + "_" + thisRef + "/" + thisAlt;
    		    						}
    		    						else
    		    						{
    		    							// SNV //
    		    							vepKey = chrom + "_" + position + "_" + ref + "/" + alt;
            		    					String vcfKey1 = chrom + "\t" + position + "\t" + id + "\t" + ref + "\t" + altAlleles;
            		    					String vcfKey2 = chrom + "\t" + position + "\t" + "." + "\t" + ref + "\t" + altAlleles;
            		    					if(vepAnnot.containsKey(id))
            		    					{
            		    						// If RS ID provided, correct the issue //
            		    						vepKey = id;
            		    					}
            		    					else if(vepAnnot.containsKey(vcfKey1))
            		    					{
            		    						// If We parsed a VCF formatted VEP file //
            		    						vepKey = vcfKey1;
            		    					}
            		    					else if(vepAnnot.containsKey(vcfKey2))
            		    					{
            		    						// If We parsed a VCF formatted VEP file //
            		    						vepKey = vcfKey2;
            		    					}
            		    					else
            		    					{
            		    						if(vepAnnot.containsKey(vepKey))
            		    						{
            		    							// We got the annotation, so good to go. //
            		    						}
            		    						else if(params.containsKey("verbose"))
            		    						{
            		    							System.err.println("Got nothing for " + vepKey);
                		    						System.err.println("Got nothing for " + vcfKey1);
                		    						System.err.println("Got nothing for " + vcfKey2);
            		    						}

            		    					}
    		    						}


    		    						// Now that we have VEP key, proceed for this alt only //

        		    					HashMap<String, String> genotypesBySample = new HashMap();
        		    					String sampleGenotypes = "";

        		    					for(int colCounter = 9; colCounter < lineContents.length; colCounter++)
        		    					{
        		    						// Get sample name, status, and gender //
        		    						if(samplesByColumn.containsKey(colCounter))
        		    						{
        		    							String sampleName = samplesByColumn.get(colCounter);
        		    							String sampleStatus = "case";
        		    							if(controlSamples.containsKey(sampleName))
        		    								sampleStatus = "control";
        		    							String sampleGender = "female";
        		    							if(maleSamples.containsKey(sampleName))
        		    								sampleGender = "male";

        		    							HashMap<String, String> sampleGenotype = MendelScan.parseGenotypeField(format, lineContents[colCounter]);

        		    							// Parse out sample genotype, depth, ref reads, var reads //
        		    							String sampleGT = ".";
        		    							if(sampleGenotype.containsKey("GT"))
        		    								sampleGT = sampleGenotype.get("GT");
        		    							String sampleDP = ".";
        		    							if(sampleGenotype.containsKey("DP"))
        		    								sampleDP = sampleGenotype.get("DP");
        		    							String sampleAD = ".";
        		    							if(sampleGenotype.containsKey("AD"))
        		    								sampleAD = sampleGenotype.get("AD");

        		    							String sampleLine = sampleStatus + "\t" + sampleGender + "\t";
        		    							sampleLine += sampleGT + "\t" + sampleDP + "\t" + sampleAD;
        		    							genotypesBySample.put(sampleName, sampleLine);

        		    							if(sampleGenotypes.length() > 0)
        		    							{
        		    								sampleGenotypes += "\t";
        		    							}

        		    							sampleGenotypes += lineContents[colCounter];

        		    						}

        		    					}

        		    					String segStatus = MendelScan.getSegregationStatus(genotypesBySample, minDepth);
        		    					double segScore = getSegregationScore(chrom, inheritanceModel, segStatus, params);

//        		    					System.out.println(vepKey + "\t" + segStatus + "\t" + segScore);

        		    					// Parse out values from info field, which should include dbSNP values //

        		    					HashMap<String, String> infoValues = MendelScan.parseInfoField(info);

        		    					// Save dbSNP ID if necessary //
        		    					if(!id.equals("."))
        		    					{
        		    						infoValues.put("RSnumber", id);
        		    					}

        		    					// Determine dbSNP status of variant //

        		    					String dbsnpStatus = MendelScan.getDbsnpStatus(infoValues);

        		    					// Assign pop score for dbSNP status //

        		    					double popScore = getPopulationScore(params, dbsnpStatus);


        	    						if(vepAnnot.containsKey(vepKey))
        	    						{
        	    							numVariantsAnnot++;

        	    							// Get the top annotation per HUGO gene //
        	    							HashMap<String, String> topByGene = topAnnotationByGene(vepAnnot.get(vepKey));

        	    							// Go through each HUGO gene, considering its top annotation //

        	    							for (String hugoGene : topByGene.keySet())
        	    							{
        	    								try {
        		    								String topAnnot = topByGene.get(hugoGene);
        		    								String[] topAnnotContents = topAnnot.split("\t");
        		    								// Note, the topAnnotScore is the integer rank of the VEP annotation //
        		    								// The annotScore double is the variant score between 0 and 1 //
        		    								int topAnnotScore = Integer.parseInt(topAnnotContents[10]);
        		    								double annotScore = getAnnotationScore(params, topAnnotScore);

        		    								// Determine the expression score, default 0.50 but otherwise rank of gene //
        		    								double expressionScore = 0.50;
        		    								if(geneRank.containsKey(hugoGene))
        		    									expressionScore = geneRank.get(hugoGene);
        		    								else if(!geneRankingProvided)
        		    									expressionScore = 1.00;		// If no gene rank provided, set all to 1. //

        		    								// Calculate the overall score //
        		    								double overallScore = segScore * popScore * annotScore * expressionScore;

        		    								// Output line //
//        		    								String outLine = chrom + "\t" + position + "\t" + ref + "\t" + altAlleles + "\t";
        		    								String outLine = chrom + "\t" + position + "\t" + ref + "\t" + alt + "\t";
        		    								// Output scores //
        		    								outLine += dfScoreLong.format(overallScore) + "\t" + dfScore.format(segScore) + "\t" + dfScore.format(popScore) + "\t" + dfScore.format(annotScore) + "\t" + dfScore.format(expressionScore) + "\t";
        		    								// Output pop info //
        		    								outLine += dbsnpStatus + "\t" + id + "\t" + info + "\t";
        		    								// Output segregation info //
        		    								outLine += segStatus + "\t";
        		    								// Append Annotation //

        		    								String ensGene = topAnnotContents[0];
        		    								String varClass = topAnnotContents[2];
        		    								String canonical = topAnnotContents[3];
        		    								String polyphen = topAnnotContents[4];
        		    								String sift = topAnnotContents[5];
        		    								String condel = topAnnotContents[6];
        		    								String txPos = topAnnotContents[7];
        		    								String aaPos = topAnnotContents[8];
        		    								String aaChange = topAnnotContents[9];

        		    								if(canonical.equals("CANONICAL"))
        		    									canonical = "YES";

        		    								outLine += hugoGene + "\t" + ensGene + "\t" + canonical + "\t";
        		    								outLine += varClass + "\t" + txPos + "\t" + aaPos + "\t" + aaChange + "\t";
        		    								outLine += polyphen + "\t" + sift + "\t" + condel + "\t";
        		    								outLine += sampleGenotypes;


        		    		    					if(stats.containsKey("variants_" + dbsnpStatus))
        		    		    					{
        		    		    						stats.put("variants_" + dbsnpStatus, stats.get("variants_" + dbsnpStatus) + 1);
        		    		    					}
        		    		    					else
        		    		    					{
        		    		    						stats.put("variants_" + dbsnpStatus, 1);
        		    		    					}

        		    								if(params.containsKey("output-file") && !outFileHandle.equals(null))
        		    									outFileHandle.println(outLine);

        		    								// Build new VCF line //

        		    								// Build new info field //
        		    								String newInfo = info;
        		    								if(newInfo.length() > 0)
        		    								{
        		    									if(newInfo.equals("."))
        		    									{
        		    										newInfo = "";
        		    									}
        		    									else
        		    									{
        		    										newInfo += ";";
        		    									}
        		    								}
        		    								newInfo += "SCORE=" + dfScore.format(overallScore) + ";";
        		    								newInfo += "SEGSCORE=" + dfScore.format(segScore) + ";";
        		    								newInfo += "POPSCORE=" + dfScore.format(popScore) + ";";
        		    								newInfo += "ANNOSCORE=" + dfScore.format(annotScore) + ";";
        		    								newInfo += "EXPRSCORE=" + dfScore.format(expressionScore) + ";";
        		    								newInfo += "VARGENE=" + hugoGene + ";VARCLASS=" + varClass;

        		    								String newVCFline = chrom + "\t" + position + "\t" + id + "\t" + ref + "\t" + altAlleles + "\t";
        		    		    					newVCFline += qual + "\t" + filter + "\t" + newInfo + "\t" + format;
        		    		    					for(int colCounter = 9; colCounter < lineContents.length; colCounter++)
        		    		    					{
        		    		    						newVCFline += "\t" + lineContents[colCounter];
        		    		    					}

        		    		    					if(params.containsKey("output-vcf"))
        		    		    						outVCFHandle.println(newVCFline);
        	    								}
        	    								catch(Exception e)
        	    								{
        	    									System.err.println("Warning: Exception thrown while processing " + hugoGene + " annotation for " + vepKey + " : " + e.getMessage());
        	    									e.printStackTrace(System.err);
        	    								}


        	    							} // End foreach hugoGene in topByGene loop


        	    						}
        	    						else
        	    						{
        	    							if(!errFileHandle.equals(null))
        	    							{
        	    								errFileHandle.println("NO_VEP_ANNOTATION" + "\t" + line);
        	    							}

        	    							if(params.containsKey("verbose"))
        	    								System.err.println("Warning: No VEP info for " + vepKey);
        	    						}

    		    					}

    	    					}

    						}
    						catch(Exception e)
    						{
    							System.err.println("Parsing exception thrown for line: " + line);
    							System.err.println(e.getMessage());
    						}



    					}
	    			}

	    		}

	    		in.close();
	    	}
	    	// Insufficient input was provided, so print usage //
	    	else
	    	{
				 System.err.println("Please provide an input file!\n" + usage);
				 System.exit(10);
	    	}


	    	System.err.println(numSamples + " samples in VCF (" + numSamplesCase + " affected, " + numSamplesControl + " unaffected, " + numSamplesMale + " male)");

	    	System.err.println(numVariants + " variants in VCF file");
	    	System.err.println(numVariantsAnnot + " matched with VEP annotation");

	    	for (String statKey : stats.keySet())
	    	{
	    		System.err.println(stats.get(statKey) + "\t" + statKey);
	    	}

	    }
	    catch(Exception e)
	    {
	    	System.err.println("Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(11);
	    }

	}	// # End PrioritizeVCF main subroutine






	/**
	 * Returns the segregation score based on segregation patterns and expected mode of inheritance
	 *
	 * @param	params	Input parameters, for determining user settings
	 * @param	status	Type of file ("positions" or "regions")
	 * @return			Double of resulting population score for this variant
	 */
	static double getSegregationScore(String chrom, String inheritanceMode, String segStatus, HashMap<String, String> params)
	{
		double segScore = 0.000;

		// Set default scores for segregation patterns //

		double scoreCaseRef = 1.00;
		double scoreCaseHet = 1.00;
		double scoreCaseHom = 1.00;
		double scoreCaseMissing = 1.00;
		double scoreControlHet = 1.00;
		double scoreControlHom = 1.00;
		double scoreControlMissing = 1.00;

		if(inheritanceMode.equals("recessive"))
		{
			// Recessive inheritance assumptions //
			scoreCaseRef = 0.10;
			scoreCaseHet = 0.50;	// For non-compound hets, let's penalize a het case //
			scoreControlHet = 1.00;	// No penalty as control could be carrier
			scoreControlHom = 0.50;
			scoreCaseMissing = 0.60;
			scoreControlMissing = 0.80;
		}
		else
		{
			// Dominant and X-linked inheritance assumptions //
			scoreCaseRef = 0.20;
			scoreCaseHom = 0.80;
			scoreCaseMissing = 0.80;
			scoreControlMissing = 0.50;
			scoreCaseHet = 1.00;
			scoreControlHet = 0.10;
			scoreControlHom = 0.01;

		}


		// Try to get the user's parameter changes //

		try {
			if(params.containsKey("seg-score-case-ref"))
				scoreCaseRef = Double.parseDouble(params.get("seg-score-case-ref"));

			if(params.containsKey("seg-score-case-het"))
				scoreCaseHet = Double.parseDouble(params.get("seg-score-case-het"));

			if(params.containsKey("seg-score-case-hom"))
				scoreCaseHom = Double.parseDouble(params.get("seg-score-case-hom"));

			if(params.containsKey("seg-score-control-het"))
				scoreControlHet = Double.parseDouble(params.get("seg-score-control-het"));

			if(params.containsKey("seg-score-control-hom"))
				scoreControlHom = Double.parseDouble(params.get("seg-score-control-hom"));

		}
		catch(Exception e)
		{
			System.err.println("Warning: Exception thrown while parsing segregation score parameters: " + e.getMessage());
		}

		// We expect all affected to be heterozygous, all controls to be reference
		String[] segContents = segStatus.split("\t");
		if(segContents.length>= 8)
		{
			// Start assuming perfect segregation //
			segScore = 1.000;
			int casesCalled = Integer.parseInt(segContents[0]);
			int casesRef = Integer.parseInt(segContents[1]);
			int casesHet = Integer.parseInt(segContents[2]);
			int casesHom = Integer.parseInt(segContents[3]);
			int casesMissing = Integer.parseInt(segContents[4]);
			int controlsCalled = Integer.parseInt(segContents[5]);
			int controlsRef = Integer.parseInt(segContents[6]);
			int controlsHet = Integer.parseInt(segContents[7]);
			int controlsHom = Integer.parseInt(segContents[8]);
			int controlsMissing = Integer.parseInt(segContents[9]);
			int controlsVariant = controlsHet + controlsHom;

			if(inheritanceMode.equals("dominant"))
			{
				// Assume 50% sensitivity to detect heterozygotes; reduce score for affecteds called Ref //

				if(casesRef > 0)
				{
					for(int caseCounter = 0; caseCounter < casesRef; caseCounter++)
					{
						segScore = segScore * scoreCaseRef;
					}
				}

				// Reduce score for affecteds homozygous for a rare dominant-acting disease //

				if(casesHom > 0 && !chrom.equals("X")&& !chrom.equals("chrX") && !chrom.equals("Y")&& !chrom.equals("chrY"))
				{
					for(int caseCounter = 0; caseCounter < casesHom; caseCounter++)
					{
						segScore = segScore * scoreCaseHom;
					}
				}

				if(controlsHet > 0)
				{
					for(int controlCounter = 0; controlCounter < controlsHet; controlCounter++)
					{
						segScore = segScore * scoreControlHet;
					}
				}

				if(controlsHom > 0)
				{
					for(int controlCounter = 0; controlCounter < controlsHom; controlCounter++)
					{
						segScore = segScore * scoreControlHom;
					}
				}

			}
			else if(inheritanceMode.equals("recessive"))
			{
				// Affecteds should be homozygous-variant; unaffecteds could be ref or het //
				// Assume 50% sensitivity to detect heterozygotes; reduce score dramatically for affecteds called Ref //

				if(casesRef > 0)
				{
					for(int caseCounter = 0; caseCounter < casesRef; caseCounter++)
					{
						segScore = segScore * scoreCaseRef;
					}
				}

				// Cases should be homozygous, so penalize if heterozygous. Compound hets done separately //

				if(casesHet > 0)
				{
					// Assume that cases called het are mis-called //
					for(int caseCounter = 0; caseCounter < casesHet; caseCounter++)
					{
						segScore = segScore * scoreCaseHet;
					}
				}

				// Controls should not be homozygous-variant; mis-calling is possible

				if(controlsHom > 0)
				{
					for(int controlCounter = 0; controlCounter < controlsHom; controlCounter++)
					{
						segScore = segScore * scoreControlHom;
					}
				}

			}
			else
			{
				System.err.println("Unrecognized inheritance model: " + inheritanceMode);
				System.exit(0);
			}


			// Regardless of inheritance, apply penalties for missing data //
			if(controlsMissing > 0)
			{
				for(int controlCounter =0; controlCounter < controlsMissing; controlCounter++)
				{
					segScore = segScore * scoreControlMissing;
				}
			}

			if(casesMissing > 0)
			{
				for(int caseCounter =0; caseCounter < casesMissing; caseCounter++)
				{
					segScore = segScore * scoreCaseMissing;
				}
			}

		}




		return(segScore);
	}

	/**
	 * Returns the population score based on dbSNP information and user-specified thresholds
	 *
	 * @param	info	HashMap of dbSNP information from ID and INFO columns
	 * @return			String with one of several possible dbSNP statuses
	 */
	static HashMap topAnnotationByGene(String annot)
	{
		HashMap<String, String> annotByGene = new HashMap<String, String>();
		HashMap<String, String> isCanonicalByGene = new HashMap<String, String>();
		HashMap<String, Integer> topScoreByGene = new HashMap<String, Integer>();

		try {
			String[] annotLines = annot.split("\n");

			for(int lineCounter = 0; lineCounter < annotLines.length; lineCounter++)
			{
				String annotLine = annotLines[lineCounter];
				String[] annotContents = annotLine.split("\t");
				String hugoGene = annotContents[1];
				String isCanonical = annotContents[3];
				int annotScore = Integer.parseInt(annotContents[10]);

				if(annotByGene.containsKey(hugoGene))
				{
					// This is the new top gene if it has same canonical status but higher score //
					if(isCanonical.equals(isCanonicalByGene.get(hugoGene)) || isCanonical.equals("YES"))
					{
						if(annotScore > topScoreByGene.get(hugoGene))
						{
							// Replace with this annot //
							annotByGene.put(hugoGene, annotLine);
							topScoreByGene.put(hugoGene, annotScore);
							isCanonicalByGene.put(hugoGene, isCanonical);
						}
					}

				}
				else
				{
					annotByGene.put(hugoGene, annotLine);
					topScoreByGene.put(hugoGene, annotScore);
					isCanonicalByGene.put(hugoGene, isCanonical);
				}
			}

		}
		catch (Exception e)
		{
			System.err.println("Warning: Exception thrown while determining top annotation: " + e.getMessage());
		}


		return(annotByGene);
	}

	/**
	 * Returns the annotation score based on VEP class user-specified thresholds
	 *
	 * @param	params	Input parameters, for determining user settings
	 * @param	status	Type of file ("positions" or "regions")
	 * @return			Double of resulting population score for this variant
	 */
	static double getAnnotationScore(HashMap<String, String> params, int vepScore)
	{
		double annotScore = 0.01;

		double annoScore1 = 0.01;	// Intergenic
		double annoScore2 = 0.01;	// Intronic
		double annoScore3 = 0.01;	// Downstream
		double annoScore4 = 0.01;	// Upstream
		double annoScore5 = 0.01;	// UTR
		double annoScore6 = 0.01;	// Noncoding gene
		double annoScore7 = 0.01;	// Mature miRNA
		double annoScore8 = 0.05;	// Synonymous or coding-unknown or partial-codon
		double annoScore9 = 0.20;	// Splice site
		double annoScore10 = 1.00;	// Nonstop
		double annoScore11 = 0.80;	// Missense + 0 damaging
		double annoScore12 = 0.95;	// Missense + 1 damaging
		double annoScore13 = 0.95;	// Missense + 2 damaging
		double annoScore14 = 0.95;	// Missense + 3 damaging
		double annoScore15 = 1.00;	// Essential splice site
		double annoScore16 = 1.00;	// Frameshift
		double annoScore17 = 1.00;	// Nonsense

		try {
			if(params.containsKey("anno-score-1"))
				annoScore1 = Double.parseDouble(params.get("anno-score-1"));
			if(params.containsKey("anno-score-2"))
				annoScore2 = Double.parseDouble(params.get("anno-score-2"));
			if(params.containsKey("anno-score-3"))
				annoScore3 = Double.parseDouble(params.get("anno-score-3"));
			if(params.containsKey("anno-score-4"))
				annoScore4 = Double.parseDouble(params.get("anno-score-4"));
			if(params.containsKey("anno-score-5"))
				annoScore5 = Double.parseDouble(params.get("anno-score-5"));
			if(params.containsKey("anno-score-6"))
				annoScore6 = Double.parseDouble(params.get("anno-score-6"));
			if(params.containsKey("anno-score-7"))
				annoScore7 = Double.parseDouble(params.get("anno-score-7"));
			if(params.containsKey("anno-score-8"))
				annoScore8 = Double.parseDouble(params.get("anno-score-8"));
			if(params.containsKey("anno-score-9"))
				annoScore9 = Double.parseDouble(params.get("anno-score-9"));
			if(params.containsKey("anno-score-10"))
				annoScore10 = Double.parseDouble(params.get("anno-score-10"));
			if(params.containsKey("anno-score-11"))
				annoScore11 = Double.parseDouble(params.get("anno-score-11"));
			if(params.containsKey("anno-score-12"))
				annoScore12 = Double.parseDouble(params.get("anno-score-12"));
			if(params.containsKey("anno-score-13"))
				annoScore13 = Double.parseDouble(params.get("anno-score-13"));
			if(params.containsKey("anno-score-14"))
				annoScore14 = Double.parseDouble(params.get("anno-score-14"));
			if(params.containsKey("anno-score-15"))
				annoScore15 = Double.parseDouble(params.get("anno-score-15"));
			if(params.containsKey("anno-score-16"))
				annoScore16 = Double.parseDouble(params.get("anno-score-16"));
			if(params.containsKey("anno-score-17"))
				annoScore17 = Double.parseDouble(params.get("anno-score-17"));

		}
		catch(Exception e)
		{
			System.err.println("Warning: Exception thrown while parsing annotation score parameters: " + e.getMessage());
		}

		try {
			if(vepScore == 1)
				annotScore = annoScore1;
			else if(vepScore == 2)
				annotScore = annoScore2;
			else if(vepScore == 3)
				annotScore = annoScore3;
			else if(vepScore == 4)
				annotScore = annoScore4;
			else if(vepScore == 5)
				annotScore = annoScore5;
			else if(vepScore == 6)
				annotScore = annoScore6;
			else if(vepScore == 7)
				annotScore = annoScore7;
			else if(vepScore == 8)
				annotScore = annoScore8;
			else if(vepScore == 9)
				annotScore = annoScore9;
			else if(vepScore == 10)
				annotScore = annoScore10;
			else if(vepScore == 11)
				annotScore = annoScore11;
			else if(vepScore == 12)
				annotScore = annoScore12;
			else if(vepScore == 13)
				annotScore = annoScore13;
			else if(vepScore == 14)
				annotScore = annoScore14;
			else if(vepScore == 15)
				annotScore = annoScore15;
			else if(vepScore == 16)
				annotScore = annoScore16;
			else if(vepScore == 17)
				annotScore = annoScore17;

		}
		catch(Exception e)
		{
			System.err.println("Warning: issue when assigning annotation score: " + e.getMessage());
		}

		return(annotScore);
	}


	/**
	 * Returns the population score based on dbSNP information and user-specified thresholds
	 *
	 * @param	params	Input parameters, for determining user settings
	 * @param	status	Type of file ("positions" or "regions")
	 * @return			Double of resulting population score for this variant
	 */
	static double getPopulationScore(HashMap<String, String> params, String status)
	{
		double popScore = 0.00;

		// Establish default parameters //
		double scoreNovel = 1.00;
		double scoreMutation = 0.95;
		double scoreKnown = 0.60;
		double scoreRare = 0.20;
		double scoreUncommon = 0.02;
		double scoreCommon = 0.01;

		try {
			if(params.containsKey("pop-score-novel"))
				scoreNovel = Double.parseDouble(params.get("pop-score-novel"));
			if(params.containsKey("pop-score-mutation"))
				scoreMutation = Double.parseDouble(params.get("pop-score-mutation"));
			if(params.containsKey("pop-score-known"))
				scoreKnown = Double.parseDouble(params.get("pop-score-known"));
			if(params.containsKey("pop-score-rare"))
				scoreRare = Double.parseDouble(params.get("pop-score-rare"));
			if(params.containsKey("pop-score-uncommon"))
				scoreUncommon = Double.parseDouble(params.get("pop-score-uncommon"));
			if(params.containsKey("pop-score-common"))
				scoreCommon = Double.parseDouble(params.get("pop-score-common"));

		}
		catch(Exception e)
		{
			System.err.println("Warning: Exception thrown while parsing population score parameters: " + e.getMessage());
		}


		try {
			if(status.equals("novel"))
			{
				popScore = scoreNovel;
			}
			else if(status.equals("mutation"))
			{
				popScore = scoreMutation;
			}
			else if(status.equals("known"))
			{
				popScore = scoreKnown;
			}
			else if(status.equals("rare"))
			{
				popScore = scoreRare;
			}
			else if(status.equals("uncommon"))
			{
				popScore = scoreUncommon;
			}
			else if(status.equals("common"))
			{
				popScore = scoreCommon;
			}
		}
		catch(Exception e)
		{
			System.err.println("Warning: Exception thrown while calculating population score: " + e.getMessage());
		}

		return(popScore);
	}

} // End Class Definition

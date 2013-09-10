/**
 * @(#)CallMpileup.java
 *
 * Copyright (c) 2009-2010 Daniel C. Koboldt and Washington University in St. Louis
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
		String usage = "USAGE: java -jar MendelScan.jar prioritize [VCF] OPTIONS\n" +
		"\tOPTIONS:\n" +
		"\t--vep-file\tVariant annotation in VEP format\n" +
		"\t--output-file\tOutput file to contain human-friendly scored variants\n" +
		"\t--output-vcf\tOutput file to contain scored variants in VCF format\n" +
		"\t--expression-file\tA list of gene expression values for tissue of interest\n" +
		"\t--inheritance\tPresumed model of inheritance: dominant, recessive, x-linked [dominant]\n";

		String vepFile = null;
		String pedFile = null;
		String geneFile = null;
		String outFile = null;
		String outVCF = null;
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
		    				maleSamples.put(sample, true);

		    			if(status.equals("1"))
		    				controlSamples.put(sample, true);
		    		}
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Warning: Exception thrown while loading pedigree information: " + e.getMessage());
		    		System.err.println(e.getLocalizedMessage());
		    		e.printStackTrace(System.err);
		    	}
	    	}



	    	HashMap<String, Double> geneRank = new HashMap();

	    	if(params.containsKey("gene-file"))
	    	{
		    	System.err.println("Loading gene expression information from " + geneFile + "...");
		    	try {
			    	geneRank = MendelScan.loadGeneExpression(geneFile);
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Warning: Exception thrown while loading gene expression information: " + e.getMessage());
		    		System.err.println(e.getLocalizedMessage());
		    		e.printStackTrace(System.err);
		    	}
	    	}


	    	HashMap<String, String> vepAnnot = null;
	    	if(params.containsKey("vep-file"))
	    	{
		    	System.err.println("Loading VEP from " + vepFile + "...");
	    		vepAnnot = MendelScan.loadVEP(vepFile);
	    		System.err.println("Loaded VEP for " + vepAnnot.size() + " variants");
	    	}


	    	// Declare file-parsing variables //

	    	String line;

	    	long numSamples = 0;
	    	long numSamplesMale = 0;
	    	long numSamplesCase = 0;
	    	long numSamplesControl = 0;
	    	long numVariants = 0;
	    	long numVariantsAnnot = 0;
	    	HashMap<String, Integer> stats = new HashMap();

			// FOrmat scores for printing //
			DecimalFormat dfScore = new DecimalFormat("#.000");


	    	// Save the sample names by column number //

	    	HashMap<Integer, String> samplesByColumn = new HashMap();


	    	// Try to proceed with the input file //

	    	if(in != null && in.ready())
	    	{
	    		// Open the output files //
	    		PrintStream outFileHandle = null;
	    		PrintStream outVCFHandle = null;
	    		if(params.containsKey("output-file"))
	    		{
	    			outFileHandle = new PrintStream( new FileOutputStream(outFile) );
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
	    					outVCFHandle.println(line);

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
	    						for(int colCounter = 9; colCounter < lineContents.length; colCounter++)
	    						{
	    							numSamples++;
	    							String sample = lineContents[colCounter];

	    							// Save sample name //
	    							samplesByColumn.put(colCounter, sample);

	    							if(controlSamples.containsKey(sample))
	    								numSamplesControl++;
	    							else
	    								numSamplesCase++;

	    							if(maleSamples.containsKey(sample))
	    								numSamplesMale++;
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

	    					// Parse the VCF main file //
	    					String chrom = lineContents[0];
	    					String position = lineContents[1];
	    					String id = lineContents[2];
	    					String ref = lineContents[3];
	    					String altAlleles = lineContents[4];
	    					String filter = lineContents[6];
	    					String info = lineContents[7];
	    					String format = lineContents[8];

	    					if(filter.equals("PASS") || filter.equals("."))
	    					{
	    						numVariants++;

		    					String[] alts = altAlleles.split(",");
		    					HashMap<String, String> genotypesBySample = new HashMap();

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

		    						}

		    					}

		    					String segStatus = getSegregationStatus(genotypesBySample, minDepth);
		    					double segScore = getSegregationScore(chrom, "dominant", segStatus);
		    					String vepKey = chrom + "_" + position + "_" + ref + "/" + altAlleles.replace(",", "/");

		    					System.out.println(vepKey + "\t" + segStatus + "\t" + segScore);

		    					// Parse out values from info field, which should include dbSNP values //

		    					HashMap<String, String> infoValues = MendelScan.parseInfoField(info);

		    					// Save dbSNP ID if necessary //
		    					if(!id.equals("."))
		    					{
		    						infoValues.put("RSnumber", id);
		    					}

		    					// Determine dbSNP status of variant //

		    					String dbsnpStatus = getDbsnpStatus(infoValues);

		    					// Assign pop score for dbSNP status //

		    					double popScore = getPopulationScore(params, dbsnpStatus);

		    					if(stats.containsKey(dbsnpStatus))
		    					{
		    						stats.put(dbsnpStatus, stats.get(dbsnpStatus) + 1);
		    					}
		    					else
		    					{
		    						stats.put(dbsnpStatus, 1);
		    					}


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

		    								// Calculate the overall score //
		    								double overallScore = segScore * popScore * annotScore * expressionScore;


		    								// Output line //
		    								String outLine = chrom + "\t" + position + "\t" + ref + "\t" + altAlleles + "\t";
		    								// Output scores //
		    								outLine += dfScore.format(overallScore) + "\t" + dfScore.format(segScore) + "\t" + dfScore.format(popScore) + "\t" + dfScore.format(annotScore) + "\t" + dfScore.format(expressionScore) + "\t";
		    								// Output pop info //
		    								outLine += dbsnpStatus + "\t" + id + "\t" + info + "\t";
		    								// Output segregation info //
		    								outLine += segStatus + "\t";
		    								// Append Annotation //
		    								outLine += hugoGene + "\t" + topAnnotContents[2] + "\t" + topAnnot.replace("\t", ";") + "\t";

		    								if(params.containsKey("output-file"))
		    									outFileHandle.println(outLine);
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
	    							System.err.println("Warning: No VEP info for " + vepKey);
	    						}


		    					// Go through VCF one alternative allele at a time //

		    					for(int altCounter = 0; altCounter < alts.length; altCounter++)
		    					{
		    						String alt = alts[altCounter];
		    					}
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
	    		System.err.println(stats.get(statKey) + "\thad VEP code " + statKey);
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
	 * Returns the population score based on dbSNP information and user-specified thresholds
	 *
	 * @param	info	HashMap of dbSNP information from ID and INFO columns
	 * @return			String with one of several possible dbSNP statuses
	 */
	static String getSegregationStatus(HashMap<String, String> genotypes, Integer minDepth)
	{
		String segStatus = "unknown";

		int casesCalled = 0;
		int casesRef = 0;
		int casesHet = 0;
		int casesHom = 0;
		int casesMissing = 0;
		int controlsCalled = 0;
		int controlsRef = 0;
		int controlsHet = 0;
		int controlsHom = 0;
		int controlsMissing = 0;

		try {
			for (String sample : genotypes.keySet())
			{
				try {
					String[] sampleContents = genotypes.get(sample).split("\t");
					String status = sampleContents[0];
					String gender = sampleContents[1];
					String gt = sampleContents[2];

					if(!gt.equals(".") && !gt.equals("./."))
					{
						// Obtain sequence depth and allele depth //
						Integer depth = 0;
						if(sampleContents[3].length() > 0 && !sampleContents[3].equals("."))
						{
							depth = Integer.parseInt(sampleContents[3]);
						}

						// Check to see if this variant matches the expectation //
						if(status.equals("case"))
						{
							casesCalled++;
							// CASE //
							if(MendelScan.isHeterozygous(gt))
							{
								casesHet++;
							}
							else
							{
								if(MendelScan.isHomozygous(gt))
								{
									casesHom++;
								}
								else if(MendelScan.isReference(gt))
								{
									// Determine if we have sufficient depth and VAF is less than 5% //
									if(depth >= minDepth)
									{
										if(sampleContents[4].length() > 0 && !sampleContents[4].equals("."))
										{
											Integer varDepth = 0;
											// Get all alternate allele depths //
											String[] readCounts = sampleContents[4].split(",");
											for(int colCounter = 1; colCounter < readCounts.length; colCounter++)
											{
												int thisDepth = Integer.parseInt(readCounts[colCounter]);
												varDepth += thisDepth;
											}

											// Determine the variant allele frequency //

											if(depth > 0)
											{
												double varFreq = (double) varDepth / (double) depth;

												// Only penalize if there's very little evidence for variant //

												if(varFreq < 0.05)
												{
													casesRef++;
												}
											}

										}
										else
										{
											casesRef++;
										}


									}



								}


							}
						}
						else
						{
							// CONTROL //
							controlsCalled++;
							if(MendelScan.isHeterozygous(gt))
							{
								controlsHet++;
							}
							else if(MendelScan.isHomozygous(gt))
							{
								controlsHom++;
							}
							else if(MendelScan.isReference(gt))
							{
								if(depth >= minDepth)
									controlsRef++;
							}

						}

					}
				}
				catch(Exception e)
				{
					System.err.println("Exception thrown while trying to parse data for " + sample + ": " + genotypes.get(sample));
				}

			}

		}
		catch(Exception e)
		{
			System.err.println("Exception thrown while calculating segregation score: " + e.getMessage());
		}

		segStatus = casesCalled + "\t" + casesRef + "\t" + casesHet + "\t" + casesHom;
		segStatus += "\t";
		segStatus += controlsCalled + "\t" + controlsRef + "\t" + controlsHet + "\t" + controlsHom;

		return(segStatus);
	}


	/**
	 * Returns the segregation score based on segregation patterns and expected mode of inheritance
	 *
	 * @param	params	Input parameters, for determining user settings
	 * @param	status	Type of file ("positions" or "regions")
	 * @return			Double of resulting population score for this variant
	 */
	static double getSegregationScore(String chrom, String inheritanceMode, String segStatus)
	{
		double segScore = 0.000;

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
			int controlsCalled = Integer.parseInt(segContents[4]);
			int controlsRef = Integer.parseInt(segContents[5]);
			int controlsHet = Integer.parseInt(segContents[6]);
			int controlsHom = Integer.parseInt(segContents[7]);
			int controlsVariant = controlsHet + controlsHom;

			if(inheritanceMode.equals("dominant"))
			{
				// Assume 50% sensitivity to detect heterozygotes; reduce score for affecteds called Ref //

				if(casesRef > 0)
				{
					for(int caseCounter = 0; caseCounter < casesRef; caseCounter++)
					{
						segScore = segScore * 0.50;
					}
				}

				// Reduce score for affecteds homozygous for a rare dominant-acting disease //

				if(casesHom > 0 && !chrom.equals("X")&& !chrom.equals("chrX") && !chrom.equals("Y")&& !chrom.equals("chrY"))
				{
					for(int caseCounter = 0; caseCounter < casesHom; caseCounter++)
					{
						segScore = segScore * 0.80;
					}
				}

				if(controlsVariant > 0)
				{
					for(int controlCounter = 0; controlCounter < controlsVariant; controlCounter++)
					{
						segScore = segScore * 0.10;
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
						segScore = segScore * 0.10;
					}
				}

				// Cases should be homozygous, so penalize if heterozygous. Compound hets done separately //

				if(casesHet > 0)
				{
					// Assume that cases called het are mis-called //
					for(int caseCounter = 0; caseCounter < casesHet; caseCounter++)
					{
						segScore = segScore * 0.50;
					}
				}

				// Controls should not be homozygous-variant; mis-calling is possible

				if(controlsHom > 0)
				{
					for(int controlCounter = 0; controlCounter < controlsHom; controlCounter++)
					{
						segScore = segScore * 0.50;
					}
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
		double annotScore = 0.00;

		try {
			if(vepScore >= 15)	// Nonsense, frameshift, splice site //
			{
				annotScore = 1.00;
			}
			else if(vepScore == 14) // Missense called damaging by 3/3 algorithms
			{
				annotScore = 0.95;
			}
			else if(vepScore == 14) // Missense called damaging by 2/3 algorithms
			{
				annotScore = 0.95;
			}
			else if(vepScore == 14) // Missense called damaging by 1/3 algorithms
			{
				annotScore = 0.95;
			}
			else if(vepScore == 14) // Missense called damaging by 0/3 algorithms
			{
				annotScore = 0.80;
			}
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
	 * @param	info	HashMap of dbSNP information from ID and INFO columns
	 * @return			String with one of several possible dbSNP statuses
	 */
	static String getDbsnpStatus(HashMap<String, String> info)
	{
		String status = "unknown";

		try {
			if(info.containsKey("G5") || info.containsKey("G5A"))
			{
				status = "common";
			}
			else if(info.containsKey("GMAF"))
			{
				Double gmaf = Double.parseDouble(info.get("GMAF"));
				if(gmaf >= 0.05)
				{
					status = "common";
				}
				else if(gmaf >= 0.01)
				{
					status = "uncommon";
				}
				else
				{
					status = "rare";
				}
			}
			else if(info.containsKey("KGPilot123"))
			{
				status = "rare";
			}
			else if(info.containsKey("RSnumber"))
			{
				status = "known";
			}
			else
			{
				status = "novel";
			}

			// Also check for dbSNP flags of mutations. These override novel/known/rare variant status //

			if(info.containsKey("MUT") || info.containsKey("CLN") || info.containsKey("PM"))
			{
				if(status.equals("novel") || status.equals("known") || status.equals("rare"))
				{
					status = "mutation";
				}
			}

		}
		catch(Exception e)
		{
			System.err.println("Warning: Exception thrown while determining dbSNP status: " + e.getMessage());
		}

		return(status);
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

		try {
			if(status.equals("novel"))
			{
				popScore = 1.00;
			}
			else if(status.equals("mutation"))
			{
				popScore = 0.95;
			}
			else if(status.equals("known"))
			{
				popScore = 0.60;
			}
			else if(status.equals("rare"))
			{
				popScore = 0.20;
			}
			else if(status.equals("uncommon"))
			{
				popScore = 0.02;
			}
			else if(status.equals("common"))
			{
				popScore = 0.01;
			}
		}
		catch(Exception e)
		{
			System.err.println("Warning: Exception thrown while calculating population score: " + e.getMessage());
		}

		return(popScore);
	}

} // End Class Definition

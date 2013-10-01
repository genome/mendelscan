/**
 * @(#)MendelScan.java
 *
 * Copyright (c) 2013 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.mendelscan;

//Import required packages //

import java.io.*;
import java.util.*;
import java.text.*;

/**
 * A set of tools for variant detection in next-generation sequence data.
 *
 * @version	2.3
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 * <BR>
 * <pre>
 * COMMANDS
 * score [vcf file] OPTIONS
 * 			Prioritize candidate variants in a VCF using segregation, population, annotation, and expression information
 * 			Input: 	VCF File
 * 					VEP annotation file
 * 			Output: Scored VCF output file
 *

 *
 * </pre>
 *
 *
 */
public class MendelScan {

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		String usage = "MendelScanScan v1.1\n\nUSAGE: java -jar MendelScan.jar [COMMAND] [OPTIONS] \n\n";
		usage = usage + "COMMANDS:\n" +
				"\tscore\t\tPrioritize variants in a VCF based on segregation, annotaton, population, and expression\n" +
				"\trhro\t\tPerform rare heterozygote rule out (RHRO) linkage analysis\n" +
				"\tsibd\t\tPerform shared identity-by-descent (SIBD) linkage analysis\n" +
				"\n";


		if(args.length > 0)
		{
			HashMap<String, String> params = getParams(args);

			// Proceed based on user's subcommand //
			if(args[0].equals("score"))
			{
				score(args, params);
			}
			else if(args[0].equals("rhro"))
			{
				rhro(args, params);
			}
			else if(args[0].equals("sibd"))
			{
				sibd(args, params);
			}
			else if(params.containsKey("help") || params.containsKey("h"))
			{
				// Print usage if -h or --help invoked //
				System.err.println(usage);
				return;
			}
			else
			{
				System.err.println("Subcommand " + args[0] + " not recognized!");
				System.err.println(usage);
			}
		}
		else
		{
			System.err.println(usage);
		}
	}


	/**
	 * Prioritize variants in a VCF
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	params			HashMap of parameters
	 */
	public static void score(String[] args, HashMap<String, String> params)
	{
		PrioritizeVCF myVCF = new PrioritizeVCF(args, params);
	}

	/**
	 * Performs rare heterozygote rule out (RHRO) analysis on a VCF
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	params			HashMap of parameters
	 */
	public static void rhro(String[] args, HashMap<String, String> params)
	{
		RareHetRuleOut myRHRO = new RareHetRuleOut(args, params);
	}

	/**
	 * Performs shared identity-by-descent (SIBD) analysis
	 *
	 * @param	args			Command-line arguments and parameters
	 * @param	params			HashMap of parameters
	 */
	public static void sibd(String[] args, HashMap<String, String> params)
	{
		SharedIBD mySIBD = new SharedIBD(args, params);
	}


	/**
	 * Parses and verifies any command-line parameters
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static HashMap getParams(String[] args)
	{
		HashMap<String, String> params = new HashMap<String, String>();

		// Parse out command line arguments //

		String arg = "";
		String value = "";
		int i = 0, j = 0;

		// Go through each argument in the command line //

		while (i < args.length)
		{
			j = i + 1;
			arg = args[i];

			// If the argument starts with a hyphen, make use of it //

			if (arg.startsWith("-"))
			{
				// Remove leading hyphens //
				while(arg.startsWith("-"))
				{
					arg = arg.replaceFirst("-", "");
				}

				// Parse out parameters followed by values //

				if (i < args.length && j < args.length && !args[j].startsWith("-"))
				{
					value = args[j];
					params.put(arg, value);
				}

				// Set other parameters to true //

				else
				{
					params.put(arg, "true");
				}
			}

			i++;
		}

		return(params);
	}


	/**
	 * Gets the infile from command line or input buffer
	 *
	 * @param	args	Command-line arguments
	 * @return			HashMap of parameter names and their values
	 */
	static BufferedReader getInfile(String[] args)
	{
		BufferedReader in = null;

	    try
	    {
	    	// Declare file-parsing variables //

	    	String line;

	    	// Check for file on command line //

	    	if(args.length > 1 && !args[1].startsWith("-"))
	    	{
	    		File infile = new File(args[1]);
	    		if(infile.exists())
	    		{
	    			// Parse the infile //
	    			System.err.println("Reading input from " + args[1]);
	    			in = new BufferedReader(new FileReader(args[1]));
	    		}
	    		else
	    		{
    				System.err.println("File not found: " + args[1] + "\n");
    				System.exit(10);
	    		}
	    	}

	    	// If no file from command line was parsed, try for piped input //

	    	if(in == null)
	    	{
		    	// Check the input stream //
		    	InputStreamReader instream = new InputStreamReader(System.in);
		    	Thread.sleep(1000);

		    	int num_naps = 0;

	    		while(!instream.ready())
	    		{
	    			System.err.println("Input stream not ready, waiting for 5 seconds...");
	    			Thread.sleep(5000);
	    			num_naps++;

	    			if(num_naps >= 100)
	    			{
	    				System.err.println("ERROR: Gave up waiting after 500 seconds...\n");
	    				System.exit(10);
	    			}
	    		}

		    	// If we have piped input, proceed with it //

		    	if(instream.ready())
		    	{
		    		System.err.println("Reading input from STDIN");
			    	in = new BufferedReader(instream);
		    	}
	    	}
	    }
	    catch(Exception e)
	    {
	    	System.err.println("ERROR: Unable to open input stream\n");
	    	System.exit(10);
	    }

		return(in);
	}


	/**
	 * Load PED sample information from a file
	 *
	 * @param	fileName	Name of the input file

	 * @return	annot		A hashmap of PED information by sample
	 */
	static HashMap loadPED(String fileName)
	{
		HashMap<String, String> ped = new HashMap<String, String>();

		BufferedReader in = null;

		// Open the infile //
		try
		{
			File infile = new File(fileName);
			in = new BufferedReader(new FileReader(infile));

			if(in != null && in.ready())
			{
				String line = "";
	    		while ((line = in.readLine()) != null)
	    		{
	    			String[] lineContents = line.split("\t");

	    			try {
	    				if(line.startsWith("#"))
	    				{
	    					// Skip comment lines in PED file //
	    				}
	    				else
	    				{
		    				String familyID = lineContents[0];
		    				String individualID = lineContents[1];
		    				String paternalID = lineContents[2];
		    				String maternalID = lineContents[3];
		    				String sex = lineContents[4];
		    				String status = lineContents[5];

		    				String pedLine = familyID + "\t" + paternalID + "\t" + maternalID + "\t" + sex + "\t" + status;

		    				ped.put(individualID, pedLine);

	    				}

	    			}
	    			catch(Exception e)
	    			{
	    				System.err.println("Error parsing PED line, so skipping: " + line);
	    			}
	    		}

	    		in.close();
			}
		}
		catch(Exception e)
		{
	    	System.err.println("ERROR: Unable to open VEP file " + fileName + " for reading\n");
	    	System.exit(10);
		}


		return(ped);
	}

	/**
	 * Load VEP annotation from a file
	 *
	 * @param	fileName	Name of the input file

	 * @return	annot		A hashmap of annotations by variant
	 */
	static HashMap loadVEP(String fileName)
	{
		HashMap<String, String> vep = new HashMap<String, String>();

		BufferedReader in = null;

		// Open the infile //
		try
		{
			File infile = new File(fileName);
			in = new BufferedReader(new FileReader(infile));

			if(in != null && in.ready())
			{
				String line = "";
	    		while ((line = in.readLine()) != null)
	    		{
	    			String[] lineContents = line.split("\t");

	    			if(line.startsWith("#"))
	    			{
	    				// VEP header line //
	    			}
	    			else
	    			{
	    				if(lineContents.length < 13)
	    				{
	    					System.err.println("Warning: VEP line had " + lineContents.length + " elements when 13+ expected; skipping");
	    				}
	    				else
	    				{
		    				String varName = lineContents[0];
		    				String varLocation = lineContents[1];
		    				String varAllele = lineContents[2];
		    				String ensGene = lineContents[3];
		    				String consequence = lineContents[6];
		    				String txPos = lineContents[7];
		    				String aaPos = lineContents[9];
		    				String aaChange = lineContents[10];
		    				String extra = lineContents[13];

		    				String hugoGene = "-";
		    				String canonical = "NO";
		    				String polyphen = "-";
		    				String sift = "-";
		    				String condel = "-";

		    				// Parse out the relevant extra information //
		    				try {
		    					if(extra.length() > 0)
		    					{
				    				String[] extraContents = extra.split(";");
				    				for (int fieldCounter = 0; fieldCounter < extraContents.length; fieldCounter++)
				    				{
				    					String[] fieldContents = extraContents[fieldCounter].split("=");
				    					if(fieldContents.length > 1)
				    					{
					    					String fieldName = fieldContents[0];
					    					String fieldValue = fieldContents[1];

					    					if(fieldName.equals("HGNC"))
					    						hugoGene = fieldValue;
					    					else if(fieldName.equals("CANONICAL"))
					    						canonical = fieldValue;
					    					else if(fieldName.equals("PolyPhen"))
					    						polyphen = fieldValue;
					    					else if(fieldName.equals("SIFT"))
					    						sift = fieldValue;
					    					else if(fieldName.equals("Condel"))
					    						condel = fieldValue;
				    					}

				    				}

		    					}

		    				}
		    				catch(Exception e)
		    				{
		    					System.err.println("Exception thrown while parsing VEP extra info for " + line);
		    					e.printStackTrace(System.err);
		    				}

		    				if(canonical.equals("YES"))
		    					canonical = "CANONICAL";

		    				String vepAnnot = ensGene;
		    				vepAnnot += "\t" + hugoGene;
		    				vepAnnot += "\t" + consequence;
		    				vepAnnot += "\t" + canonical;
		    				vepAnnot += "\t" + polyphen;
		    				vepAnnot += "\t" + sift;
		    				vepAnnot += "\t" + condel;
		    				//if(!txPos.equals("-"))

		    				vepAnnot += "\t" + txPos;
		    				vepAnnot += "\t" + aaPos;
		    				vepAnnot += "\t" + aaChange;

		    				int vepScore = getVEPscore(vepAnnot);
		    				vepAnnot += "\t" + vepScore;


		    				if(vep.containsKey(varName))
		    				{
		    					// If this is not canonical but the previous annotation is, keep it
		    					if(vep.get(varName).contains("CANONICAL"))
		    					{
		    						if(extra.contains("CANONICAL"))
		    						{
		    							// Add if they're both canonical which means different genes //
		    							vep.put(varName, vep.get(varName) + "\n" + vepAnnot);
		    						}
		    					}
		    					else if(extra.contains("CANONICAL"))
		    					{
		    						// THis one is canonical, previous was not, so overwrite //
		    						vep.put(varName, vepAnnot);
		    					}
		    					else
		    					{
		    						// Append if neither canonical //
		    						vep.put(varName, vep.get(varName) + "\n" + vepAnnot);
		    					}

		    				}
		    				else
		    				{
		    					vep.put(varName, vepAnnot);
		    				}
	    				}


	    			}
	    		}

	    		in.close();
			}
		}
		catch(Exception e)
		{
	    	System.err.println("ERROR: Unable to open VEP file " + fileName + " for reading\n");
	    	System.exit(10);
		}


		return(vep);
	}


	/**
	 * Assign numerical value to a VEP annotation class
	 *
	 * @param	fileName	Name of the input file

	 * @return	annot		A hashmap of annotations by variant
	 */
	static int getVEPscore (String vepAnnot)
	{
		int score = 0;
		String[] vepContents = vepAnnot.split("\t");
		String consequence = vepContents[2];
		String polyphen = vepContents[4];
		String sift = vepContents[5];
		String condel = vepContents[6];

		if(consequence.contains("STOP_GAINED"))
		{
			score = 17;
		}
		else if(consequence.contains("FRAMESHIFT_CODING"))
		{
			score = 16;
		}
		else if(consequence.contains("ESSENTIAL_SPLICE_SITE"))
		{
			score = 15;
		}
		else if(consequence.contains("NON_SYNONYMOUS_CODING"))
		{
			int numSayDeleterious = 0;
			int numSayNeutral = 0;

			if(polyphen.contains("damaging"))
				numSayDeleterious++;

			if(sift.contains("deleterious"))
				numSayDeleterious++;

			if(condel.contains("deleterious"))
				numSayDeleterious++;

			// Score is 11 plus the number of algorithms calling it deleterious, so max 14 //
			score = 11 + numSayDeleterious;
		}
		else if(consequence.contains("STOP_LOST"))
		{
			score = 10;
		}
		else if(consequence.contains("SPLICE_SITE"))
		{
			score = 9;
		}
		else if(consequence.contains("SYNONYMOUS_CODING") || consequence.contains("CODING_UNKNOWN") || consequence.contains("PARTIAL_CODON") || consequence.contains("COMPLEX_INDEL"))
		{
			score = 8;
		}
		else if(consequence.contains("WITHIN_MATURE_miRNA"))
		{
			score = 7;
		}
		else if(consequence.contains("WITHIN_NON_CODING_GENE"))
		{
			score = 6;
		}
		else if(consequence.contains("UTR"))
		{
			score = 5;
		}
		else if(consequence.contains("UPSTREAM"))
		{
			score = 4;
		}
		else if(consequence.contains("DOWNSTREAM"))
		{
			score = 3;
		}
		else if(consequence.contains("INTRONIC"))
		{
			score = 2;
		}
		else if(consequence.contains("INTERGENIC"))
		{
			score = 1;
		}

		return(score);
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
	 * Determines the numbers of cases/controls ref/het/hom for a variant
	 *
	 * @param	genotypes	Gender, disease status, and genotype of each sample
	 * @param	minDepth	Integer of minimum depth threshold for confident genotype calling
	 * @return	segStatus	A string with counts of cases/controls ref/het/hom at variant position
	 */
	static String getSegregationStatus(HashMap<String, String> genotypes, Integer minDepth)
	{
		String segStatus = "unknown";
		String altFreqs = "";

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

		// FOrmat scores for printing //
		DecimalFormat dfFreq = new DecimalFormat("0.000");

		try {
			for (String sample : genotypes.keySet())
			{
				try {
					String[] sampleContents = genotypes.get(sample).split("\t");
					String status = sampleContents[0];
					String gender = sampleContents[1];
					String gt = sampleContents[2];
					String sampleDP = sampleContents[3];
					String sampleAD = sampleContents[4];

					int sampleReads1 = 0;
					int sampleReads2 = 0;

					if(!gt.equals(".") && !gt.equals("./."))
					{
						// Obtain sequence depth and allele depth //
						Integer depth = 0;

						try {
							if(sampleContents[3].length() > 0 && !sampleContents[3].equals("."))
							{
								depth = Integer.parseInt(sampleContents[3]);
							}

							String[] adContents = sampleAD.split(",");
							for(int aCounter = 0; aCounter < adContents.length; aCounter++)
							{
								int thisReads2 = Integer.parseInt(adContents[aCounter]);
								if(thisReads2 > sampleReads2)
									sampleReads2 = thisReads2;
							}
						}
						catch (Exception e)
						{

						}

						double sampleVAF = 0.00;

						if(depth > 0)
						{
							sampleReads1 = depth - sampleReads2;
							sampleVAF = (double) sampleReads2 / (double) depth;
							altFreqs += "\t" + dfFreq.format(sampleVAF);
						}
						else
						{
							altFreqs += "\tNA";
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
										// If SampleVAF is appreciable, don't count as ref //
										if(sampleVAF >= 0.05)
										{
											// Just don't count it, or change it to het? //
											if(sampleVAF >= 0.10)
											{
												casesHet++;
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
					else
					{
						// Genotype was filtered, so altfreqs not calculated //
						altFreqs += "\tNA";
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
		segStatus += controlsCalled + "\t" + controlsRef + "\t" + controlsHet + "\t" + controlsHom;// + altFreqs;

		return(segStatus);
	}



	/**
	 * Load list of genes ranked by gene expression in tissue(s) of interest
	 *
	 * @param	fileName	Name of the input file

	 * @return	genes		A hashmap of genes ranked by expression in tissue(s) of interest
	 */
	static HashMap loadGeneExpression(String fileName)
	{
		HashMap<String, Double> genes = new HashMap<String, Double>();

		BufferedReader in = null;

		// Open the infile //
		try
		{
			File infile = new File(fileName);
			in = new BufferedReader(new FileReader(infile));

			if(in != null && in.ready())
			{
				String line = "";
				Integer lineCounter = 0;
	    		while ((line = in.readLine()) != null)
	    		{
	    			lineCounter++;
	    			String[] lineContents = line.split("\t");

	    			try {
	    				String gene = lineContents[0];
	    				genes.put(gene, (double) lineCounter);
	    			}
	    			catch(Exception e)
	    			{
	    				System.err.println("Error parsing PED line, so skipping: " + line);
	    			}
	    		}

	    		in.close();

	    		// Now that we have the total number of lines, go back through and compute the rank of each gene //
	    		for (String gene : genes.keySet())
	    		{
	    			double pctRank = 1.00 - (genes.get(gene) / (double) lineCounter);
	    			genes.put(gene, pctRank);

	    		}


			}
		}
		catch(Exception e)
		{
	    	System.err.println("ERROR: Unable to open VEP file " + fileName + " for reading\n");
	    	System.exit(10);
		}


		return(genes);
	}


	/**
	 * Load a BED file of centromeres by chromosome
	 *
	 * @param	fileName	Name of the input file

	 * @return	genes		A hashmap of centromeres by chromosome
	 */
	static HashMap loadCentromeres(String fileName)
	{
		HashMap<String, String> regions = new HashMap<String, String>();

		BufferedReader in = null;

		// Open the infile //
		try
		{
			File infile = new File(fileName);
			in = new BufferedReader(new FileReader(infile));

			if(in != null && in.ready())
			{
				String line = "";
				Integer lineCounter = 0;
	    		while ((line = in.readLine()) != null)
	    		{
	    			lineCounter++;
	    			String[] lineContents = line.split("\t");

	    			try {
	    				if(lineContents.length >= 3)
	    				{
		    				String chrom = lineContents[0];
		    				Integer chrStart = Integer.parseInt(lineContents[1]);
		    				Integer chrStop = Integer.parseInt(lineContents[2]);
		    				regions.put(chrom, chrStart + "\t" + chrStop);
	    				}

	    			}
	    			catch(Exception e)
	    			{
	    				System.err.println("Error parsing PED line, so skipping: " + line);
	    			}
	    		}

	    		in.close();


			}
		}
		catch(Exception e)
		{
	    	System.err.println("ERROR: Unable to open VEP file " + fileName + " for reading\n");
	    	System.exit(10);
		}


		return(regions);
	}



	/**
	 * Parse the INFO field of a VCF file
	 *
	 * @param	infoField	The INFO field contents

	 * @return	info		A hashmap of info field contents
	 */
	static HashMap parseInfoField (String infoField)
	{
		HashMap<String, String> info = new HashMap<String, String>();

		try {
			// Parse out relevant information //
			if(infoField.length() > 0)
			{
				String[] infoContents = infoField.split(";");

				if(infoContents.length > 0)
				{
					for(int valueCounter = 0; valueCounter < infoContents.length; valueCounter++)
					{
						String thisInfo = infoContents[valueCounter];
						if(thisInfo.contains("="))
						{
							String[] thisContents = thisInfo.split("=");
							if(thisContents.length > 1)
							{
								String fieldName = thisContents[0];
								String fieldValue = thisContents[1];
								info.put(fieldName, fieldValue);
							}
						}
						else
						{
							info.put(thisInfo, "1");
						}
					}
				}
			}
		}
		catch (Exception e)
		{
			System.err.println("Warning: Issue parsing VCF info field " + infoField + "\n" + e.getMessage());
		}

		return(info);
	}

	/**
	 * Parse a VCF genotype into a hash of fieldname-fieldvalue pairs
	 *
	 * @param	format			The format field of the VCF, e.g. GT:GQ:DP
	 * @param	genotypeField	The actual VCF genotype, e.g. 0/1:20:32

	 * @return	info		A hashmap of genotype information by column name
	 */
	static HashMap parseGenotypeField (String format, String genotypeField)
	{
		HashMap<String, String> info = new HashMap<String, String>();

		try {
			// Parse out relevant information //
			if(genotypeField.length() > 0)
			{
				String[] genotypeContents = genotypeField.split(":");
				String[] fieldNames = format.split(":");

				if(genotypeContents.length > 0 && fieldNames.length > 0)
				{
					for(int valueCounter = 0; valueCounter < genotypeContents.length; valueCounter++)
					{
						if(valueCounter < fieldNames.length)
						{
							String fieldName = fieldNames[valueCounter];
							String fieldValue = genotypeContents[valueCounter];
							info.put(fieldName, fieldValue);
						}

					}
				}
			}
		}
		catch (Exception e)
		{
			System.err.println("Warning: Issue parsing VCF info field " + genotypeField + "\n" + e.getMessage());
		}

		return(info);
	}


	/**
	 * Return true if a VCF genotype is heterozygous
	 *
	 * @param	genotype	The actual VCF genotype, e.g. 0/1
	 * @return	isHet		True if heterozygous, false if missing, ref, or homozygous
	 */
	static boolean isHeterozygous (String gt)
	{
		try {
			if(gt.contains("/"))
			{
				String[] gtContents = gt.split("/");
				if(gtContents.length > 1)
				{
					String a1 = gtContents[0];
					String a2 = gtContents[1];
					if(!a1.equals(a2))
					{
						return true;
					}
				}
			}
			else if(gt.contains("|"))
			{
				String[] gtContents = gt.split("|");
				if(gtContents.length > 1)
				{
					String a1 = gtContents[0];
					String a2 = gtContents[1];
					if(!a1.equals(a2))
					{
						return true;
					}
				}
			}
		}
		catch(Exception e)
		{

		}


		return false;
	}

	/**
	 * Return true if a VCF genotype is reference
	 *
	 * @param	genotype	The actual VCF genotype, e.g. 0/1
	 * @return	isHet		True if heterozygous, false if missing, ref, or homozygous
	 */
	static boolean isReference (String gt)
	{
		try {
			if(gt.contains("/"))
			{
				String[] gtContents = gt.split("/");
				if(gtContents.length > 1)
				{
					String a1 = gtContents[0];
					String a2 = gtContents[1];
					if(a1.equals(a2) && a1.equals("0"))
					{
						return true;
					}
				}
			}
			else if(gt.contains("|"))
			{
				String[] gtContents = gt.split("|");
				if(gtContents.length > 1)
				{
					String a1 = gtContents[0];
					String a2 = gtContents[1];
					if(a1.equals(a2) && a1.equals("0"))
					{
						return true;
					}
				}
			}
		}
		catch(Exception e)
		{

		}


		return false;
	}

	/**
	 * Return true if a VCF genotype is homozygous-variant
	 *
	 * @param	genotype	The actual VCF genotype, e.g. 0/1
	 * @return	isHet		True if heterozygous, false if missing, ref, or homozygous
	 */
	static boolean isHomozygous (String gt)
	{
		try {
			if(gt.contains("/"))
			{
				String[] gtContents = gt.split("/");
				if(gtContents.length > 1)
				{
					String a1 = gtContents[0];
					String a2 = gtContents[1];
					if(a1.equals(a2) && !a1.equals("0"))
					{
						return true;
					}
				}
			}
			else if(gt.contains("|"))
			{
				String[] gtContents = gt.split("|");
				if(gtContents.length > 1)
				{
					String a1 = gtContents[0];
					String a2 = gtContents[1];
					if(a1.equals(a2) && !a1.equals("0"))
					{
						return true;
					}
				}
			}
		}
		catch(Exception e)
		{

		}


		return false;
	}

}

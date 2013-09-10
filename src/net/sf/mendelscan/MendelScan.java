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
	    				String familyID = lineContents[0];
	    				String individualID = lineContents[1];
	    				String paternalID = lineContents[2];
	    				String maternalID = lineContents[3];
	    				String sex = lineContents[4];
	    				String status = lineContents[5];

	    				String pedLine = familyID + "\t" + paternalID + "\t" + maternalID + "\t" + sex + "\t" + status;

	    				ped.put(individualID, pedLine);
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
		else if(consequence.contains("SYNONYMOUS_CODING"))
		{
			score = 9;
		}
		else if(consequence.contains("CODING_UNKNOWN") || consequence.contains("PARTIAL_CODON") || consequence.contains("COMPLEX_INDEL"))
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

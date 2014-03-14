/**
 * @(#)RareHetRuleOut.java
 *
 * Copyright (c) 2013 Daniel C. Koboldt and Washington University in St. Louis
 *
 * COPYRIGHT
 */

package net.sf.mendelscan;

import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.TreeMap;

/**
 * A class for performing rare heterozygote rule out analysis to map dominant disease genes
 *
 * @version	1.1
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */

public class RareHetRuleOut {
	public RareHetRuleOut(String[] args, HashMap<String, String> params)
	{
		String usage = "USAGE: java -jar MendelScan.jar rhro [VCF] OPTIONS\n" +
		"\tOPTIONS:\n" +
		"\t--ped-file\tPedigree file in 6-column tab-delimited format\n" +
		"\t--centromere-file\tA tab-delimited, BED-like file indicating centromere locations by chromosome" +
		"\t--output-file\tOutput file to contain informative variants\n" +
		"\t--output-windows\tOutput file to contain RHRO windows. Otherwise they print to STDOUT\n" +
		"\t--inheritance\tPresumed model of inheritance: dominant, recessive, x-linked [dominant]\n";


		String pedFile = null;
		String outFile = null;
		String outWindows = null;
		String centromereFile = null;
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
			if(params.containsKey("ped-file"))
				pedFile = params.get("ped-file");

			if(params.containsKey("output-file"))
				outFile = params.get("output-file");

			if(params.containsKey("output-windows"))
				outWindows = params.get("output-windows");

			if(params.containsKey("inheritance"))
				inheritanceModel = params.get("inheritance");

			if(params.containsKey("centromere-file"))
				centromereFile = params.get("centromere-file");
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
	    	HashMap<String, String> centromeresByChrom = new HashMap();

	    	// Parse out sample information //
	    	HashMap<String, Boolean> controlSamples = new HashMap();
	    	HashMap<String, Boolean> maleSamples = new HashMap();

	    	if(params.containsKey("centromere-file"))
	    	{
	    		centromeresByChrom = MendelScan.loadCentromeres(centromereFile);
	    	}


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

	    	// Declare file-parsing variables //

	    	String line;

	    	long numSamples = 0;
	    	long numSamplesMale = 0;
	    	long numSamplesCase = 0;
	    	long numSamplesControl = 0;
	    	long numVariants = 0;

	    	TreeMap<String, Integer> stats = new TreeMap();

			// FOrmat scores for printing //
			DecimalFormat dfScore = new DecimalFormat("0.000");


	    	// Save the sample names by column number //

	    	HashMap<Integer, String> samplesByColumn = new HashMap();

	    	// Try to proceed with the input file //

	    	if(in != null && in.ready())
	    	{
	    		// Open the output files //
	    		PrintStream outFileHandle = null;
	    		PrintStream outWindowsHandle = null;
	    		if(params.containsKey("output-file"))
	    		{
	    			outFileHandle = new PrintStream( new FileOutputStream(outFile) );
	    		}

	    		if(params.containsKey("output-windows"))
	    		{
	    			outWindowsHandle = new PrintStream( new FileOutputStream(outWindows) );
	    		}

	    		// Centromere parameters //

	    		String centroChrom = "";
	    		Integer centroStart = 0;
	    		Integer centroStop = 0;

	    		// Windowing parameters //

	    		String currentChrom = "";
	    		String windowChrom = "";
	    		Integer windowStart = 0;
	    		Integer windowStop = 0;
	    		int windowMarkers = 0;

	    		// Parse the infile line by line //

	    		while ((line = in.readLine()) != null)
	    		{
	    			String[] lineContents = line.split("\t");

	    			if(line.startsWith("#"))
	    			{
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
    						try
    						{
    	    					// Parse the VCF main file //
    	    					String chrom = lineContents[0];
    	    					Integer position = Integer.parseInt(lineContents[1]);
    	    					String id = lineContents[2];
    	    					String ref = lineContents[3];
    	    					String altAlleles = lineContents[4];
    	    					String qual = lineContents[5];
    	    					String filter = lineContents[6];
    	    					String info = lineContents[7];
    	    					String format = lineContents[8];


    	    					if(!chrom.equals(currentChrom))
    	    					{
    	    						if(windowMarkers > 0)
    	    						{
			    						String reason = "newchrom";
		    							String windowLine = windowChrom + "\t" + windowStart + "\t" + windowStop + "\t" + windowMarkers + "\t" + reason;
		    							if(params.containsKey("output-windows"))
		    							{
		    								outWindowsHandle.println(windowLine);
		    							}
		    							else
		    							{
		    								System.out.println(windowLine);
		    							}

		    							// Count the window //

		    							if(stats.containsKey("windows"))
	    		    					{
	    		    						stats.put("windows", stats.get("windows") + 1);
	    		    					}
	    		    					else
	    		    					{
	    		    						stats.put("windows", 1);
	    		    					}
    	    						}
    	    						// Reset everything to zero if we're not currently in a window //
    	    						windowChrom = chrom;
    	    						windowStart = 1;
    	    						windowStop = 0;
    	    						windowMarkers = 0;
    	    					}

    	    					// If the centromere chromosome has changed, try to get it for this chrom //
    	    					if((centroChrom.length() == 0) || (centroChrom.length() > 0 && !centroChrom.equals(chrom)))
    	    					{
    	    						try {
        	    						// Try to get centromere coordinates //
    		    						if(centromeresByChrom.containsKey(chrom))
    		    						{
    		    							String[] centroContents = centromeresByChrom.get(chrom).split("\t");
    		    							centroStart = Integer.parseInt(centroContents[0]);
    		    							centroStop = Integer.parseInt(centroContents[1]);
    		    							centroChrom = chrom;
    		    							System.err.println("Centromere for " + chrom + " set to " + centroStart + "-" + centroStop);
    		    						}
    	    						}
    	    						catch (Exception e)
    	    						{
    	    							System.err.println("Warning: issue when trying to parse centromere: " + e.getMessage());
    	    						}

    	    					}




		    					boolean inCentromere = false;
		    					if(position > centroStart && position < centroStop)
		    						inCentromere = true;


    	    					if(!inCentromere && (filter.equals("PASS") || filter.equals(".")))
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

    		    							String sampleFT = ".";
    		    							if(sampleGenotype.containsKey("FT"))
    		    							{
    		    								sampleFT = sampleGenotype.get("FT");
    		    								// If filter status is non-passing, set it to missing //
    		    								if(!sampleFT.equals(".") && !sampleFT.equals("PASS"))
    		    									sampleGT = ".";
    		    							}
    		    							String sampleLine = sampleStatus + "\t" + sampleGender + "\t";
    		    							sampleLine += sampleGT + "\t" + sampleDP + "\t" + sampleAD;
    		    							genotypesBySample.put(sampleName, sampleLine);

    		    						}

    		    					}

    		    					// Parse out values from info field, which should include dbSNP values //

    		    					HashMap<String, String> infoValues = MendelScan.parseInfoField(info);

    		    					// Save dbSNP ID if necessary //
    		    					if(!id.equals("."))
    		    					{
    		    						infoValues.put("RSnumber", id);
    		    					}

    		    					// Determine dbSNP status of variant //

    		    					String dbsnpStatus = MendelScan.getDbsnpStatus(infoValues);

    		    					if(stats.containsKey("dbSNP " + dbsnpStatus))
    		    					{
    		    						stats.put("dbSNP " + dbsnpStatus, stats.get("dbSNP " + dbsnpStatus) + 1);
    		    					}
    		    					else
    		    					{
    		    						stats.put("dbSNP " + dbsnpStatus, 1);
    		    					}

    		    					// Get Segregation Status //
    		    					String segStatus = MendelScan.getSegregationStatus(genotypesBySample, minDepth);
    		    					String[] segContents = segStatus.split("\t");
    		    					int casesCalled = Integer.parseInt(segContents[0]);
    		    					int casesRef = Integer.parseInt(segContents[1]);
    		    					int casesHet = Integer.parseInt(segContents[2]);
    		    					int casesHom = Integer.parseInt(segContents[3]);
    		    					int controlsCalled = Integer.parseInt(segContents[4]);
    		    					int controlsRef = Integer.parseInt(segContents[5]);
    		    					int controlsHet = Integer.parseInt(segContents[6]);
    		    					int controlsHom = Integer.parseInt(segContents[7]);
    		    					int controlsVariant = controlsHet + controlsHom;

    		    					String markerType = "";
	    							double pctCasesHet = 0.00;
	    							if(casesCalled > 0)
	    								pctCasesHet = (double) casesHet / (double) casesCalled;

    		    					// Part 1: Homozygous difference bewteen affected //
    		    					if(casesRef > 0 && casesHom > 0)
    		    					{
    		    						markerType = "RuleOut";
//    		    						System.err.println(segStatus + "\tRuleOut");
    		    					}
    		    					// Part 2: Rare Heterozygote //
    		    					else if(!dbsnpStatus.equals("common") && !dbsnpStatus.equals("uncommon"))
    		    					{
    		    						// If rare variant not found in controls... //
    		    						if(controlsVariant == 0 && casesCalled > 0 && casesHet > 0)
    		    						{
    		    							if(pctCasesHet > 0.90)
        		    							markerType = "SharedHet";
    		    							else
    		    								markerType = "RareHet";
//    		    							System.err.println(segStatus + "\tSharedHet\t" + pctCasesHet);
    		    						}
    		    					}

    		    					if(markerType.length() > 0)
    		    					{
    		    						// Count it //

        		    					if(stats.containsKey("markers " + markerType))
        		    					{
        		    						stats.put("markers " + markerType, stats.get("markers " + markerType) + 1);
        		    					}
        		    					else
        		    					{
        		    						stats.put("markers " + markerType, 1);
        		    					}

    		    						// Print it //
    		    						String markerLine = chrom + "\t" + position + "\t" + id + "\t" + dbsnpStatus + "\t" + markerType + "\t" + pctCasesHet + "\t" + segStatus;
    		    						if(params.containsKey("output-file"))
    		    						{
    		    							outFileHandle.println(markerLine);
    		    						}
    		    						else
    		    						{
    		    							System.out.println(markerLine);
    		    						}
    		    					}
    		    					else
    		    					{
    		    						// Not an informative marker, so count it //
        		    					if(stats.containsKey("markers not informative"))
        		    					{
        		    						stats.put("markers not informative", stats.get("markers not informative") + 1);
        		    					}
        		    					else
        		    					{
        		    						stats.put("markers not informative", 1);
        		    					}
    		    					}



    		    					// If we have a window going, and either centromere or new chromosome was reached, end it //
    		    					if(windowMarkers > 0 && windowStop < centroStop && position > centroStart)
    		    					{
    		    							String reason = "centromere";
    		    							// Adjust window to position just before centromere //
    		    							windowStop = centroStart - 1;
    		    							// windowStop = position - 1; We don't update window stop because now we're on a new chromosome
    		    							String windowLine = windowChrom + "\t" + windowStart + "\t" + windowStop + "\t" + windowMarkers + "\t" + reason;
    		    							if(params.containsKey("output-windows"))
    		    							{
    		    								outWindowsHandle.println(windowLine);
    		    							}
    		    							else
    		    							{
    		    								System.out.println(windowLine);
    		    							}

    		    							// Count the window //

    		    							if(stats.containsKey("windows"))
            		    					{
            		    						stats.put("windows", stats.get("windows") + 1);
            		    					}
            		    					else
            		    					{
            		    						stats.put("windows", 1);
            		    					}

    			    						// Reset window //
    			    						windowMarkers = 0;
    			    						windowChrom = "";
    			    						windowStop = 0;
    			    						// Set next position as earliest possible window start //
    			    						windowStart = centroStop + 1;

    		    					}


    		    					// If an informative marker type was obtained, use it for windowing //

    		    					if(markerType.equals("RuleOut"))
    		    					{
    		    						// If we had a region going, terminate it here //
    		    						if(windowMarkers > 0)
    		    						{
    		    							windowStop = position - 1;

    		    							// Count the window //

    		    							if(stats.containsKey("windows"))
            		    					{
            		    						stats.put("windows", stats.get("windows") + 1);
            		    					}
            		    					else
            		    					{
            		    						stats.put("windows", 1);
            		    					}

    		    							// Print the window //

    		    							String windowLine = windowChrom + "\t" + windowStart + "\t" + windowStop + "\t" + windowMarkers + "\truleout";
    		    							if(params.containsKey("output-windows"))
    		    							{
    		    								outWindowsHandle.println(windowLine);
    		    							}
    		    							else
    		    							{
    		    								System.out.println(windowLine);
    		    							}
    		    						}

    		    						// Reset window //
    		    						windowMarkers = 0;
    		    						windowChrom = chrom;
    		    						windowStop = 0;
    		    						// Set next position as earliest possible window start //
    		    						windowStart = position + 1;
    		    					}
    		    					else if(markerType.equals("SharedHet"))
    		    					{
    		    						// Start or continue a window //
    		    						if(windowMarkers > 0)
    		    						{
    		    							// Continue a window //
    		    							windowStop = position;
    		    							windowMarkers++;
    		    						}
    		    						else
    		    						{
    		    							// If we don't have a start position, set it now //
    		    							if(!windowChrom.equals(chrom))
    		    								windowStart = 1;
    		    							else if(windowStart == 0)	// Start it at first base of chromosome if no reason not to //
    		    								windowStart = 1;

    		    							windowChrom = chrom;

    		    							windowStop = position;
    		    							windowMarkers = 1;
    		    						}
    		    					}
    		    					else
    		    					{
    		    						// Neither shared het nor rule out, but if we have a window going, move it forward
    		    						if(windowMarkers > 0)
    		    							windowStop = position;
    		    						// Otherwise if we have a different chromosome //
    		    						else if(!windowChrom.equals(chrom))
    		    							windowChrom = chrom;
    		    					}


    	    					}

    	    					currentChrom = chrom;

    						}
    						catch(Exception e)
    						{
    							System.err.println("Parsing exception thrown for line: " + line);
    							System.err.println(e.getMessage());
    						}



    					}
	    			}

	    		}

				// If we had a region going, print it out //
				if(windowMarkers > 0)
				{
					// Count the window //

					if(stats.containsKey("windows"))
					{
						stats.put("windows", stats.get("windows") + 1);
					}
					else
					{
						stats.put("windows", 1);
					}

					// Print the window //
					String windowLine = windowChrom + "\t" + windowStart + "\t" + windowStop + "\t" + windowMarkers + "\tend";
					if(params.containsKey("output-windows"))
					{
						outWindowsHandle.println(windowLine);
					}
					else
					{
						System.out.println(windowLine);
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


	    	// Print Summary Statistics //

	    	System.err.println(numSamples + " samples in VCF (" + numSamplesCase + " affected, " + numSamplesControl + " unaffected, " + numSamplesMale + " male)");

	    	System.err.println(numVariants + " variants in VCF file");

	    	for (String statKey : stats.keySet())
	    	{
	    		System.err.println(stats.get(statKey) + "\t" + statKey);
	    	}

	    	//TreeMap<String, Object> statMap = new TreeMap<String, Object>();



	    }
	    catch(Exception e)
	    {
	    	System.err.println("Exception: " + e.getLocalizedMessage());
	    	e.printStackTrace(System.err);
	    	System.exit(11);
	    }
	}

}

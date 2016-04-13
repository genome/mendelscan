/**
 * @(#)SharedIBD.java
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
import java.util.HashMap;
import java.util.TreeMap;

/**
 * A class for performing shared identity-by-descent analysis to map dominant disease genes
 *
 * @version	1.1
 *
 * @author Daniel C. Koboldt <dkoboldt@genome.wustl.edu>
 *
 */

public class SharedIBD {
	public SharedIBD(String[] args, HashMap<String, String> params)
	{
		String usage = "USAGE: java -jar MendelScan.jar sibd [FIBD] OPTIONS\n" +
		"\tFIBD: The uncompressed FastIBD output file from BEAGLE\n" +
		"\tOPTIONS:\n" +
		"\t--ped-file\tPedigree file in 6-column tab-delimited format\n" +
		"\t--markers-file\tMarkers file in BEAGLE format\n" +
		"\t--centromere-file\tA tab-delimited, BED-like file indicating centromere locations by chromosome" +
		"\t--output-file\tOutput file to contain IBD markers with chromosomal coordinates\n" +
		"\t--output-windows\tOutput file to contain RHRO windows. Otherwise they print to STDOUT\n" +
		"\t--ibd-score-threshold\tMaximum Beagle FastIBD score below which segments will be used [10e-10]\n" +
		"\t--window-resolution\tWindow size in base pairs to use for SIBD region binning [100000]\n" +
		"\t--inheritance\tPresumed model of inheritance: dominant, recessive, x-linked [dominant]\n";

		String pedFile = null;
		String markersFile = null;
		String outFile = null;
		String outWindows = null;
		String centromereFile = null;
		String inheritanceModel = "dominant";
		String ibdScoreThreshold = "10E-10";
		Integer minDepth = 20;
		Integer windowResolution = 100000;

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

			if(params.containsKey("markers-file"))
				markersFile = params.get("markers-file");

			if(params.containsKey("output-file"))
				outFile = params.get("output-file");

			if(params.containsKey("output-windows"))
				outWindows = params.get("output-windows");

			if(params.containsKey("inheritance"))
				inheritanceModel = params.get("inheritance");

			if(params.containsKey("ibd-score-threshold"))
			{
				ibdScoreThreshold = params.get("ibd-score-threshold");
				try {
					System.err.println("IBD score threshold: " + Double.valueOf(ibdScoreThreshold));
				}
				catch (Exception e)
				{
					System.err.println("Warning: invalid IBD score threshold provided; please use scientific notation e.g. 10E-10\n");
					return;
				}
			}

			if(params.containsKey("window-resolution"))
			{
				try {
					windowResolution = Integer.parseInt(params.get("window-resolution"));
				}
				catch(Exception e)
				{
					System.err.println("Error: unable to parse an integer window resolution from the value you provided. Provide something like: 100000\n");
					return;
				}
			}

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

	    	TreeMap<String, Integer> stats = new TreeMap();

	    	HashMap<String, String> pedInfo = null;
	    	HashMap<String, String> centromeresByChrom = new HashMap();

	    	HashMap<String, String> pairSegments = new HashMap();

	    	// Parse out sample information //
	    	HashMap<String, Boolean> affectedSamples = new HashMap();
	    	HashMap<String, Boolean> controlSamples = new HashMap();
	    	HashMap<String, Boolean> maleSamples = new HashMap();

	    	if(params.containsKey("centromere-file"))
	    	{
	    		centromeresByChrom = MendelScan.loadCentromeres(centromereFile);
	    	}

	    	// Load the pedigree file //

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
//		    			String familyID = pedLineContents[0];
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
		    				affectedSamples.put(sample, true);
		    			}
		    		}

		    		System.err.println(numMales + " males, " + numCases + " cases, " + numControls + " controls");
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Warning: Exception thrown while loading pedigree information: " + e.getMessage());
		    		System.err.println(e.getLocalizedMessage());
		    		e.printStackTrace(System.err);
		    		return;
		    	}
	    	}
	    	else
	    	{
	    		System.err.println("ERROR: Please provide a pedigree file!\n" + usage);
	    		return;
	    	}

	    	// Determine number of affected pairs //
	    	HashMap<String, Boolean> affectedPairs = new HashMap();
	    	int numAffectedPairs = 0;

	    	for (String sample1 : affectedSamples.keySet())
	    	{
	    		for (String sample2 : affectedSamples.keySet())
	    		{
	    			if(!sample1.equals(sample2))
	    			{
	    				String key1 = sample1 + "\t" + sample2;
	    				String key2 = sample2 + "\t" + sample1;
	    				if(affectedPairs.containsKey(key1) || affectedPairs.containsKey(key2))
	    				{
	    					// This pair already counted, so ignore it //
	    				}
	    				else
	    				{
	    					// Count affected pairs //
	    					numAffectedPairs++;
	    					affectedPairs.put(key1, true);
	    				}
	    			}
	    		}
	    	}

	    	stats.put("affected pairs", numAffectedPairs);


	    	// Load the markers file //
	    	String[] markers = null;
	    	if(params.containsKey("markers-file"))
	    	{
		    	System.err.println("Loading markers from " + markersFile + "...");
		    	try {
		    		markers = loadMarkers(markersFile);
		    		stats.put("markers in BEAGLE markers file", markers.length);
		    	}
		    	catch(Exception e)
		    	{
		    		System.err.println("Exception while loading markers: " + e.getMessage());
		    		return;
		    	}
	    	}
	    	else
	    	{
	    		System.err.println("ERROR: Please provide a markers file!\n" + usage);
	    		return;
	    	}

    		// Open the output files //
    		PrintStream outFileHandle = null;
    		PrintStream outWindowsHandle = null;
    		if(params.containsKey("output-file"))
    		{
    			outFileHandle = new PrintStream( new FileOutputStream(outFile) );
    			outFileHandle.println("chrom\tchr_start\tchr_stop\tsample1\tsample2\tindex_start\tindex_stop\tp_value\tnum_markers");
    		}

    		if(params.containsKey("output-windows"))
    		{
    			outWindowsHandle = new PrintStream( new FileOutputStream(outWindows) );
    			outWindowsHandle.println("chrom\tchr_start\tchr_stop\ttotal_pairs\tpairs_ibd");
    		}

    		// Centromere parameters //

    		String centroChrom = "";
    		Integer centroStart = 0;
    		Integer centroStop = 0;

    		// Parse the infile of FIBD results line by line //
	    	String line = "";

    		while ((line = in.readLine()) != null)
    		{
    			String[] lineContents = line.split("\t");
    			if(lineContents.length >= 5)
    			{
    				try {
    					String sample1 = lineContents[0];
    					String sample2 = lineContents[1];

    					// Proceed only with lines for non-control samples //

    					if(affectedSamples.containsKey(sample1) && affectedSamples.containsKey(sample2))
    					{
        					Integer index1 = Integer.parseInt(lineContents[2]);
        					Integer index2 = Integer.parseInt(lineContents[3]);
        					Double pvalue = Double.valueOf(lineContents[4]);

        					if(pvalue < Double.valueOf(ibdScoreThreshold))
        					{
        						// IBD segment significant //
        						if(stats.containsKey("IBD segments significant"))
        						{
        							stats.put("IBD segments significant", stats.get("IBD segments significant") + 1);
        						}
        						else
        						{
        							stats.put("IBD segments significant", 1);
        						}

            					Integer numMarkers = index2 - index1;

            					if(index1 >= markers.length)
            						index1 = markers.length - 1;

            					if(index2 >= markers.length)
            						index2 = markers.length - 1;

            					String markerLine1 = markers[index1];
            					String markerLine2 = markers[index2];

            					String[] markerLineContents1 = markerLine1.split("\t");
            					String[] markerLineContents2 = markerLine2.split("\t");

            					String[] chrPos1 = markerLineContents1[0].split(":");
            					String chrom = chrPos1[0];
            					String chrStart = chrPos1[1];

            					String[] chrPos2 = markerLineContents2[0].split(":");
            					String chrom2 = chrPos2[0];
            					String chrStop = chrPos2[1];

            					String pairKey = sample1 + "\t" + sample2;
            					if(pairSegments.containsKey(pairKey))
            					{
            						pairSegments.put(pairKey, pairSegments.get(pairKey) + "\n" + chrom + "\t" + chrStart + "\t" + chrStop + "\t" + numMarkers);
            					}
            					else
            					{
            						pairSegments.put(pairKey, chrom + "\t" + chrStart + "\t" + chrStop + "\t" + numMarkers);
            					}

            					// Print the line with chromosomal coordinates //

            					String outLine = chrom + "\t" + chrStart + "\t" + chrStop + "\t" + line + "\t" + numMarkers;
            					if(params.containsKey("output-file"))
            					{
            						outFileHandle.println(outLine);
            					}
            					else
            					{
                					System.out.println(outLine);
            					}

        					}
        					else
        					{
        						// IBD segment not significant //
        						if(stats.containsKey("IBD segments not significant"))
        						{
        							stats.put("IBD segments not significant", stats.get("IBD segments not significant") + 1);
        						}
        						else
        						{
        							stats.put("IBD segments not significant", 1);
        						}
        					}



    					}
    					else
    					{
    						System.err.println("Samples not affected: " + sample1 + " " + sample2);
    					}




    				}
    				catch(Exception e)
    				{
    					System.err.println("Warning: Exception thrown while parsing this line: " + line + "\n in the FIBD file: " + e.getMessage() + "\n" + e.getLocalizedMessage());
    					e.printStackTrace(System.err);
    				}
    			}
    			else
    			{
    				System.err.println("Warning: line in FIBD file has less than 5 columns and is ignored: " + line);
    			}
    		}

    		in.close();


    		// GO through and do the windowing //
    		String windowChrom = "";
    		TreeMap<Integer, Integer> windowIBDpairs = new TreeMap();
    		HashMap<Integer, Boolean> positionCounted = new HashMap();

    		for (String pairKey : pairSegments.keySet())
    		{
    			String[] pairLines = pairSegments.get(pairKey).split("\n");
    			for(int lineCounter = 0; lineCounter < pairLines.length; lineCounter++)
    			{
    				String pairLine = pairLines[lineCounter];
    				String[] lineContents = pairLine.split("\t");

    				String chrom = lineContents[0];
    				windowChrom = chrom;
    				Integer chrStart = Integer.parseInt(lineContents[1]);
    				Integer chrStop = Integer.parseInt(lineContents[2]);

    				Integer windowStart = chrStart / windowResolution;
    				Integer windowStop = chrStop / windowResolution;

//    				System.err.println(pairKey + "\t" + pairLine + "\t" + windowStart + "-" + windowStop);

    				for(Integer windowPos = windowStart; windowPos <= windowStop; windowPos++)
    				{
    					if(positionCounted.containsKey(windowPos))
    					{
    						// Already counted this position for this sample pair, so skip it //
    					}
    					else
    					{
    						// Add to the total or start total at one //
    						if(windowIBDpairs.containsKey(windowPos))
    						{
    							windowIBDpairs.put(windowPos, windowIBDpairs.get(windowPos) + 1);
    						}
    						else
    						{
    							windowIBDpairs.put(windowPos, 1);
    						}
    					}
    				}
    			}
    		}

    		// Now that we've tallied up the windows, let's go through them one at a time and print the numbers //

    		for (Integer windowPos : windowIBDpairs.keySet())
    		{
    			Integer windowStart = windowPos * windowResolution;
    			Integer windowStop = windowStart + windowResolution - 1;
    			Integer ibdPairs = windowIBDpairs.get(windowPos);

    			String outLine = windowChrom + "\t" + windowStart + "\t" + windowStop + "\t" + stats.get("affected pairs") + "\t" + ibdPairs;

    			if(params.containsKey("output-windows"))
    			{
    				outWindowsHandle.println(outLine);
    			}
    			else
    			{
    				System.err.println(outLine);
    			}

    		}


	    	// Print summary statistics //

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


	}

	/**
	 * Load markers from the BEAGLE markers file
	 *
	 * @param	fileName	Name of the input file

	 * @return	annot		A hashmap of PED information by sample
	 */
	static String[] loadMarkers(String fileName)
	{
		String lines = "";


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
	    			if(lines.length() > 0)
	    			{
	    				lines += "\n" + line;
	    			}
	    			else
	    			{
	    				lines = line;
	    			}

	    		}

	    		in.close();

	    		String[] markers = lines.split("\n");
	    		return(markers);
			}
		}
		catch(Exception e)
		{
	    	System.err.println("ERROR: Exception thrown while parsing markers file " + fileName + "\n");
	    	System.exit(10);
		}

		// IF we haven't returned markers here there's an issue, so return empty array //

		String[] markers = new String[1];
		return(markers);

	}

}

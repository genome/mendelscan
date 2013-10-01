## Documentation ##

MendelScan is a command-line program with multiple subcommands (e.g. score, rhro, and sibd). Each subcommand has a unique set of inputs and outputs. For the list of available subcommands, enter:
	java -jar MendelScan.jar --help

## Available Subcommands ##

These subcommands are currently supported:

	java -jar MendelScan.jar score 		Prioritize a VCF 
	java -jar MendelScan.jar rhro		Perform RHRO analysis
	java -jar MendelScan.jar sibd		Perform SIBD analysis
	
For detailed usage information, enter the subcommand followed by -h or --help, e.g.:

	java -jar MendelScan.jar score -h

For those familiar with Java, the auto-generated Javadoc documentation may be useful as well.

## How to Cite ##

A manuscript describing MendelScan is currently under review. In the interim, please cite MendelScan by noting the version and citing http://gmt.genome.wustl.edu/mendelscan.

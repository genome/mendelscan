## Documentation ##

MendelScan is a command-line program with multiple subcommands (e.g. score, rhro, and sibd). Each subcommand has a unique set of inputs and outputs. For the list of available subcommands, enter:

<pre class="terminal">
java -jar MendelScan.jar --help
</pre>

## Available Subcommands ##

These subcommands are currently supported:

<pre class="terminal">
java -jar MendelScan.jar score &#35; Prioritize a VCF
java -jar MendelScan.jar rhro  &#35; Perform RHRO analysis
java -jar MendelScan.jar sibd  &#35; Perform SIBD analysis
</pr>

For detailed usage information, enter the subcommand followed by -h or --help, e.g.:

<pre class="terminal">
	java -jar MendelScan.jar score -h
</pre>

For those familiar with Java, the auto-generated Javadoc documentation may be useful as well.

## How to Cite ##

A manuscript describing MendelScan is currently under review. In the interim, please cite MendelScan by noting the version and citing http://gmt.genome.wustl.edu/mendelscan.

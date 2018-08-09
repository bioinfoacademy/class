#!/usr/bin/perl -w
#This script gets gets the pubmed abstracts for a given term
use LWP::Simple;

#search term to find
$search_term = "breast cancer";

#replace space with +
$search_term =~ s/\s/+/g;

#print $search_term;

#maximum number of results to retrieve
$retmax = 10;

#base url
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';

#set the database to search
$db_name = 'pmc';

#generate the query url
$url = $base."esearch.fcgi?db=$db_name&retmax=$retmax&term=$search_term";

#print $url;

#Submit the search and retrieve the XML based results
$esearch_result=get($url);

#print "$esearch_result";

#extract paper IDs using match regex . anything * anytime, globally
@ids = ($esearch_result =~ m|.*<Id>(.*)</Id>.*|g);

#print join(',',@ids);

#delete old abstract.txt file
unlink "abstracts.txt";

#loop through all the ids
foreach $id (@ids)
	{
	#print "$id\n";
	#get abstract for each pubmed id
	$fetchurl = $base."efetch.fcgi?db=pubmed&id=$id&retmode=text&rettype=abstract";
	#print get($fetchurl);
	#open a file for appending the output
	open(OUTFILE,'>>','abstracts.txt');
	#get the results and print to the filehandle
	print OUTFILE get($fetchurl);
	#close file
	close OUTFILE;
	}


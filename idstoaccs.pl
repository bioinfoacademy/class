#!/usr/local/bin/perl -w
use LWP::Simple;

#query list of ids
$gi_list = '24475906,224465210,50978625,9507198';

#assemble the URL
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "efetch.fcgi?db=nucleotide&id=$gi_list&rettype=acc";

#post the URL query  to NCBI
$output = get($url);

#show the output in the screen
print "$output";

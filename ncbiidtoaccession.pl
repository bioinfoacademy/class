use LWP::Simple;

#query
$gi_list = '24475906,224465210,50978625,9507198';

#url
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$url = $base . "efetch.fcgi?db=nucleotide&id=$gi_list&rettype=acc";

#post
$output = get($url);

#print
print "$output";



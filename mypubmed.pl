use LWP::Simple;

unlink "abstract.txt";

$search_term = "coronavirus SARS";
$search_term =~ s/\s/+/g;
$retmax = 10;
$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
$db_name = 'pmc';
$url = $base."esearch.fcgi?db=$db_name&retmax=$retmax&term=$search_term";
$esearch_result=get($url);

@ids=($esearch_result =~ m|.*<Id>(.*)</Id>.*|g);

foreach $id (@ids)
    {
    print "$id\n";
    $fetchedurl = $base."efetch.fcgi?db=pubmed&id=$id&retmode=text&rettype=abstract";
    #print $fetchedurl,"\n";
    open(OUTFILE,'>>','abstract.txt');
    print OUTFILE get($fetchedurl);
    sleep 1;
    close OUTFILE;
    }

    

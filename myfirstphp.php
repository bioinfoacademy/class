<?php

#Author - Vijay Nagarajan
print "Hello world\n";

#variable
$DNA = 'AGCTAGCTGACTATACGTGACTAGCTATCG';
print $DNA."\n";

$DNA2 = 'ACCGATGCTAGCTAGCTATGCTA';
print "here are the original dna fragments:\n\n";
print $DNA."\n";
print $DNA2."\n";
$DNA3 = "$DNA$DNA2";
print "here is my concatenation:\n\n";
print $DNA3."\n";

$RNA = $DNA3;
$RNA = str_replace("T","U",$RNA);
print "here is my rna:\t";
print "$RNA\n";

?>


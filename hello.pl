#!/usr/bin/perl

#Author - Vijay Nagarajan PhD
#print a string - hello vijay
print "Hello good job \n";

#variable
$DNA = 'ACGGGAGGACGGGAAAATTACTACGGCATTAGC';
print $DNA,"\n";

#string concatenation
$DNA2 = 'ATAGTGCCGTGAGAGTGATGTAGTA';
print "Here are the original two DNA fragments:\n\n";
print $DNA,"\n";
print $DNA2,"\n\n";
$DNA3 = "$DNA$DNA2";
print "Here is the concatenation of the first two fragments (version 1):\n\n";
print $DNA3,"\n";

#Transcribe DNA to RNA
$RNA = $DNA3;
$RNA =~ s/T/U/g;
print "Here is the result of transcribing the DNA to RNA:\n\n";
print "$RNA\n";




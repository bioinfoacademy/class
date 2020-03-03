$DNA = 'AGCTAGCTATCGTACATCGATCGATGCTGATCGATGCATCG';
$DNA2 = 'GCATCGATGCTAGCTAGTGCACTA';

print "Here are the original two DNA fragments:\n\n";
print $DNA,"\n";
print $DNA2,"\n\n";

$DNA3 = "$DNA$DNA2";
print "Here is the concatenation of the first two fragments:\t";
print $DNA3,"\n";

#we are going to transcribe now
$RNA = $DNA;
$RNA =~ s/T/U/g;
print "Here is my mrna:\n\n";
print "$RNA\n";



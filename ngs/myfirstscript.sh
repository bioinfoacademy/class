#display text to user and wait for input
echo -n "Enter your file name with full path: "

#read user input in to the variable myinputfile
read myinputfile

#this is a loop from 1 to 22
for i in {1..22}
#this block does whatever i want with the loop values
do
#i am just echoing the word chr
    echo -n "chr"
    echo -n $i
    echo -n ":"
#i am grepping for all the peaks that are in a particular chr and count them
    grep -w chr$i $myinputfile | wc -l
done


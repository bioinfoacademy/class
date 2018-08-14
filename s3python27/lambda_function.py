from __future__ import print_function

import json
import urllib
import boto3

print('Loading function')

s3 = boto3.client('s3')
s4 = boto3.resource("s3")

def lambda_handler(event, context):
    #open a file with a list of numbers, calculate total number of items, sum and average and write it to a new file
    # Get the object from the event and show its content type
    bucket = event['Records'][0]['s3']['bucket']['name']
    key = urllib.unquote_plus(event['Records'][0]['s3']['object']['key'].encode('utf8'))
    print(bucket + key)
    
    #start an empty array
    numberlist=[]
    
    #get the file contents
    try:
        response = s3.get_object(Bucket=bucket, Key=key)
        #extract data from the s3 object
        s3data = response['Body'].read()
        #convert object to string
        stringdata = str(s3data)
        #create an array, splitting string by new line
        numberlist = stringdata.split("\n") 
        #convert string to int in array
        numberlist = map(int, numberlist)
        
        print(numberlist)
        #sum of array
        nsum=sum(numberlist)
        #average of list
        nave=round(nsum/len(numberlist),2)
        #lenght of list
        nlen=len(numberlist)
        
        #write the output to string
        stringg = "Total number of items :"+str(nlen)+"\n"+"Sum is :"+str(nsum)+"\n"+"Average is :"+str(nave)
        #encode string as ascii
        encoded_string = stringg.encode("utf-8")
        file_name = "sum.txt"
        #write string to file and put object to the bucket
        s4.Bucket(bucket).put_object(Key=file_name, Body=encoded_string)


    except Exception as e:
        print(e)
        print('Error getting object {} from bucket {}. Make sure they exist and your bucket is in the same region as this function.'.format(key, bucket))
        raise e
#!/bin/bash
#put this in the same directory as the txt files which contain lists of PFNs

#the elements in both arrays should have similar names
txtFiles=('testTTBar.txt')
outputTxtFiles=('testTTBar_tailored.txt')

for r in ${!txtFiles[*]}
do
	#count the number of lines in the file
	f=$(wc -l < ${txtFiles[$r]})
	l=$(($f-1))
	#echo $l
	#echo $f
	
	#replace /eos with file:/eos
	eval "sed 's&/eos&file:/eos&g' ${txtFiles[$r]} > temp.txt"
	
	#append a comma to the end of every line except the last line
	eval "sed '1,${l}s&root&root,&g' temp.txt > tempTwo.txt"

	#get rid of all new line characters, so that all files are listed on one line
	#sed works line by line, so it is cumbersome to use sed to replace newline characters
	eval "tr -d '\n' < tempTwo.txt > ${outputTxtFiles[$r]}"

	rm temp.txt tempTwo.txt

done


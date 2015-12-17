#!/bin/bash

cat > tmp/datasetToMarkdown.awk <<EOF
(NR==1){
  print "| Dataset                                      | Dataset name  | cross section | filter efficiency | number of events | new cross section | new cross section error | "
  print "| -------------------------------------------- | ------------- | ------------- | ----------------- | ---------------- | ----------------- | ----------------------- | "
}


(NF!=0){
  if(match(\$1,"#")){
    print "|                                              |               |               |                   |                  |   |  |"
  }else{
    print "| ", \$1, " | ", \$2, " | ", \$3, " | ", \$4, " | ", \$5, " | ", \$6, " | ", \$7, " |"
  }
}



EOF

awk -f tmp/datasetToMarkdown.awk configs/datasets.dat

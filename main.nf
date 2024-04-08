$HOSTNAME = ""
params.outdir = 'results'  

//* params.nproc =  20  //* @input @description:"How many processes to use for each step. Default 1"
params.mate="single"
//* params.projectDir = "${projectDir}"  //* @input @description:"How many processes to use for each step. Default 1"
params.metadata.metadata = "${params.projectDir}/tools.json"



if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 

if (params.reads){
Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_0_reads_g_15}
 } else {  
	g_0_reads_g_15 = Channel.empty()
 }

Channel.value(params.mate).into{g_1_mate_g_15;g_1_mate_g_16}


process unizp {

input:
 set val(name),file(reads) from g_0_reads_g_15
 val mate from g_1_mate_g_15

output:
 set val(name),file("*.fastq")  into g_15_reads0_g_16

script:

readArray = reads.toString().split(' ')	
R1 = readArray[0]
R2 = readArray[1]

"""
case "$R1" in
*.gz | *.tgz ) 
        gunzip -c $R1 > R1.fastq
        ;;
*)
        cp $R1 ./R1.fastq
        echo "$R1 not gzipped"
        ;;
esac

case "$R2" in
*.gz | *.tgz ) 
        gunzip -c $R2 > R2.fastq
        ;;
*)
        cp $R2 ./R2.fastq
        echo "$R2 not gzipped"
        ;;
esac
"""
}


process fatsq_to_fasta {

input:
 set val(name),  file(reads) from g_15_reads0_g_16
 val mate from g_1_mate_g_16

output:
 set val(name),  file("*fasta")  into g_16_airr_fasta_file0_g_10

#shell example: 

#!/bin/sh 



script:
	
readArray = reads.toString().split(' ')	

if(mate=="pair"){
	R1 = readArray.grep(~/.*R1.*/)[0]
	R2 = readArray.grep(~/.*R2.*/)[0]
	
	R1n = R1.replace('.fastq','')
	R2n = R2.replace('.fastq','')
	
	"""
	 awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}'  ${R1n}.fastq > ${R1n}.fasta
	 awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}'  ${R2n}.fastq > ${R2n}.fasta
	"""
	
}else{

	"""
	 awk 'NR%4==1{printf ">%s\n", substr($0,2)}NR%4==2{print}' ${name}.fastq > ${name}.fasta
	"""
}
}


process vdjbase_input {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${chain}$/) "reads/$filename"}
input:
 set val(name),file(reads) from g_16_airr_fasta_file0_g_10

output:
 file "${chain}"  into g_10_germlineDb00

script:
chain = params.vdjbase_input.chain



'''
#!/bin/sh 
mkdir ${chain}
mv ${reads} ${chain}
'''

}


process metadata {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.json$/) "metadata/$filename"}

output:
 file "*.json"  into g_12_jsonFile00

script:
metadata = params.metadata.metadata
"""
#!/usr/bin/env Rscript

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite")
}
library(jsonlite)

data <- read_json("${metadata}") 

versions <- lapply(1:length(data), function(i){
	
	docker <- data[i]
	tool <- names(data)[i]
	
	if(grepl("Custom", docker)){
		ver <- "0.0"
	}else{
		ver <- system(paste0(tool," --version"), intern = TRUE)
		ver <- gsub(paste0(tool,": "), "", ver)
	}
	ver
	
})

names(versions) <- names(data)

json_data <- list(
  sample = list(
    data_processing = list(
      preprocessing = list(
        software_versions = versions
	   )
	 )
  )
)

# Convert to JSON string without enclosing scalar values in arrays
json_string <- toJSON(json_data, pretty = TRUE, auto_unbox = TRUE)
print(json_string)
# Write the JSON string to a file
writeLines(json_string, "pre_processed_metadata.json")
"""

}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}

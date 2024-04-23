$HOSTNAME = ""
params.outdir = 'results'  

//* params.nproc =  20  //* @input @description:"How many processes to use for each step. Default 1"
params.mate="single"
params.projectDir="${projectDir}"

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

if(mate=="pair"){
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
}else{
	"""
	case "$reads" in
	*.gz | *.tgz ) 
	        gunzip -c $reads > R1.fastq
	        ;;
	*)
	        cp $reads ./R1.fastq
	        echo "$reads not gzipped"
	        ;;
	esac
	"""
}
}


process fatsq_to_fasta {

input:
 set val(name),  file(reads) from g_15_reads0_g_16
 val mate from g_1_mate_g_16

output:
 set val(name),  file("*fasta")  into g_16_airr_fasta_file0_g_10

script:

readArray = reads.toString().split(' ')	
if(mate=="pair"){
	R1 = readArray[0]
	R2 = readArray[0]
	
	R1n = R1.replace('.fastq','')
	R2n = R2.replace('.fastq','')
	
	"""
	 awk 'NR%4==1{printf ">%s\\n", substr(\$0,2)}NR%4==2{print}'  ${R1n}.fastq > ${R1n}.fasta
	 awk 'NR%4==1{printf ">%s\\n", substr(\$0,2)}NR%4==2{print}'  ${R2n}.fastq > ${R2n}.fasta
	"""
	
}else{

	"""
	 awk 'NR%4==1{printf ">%s\\n", substr(\$0,2)}NR%4==2{print}' ${reads} > ${name}.fasta
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

"""
mkdir ${chain}
mv ${reads} ${chain}
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

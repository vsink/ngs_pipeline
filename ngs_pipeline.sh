#!/bin/bash
# This script depends on: samtools, bwa, 
# It was tested under Linux (Arch Linux). Supporting non-Linux systems is unknown.
# It is used on 'as-is' basis. Authors are not liable for any damages arising from its use.
# Color text table

##################
# system options #
#################

bin="/usr/bin"
use_samtools_version="new" #new
samtools_old="/home/slava/work/tools/samtools-0.1.19/samtools"
samtools_new="/usr/bin/samtools"
bwa="/bin/bwa"
bowtie="/bin/bowtie"
bcftools="/home/slava/work/tools/samtools-0.1.19/bcftools/bcftools"

###############################################################################
bamsort=1
###############################################################################

c_std="[0;39m"
c_h_std="[1;37m"
c_pink="[0;35m"
c_h_pink="[1;35m"
c_red="[0;31m"
c_h_red="[1;31m"
c_cayan="[0;36m"
c_h_cayan="[1;36m"
c_yellow="[1;33m"
c_green="[0;32m"
c_h_green="[1;32m"
c_blue="[0;34m"
c_h_blue="[1;34m"
c_l_blue="[1;96m"

workdir=`pwd`


mapping="no"
sam2bam="no"
map_format="fastq"
snp_calling="no"
makefq="no"
enable_scripts="no"
overwrite_files="no"
delete_reads="no"

bam_unsort_ext=".unsorted.bam"
bam_sort_ext=".sort"
delete_unsorted_bam="yes"

mapper="bwa"
map_algorithm="mem" #map, samse, sampe
paired="no"
paired_ext1="_1"
paired_ext2="_2"


snp_caller="samtools"
snp_caller_ext=".vcf"


function help_me {
#echo "${c_h_blue}"
#echo "  _   _  _____  _____                _       _    "
#echo " | \ | |/ ____|/ ____|              (_)     | |   "
#echo " |  \| | |  __| (___   ___  ___ _ __ _ _ __ | |_  "
#echo " | .   | | |_ |\___ \ / __|/ __| '__| | '_ \| __| "
#echo " | |\  | |__| |____) |\__ \ (__| |  | | |_) | |_  "
#echo " |_| \_|\_____|_____/ |___/\___|_|  |_| .__/ \__| "
#echo "                  ______              | |         "
#echo "                 |______|             |_|         "
#echo "${c_std}"
#echo "version : 0.6b (build 150316)"

printf "Help\nKeys:\n-m\t\t: Map to reference\n-c\t\t: Sam2Bam\n-s\t\t: SNPcalling\n-p\t\t: Paired Reads\n-r FILE.FASTA\t: Reference File\n\n"
exit
}




while getopts "mcshf:r:hp" opt
do
case $opt in
#i) pileup_file=$OPTARG;;
m) mapping="yes";;
c) sam2bam="yes";;
s) snp_calling="yes";;
p) paired="yes";;
f) map_format=$OPTARG;;
r) reference=$OPTARG;;
h) help_me;;
esac
done



# echo $pileup_file
#reading options from pileup file
# if [ -f ${workdir}/$pileup_file ]
# then	
# 	source $pileup_file
# 	printf "Pileup file\t[OK]\n" 
# else
#     echo "${c_red}ERROR: Pileup file is not found. This is a critical error! Script will be aborted!${c_std}"
#     exit
# fi

#check system settings 

#check bin directory
if [ ! -d "$bin" ]; then
  # Control will enter here if $DIRECTORY exists.
  echo "${c_red}Ошибка: $bin папка на найдена. Программа будет остановлена!${c_std}"
  exit
fi


#check samtools 
if [[ -f $samtools_old && $use_samtools_version == "old" ]];then
  # Control will enter here if $DIRECTORY exists.
samtols_exist=1
printf "Samtools\t[OK] Старая версия!\n"
elif [[ -f $samtools_new && $use_samtools_version == "new" ]];then
  # Control will enter here if $DIRECTORY exists.
samtols_exist=1
printf "Samtools\t[OK] Новая версия!\n"

else
printf "Samtools не найдена\t\t[ОШИБКА]\n"
samtols_exist=0    
fi

#check bwa 
if [ -f $bwa ]; then
  # Control will enter here if $DIRECTORY exists.
bwa_exist=1    
printf "Bwa\t\t[OK]\n"
else
printf "Bwa не найдена\t\t[ОШИБКА]\n"
bwa_exist=0
fi

#check bowtie 
if [ -f $bowtie ]; then
  # Control will enter here if $DIRECTORY exists.
bowtie_exist=1    
printf "Bowtie\t\t[OK]\n"
else
printf "Bowtie не найдена\t\t[ОШИБКА]\n"
bowtie_exist=0
fi

#check bowtie 
if [ -f $bcftools ]; then
  # Control will enter here if $DIRECTORY exists.
bcftools_exist=1    
printf "bcftools\t[OK]\n"
else
printf "bcftools не найдена\t\t[ОШИБКА]\n"
bcftools_exist=0    
fi

# if [ -f $freebayes ]; then
#   # Control will enter here if $DIRECTORY exists.
# freebayes_exist=1    
# echo "Freebayes...[OK]"
# else
# echo "Freebayes not found"
# freebayes_exist=1    =0
# fi


#check java
# if [ -f ${bin}/java ]; then
#   # Control will enter here if $DIRECTORY exists.
# java_exists=1
# echo "Java...[OK]"    
# else
# echo "Java not found"
# java_exist=0
# fi

#check python2
# if [ -f ${bin}/python2 ]; then
#   # Control will enter here if $DIRECTORY exists.
# python_exists=1
# echo "Python2...[OK]"
# else
# echo "Python2 not found"
# python_exists=0
# fi

#check snpEff
#if [ -f $snp_eff ]; then
  # Control will enter here if $DIRECTORY exists.
#snpeff_exists=1
#echo "snpEff...[OK]"    
#else
#echo "snpEff not found"
#snpeff_exists=0
#fi

#check reference file
if [ -f ${workdir}/$reference ]
then
    printf "Референс\t[OK]\n" 
    reference_base=${reference%.*}   
else
    echo "${c_red}ОШИБКА: Референсный геном не найден. Программа будет остановлена!${c_std}"
    exit
fi



if [[ ! -f "${reference}.amb" && $bwa_exist == 1 ]];then
echo "Bwa словарь не найден и будет создан"
$bwa index -a bwtsw ${workdir}/$reference
fi

if [[ ! -f "${reference}.1.ebwt" && $bowtie_exist == 1 ]];then
echo "Bowtie словарь не найден и будет создан"
bowtie-build $reference $reference
fi

#screen pileup steps
printf "\n*** Пожалуйста, внимательно прочитайте нижеследующее!!! *** \n\n"

# check mapping step
if  [[ $mapping == "yes" ]]; then
#echo "# You choose $mapper as mapper"
map=1
else
map=0
fi

#check if map format is fastq and mapping is enabled
if  [[ $map_format == "fastq" && $map == 1 ]]; then
echo "-> Вы хотите картировать fastq файлы при помощи программы $mapper"
fastq=1
#sff=0
fi

#check if map format is sff and mapping is enabled
if  [[ $map_format == "sff" && $map == 1 ]]; then
#(The first sff files will be converted to fastq format  and after that analysis will continue with converted files. In working directory python script (sff2fastq.py) automaticaly will be created)"
echo "-> Вы хотите картировать sff файлы при помощи программы $mapper" 
#fastq=0
sff=1
fi

if ([[ $paired == "yes" ]]);then
printf "${c_h_red}Внимание! Вы указали в настройках использование paired reads!!! ${c_std}\n"
fi

#if sff2fastq.py file not found will create it. It depends on biopython!!
if  ([[ $sff == 1 && ! -f ${workdir}/sff2fastq.py ]]) && ([[ $sff == 1 &&  $map == 1 ]]); then
#echo "# sff2fastq.py doesnt exist and was created!"
cat > sff2fastq.py << EOL
#!/usr/bin/python
### Generates a fastq file from a 454 SFF file
import sys
from Bio import SeqIO
SFFfile = sys.argv[1]

SeqIO.convert(SFFfile, "sff-trim", SFFfile[:-4]+".fastq", "fastq")

EOL
fi

# check bam2 sam conversion step
if  [[ $sam2bam == "yes" ]]; then
echo "-> Вы хотите конвертировать SAM -> BAM формат"
bam=1
fi

# # check if you want to sort bam files
# if  [[ $bam_sort == "yes" && $bam == 1 ]]; then
# echo "-> After that you want to sort BAM file"
# bamsort=1
# fi

# check if you want to sort bam files
if  [[ $makefq == "yes" ]]; then
fqenable=1
echo "-> You want to make fq (consensus) file"
fi


# check if you want to sort bam files
if  ([[ $snp_calling == "yes" && $samtols_exist == 1 ]]); then
snpcall=1
echo "-> Вы желаете определить SNP при помощи $snp_caller"
fi

# check if you want to annotate vcf
#if  ([[ $use_snpeff == "yes" && $snpeff_exists == 1 ]]); then
#echo "-> You want to annotate vcf files by SNPMiner2?"
#annotate=1
#fi

# check if you want to annotate vcf
if  ([[ $enable_scripts == "yes" ]]); then
echo "-> You want to run scripts after analysis."
script=1
fi

# Procedures


function sff_conversion {
sff_file_exists=`ls -1 *.sff 2>/dev/null | wc -l`
if [[ ($sff_file_exists != 0) &&  ( $sff == 1) ]];then
	printf "\n"
	echo "${c_h_cayan}       **** SFF2FASTQ CONVERSION ****       ${c_std}"
	printf "\n"
    for i in ${workdir}/*.sff
		do        
		local sff_file=`basename -s .sff $i`		
		printf "${c_h_cayan}${sff_file^}....RUN${c_std}\n"		
		${bin}/python2 sff2fastq.py $sff_file.sff		
		printf "${c_h_cayan}${sff_file^}....DONE${c_std}\n"
		done
    
    fi

}

# Mapping

function bwa_mapping {
#Check files 
fastq_file_exists=`ls -1 *.fastq 2>/dev/null | wc -l`
    if [[ $fastq_file_exists != 0 ]];then
	#Run mapping		
	printf "${c_h_blue}\n\t<<<КАРТИРОВАНИЕ С ПОМОЩЬЮ BWA>>>${c_std}\n\n"
	
	if [[ $map_algorithm == "samse" ]];then
	for i in ${workdir}/*.fastq
		do        
		local file=`basename -s .fastq $i`		
		printf "${c_h_blue}${file}:\n\t\t-->${c_std}\n"
		${bin}/bwa aln $reference $file.fastq  > $file.sai;${bin}/bwa samse $reference -n 1 $file.sai $file.fastq  > $file.sam;rm $file.sai

		printf "${c_h_blue}\n\t<--${c_std}\n"
		done
	#elif [[ $map_algorithm == "sampe" ]];then
	
	elif [[ $map_algorithm == "mem" ]];then
	if ([[ $paired == "yes" ]]);then
	r1=()
	cnt=0
	cnt1=0
	for i in $(ls *$paired_ext1.fastq | uniq)
	do
	r1[$cnt]=$i
	cnt=$((cnt+1))
	done
	for i in $(ls *$paired_ext2.fastq | uniq)
	do
	local file=`basename -s .fastq $i|cut -d"_" -f1`
	printf "${c_h_blue}${file^}:\n\t-->${c_std}\n"
	echo ${r1[cnt1]} $i
	${bin}/bwa mem -t 2 -v 0 $reference ${r1[cnt1]} $i  > $file.sam 
	printf "${c_h_blue}\n\t<--${c_std}\n"
	cnt1=$((cnt1+1))
	done
	else
	for i in ${workdir}/*.fastq
		do        
		local file=`basename -s .fastq $i`	
		printf "${c_h_blue}${file^}:\n\t-->${c_std}\n"
		${bin}/bwa mem $reference $file.fastq  > $file.sam
		printf "${c_h_blue}${file^}:\n\t<--${c_std}\n"
		if [ $delete_reads == "yes" ]; then
			rm $file.fastq
		fi	
		done
			
	fi	
	else
	
	printf "Error! Not found any correct algorithm for BWA mapping! The programm will be stopped."
	
	fi
	#echo "${c_h_blue}Mapping of ${FILE}${c_std} ${c_h_red}[STARTED]${c_std}"
#	bwa mem " + @ref_file + " " + file + " > " + loc_sam_output
	#else
	#echo "No any fastq file is found"	
    fi    
} 

function bowtie_mapping {
#Check files 
fastq_file_exists=`ls -1 *.fastq 2>/dev/null | wc -l`
    if [[ $fastq_file_exists != 0 ]];then
	#Run mapping		
	printf "${c_h_blue}\n\t<<<КАРТИРОВАНИЕ С ПОМОЩЬЮ BOWTIE>>>\n\n${c_std}"
	
	for i in ${workdir}/*.fastq
		do        
		local file=`basename -s .fastq $i`		
		printf "${c_h_blue}${file^}:\n\t-->${c_std}\n"		
		${bin}/bowtie $reference -t -q ${file}.fastq  -S ${file}.sam		
		printf "${c_h_blue}\n\t<--${c_std}\n"
		if [ $delete_reads == "yes" ]; then
			rm $file.fastq
		fi	
		done
	#echo "${c_h_blue}Mapping of ${FILE}${c_std} ${c_h_red}[STARTED]${c_std}"
#	bwa mem " + @ref_file + " " + file + " > " + loc_sam_output
	#else
	#echo "No any fastq file is found"	
	
    fi    
} 

function sam_conversion {	
	sam_file_exists=`ls -1 *.sam 2>/dev/null | wc -l`
	  if [[ $sam_file_exists != 0 ]]; then		  
		  printf "${c_yellow}\n\t<<<КОНВЕРТИРУЕМ SAM-->BAM>>>\n\n${c_std}"
		  #sam2bam conversion
		  for i in ${workdir}/*.sam
		  do        
		  local file=`basename -s .sam $i`  
		  temp_file=$(mktemp)
		  
		  if [[ -f ${workdir}/$file$bam_sort_ext.bam && $overwrite_files == "no" ]]; then
                                printf "${c_yellow}$file$bam_sort_ext.bam...[ПРОПУЩЕН]${c_std}\n"
                                continue
                                fi
		  printf "${c_yellow}${file^}:\n\t-->${c_std}\n"
		  if [ $use_samtools_version == "old" ];then
		  echo -n "КОНВЕРТАЦИЯ..."
		  $samtools_old import $reference $file.sam $file$bam_unsort_ext				  
		  printf "OK\n"
		  bam_sorting
		  elif [ $use_samtools_version == "new" ];then
		  echo -n "КОНВЕРТАЦИЯ..."
		  $samtools_new fixmate -O bam $file.sam $file$bam_unsort_ext
		  printf "OK\n"
		  bam_sorting
		  fi
		  printf "${c_yellow}\n\t<--${c_std}\n"
		  if [[ $delete_sam == "yes" &&  -n $delete_sam ]]; then
			rm $file.sam
		  fi					  
		  
		  done
	
	fi	
}


function bam_sorting {
unsorted_bam_exists=`ls -1 *$bam_unsort_ext 2>/dev/null | wc -l`
		if [[ $unsorted_bam_exists != 0 ]];then		
		#bam sort
		
# 		printf "${c_yellow}\n\t<<<BAM SORTING>>>\n\n${c_std}"
		
		
			for i in ${workdir}/*$bam_unsort_ext
			do        			
				
				local bam_file=`basename -s $bam_unsort_ext $i`
				if [[ -f ${workdir}/$bam_file$bam_sort_ext.bam && $overwrite_files == "no" ]]; then
#                                 printf "${c_yellow}$bam_file$bam_sort_ext.bam\t[FOUND][PASSED]${c_std}"
                                continue
                                fi
# 				printf "${c_yellow}${bam_file^}:\n\t-->${c_std}\n"
				if [ $use_samtools_version == "old" ];then                                    
					$samtools_old sort $bam_file$bam_unsort_ext $bam_file$bam_sort_ext
					if [ $delete_unsorted_bam == "yes" ]; then
						rm $bam_file$bam_unsort_ext
					fi				
				elif [ $use_samtools_version == "new" ];then				
				echo -n "СОРТИРОВКА..."
				$samtools_new sort -O bam -o $bam_file$bam_sort_ext.bam -T  /tmp/$bam_file $bam_file$bam_unsort_ext				
				printf "OK\n"				
				echo -n "ПОСТРОЕНИЕ ИНДЕКСА..."
				$samtools_new index $bam_file$bam_sort_ext.bam
				printf "OK\n"				
# 					if [ $delete_unsorted_bam == "yes" ]; then
						rm $bam_file$bam_unsort_ext
				
# 					fi					  
				fi			
				#mv $ ${bam_file}_sort ${bam_file}1.bam
# 				printf "${c_yellow}\t<--${c_std}\n"
			done
		fi

}

function snp_calling {
sorted_bam_exists=`ls -1 *$bam_sort_ext.bam 2>/dev/null | wc -l`

	if [[ $snp_caller == "samtools"  && $sorted_bam_exists != 0 ]];then
		#snp_calling
		
		printf "${c_h_green}\n\t<<<SNP CALLING>>>\n\n${c_std}"
		
			for i in ${workdir}/*$bam_sort_ext.bam
			do        
			
				local bam_file=`basename -s $bam_sort_ext.bam $i`   						
				
			#samtools mpileup -f $reference  $bam_file$bam_sort_ext.bam > $bam_file$pileup_ext
			if [[ -f ${workdir}/$bam_file$snp_caller_ext && $overwrite_files == "no" ]]; then
			printf "${c_h_green}$bam_file$snp_caller_ext...[ПРОПУЩЕН]${c_std}\n"
                        continue
                        else
                        printf "${c_h_green}${bam_file^}:\n\t-->${c_std}\n"
                        fi
			if [ $use_samtools_version == "old" ];then
			$samtools_old mpileup -uf $reference $bam_file$bam_sort_ext.bam | $bcftools view -vcg - > $bam_file$snp_caller_ext		
			elif [ $use_samtools_version == "new" ];then				
			$samtools_new mpileup -ugf $reference $bam_file$bam_sort_ext.bam | bcftools call -vmo $bam_file$snp_caller_ext
			#tabix -p vcf $bam_file.vcf.gz
			fi
				#printf "${c_h_green}${bam_file^}\t<--${c_std}\n"
				printf "${c_h_green}\t<--${c_std}\n"
			done
		
		fi
	
	#${BIN_DIR}samtools mpileup -f ${REF_DIR}${REF_FILE} ${FILE}_sort.bam > ${FILE}.pileup
	#fi
    #echo "${c_h_blue}SNP CALLING OF ${FILE}${c_std} BY SAMTOOLS ${c_h_red}[STARTED]${c_std}"
    #${BIN_DIR}samtools mpileup -uf ${REF_DIR}${REF_FILE} ${FILE}_sort.bam | bcftools view -vcg - > ${FILE}.st.vcf
#	fi


}


function making_fq {
sorted_bam_exists=`ls -1 *$bam_sort_ext.bam 2>/dev/null | wc -l`


	if [[ $sorted_bam_exists != 0 ]];then
		#snp_calling
		printf "\n"
		echo "${c_l_blue}       **** CONSENSUS (FQ) CREATING ****       ${c_std}"
		printf "\n"
			for i in ${workdir}/*$bam_sort_ext.bam
			do        			
			#samtools mpileup -uf ref.fa readset_ref_bwa.bam | bcftools view -cg - | vcfutils.pl vcf2fq > readset_ref_bwa_cons.fq
				local bam_file=`basename -s $bam_sort_ext.bam $i`    				
				if [[ -f ${workdir}/$bam_file.fq && $overwrite_files == "no" ]]; then
			printf "${c_l_blue}$bam_file.fq....PASSED${c_std}\n"
                        continue
                        elif [[ -f ${workdir}/$bam_file.fasta && $overwrite_files == "no" ]]; then
			echo "$bam_file.fasta....PASSED"
                        continue
                        fi
				printf "${c_l_blue}${bam_file^}:\n\t-->${c_std}\n"
			#samtools mpileup -f $reference  $bam_file$bam_sort_ext.bam > $bam_file$pileup_ext
			$samtools_old mpileup -uf $reference $bam_file$bam_sort_ext.bam | $bcftools view -cg - | vcfutils.pl vcf2fq > $bam_file.fq
				printf "${c_l_blue}\n\t<--${c_std}\n"
			done		
				
		fi
	
	#${BIN_DIR}samtools mpileup -f ${REF_DIR}${REF_FILE} ${FILE}_sort.bam > ${FILE}.pileup
	#fi
    #echo "${c_h_blue}SNP CALLING OF ${FILE}${c_std} BY SAMTOOLS ${c_h_red}[STARTED]${c_std}"
    #${BIN_DIR}samtools mpileup -uf ${REF_DIR}${REF_FILE} ${FILE}_sort.bam | bcftools view -vcg - > ${FILE}.st.vcf
#	fi


}


# function snpeff_annotation {
# vcf_exists=`ls -1 *$snp_caller_ext 2>/dev/null | wc -l`
# 	if [[ $vcf_exists != 0 ]];then
# 	printf "\n"
# 	echo "${c_pink}       **** SNPEF ANNOTATION ****       ${c_std}"
# 	printf "\n"
# 			for i in ${workdir}/*$snp_caller_ext
# 			do        
# 				#local vcf_file=`basename -s $snp_caller_ext $i`  
# 				#local vcf_file= ${i%.*}
# 				
# 				t=${i#*.}
# 				if [ ! ${t%$snpeff_output_ext} == "eff.vcf" ]; then				
# 				local vcf_file=`basename -s $snp_caller_ext $i` 	
# 				else
# 				continue
# 				fi
# 				
# 				if [[ -f ${workdir}/$vcf_file$snpeff_output_ext && $overwrite_files == "no" ]]; then
# 			printf "${c_pink}$vcf_file$snpeff_output_ext....PASSED${c_std}\n"
#                         continue
#                         else
#                         printf "${c_pink}${vcf_file^}....RUN${c_std}\n"
#                         fi
# 				
# 			#samtools mpileup -f $reference  $bam_file$bam_sort_ext.bam > $bam_file$pileup_ext
# 			#samtools mpileup -uf $reference $bam_file$bam_sort_ext.bam | bcftools view -vcg - > $bam_file$snp_caller_ext	
# 			if [[ $snpeff_rename_chromosome == "yes" ]];then
# 			sed -i "s/"$snpeff_original_name"/"$snpeff_new_name"/g"  $vcf_file$snp_caller_ext	
# 			printf "%s replaced by: %s\n" $snpeff_original_name $snpeff_new_name
# 			fi
# 			${bin}/java -jar  $snpeff -c $snpeff_config $snpeff_genome_name -v $vcf_file$snp_caller_ext $snpeff_filter > $vcf_file$snpeff_output_ext
# 				printf "${c_pink}${vcf_file^}....DONE${c_std}\n"
# 			done
# 	
# 	fi
# 
# 
# 
# 
# }

function scripts_executing {
n=$script_count
str="cmd_script"
for ((i=1; i<=n; i++)) {
   eval \$$str$i
   
}
}


GLOBAL_TIME_START=$(date +"%s") # variable of time of starting analysis

#prompt to start analysis
read -p "Начинаем анализ? [yn]" answer
if [[ $answer = y ]];then
  # run the command
  #echo "Analysis started"  
  if [[ $sff == 1 ]];then
	sff_conversion
  fi
  if [[ $map == 1 && $mapper == 'bwa' ]];then
	bwa_mapping	
  fi
  if [[ $map == 1 && $mapper == 'bowtie' ]];then
	bowtie_mapping	
  fi
  if [[ $bam == 1 ]];then
	sam_conversion
  fi 
#   if [[ $bamsort == 1 ]];then
# 	bam_sorting
#   fi
  
  if [[ $fqenable == 1 ]];then
	making_fq
  fi
  
  if [[ $snpcall == 1 ]];then
	snp_calling
  fi
#   if [[ $annotate == 1 ]];then
# 	snpeff_annotation
#   fi  
  if [[ $script == 1 ]];then
	scripts_executing
  fi
  
fi


GLOBAL_TIME_END=$(date +"%s")
GLOBAL_TIME_DIFF=$(($GLOBAL_TIME_END-$GLOBAL_TIME_START))
printf "\n"
echo -e "-----------------------------------\n"$(($GLOBAL_TIME_DIFF / 60)) minutes and $(($GLOBAL_TIME_DIFF % 60)) seconds elapsed."\n-----------------------------------\n"
echo -e ${c_h_std}'(c)Viacheslav V. Sinkov, 2014-2016 (vsinkov_at_gmail.com)'${c_std}
printf "\n"

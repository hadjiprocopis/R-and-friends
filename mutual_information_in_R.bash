#!/bin/bash

# program by Andreas Hadjiprocopis
# livantes at soi.city.ac.uk
# Institute of Cancer Research
# London, 2010
# Free for everyone to copy, use, modify / quote author

appname=`basename "$0"`
bindir=$( cd "$( dirname "$0" )" && pwd )
now=`date`
rand=$$
param_ifs=':'
seed=$$
png_width=1024
png_height=768
png_fontsize=18

#extra
num_bins=100

function display_help_exit {
	echo "Usage : $1 -i input_file -o output_file [-c col_num ... -c ... | -C col_num ... -C ...] [-z col_name ... -z ... | -Z col_name ... -Z ...] [-1] [-e input_data_separator] [-E output_data_separator] [-R R_output_script_name] [-2 seed] [-3 png_width:png_height:png_fontsize]"
	echo "other params: [-b num_bins] "
	echo "This script will calculate the mutual information between all possible pairs of the specified columns-to-do of the input file. Unfortunately, the number of elements in each column must be the same."
	echo "The way to do that is 1) standardize (i.e. zero mean, unit stdev) each column, 2) put all the standardize columns together, 3) calculate the histogram of all the data together, with specified number of breaks (or bins) using -b (default is ${num_bins}), 4) label each value on each column according to which bin in the histogram it falls, e.g. if it falls to first bin, its label is 1, second bin, 2 etc., 5) calculate the mutual information between each pair of columns of labels"
	echo "The output file can then be plotten using the 'heatmap_in_R.bash -i inp -o out' script"
	echo "fixed params"
	echo " -i inputFile"
	echo "	defines the input file."
	echo " -o outputFile"
	echo "	defines the output file."
	echo " -c col_num"
	echo "	choose the column number you wish to operate as a comma-separated list of column numbers (starting from 1) - repeat -c options if necessary, all columns will be processed unless this option is specified"
	echo " -C col_num"
	echo "	choose the column number you DO NOT wish to operate as a comma-separated list of column numbers (starting from 1) - repeate -C options if necessary, all columns will be processed unless this option is specified"
	echo " -z col_name"
	echo "	choose the column name you wish to operate - repeat -z options if necessary, all columns will be processed unless this option is specified. Use unique separator (-1), default is '${param_ifs}' or repeat -z options"
	echo " -Z col_nam"
	echo "	choose the column number you DO NOT wish to operate as a list of column names or repeate -Z options. Use unique separator (-1), default is '${param_ifs}' or repeat -Z options"
	echo " -1"
	echo "	include columns not selected for processing to the output, default is not to print un-processed columns"
	echo " -2 random_number_generator_seed"
	echo "	useful to repeat experiments, otherwise seed is selected at random using the process id"
	echo " -3 png_width:png_height:png_fontsize : the size of the output image if any in pixels and font size, default: ${png_width}:${png_height}:${png_fontsize}"
# standard features
	echo " -h"
	echo "	optional, that the first row is not numerical but it holds column names, these names will be copied to the output. If the script fails, this could be the reason (i.e. calculations with column names)."
	echo " -e sep"
	echo "	      optional, 'sep' is the separator between columns in the input file. If the script fails, this could be the reason (i.e. spaces are used as col separator but in actual terms they are tab separated)."
	echo " -R filename"
	echo "	      specify the name of the file to contain the R script to execute (because this bash file"
	echo "	      creates an R script which then runs. If this option is not used, then the output script name will be deleted."
	echo "	Otherwise, the R script will not be deleted and you can  re-run it again using R < script_name"
	echo ""
	echo "Other options:"
	echo "program by Andreas Hadjiprocopis"
	echo "andreashad2 at gmail.com"
	echo "Institute of Cancer Research"
	echo "London, 2010"
	echo "Free for everyone to copy, use, modify / quote author"
	echo ""
	exit 1
}

function column_names_to_column_numbers {
	local   _inputfile="$1"; shift # must have header
	local   _ifs="$1"; shift
	local   _column_names=${@}
	local   _rand=$$
	local   _i _ret="" _ifs2=""
	if [ "${_ifs}" != "" ]; then _ifs2="--ifs '${_ifs}'"; fi
	# dont echo anything apart from the last one before return (unless > /dev/stderr)
#       echo extract_columns.pl --input "${_inputfile}" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit > /dev/stderr
cat << EOC | sh
	extract_columns.pl --input "${_inputfile}" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit > /dev/null
EOC
	if [ $? -ne 0 ]; then
		echo "$appname : ${FUNCNAME}, line ${LINENO} : call to the following command has failed"
		echo extract_columns.pl --input \"${_inputfile}\" --output /tmp/${_rand}.cn2cn --include_column_names_in_output --have_column_names_header ${_ifs2} --print_columns_and_exit
		return 1
	fi
	oldIFS="${IFS}"; IFS=${param_ifs}
	for _a_col_name in ${_column_names[*]}; do
#	      echo "checking for ${_a_col_name} (/tmp/${_rand}.cn2cn)" > /dev/stderr
		_a_col=`awk -F '\t' -v n="${_a_col_name}" '{if($4==n){print ($2+1); exit 2}}' /tmp/${_rand}.cn2cn`
		if [ $? -ne 2 ]; then
#		       echo "$appname : ${FUNCNAME}, line ${LINENO} : could not find column name '${_a_col_name}' in the header of the input file '${_inputfile}'."
#		       return 1
			continue # baby
		fi
#	       echo "found ${_a_col}" > /dev/stderr
		_ret+=",${_a_col}"
	done
	IFS="${oldIFS}"
	echo "${_ret#,}"
	rm -f /tmp/${_rand}.cn2cn
	return 0
}
inputFile=""
outputFile=""
# for standard arguments
separator=''; Rseparator=''
separatorOut='\t'
header="FALSE"
Rscript_name="${rand}.R"
outputPrecision="5"
doColumns=""
dontDoColumns=""
doColumnsNames=""
dontDoColumnsNames=""
do_only_selected_columns=0
include_discarded_columns_in_output=0
delete_Rscript_when_done=1

if [ "$*" = "" ]; then
	display_help_exit $appname
fi

while getopts "c:C:z:Z:i:o:R:Hhe:E:p:12:3:b:" OPTION
do
	case $OPTION in
		c)
			if [ "${doColumns}" == "" ]; then doColumns="${OPTARG}"; else doColumns="${doColumns},$OPTARG"; fi
			do_only_selected_columns=1
			;;
		C)
			if [ "${dontDoColumns}" == "" ]; then dontDoColumns="${OPTARG}"; else dontDoColumns="${dontDoColumns},$OPTARG"; fi
			do_only_selected_columns=1
			;;
		z)
			if [ "${doColumnsNames}" == "" ]; then doColumnsNames="${OPTARG}"; else doColumnsNames="${doColumnsNames}${param_ifs}$OPTARG"; fi
			do_only_selected_columns=1
			;;
		Z)
			if [ "${dontDoColumnsNames}" == "" ]; then dontDoColumnsNames="${OPTARG}"; else dontDoColumnsNames="${dontDoColumnsNames}${param_ifs}$OPTARG"; fi
			do_only_selected_columns=1
			;;
		1)
			include_discarded_columns_in_output=1
			;;
		i)
			if [ "$inputFile" != "" ]; then
				echo "$appname : you can specify only one input file"
				exit 1
			fi
			inputFile="$OPTARG"
			;;
		o)
			if [ "$outputFile" != "" ]; then
				echo "$appname : you can specify only one output folder"
				exit 1
			fi
			outputFile=$OPTARG
			;;
		# standard arguments
		2)
			seed=$OPTARG
			;;
		3)
			png_width=`echo "${OPTARG}" | cut -d':' -f1`
			png_height=`echo "${OPTARG}" | cut -d':' -f2`
			png_fontsize=`echo "${OPTARG}" | cut -d':' -f3`
			;;
		R)
			delete_Rscript_when_done=0
			Rscript_name="$OPTARG"
			;;
		H)
			display_help_exit "$appname"
			;;
		h)
			# in case there is a header (column names at the first row of file)
			header="TRUE"
			;;
		e)
			separator="$OPTARG"; Rseparator=",sep='${separator}'"
			;;
		E)
			separatorOut=$OPTARG
			;;
		p)
			outputPrecision=$OPTARG
			;;
		# other options
		b)
			num_bins=$OPTARG
			;;
	esac
done

if [ ! -f "$inputFile" ]; then
	echo "$0 : input file $inputFile does not exist."
	exit 1
fi
if [ "$Rscript_name" = "" ]; then
	echo "$0 : script file name is not defined (use the -R option to do it)."
	exit 1
fi
if [ "$outputFile" = "" ]; then
	echo "$0 : you need to specify an output file using the '-o' option."
	exit 1
fi
if [[ ("${doColumns}" != "" || "${doColumnsNames}" != "") && ("${dontDoColumns}" != "" || "${dontDoColumnsNames}" != "") ]]; then
	echo "$0 : either -c or -C options but not both can be used (either specify columns to do, or columns not to do - both makes no sense)"
	exit 1
fi

if [ "${doColumnsNames}" != "" ]; then
	crap=$(column_names_to_column_numbers "${inputFile}" "${separator}" ${doColumnsNames})
	if [ $? -ne 0 ]; then echo ${crap}; exit 1; fi
	if [ "${crap}" != "" ]; then
		doColumns+="${crap}"
		echo "$appname : extra column numbers to do via their column names : ${crap}"
	fi
fi
if [ "${dontDoColumnsNames}" != "" ]; then
	crap=$(column_names_to_column_numbers "${inputFile}" "${separator}" ${dontDoColumnsNames})
	if [ $? -ne 0 ]; then echo ${crap}; exit 1; fi
	if [ "${crap}" != "" ]; then
		dontDoColumns+="${crap}"
		echo "$appname : extra column numbers not to do via their column names : ${crap}"
	fi
fi
doColumns=${doColumns/X,/}; doColumns=${doColumns##,}; dontDoColumns=${dontDoColumns/X,/}; dontDoColumns=${dontDoColumns##,};

cat << EOC > "$Rscript_name"
set.seed($seed); write("$appname : starting R, seed is $seed", stderr());

# R script by $appname on $now
write("$appname : reading data from '$inputFile'", stderr()); gotError<-FALSE
# fill=T allows for uneven number of rows per column
tryCatch({oridata = read.table("$inputFile", fill=T, header=$header $Rseparator, check.names=FALSE)}, error=function(ex){write(paste("$appname : error reading file $inputFile\n",ex,sep=""),stderr()); gotError<-TRUE});
if( gotError == TRUE ){ quit(status=1) }
write("$appname : done.", stderr())
actual_num_cols=ncol(oridata)
if( "${dontDoColumns}" == "" ){
	if( "${doColumns}" == "" ){
		# no columns specified, do all of them
		columns_to_do = seq(from=1, to=actual_num_cols)
		must_do_this_column=rep(TRUE, actual_num_cols)
	} else { columns_to_do = c($doColumns); must_do_this_column=rep(FALSE, actual_num_cols); must_do_this_column[columns_to_do]=TRUE }
} else {
	must_do_this_column=rep(FALSE, actual_num_cols)
	crap=c($dontDoColumns); num_crap = length(crap); columns_to_do=c();
	for(i in 1:actual_num_cols){
		found = 0;
		for(j in 1:num_crap){ if( crap[j] == i ){ found=1; break; } }
		if( found == 0 ){ columns_to_do=append(columns_to_do, i) }
	}
	must_do_this_column[columns_to_do]=TRUE
}
num_columns_to_do = length(columns_to_do); num_rows = nrow(oridata)
if( num_rows==0 ){ write("$appname : no rows left to do", stderr()) }
if( num_columns_to_do==0 ){ write("$appname : no columns left to do", stderr()) }
if( (num_columns_to_do==0) || (num_rows==0) ){ quit(status=0) }

source("${bindir}/library_in_R.R")
library(infotheo)

# this is where the actual program starts
data_out=NULL
labeled_data=list(); standardised_data=list()
mutual_information=matrix(ncol=num_columns_to_do,nrow=num_columns_to_do)
# make a histogram OF ALL THE DATA after standardisation -- outliers are assumed to be removed before running this program
for(in_col_id1 in 1:num_columns_to_do){
#	stdev = sd(oridata[,in_col_id1])
#	if( stdev <= ${min_stdev_of_column_to_include} )
	a_col_id = columns_to_do[in_col_id1]
	tmpdat = oridata[,a_col_id]; tmpdat=tmpdat[!is.na(tmpdat)]
	standardised_data[[in_col_id1]] = (tmpdat-mean(tmpdat)) / sd(tmpdat)
}
total_data=unlist(standardised_data)
total_histogram=hist(total_data, breaks=${num_bins}, plot=F); total_data=NULL
for(in_col_id1 in 1:num_columns_to_do){
	labeled_data[[in_col_id1]] = label_data(total_histogram\$breaks, standardised_data[[in_col_id1]])
}
cat(paste("$appname : calculating ", num_columns_to_do, " columns :", sep=''), file=stderr())
for(in_col_id1 in 1:num_columns_to_do){
	a_col1 = columns_to_do[in_col_id1]
	cat(" ", a_col1, file=stderr())
	for(in_col_id2 in 1:num_columns_to_do){
		if( in_col_id1 == in_col_id2 ){ next }
		a_col2 = columns_to_do[in_col_id2]
		mutual_information[in_col_id1,in_col_id2] = mutinformation(labeled_data[[in_col_id1]], labeled_data[[in_col_id2]])
	}
}
write(". Done all.",stderr())

# precision and when this precision is not enough
options(digits=$outputPrecision, scipen=100)
if( "$header" == "TRUE" ){ namr=names(oridata[,columns_to_do]) } else { namr=paste("col:",columns_to_do,sep="") }; namc=c("AA", namr)
write.table(cbind(namr,mutual_information), file="$outputFile", sep="$separatorOut", col.names=namc, row.names=F, append=FALSE, quote=FALSE)
write("$appname : done with precision $outputPrecision", stderr())
EOC

echo "$appname : executing R script ($Rscript_name)"
R --vanilla < "$Rscript_name" > /tmp/rlog$rand.txt
if [ $? -ne 0 ]; then
#	cat /tmp/rlog$rand.txt; echo "$appname : error executing R script ($Rscript_name)";
	echo "Command line:"
	echo "$0 ${*}"
	exit 1
fi
echo "Command line:"
echo "$0 ${*}"
rm -f /tmp/rlog$rand.txt /tmp/$rand.* 
if [ "${delete_Rscript_when_done}" -eq 1 ]; then rm -f "$Rscript_name"; fi
exit 0

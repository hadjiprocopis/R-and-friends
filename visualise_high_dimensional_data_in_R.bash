#!/bin/bash

# program by Andreas Hadjiprocopis
# livantes at soi.city.ac.uk
# Institute of Cancer Research
# London, 2010
# Free for everyone to copy, use, modify / quote author

chr() {
  printf \\$(printf '%03o' $1)
}

ord() {
  printf '%d' "'$1"
}


function display_help_exit {
	echo "Usage : $1 -i input_file -o output_file [-h] [-e input_separator] [-E output_separator] [-S seed] [-s shape_column] [-C colors_column]"
	echo "This script will read a N-dimensional (=N columns) data from input_file and cluster it"
	echo " -i input_file"
	echo "        defines the input file."
	echo " -o output_file"
	echo "		where the results go"
	echo " -S seed"
	echo "          seed to the random number generator"
	echo " -c col"
	echo "	      choose the column number you wish to fit (use -c option repeatedly for more columns)"
	echo " -h"
	echo "        optional, that the first row is not numerical but it holds column names"
	echo " -e sep"
	echo "	      optional, 'sep' is the separator between columns in the input file, default is white space"
	echo "        which is '' (i.e. nothing in between the quotes)"
	echo " -R filename"
	echo "	      specify the name of the file to contain the R script to execute (because this bash file"
	echo "	      creates an R script which then runs)"
	echo "	      Example:"
	echo "$app_name -i a.txt -o abc -h"
	echo "read column 16 from 'a.txt' and write to output basename 'abc'"
	echo "The distribution to test is 'gamma' and the number of bins when"
	echo "constructing histograms is 5"
	echo "-h means that the input file's first row is not numerical but it consists of row names "
	echo ""
	echo "program by Andreas Hadjiprocopis"
	echo "andreashad2 at gmail.com"
	echo "Institute of Cancer Research"
	echo "London, 2010"
	echo "Free for everyone to copy, use, modify / quote author"
	echo ""
	exit 1
}

input_file=""
output_file=""
col1=""
separator="\t"
separatorOut="\t"
shapes_column=0
colors_column=0
header="FALSE"
app_name="$0"
rand=$$
scaleData=""
dimensions=2
Rscript_name="${rand}.R"
limits=""
#Seed="${rand}"
Seed=$$

if [ "$*" = "" ]; then
	display_help_exit $app_name
fi

app_name=`basename $0`
now=`date`

while getopts "i:o:s:c:hR:e:E:HS:d:C:S:" OPTION
do
	case $OPTION in
		i)
			input_file="$OPTARG"
			;;
		o)
			output_file=$OPTARG
			;;
		c)
			col1="$col1,$OPTARG"
			;;
		S)
			shapes_column="$OPTARG"
			;;
		C)
			colors_column="$OPTARG"
			;;
		H)
			display_help_exit
			;;
		h)
			# in case there is a header (column names at the first row of file)
			header="TRUE"
			;;
		e)
			separator="$OPTARG"
			;;
		E)
			separatorOut="$OPTARG"
			;;
		d)
			dimensions="$OPTARG"
			if [ ${dimensions} -ne 2 ] && [ ${dimensions} -ne 3 ]; then
				echo "$0 : number of dimensions must be 2 or 3 and not $dimensions"
				exit 1
			fi
			;;
	esac
done

if [ ! -f "$input_file" ]; then
	echo "$0 : input file $input_file does not exist."
	exit 1
fi
if [ "$Rscript_name" = "" ]; then
	echo "$0 : script file name is not defined (use the -R option to do it)."
	exit 1
fi
if [ "$output_file" = "" ]; then
	echo "$0 : you need to specify an output file using the '-o' option."
	exit 1
fi

if [[ "$output_file" =~ \.pdf|\.PDF$ ]]; then device_open_command="pdf(file=\"${output_file}\")";
elif [[ "$output_file" =~ \.png|\.PNG$ ]]; then device_open_command="png(filename=\"${output_file}\")";
else
	echo "$0 : don't know the file format of output file ${output_file}, I only know how to write PDF and PNG and therefore output filenames must end in one of these extensions."
	exit 1
fi

col1=${col1/,/}
if [ ${Seed} = "" ]; then Seed="${rand}"; fi

# first do a pca excluding the columns defined for shape and color
if [ $header == "TRUE" ]; then sed '1d' "${input_file}" > /tmp/$rand.1; head -1 "${input_file}" > /tmp/$rand.header; else cp "${input_file}"  /tmp/$rand.1; fi
num_dimensions_aux=0
rm -f /tmp/$rand.shapes  /tmp/$rand.unique_shapes
shapes_load_command=""; shapes_use_command=""
if [ ${shapes_column} -gt 0 ]; then
	awk -v sc=${shapes_column} '{print $sc}' /tmp/$rand.1 > /tmp/$rand.shapes
	if [ "${header}" == "TRUE" ]; then header_str=`awk -v sc=${shapes_column} '{print $sc}' /tmp/$rand.header`; else header_str="shapes"; fi
	sort -u /tmp/$rand.shapes > /tmp/$rand.unique_shapes
	ns=`awk 'END{print NR}' /tmp/$rand.unique_shapes`
	all_shapes_to_str=""; all_shapes_from_str=""
	if [ ${ns} -gt 18 ]; then
		echo "$0 : R supports only a maximum of 18 different shapes in plottingyour shape column (# ${shapes_column}) has $ns distinct values, will resort to the ASCII table."
		rs=0; declare -a chars=(`perl -e 'print join(" ", (a..z))." ".join(" ", (A..Z))." ".join(" ", (0..9))."\n";'`)
		if [ ${ns} -gt ${#chars[@]} ]; then
			echo "$0 : the ASCII table allows only ${#chars[@]} different readable characters, your shape column (# ${shapes_column}) has $ns distinct values, can not do that."
			exit 1
		else
			while read line; do
				gsed -e "s/\<$line\>/__${chars[rs]}__/g" /tmp/$rand.shapes > /tmp/$rand.2 && mv /tmp/$rand.2 /tmp/$rand.shapes
				all_shapes_to_str+="'""${chars[rs]},""'"; all_shapes_from_str+="$line,"
				let rs++
			done < /tmp/$rand.unique_shapes
			shapes_load_command="shapes=scan(\"/tmp/$rand.shapes\", what='character')"
		fi			
	else
		rs=1
		while read line; do
			gsed -e "s/\<$line\>/__${rs}__/g" /tmp/$rand.shapes > /tmp/$rand.2 && mv /tmp/$rand.2 /tmp/$rand.shapes
			all_shapes_to_str+="$rs,"; all_shapes_from_str+="$line,"
			let rs++
		done < /tmp/$rand.unique_shapes
		shapes_load_command="shapes=scan(\"/tmp/$rand.shapes\")"
	fi
	sed -e 's/__//g' /tmp/$rand.shapes > /tmp/$rand.2 && mv /tmp/$rand.2 /tmp/$rand.shapes
#	echo "__Shapes__" > /tmp/$rand.2; cat  /tmp/$rand.2 /tmp/$rand.shapes > /tmp/$rand.3 ; mv /tmp/$rand.3 /tmp/$rand.shapes
#	shapes_load_command="shapes=read.table(\"/tmp/$rand.shapes\", header=F)"
	shapes_use_command=",pch=shapes"
	all_shapes_to_str=${all_shapes_to_str/%,/}; all_shapes_from_str=${all_shapes_from_str/%,/}
	shapes_legend_command="legend(\"topleft\",,title=""'""${header_str}""'"",c(${all_shapes_from_str}),pch=c(${all_shapes_to_str}))"
	let num_dimensions_aux++;
fi
rm -f /tmp/$rand.colors  /tmp/$rand.unique_colors
colors_load_command=""; colors_use_command=""
if [ ${colors_column} -gt 0 ]; then
	awk -v sc=${colors_column} '{print $sc}' /tmp/$rand.1 > /tmp/$rand.colors
	if [ "${header}" == "TRUE" ]; then header_str=`awk -v sc=${colors_column} '{print $sc}' /tmp/$rand.header`; else header_str="colors"; fi
	sort -u /tmp/$rand.colors > /tmp/$rand.unique_colors
	ns=`awk 'END{print NR}' /tmp/$rand.unique_colors`
	all_colors_to_str=""; all_colors_from_str=""
	# R supports 657 colors, first is white so omit that, but anyway is quite stupid
	if [ ${ns} -gt 656 ]; then
		echo "$0 : R supports only a maximum of 657-1 different colors in plotting your color column (# ${colors_column}) has $ns distinct values, can not do that."
		exit 1
	else
		rs=1
		while read line; do
			gsed -e "s/\<$line\>/__${rs}__/g" /tmp/$rand.colors > /tmp/$rand.2 && mv /tmp/$rand.2 /tmp/$rand.colors
			all_colors_to_str+="${rs},"; all_colors_from_str+="${line},"
			let rs++
		done < /tmp/$rand.unique_colors
		sed -e 's/__//g'  /tmp/$rand.colors > /tmp/$rand.2 && mv /tmp/$rand.2 /tmp/$rand.colors
	fi	
#	echo "__Colors__" > /tmp/$rand.2; cat  /tmp/$rand.2 /tmp/$rand.colors > /tmp/$rand.3 ; mv /tmp/$rand.3 /tmp/$rand.colors
	colors_load_command="colors=read.table(\"/tmp/$rand.colors\", header=F)"
	colors_use_command=",col=colors[,1]"
	all_colors_to_str=${all_colors_to_str/%,/}; all_colors_from_str=${all_colors_from_str/%,/}
	colors_legend_command="legend(\"topright\",,title=""'""${header_str}""'"",c(${all_colors_from_str}),c(${all_colors_to_str}))"
	let num_dimensions_aux++;
fi
awk -v sc=${shapes_column} -v cc=${colors_column} '{for(i=1;i<=NF;i++){if((i==sc)||(i==cc))continue; printf $i" "}printf "\n"}' /tmp/$rand.1 > /tmp/$rand.2
awk -v sc=${shapes_column} -v cc=${colors_column} '{for(i=1;i<=NF;i++){if((i==sc)||(i==cc))continue; printf $i" "}printf "\n"}' /tmp/$rand.header > /tmp/$rand.small_header
num_dims=`awk '{print NF; exit}' /tmp/$rand.2`; num_dims=$((num_dims-num_dimensions_aux))
if [ "${num_dims}" -gt 2 ]; then
	echo "$0 : executing : " pca_in_R.bash -i /tmp/$rand.2 -o /tmp/$rand.pca -l /tmp/$rand.pca.loadings
	pca_in_R.bash -i /tmp/$rand.2 -o /tmp/$rand.pca -l /tmp/$rand.pca.loadings
	if [ $? -ne 0 ]; then
		echo "$0 : the following command has failed"
		echo pca_in_R.bash -i /tmp/$rand.2 -o /tmp/$rand.pca -l /tmp/$rand.pca.loadings
		exit 1
	fi
	awk -v d=${dimensions} -v sep="${separator}" '{for(i=1;i<=d;i++)printf $i""sep; printf "\n"}' /tmp/$rand.pca > /tmp/$rand.4
else cp /tmp/$rand.2 /tmp/$rand.4; fi

cat << EOC > "$Rscript_name"
# R script by $app_name on $now
library(cluster)
write(paste("$app_name : random number generator seed : ", $Seed, sep=""), file="$output_file", append=FALSE);
set.seed($Seed)

tryCatch({ data = read.table("/tmp/$rand.4", sep="$separator", header=F)}, error=function(ex){write(paste("$app_name : error reading file $input_file\n",ex,sep=""),stderr())});
dataAllDimensions = as.matrix(data)
dataSelectedColumn = as.matrix(data[,$col1])
num_rows_data = length(dataSelectedColumn)
dataIn = dataSelectedColumn

${shapes_load_command}
${colors_load_command}

${device_open_command}
plot(dataIn[,1]~dataIn[,2] ${colors_use_command} ${shapes_use_command})
${colors_legend_command}
${shapes_legend_command}
dev.off()

write("$app_name : image written to file $output_file", stderr())

quit(status=0)
EOC

echo "$app_name : executing R script ($Rscript_name)"
R --vanilla < "$Rscript_name" > /tmp/rlog$rand.txt
if [ $? -ne 0 ]; then
        cat /tmp/rlog$rand.txt; echo "$appname : error executing R script ($Rscript_name)"
        echo "Command line:"
        echo "$0 ${*}"
        exit 1
fi      
echo "Command line:"
echo "$0 ${*}"
rm -f /tmp/$rand.* /tmp/rlog"$rand".txt "${Rscript_name}"
exit 0

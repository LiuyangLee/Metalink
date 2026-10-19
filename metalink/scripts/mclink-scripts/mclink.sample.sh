#!/bin/bash
:<<!EOF!
 * @Date: 2023-06-14 16:58:11
 * @LastEditors: LiuyangLi
 * @LastEditTime: 2024-09-21 23:00:54
 * @FilePath: /liliuyang/software/mclink/mclink.sh
 * @Description:
!EOF!

help_message () {
	echo ""
	echo "Usage: mclink.sh [options] -i RelativeAbundance_input_file [-f SumCoverage_input_file] -k KO_pathway_info -s script_DIR -p TRUE -o output_DIR"
	echo "Note1: Make sure to provide directory."
	echo ""
	echo "Options:"
	echo ""
	echo "	-i STR          RelativeAbundance_input_file: primary input matrix (RA scale); the Keep filter mask is derived from it"
	echo "  -f STR          SumCoverage_input_file: optional second abundance matrix (SumCoverage scale); filtered with the RA-derived Keep mask and processed under <output_DIR>/SumCoverage"
	echo "	-k STR          file for KO_pathway_info"
	echo "	-s STR          directory for script"
	echo "	-p STR          pathway_split TRUE or FALSE"
	echo "  -n INT          minsample: minimum number of samples in which a KO must be detected (default: 100)"
	echo "  -r NUM          prevalence: minimum detection rate (0-1) of a KO across all samples (default: 0)"
	echo "  -a NUM          abundance: minimum mean relative abundance of a KO across all samples (default: 0)"
	echo "	-o STR          directory for output and error files (if exist). e.g., All_modules.tsv"
	echo "  -e NUM          completeness_threshold: threshold for filtering by completeness (default: 0.5)"
	echo "  -l STR          plus_scale_method: scaling method for plus (default: mean)"
	echo "  -m STR          comma_scale_method: scaling method for comma (default: max)"
	echo "  -c INT          n_cores: number of cores for mclink module-parallel processing (default: 1)"
	echo "";}

########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################
minsample=100
prevalence=0
abundance=0
completeness_threshold=0.5
plus_scale_method=mean
comma_scale_method=max
n_cores=1
## load in params
OPTS=`getopt -o hi:k:s:p:n:r:a:o:e:l:m:c:f: --long RelativeAbundance_input_file:,SumCoverage_input_file:,KO_pathway_info:,script_DIR:,pathway_split:,minsample:,prevalence:,abundance:,output_DIR:,completeness_threshold:,plus_scale_method:,comma_scale_method:,n_cores:,help -- "$@"`
if [ $? -ne 0 ]; then help_message; exit 1; fi

# loop through input params
while true; do
		case "$1" in
				-i) RelativeAbundance_input_file=$2;       shift 2;;
				-f) SumCoverage_input_file=$2;   shift 2;;
				-k) KO_pathway_info=$2;  shift 2;;
				-s) script_DIR=$2;       shift 2;;
				-p) pathway_split=$2;    shift 2;;
				-n) minsample=$2;        shift 2;;
				-r) prevalence=$2;       shift 2;;
				-a) abundance=$2;        shift 2;;
				-o) output_DIR=$2;       shift 2;;
				-e) completeness_threshold=$2; shift 2;;
				-l) plus_scale_method=$2;  shift 2;;
				-m) comma_scale_method=$2; shift 2;;
				-c) n_cores=$2;          shift 2;;
				-h|--help) help_message; exit 0;;
				--) shift; break;;
				*) break;;
		esac
done
########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################

for var in minsample prevalence abundance; do
	eval val=\$$var
	if [ -z "$val" ]; then
		echo "Error: $var not specified."
		help_message; exit 1
	fi
done

#check if the pathway_split exists
if [ -z $pathway_split ]; then
	echo "Please input the pathway_split. Exiting..."
	help_message; exit 1
fi


#check if the RelativeAbundance_input_file exists
if [ -z ${RelativeAbundance_input_file} ]; then
	echo "Please input the RelativeAbundance_input_file."
	help_message; exit 1
else
	if [ ! -s $RelativeAbundance_input_file ]; then
			echo "File $RelativeAbundance_input_file does not exist. Exiting..."
			help_message; exit 1
	fi
fi

#check if the SumCoverage_input_file exists when provided
if [ ! -z ${SumCoverage_input_file} ]; then
	if [ ! -s $SumCoverage_input_file ]; then
			echo "File $SumCoverage_input_file does not exist or is empty. Exiting..."
			help_message; exit 1
	fi
fi


#check if the KO_pathway_info was input
if [ -z ${KO_pathway_info} ]; then
	echo "Please input the KO pathway info."
	help_message; exit 1
else
	#check if the KO_pathway_info file exists
	if [ ! -s $KO_pathway_info ]; then
			echo "$KO_pathway_info does not exist. Exiting..."
			help_message; exit 1
	fi
fi

# Checks for scripts folder
if [ -z ${script_DIR} ]; then
	echo "Please input the script_DIR."
	help_message; exit 1
else
	if [ ! -s ${script_DIR}/mclink.sample.cmd.R ]; then
			echo "The file ${script_DIR}/mclink.sample.cmd.R does not exist."
			help_message; exit 1
	fi
	#if [ ! -s ${script_DIR}/M1S4.Heatmap.cmd.R ]; then
	#		echo "The file ${script_DIR}/M1S4.Heatmap.cmd.R does not exist."
	#		help_message; exit 1
	#fi
	#if [ ! -s ${script_DIR}/M1S5.Heatmap.cmd.R ]; then
	#		echo "The file ${script_DIR}/M1S5.Heatmap.cmd.R does not exist."
	#		help_message; exit 1
	#fi
fi


#check if the output_DIR dir exists
if [ -z ${output_DIR} ]; then
	echo "Please input the output_DIR."
	help_message; exit 1
fi
########################################################################################################
########################                       MKDIR 1 !                        ########################
########################################################################################################

if [ ! -d ${output_DIR} ]; then mkdir -p ${output_DIR};
else
	echo "Warning: ${output_DIR} already exists."
fi


output_KO_Sample_by_pathway_presence=${output_DIR}/Completeness/mc_by_pathway/
output_Heatmap_by_pathway_presence=${output_DIR}/Completeness/mc_by_pathway_heatmap/

if [ ! -d ${output_KO_Sample_by_pathway_presence} ]; then mkdir -p ${output_KO_Sample_by_pathway_presence};
else
    echo "Warning: ${output_KO_Sample_by_pathway_presence} already exists."
		rm -r ${output_KO_Sample_by_pathway_presence}/*.tsv
fi

if [ ! -d ${output_Heatmap_by_pathway_presence} ]; then mkdir -p ${output_Heatmap_by_pathway_presence};
else
    echo "Warning: ${output_Heatmap_by_pathway_presence} already exists."
		rm -rf ${output_Heatmap_by_pathway_presence}/*.pdf ${output_Heatmap_by_pathway_presence}/*.html
fi

########################################################################################################
########################                       MKDIR 2 !                        ########################
########################################################################################################

output_KO_Sample_by_pathway_log_transformed=${output_DIR}/RelativeAbundance/log_scale_by_pathway/
output_KO_Sample_by_pathway_log_transformed_filtered=${output_DIR}/RelativeAbundance/log_scale_by_pathway_filtered/
#output_Heatmap_by_pathway_log_transformed=${output_DIR}/RelativeAbundance/log_scale_by_pathway_heatmap/
output_Heatmap_by_pathway_log_transformed_filtered=${output_DIR}/RelativeAbundance/log_scale_by_pathway_filtered_heatmap/

if [ ! -d ${output_KO_Sample_by_pathway_log_transformed} ]; then mkdir -p ${output_KO_Sample_by_pathway_log_transformed};
else
    echo "Warning: ${output_KO_Sample_by_pathway_log_transformed} already exists."
		rm -r ${output_KO_Sample_by_pathway_log_transformed}/*.tsv
fi

if [ ! -d ${output_KO_Sample_by_pathway_log_transformed_filtered} ]; then mkdir -p ${output_KO_Sample_by_pathway_log_transformed_filtered};
else
    echo "Warning: ${output_KO_Sample_by_pathway_log_transformed_filtered} already exists."
		rm -r ${output_KO_Sample_by_pathway_log_transformed_filtered}/*.tsv
fi

#if [ ! -d ${output_Heatmap_by_pathway_log_transformed} ]; then mkdir -p ${output_Heatmap_by_pathway_log_transformed};
#else
#    echo "Warning: ${output_Heatmap_by_pathway_log_transformed} already exists."
#		rm -rf ${output_Heatmap_by_pathway_log_transformed}/*.pdf ${output_Heatmap_by_pathway_log_transformed}/*.html
#fi

if [ ! -d ${output_Heatmap_by_pathway_log_transformed_filtered} ]; then mkdir -p ${output_Heatmap_by_pathway_log_transformed_filtered};
else
    echo "Warning: ${output_Heatmap_by_pathway_log_transformed_filtered} already exists."
		rm -rf ${output_Heatmap_by_pathway_log_transformed_filtered}/*.pdf ${output_Heatmap_by_pathway_log_transformed_filtered}/*.html
fi
########################################################################################################
########################                       MKDIR 3 !                        ########################
########################################################################################################

output_KO_Sample_by_pathway_avg_transformed=${output_DIR}/RelativeAbundance/avg_scale_by_pathway/
output_KO_Sample_by_pathway_avg_transformed_filtered=${output_DIR}/RelativeAbundance/avg_scale_by_pathway_filtered/
#output_Heatmap_by_pathway_avg_transformed=${output_DIR}/RelativeAbundance/avg_scale_by_pathway_heatmap/
output_Heatmap_by_pathway_avg_transformed_filtered=${output_DIR}/RelativeAbundance/avg_scale_by_pathway_filtered_heatmap/

if [ ! -d ${output_KO_Sample_by_pathway_avg_transformed} ]; then mkdir -p ${output_KO_Sample_by_pathway_avg_transformed};
else
    echo "Warning: ${output_KO_Sample_by_pathway_avg_transformed} already exists."
		rm -r ${output_KO_Sample_by_pathway_avg_transformed}/*.tsv
fi

if [ ! -d ${output_KO_Sample_by_pathway_avg_transformed_filtered} ]; then mkdir -p ${output_KO_Sample_by_pathway_avg_transformed_filtered};
else
    echo "Warning: ${output_KO_Sample_by_pathway_avg_transformed_filtered} already exists."
		rm -r ${output_KO_Sample_by_pathway_avg_transformed_filtered}/*.tsv
fi

#if [ ! -d ${output_Heatmap_by_pathway_avg_transformed} ]; then mkdir -p ${output_Heatmap_by_pathway_avg_transformed};
#else
#    echo "Warning: ${output_Heatmap_by_pathway_avg_transformed} already exists."
#		rm -rf ${output_Heatmap_by_pathway_avg_transformed}/*.pdf ${output_Heatmap_by_pathway_avg_transformed}/*.html
#fi

if [ ! -d ${output_Heatmap_by_pathway_avg_transformed_filtered} ]; then mkdir -p ${output_Heatmap_by_pathway_avg_transformed_filtered};
else
    echo "Warning: ${output_Heatmap_by_pathway_avg_transformed_filtered} already exists."
		rm -rf ${output_Heatmap_by_pathway_avg_transformed_filtered}/*.pdf ${output_Heatmap_by_pathway_avg_transformed_filtered}/*.html
fi

########################################################################################################
########################                       MKDIR 4 !                        ########################
########################################################################################################

output_KO_Sample_by_pathway_log_transformed_SumCoverage=${output_DIR}/SumCoverage/log_scale_by_pathway/
output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage=${output_DIR}/SumCoverage/log_scale_by_pathway_filtered/
output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage=${output_DIR}/SumCoverage/log_scale_by_pathway_filtered_heatmap/

output_KO_Sample_by_pathway_avg_transformed_SumCoverage=${output_DIR}/SumCoverage/avg_scale_by_pathway/
output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage=${output_DIR}/SumCoverage/avg_scale_by_pathway_filtered/
output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage=${output_DIR}/SumCoverage/avg_scale_by_pathway_filtered_heatmap/

if [ ! -z ${SumCoverage_input_file} ]; then

	if [ ! -d ${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage} ]; then mkdir -p ${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage};
	else
	    echo "Warning: ${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage} already exists."
			rm -r ${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage}/*.tsv
	fi

	if [ ! -d ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage} ]; then mkdir -p ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage};
	else
	    echo "Warning: ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage} already exists."
			rm -rf ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage}/*.pdf ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage}/*.html
	fi

	if [ ! -d ${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage} ]; then mkdir -p ${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage};
	else
	    echo "Warning: ${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage} already exists."
			rm -r ${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage}/*.tsv
	fi

	if [ ! -d ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage} ]; then mkdir -p ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage};
	else
	    echo "Warning: ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage} already exists."
			rm -rf ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage}/*.pdf ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage}/*.html
	fi

fi

########################################################################################################
########################                    BEGIN PIPELINE !                    ########################
########################################################################################################
echo "[$(date '+%Y-%m-%d %H:%M:%S')] mclink module-parallel processing using ${n_cores} core(s)."
set -x
set -e
sumcov_arg=""
if [ ! -z ${SumCoverage_input_file} ]; then
	sumcov_arg="--in_KO_Sample_wide_SumCoverage ${SumCoverage_input_file}"
fi

Rscript ${script_DIR}/mclink.sample.cmd.R \
		--in_KO_Sample_wide_RelativeAbundance ${RelativeAbundance_input_file} \
		--in_KO_pathway_info ${KO_pathway_info} \
		--table_feature abundance \
		--plus_scale_method ${plus_scale_method} \
		--comma_scale_method ${comma_scale_method} \
		--minsample ${minsample} \
		--prevalence ${prevalence} \
		--abundance ${abundance} \
		--out_DIR ${output_DIR} \
		--split_by_pathway ${pathway_split} \
		--n_cores ${n_cores} \
		${sumcov_arg}

Rscript ${script_DIR}/mclink.sample.cmd.R \
		--in_KO_Sample_wide_RelativeAbundance ${RelativeAbundance_input_file} \
		--in_KO_pathway_info ${KO_pathway_info} \
		--table_feature completeness \
		--plus_scale_method ${plus_scale_method} \
		--comma_scale_method max \
		--minsample ${minsample} \
		--prevalence ${prevalence} \
		--abundance ${abundance} \
		--out_DIR ${output_DIR} \
		--split_by_pathway ${pathway_split} \
		--n_cores ${n_cores}

sumcov_filtered_arg=""
if [ -s ${output_DIR}/SumCoverage/mc_abundance.tsv ]; then
	sumcov_filtered_arg="--abundance_SumCoverage ${output_DIR}/SumCoverage/mc_abundance.tsv --output_SumCoverage ${output_DIR}/SumCoverage/mc_abundance_filtered_by_completeness.tsv"
fi

if [ -s ${output_DIR}/RelativeAbundance/mc_abundance.tsv ] && [ -s ${output_DIR}/Completeness/mc_completeness.tsv ]; then
	Rscript ${script_DIR}/mclink_abundance_filtered_by_completeness.cmd.R \
		--abundance_RelativeAbundance ${output_DIR}/RelativeAbundance/mc_abundance.tsv \
		--completeness ${output_DIR}/Completeness/mc_completeness.tsv \
		--output_RelativeAbundance ${output_DIR}/RelativeAbundance/mc_abundance_filtered_by_completeness.tsv \
		--completeness-threshold ${completeness_threshold} \
		--strategy zero \
		--meta-cols "Module_Entry,Level_2,Level_3,Module_Name,Definition" \
		${sumcov_filtered_arg}
else
	echo "SKIP abundance filtering by completeness: RelativeAbundance/mc_abundance.tsv or Completeness/mc_completeness.tsv is missing/empty (mclink was skipped upstream, see the 'nothing to compute' message above). Check whether any sequences were functionally classified for this input."
fi

num_samples=$(head -1 ${RelativeAbundance_input_file} | awk -F'\t' '{print NF-1}')
echo "Detected ${num_samples} sample(s) in input."

for file in ${output_KO_Sample_by_pathway_presence}/*.tsv; do
    basename=$(basename ${file})
    echo "Processing: ${basename}"

    log_transformed_file=${output_KO_Sample_by_pathway_log_transformed}/${basename}
		log_transformed_file_filtered=${output_KO_Sample_by_pathway_log_transformed_filtered}/${basename}
    [[ -f "${log_transformed_file}" ]] || { echo "SKIP ${basename}: no log_transformed abundance file"; continue; }
		log_transformed_file_SumCoverage=${output_KO_Sample_by_pathway_log_transformed_SumCoverage}/${basename}
		log_transformed_file_filtered_SumCoverage=${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage}/${basename}
		sumcov_log_filtered_arg=""
		if [ -f "${log_transformed_file_SumCoverage}" ]; then
			sumcov_log_filtered_arg="--abundance_SumCoverage ${log_transformed_file_SumCoverage} --output_SumCoverage ${log_transformed_file_filtered_SumCoverage}"
		fi
		Rscript ${script_DIR}/mclink_abundance_filtered_by_completeness.cmd.R \
        --abundance_RelativeAbundance ${log_transformed_file} \
        --completeness ${file} \
        --output_RelativeAbundance ${log_transformed_file_filtered} \
        --completeness-threshold ${completeness_threshold} \
        --strategy zero \
        --meta-cols Module_Name \
        ${sumcov_log_filtered_arg}

    if [ ${num_samples} -ge 2 ]; then
        avg_transformed_file=${output_KO_Sample_by_pathway_avg_transformed}/${basename}
			avg_transformed_file_filtered=${output_KO_Sample_by_pathway_avg_transformed_filtered}/${basename}
        [[ -f "${avg_transformed_file}" ]] || { echo "SKIP ${basename}: no avg_transformed abundance file"; continue; }
			avg_transformed_file_SumCoverage=${output_KO_Sample_by_pathway_avg_transformed_SumCoverage}/${basename}
			avg_transformed_file_filtered_SumCoverage=${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage}/${basename}
			sumcov_avg_filtered_arg=""
			if [ -f "${avg_transformed_file_SumCoverage}" ]; then
				sumcov_avg_filtered_arg="--abundance_SumCoverage ${avg_transformed_file_SumCoverage} --output_SumCoverage ${avg_transformed_file_filtered_SumCoverage}"
			fi
			Rscript ${script_DIR}/mclink_abundance_filtered_by_completeness.cmd.R \
            --abundance_RelativeAbundance ${avg_transformed_file} \
            --completeness ${file} \
            --output_RelativeAbundance ${avg_transformed_file_filtered} \
            --completeness-threshold ${completeness_threshold} \
            --strategy zero \
            --meta-cols Module_Name \
            ${sumcov_avg_filtered_arg}
    fi
done

set +e
Rscript ${script_DIR}/mclink.M1S4.Heatmap.cmd.R \
    --user_Sample_in_DIR ${output_KO_Sample_by_pathway_presence} \
    --out_DIR ${output_Heatmap_by_pathway_presence}
Rscript ${script_DIR}/mclink.M1S4.Heatmap.cmd.R \
    --user_Sample_in_DIR ${output_KO_Sample_by_pathway_log_transformed_filtered} \
    --out_DIR ${output_Heatmap_by_pathway_log_transformed_filtered}

if [ ${num_samples} -ge 2 ]; then
    Rscript ${script_DIR}/mclink.M1S5.Heatmap.cmd.R \
        --user_Sample_in_DIR ${output_KO_Sample_by_pathway_avg_transformed_filtered} \
        --out_DIR ${output_Heatmap_by_pathway_avg_transformed_filtered}
fi


Rscript ${script_DIR}/mclink.M1S4.Heatmap.html.cmd.R \
  --mc_completeness_DIR ${output_KO_Sample_by_pathway_presence} \
  --mc_detected_KOs ${output_DIR}/Completeness/mc_detected_KOs.tsv \
  --out_DIR ${output_Heatmap_by_pathway_presence}
Rscript ${script_DIR}/mclink.M1S4.Heatmap.html.cmd.R \
  --mc_completeness_DIR ${output_KO_Sample_by_pathway_log_transformed_filtered} \
  --mc_detected_KOs ${output_DIR}/Completeness/mc_detected_KOs.tsv \
  --out_DIR ${output_Heatmap_by_pathway_log_transformed_filtered}

if [ ${num_samples} -ge 2 ]; then
    Rscript ${script_DIR}/mclink.M1S5.Heatmap.html.cmd.R \
      --mc_completeness_DIR ${output_KO_Sample_by_pathway_avg_transformed_filtered} \
      --mc_detected_KOs ${output_DIR}/Completeness/mc_detected_KOs.tsv \
      --out_DIR ${output_Heatmap_by_pathway_avg_transformed_filtered}
fi

if [ -s ${output_DIR}/SumCoverage/mc_abundance.tsv ]; then
    Rscript ${script_DIR}/mclink.M1S4.Heatmap.cmd.R \
        --user_Sample_in_DIR ${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage} \
        --out_DIR ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage}

    if [ ${num_samples} -ge 2 ]; then
        Rscript ${script_DIR}/mclink.M1S5.Heatmap.cmd.R \
            --user_Sample_in_DIR ${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage} \
            --out_DIR ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage}
    fi

    Rscript ${script_DIR}/mclink.M1S4.Heatmap.html.cmd.R \
      --mc_completeness_DIR ${output_KO_Sample_by_pathway_log_transformed_filtered_SumCoverage} \
      --mc_detected_KOs ${output_DIR}/SumCoverage/mc_detected_KOs.tsv \
      --out_DIR ${output_Heatmap_by_pathway_log_transformed_filtered_SumCoverage}

    if [ ${num_samples} -ge 2 ]; then
        Rscript ${script_DIR}/mclink.M1S5.Heatmap.html.cmd.R \
          --mc_completeness_DIR ${output_KO_Sample_by_pathway_avg_transformed_filtered_SumCoverage} \
          --mc_detected_KOs ${output_DIR}/SumCoverage/mc_detected_KOs.tsv \
          --out_DIR ${output_Heatmap_by_pathway_avg_transformed_filtered_SumCoverage}
    fi
fi
########################################################################################################
########   GROUP KO AND ABUNDANCE BY PATHWAY AND MODULE PIPELINE SUCCESSFULLY FINISHED!!!      #########
########################################################################################################
echo -e "\nGROUP KO AND ABUNDANCE BY PATHWAY AND MODULE PIPELINE SUCCESSFULLY FINISHED!!!\n"

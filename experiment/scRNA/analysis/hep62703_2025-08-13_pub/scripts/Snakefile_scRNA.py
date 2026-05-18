
##run with
#snakemake -nr --cores 6 -s Snakefile_scRNA.py

#The same directory of the Snakefile_scRNA.py must contain all Rmd files, functions.R and config.yaml

#configfile: 'config.yaml' #possiblity to read in the configfile into dict 'config'
#possiblity to validate the config file: https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html

##################################### load libraries required for Snakemake #########################

##on the cluster, start snakemake environment:
#conda init bash
#source ~/.bashrc
#conda activate  /misc/software/ngs/pipelines/snakePipes/snakePipes-2.5.1
##run the workflow with dry mode and reason:
#snakemake -nr -s Snakefile_scRNA.py

##show workflow dag
#snakemake -nr -s Snakefile_scRNA.py --dag  | dot -Tsvg > dag.svg
##force rerunning a rule
#snakemake -nr -f -s Snakefile_scRNA.py doublet_detection

import yaml
from datetime import date
from pathlib import Path
import os.path
import os
from collections.abc import MutableMapping
from snakemake.logging import logger
import time
import jsonschema

##################################### setup #########################

#get date (not required at the moment)
#TODAY = date.today().strftime("%Y-%m-%d")



############### parts to adapt when new rule is added ##########

#To add a new script to the workflow:
#1. Define new parameters in the config.yaml that influence the script
#2. Add the rule with a call to the script, inputs and outputs to this file in a similar manner as the other rules
#3. Adapt the dictionaries and all_outputs string in this section
#4. Adapt the schema to define allowed formats for the parameters
#5. Add the output file to rule ALL
#6. Add the script to the init_project.sh file so it is downloaded

#define which parameters should trigger which rule (the order of the rules should correspond to the order of execution)
rule_dict = {
  'preprocessing' : ['libs_exclude', 'filters__min_cells'],
  'metadata_specification' : [
    'metadata_spec__metadata_spec_file',
    'metadata_spec__shorten_sample_names',
    'metadata_spec__n_char_shortened'
  ],
  'doublet_detection' : [],
  'qc_report' : ['filters__min_lib_ncounts', 'filters__min_n_feats', 'filters__max_pct_mt'],
  'dim_reduction' : [
    'scParam__PCfrom',
    'scParam__PCto',
    'scParam__VarFeat',
    'integration__integrate',
    'integration__batch_info_file',
    'integration__methods'
  ], #'scParam__jack'
  'clustering' : ['scParam__cluster_on', 'scParam__dimred'],
  'find_all_markers' : ['scParam__dimred', 'scParam__clusters'],
  'plot_markers' : ['markers', 'scParam__dimred', 'scParam__clusters'],
  'plot_signatures' : ['signatures', 'scParam__dimred', 'scParam__clusters'],
  'singler' : ['singlerRefpath', 'scParam__dimred'],
  'azimuth' : [],
  'azimuth_analysis' : [],
  'scGate' : [],
  'annotation_CD4_CD8' : [],
  'annotation_T_cell_subtypes' : [],
  'set_model_variables' : [],
  'cluster_assessment_umaps' : [],
  'fit_DE_cond_nebula_CD4_CD8' : [],
  'stats_DE_cond_nebula_CD4_CD8' : [],
  'GSEA_cond_nebula_CD4_CD8_A' : [],
  'GSEA_analysis_cond_CD4_CD8_A' : [],
  'GSEA_cond_barplots_CD4_CD8_A' : [],
  'get_cell_counts' : [],
  'DA_glm_nb' : [],
  'umap_plots_manuscript' : [],
  'shinycell' : [],
  'liana' : []
}

#Define what should be written in the all_outputs.html file
all_outputs_html_content = """<pre>
Read in filtered counts output from Cellranger and do demultiplexing if necessary. Apply min_cells threshold on the features.<br>Save all samples separately as Seurat objects and combined together as Seurat object for further analyses.
<a href="Preprocessing.html">01. Preprocessing</a>
<br><br>
Add doublet scores based on scDblFinder.
<a href="Doublet_Detection.html">02. Doublet detection</a>
<br><br>
Show quality control analyses and filter out cells based on user-defined thresholds for min_lib_ncounts, min_n_feats and max_pct_mt.
<a href="QC_Report.html">03. QC report</a>
<br><br>
Perform PCA and secondary dimensionality reductions with UMAP for different PC numbers, add cell cycle scores and graph-based clusters.
<a href="Dim_Reduction.html">04. Dimensionality reduction</a>
<br><br>
Add cell cycle scores and graph-based clusters and plot them together with available meta data in the UMAP
<a href="Clustering.html">05. Clustering</a>
<br><br>
Plot user-defined marker genes as heatmap, dotplot, violins and in UMAP space.
<a href="Plot_Markers.html">06. Marker plots</a>
<br><br>
Add signature scores based on user-defined signature gene sets and plot them as heatmap and in UMAP space.
<a href="Plot_Signatures.html">07. Signature plots</a>
<br><br>
Find and visualize marker genes based on the chosen group defined in 'clusters'.<br>Perform GO enrichment analysis on the marker gene sets.
<a href="Find_All_Markers.html">08. Cluster-level marker genes</a>
<br><br>
Perform automated cell type classification with SingleR using multiple bulk-RNA reference datasets.
<a href="Singler.html">09. SingleR annotation</a>
<br><br>
Run Liana to get ligand-receptor interaction predictions from multiple methods.
<a href="Liana_interactions.html">10. Liana interactions</a>
<br><br>
Generate a shiny app for interactive exploration of the dataset.
<a href="ShinyCell.html">11. ShinyCell</a>
<br><br>
</pre>"""

##################################### custom functions ##########


def read_yaml(filename):
  #read yaml file into a dictionary, return empty dictionary when file is empty
  with open(filename, "r") as stream:
      try:
          obj = yaml.safe_load(stream)
          if obj is None:
            obj = {}
      except yaml.YAMLError as exc:
          logger.info(exc)
  return obj


def validate(obj, schema):
  # this will simply crash the program with a stacktrace if the
  # config does not validate against the schema
  schemauri=Path("./schema/").resolve().as_uri()
  schemaresolver = jsonschema.validators.RefResolver(
            base_uri=schemauri+"/",
            referrer=True)
  jsonschema.validate(instance = obj,
                      schema = {"$ref": schema},
                      resolver=schemaresolver)

  #if file is empty the type is "NoneType". Initialize an empty dict in that case
  if not isinstance(obj, dict): 
    obj = {}
  return obj


def touch_if_not_exists(filepath):
  #create an empty file if it does not exist
  file_exists = os.path.isfile(filepath) 
  if not file_exists: 
    Path(filepath).touch()


def touch_if_new(new_keys, rule_dict):
  #update the timestamp for the dummy file indicating whether the respective rule should be rerun if
  #any new_keys are in the rule dictionray
  for rule in rule_dict.keys():
     if any(x in rule_dict[rule] for x in new_keys):
      logger.info('Rerun ' + rule)
      Path(SNAKEDIR + '.run_' + rule).touch()
   
 
def reset_if_not_modified(same_keys, rule_dict):
  #reset the timestamps on dummy files to an older date than the output files 
  #if all entries in the rule dictionary for a rule are within same_keys
  for rule in rule_dict.keys():
    if all([x in same_keys for x in rule_dict[rule] ]):
      logger.info(rule + ' parameters are unchanged')
      os.system('touch -a -m -t 200001010000.00 ' + SNAKEDIR + '.run_' + rule)   


def flatten_dict(d: MutableMapping, parent_key: str = '', sep: str ='__') -> MutableMapping:
  #convert a nestted dictionary to a one-level dictionary
  #https://www.freecodecamp.org/news/how-to-flatten-a-dictionary-in-python-in-4-different-ways/
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, MutableMapping):
            items.extend(flatten_dict(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)


def nest_dict(flat):
  #convert a flat dictionary to a nested dictionary
  ##https://stackoverflow.com/questions/50607128/creating-a-nested-dictionary-from-a-flattened-dictionary
    result = {}
    for k, v in flat.items():
        _nest_dict_rec(k, v, result)
    return result


def _nest_dict_rec(k, v, out, sep='__'):
    k, *rest = k.split(sep, 1)
    if rest:
        _nest_dict_rec(rest[0], v, out.setdefault(k, {}))
    else:
        out[k] = v


def dict_compare(d1, d2):
  #compare two dictionaries and return 4 sets including keys that are added, removed, modified or the same in d1 compared to d2
  #https://stackoverflow.com/questions/4527942/comparing-two-dictionaries-and-checking-how-many-key-value-pairs-are-equal
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    added = d1_keys - d2_keys
    removed = d2_keys - d1_keys
    #modified = {o : (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    modified = set(o for o in shared_keys if d1[o] != d2[o])
    same = set(o for o in shared_keys if d1[o] == d2[o])
    return added, removed, modified, same
  

def get_new_key_values(d1, d2):
  #return a dictionary only containing entries in d1 that are new or modified compared to d2
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    added = d1_keys - d2_keys
    removed = d2_keys - d1_keys
    #modified = {o : (d1[o], d2[o]) for o in shared_keys if d1[o] != d2[o]}
    modified = {o : d1[o] for o in shared_keys if d1[o] != d2[o]}
    added = {o : d1[o] for o in added}
    new_key_values = {**modified,**added}
    return new_key_values


def get_unmodified_key_values(d1, d2):
   #return a dictionary only containing entries that are the same in d1 and d2
    d1_keys = set(d1.keys())
    d2_keys = set(d2.keys())
    shared_keys = d1_keys.intersection(d2_keys)
    unmodified_dict = {o : d1[o] for o in shared_keys if d1[o] == d2[o]}
    return unmodified_dict


def find_unused_params(new_keys, rule_dict):
  #compare the values in the rule dictionary with the new_keys to find the ones that do not
  #influence any of the rules. Return them as set.
  used_params = set(item for sublist in rule_dict.values() for item in sublist)
  unused = []
  for key in new_keys:
    if key not in used_params:
        unused.append(key)
  return(set(unused))


def initialize_snakedir(snakedirpath, rule_dict, init_dict, html_content):
  #initialize directory
  Path(snakedirpath).mkdir(parents=True, exist_ok=True)
  #initialize dummy files for each workflow rule
  for rule in rule_dict.keys():
    touch_if_not_exists(snakedirpath + '.run_' + rule)
    touch_if_not_exists(snakedirpath + '.' + rule + '.yaml')

  #initialize the "old" config, against which the config.yaml file is compared
  CONFIG_SNAKE = snakedirpath + "CONFIG_SNAKE.yaml"
  if not os.path.isfile(CONFIG_SNAKE):
    with open(CONFIG_SNAKE, 'w') as file:
      documents = yaml.dump(init_dict, file)
    
  with open('all_outputs.html', 'w') as file:
    file.write(html_content)
    
  return CONFIG_SNAKE


def write_rule_config(new_config_dict, rule, out_yaml):
  #write out the part of new_config_dict that influences a specified rule into out_yaml file
  dict_out = {}
  for i in new_config_dict.keys():
    if i in rule_dict[rule]:
      dict_out[i] = new_config_dict[i]
  with open(out_yaml, 'w') as f:
      yaml.dump(dict_out, f)

  
def merge_rule_config(snakedirpath, out_yaml, rule_dict, out_key_val_dict):
  #merge all rule-specific yaml files and the unmodified key-value pairs in a final dictionary
  #which is saved in the end of the workflow
  for rule in rule_dict.keys():
    temp_dict = read_yaml(snakedirpath + '.' + rule + '.yaml')
    #dict on the right overwrites key:value pairs on the left
    out_key_val_dict = { **out_key_val_dict, **temp_dict}
  print(out_key_val_dict)
  own_config = nest_dict(out_key_val_dict)
  with open(out_yaml, 'w') as f:
      yaml.dump(out_key_val_dict, f)


##################################### Container Setup

container_setup = read_yaml('container_setup.yaml')
container_path = os.path.join(container_setup['container_path'], container_setup['container_name'])
container_bind = '--bind ' + ','.join(container_setup['bind'])
container_command = 'singularity exec ' + container_bind + ' ' + container_path

##################################### read config and compare its entries with previous run ##########

SCRIPTS_DIR = os.getcwd() + '/'
ANALYSIS_DIR = os.environ['ANALYSIS_DIR'] + '/'
CONFIG =  'config.yaml'

#initialize hidden directory containing dummy files that trigger the rules and the old version of the config
SNAKEDIR = '.workflow_call/'
CONFIG_SNAKE = initialize_snakedir(snakedirpath = SNAKEDIR, rule_dict = rule_dict, init_dict = {}, html_content = all_outputs_html_content)

#read new config and previous config
#config_dict = read_yaml(ANALYSIS_DIR + CONFIG_DIR + "config.yaml") #gets loaded by specifying --configfile
config_dict = read_yaml(CONFIG)
try:
  validate(config_dict, "config.schema.json")
except Exception as exc:
  raise Exception(f"{CONFIG} contains errors: {exc}")

#does not need to be validated since not manually edited:
config_old_dict = read_yaml(CONFIG_SNAKE) 

#compare the old with the new config
flat_new = flatten_dict(config_dict)
flat_old = flatten_dict(config_old_dict)
added, removed, modified, same = dict_compare(flat_new, flat_old)
unused = find_unused_params(added.union(removed, modified, same), rule_dict)
added = added - unused
removed = removed - unused
modified = modified - unused
same = same - unused
logger.info('Added parameters: ' +  ', '.join(added))
logger.info('Removed parameters: ' + ', '.join(removed))
logger.info('Modified parameters: ' +  ', '.join(modified))
logger.info('Unchanged parameters: ' +  ', '.join(same))
logger.info('Parameters not used by any rules: ' +  ', '.join(unused))

# get key:value pairs for parameters that are newly added or changed
updated_key_val_dict = get_new_key_values(flat_new, flat_old)
# get key:value pairs for parameters that have not been modified and are relevant for the rules
unmodified_key_val_dict = get_unmodified_key_values(flat_new, flat_old)


#modify timestamps on the dummy files if the rule should be rerun and change to old timestamp if all parameters are the same
touch_if_new(added, rule_dict)
touch_if_new(modified, rule_dict)
reset_if_not_modified(same, rule_dict)

######################################### RULES ###########################################################

#requires 'mail' to be installed
onsuccess:
        logger.info("Workflow finished, no error")
        merge_rule_config(SNAKEDIR, CONFIG_SNAKE, rule_dict, unmodified_key_val_dict)
        #shell("mail -s 'Workflow finished' " + E_MAIL + " < {log}")

onerror:
        logger.info("An error occurred")
        merge_rule_config(SNAKEDIR, CONFIG_SNAKE, rule_dict, unmodified_key_val_dict)
        #shell("mail -s 'an error occurred' " + E_MAIL + " < {log}")


rule ALL:
    input:
        liana_done = "Liana_interactions.html",
        shinycell_done = "ShinyCell.html",
        azimuth_analysis_done = "Azimuth_analysis.html",
        cluster_assessment_umaps_done = "cluster_assessment_umaps.html",
        GSEA_analysis_cond_CD4_CD8_A_done = "GSEA_analysis_cond_CD4-CD8_A.html",
        GSEA_cond_barplots_CD4_CD8_A_done = "GSEA_cond_barplots_CD4-CD8_A.html",
        DA_glm_nb_done = "DA_glm_nb.html",
        umap_plots_manuscript_done = "umap_plots_manuscript.html"


rule preprocessing:
    input:
        config_check = SNAKEDIR + ".run_preprocessing",
        script = "Preprocessing.Rmd"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="preprocessing"
    threads: 6
    output:
        #ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_raw_merged.rds"
        "Preprocessing.html"
    run:
        print("Running Preprocessing.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, 'preprocessing', '.workflow_call/.preprocessing.yaml')

rule metadata_specification:
    input:
        config_check = SNAKEDIR + ".run_metadata_specification",
        script = "Metadata_Specification.Rmd",
        infile = "Preprocessing.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="metadata_specification"
    threads: 6
    output:
        "Metadata_Specification.html"
    run:
        print("Running Metadata_Specification.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, 'metadata_specification', '.workflow_call/.metadata_specification.yaml')

rule doublet_detection:
    input:
        config_check = SNAKEDIR + '.run_doublet_detection',
        script =  "Doublet_Detection.Rmd",
        infile = "Metadata_Specification.html" #ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_raw_merged.rds",
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="doublet_detection"
    threads: 6
    output:
        "Doublet_Detection.html"
        #ANALYSIS_DIR + config_dict['outputpath']['table'] + "/scDblFinder_predictions.csv"
    run:
        print("Running Doublet_Detection.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, 'doublet_detection', '.workflow_call/.doublet_detection.yaml')

rule qc_report:
    input:
        config_check = SNAKEDIR + '.run_qc_report',
        script =  "QC_Report.Rmd",
        infile = "Doublet_Detection.html" #ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_raw_merged.rds"
        #infile2 = 'Doublet_Detection.html'
        #infile = ANALYSIS_DIR + config_dict['outputpath']['table'] + "/scDblFinder_predictions.csv"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="qc_report"
    threads: 6
    output:
        "QC_Report.html"
        #ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_filtered.rds"
    run:
        print("Running QC_Report.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        time.sleep(2)
        write_rule_config(flat_new, 'qc_report', '.workflow_call/.qc_report.yaml')

rule dim_reduction:
    input:
        config_check = SNAKEDIR + '.run_dim_reduction',
        script =  "Dim_Reduction.Rmd",
        #harmony_file = ANALYSIS_DIR + (config_dict['harmony_csv'] if config_dict['harmony_csv'] is not None else ""),
        infile = "QC_Report.html" #ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_filtered.rds"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="dim_reduction"
    threads: 6
    output:
        "Dim_Reduction.html"
    run:
        print("Running Dim_Reduction.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        time.sleep(2)
        write_rule_config(flat_new, 'dim_reduction', '.workflow_call/.dim_reduction.yaml')

rule clustering:
    input:
        config_check = SNAKEDIR + '.run_clustering',
        script =  "Clustering.Rmd",
        infile = "Dim_Reduction.html" #ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_filtered.rds"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="clustering"
    threads: 6
    output:
        "Clustering.html"
    run:
        print("Running Clustering.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        time.sleep(2)
        write_rule_config(flat_new, 'clustering', '.workflow_call/.clustering.yaml')

rule find_all_markers:
    input:
        config_check = SNAKEDIR + '.run_find_all_markers',
        script =  "Find_All_Markers.Rmd",
        clustering_done = "Clustering.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="find_all_markers"
    threads: 6
    output:
        "Find_All_Markers.html"
    run:
        print("Running Find_All_Markers.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, 'find_all_markers', '.workflow_call/.find_all_markers.yaml')

rule plot_markers:
    input:
        config_check = SNAKEDIR + '.run_plot_markers',
        #infile = ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_filtered.rds",
        marker_file = ANALYSIS_DIR + config_dict['markers'] if config_dict['markers'] is not None else "/misc/rci/refdata/marker_dummy.txt",
        script =  "Plot_Markers.Rmd",
        find_all_markers_done = "Find_All_Markers.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="plot_markers"
    threads: 6
    output:
        "Plot_Markers.html"
    run:
        print("Running Plot_Markers.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, 'plot_markers', '.workflow_call/.plot_markers.yaml')

rule plot_signatures:
    input:
        config_check = SNAKEDIR + '.run_plot_signatures',
        signature_file = ANALYSIS_DIR + config_dict['signatures'] if config_dict['signatures'] is not None else "/misc/rci/refdata/marker_dummy.txt",
        #infile = ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_filtered.rds",
        script =  "Plot_Signatures.Rmd",
        plot_markers_done = "Plot_Markers.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="plot_signatures"
    threads: 6
    output:
        "Plot_Signatures.html"
    run:
        print("Running Plot_Signatures.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        time.sleep(2)
        write_rule_config(flat_new, 'plot_signatures', '.workflow_call/.plot_signatures.yaml')

rule singler:
    input:
        config_check = SNAKEDIR + '.run_singler',
        #infile = ANALYSIS_DIR + config_dict['outputpath']['data'] + "/SeuratObj_filtered.rds",
        script =  "Singler.Rmd",
        plot_signatures_done = "Plot_Signatures.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname="singler"
    threads: 6
    output:
        "Singler.html"
    run:
        print("Running Singler.Rmd with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        #update_analparams(updated_key_val_dict, 'singler', CONFIG_SNAKE)
        write_rule_config(flat_new, 'singler', '.workflow_call/.singler.yaml')

rule_name = "azimuth"

rule azimuth:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name.capitalize()}.Rmd",
        singler_done = "Singler.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name.capitalize()}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "azimuth_analysis"

rule azimuth_analysis:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name.capitalize()}.Rmd",
        azimuth_done = "Azimuth.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name.capitalize()}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "scGate"

rule scGate:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name}.Rmd",
        azimuth_done = "Azimuth.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "annotation"
config_id = "CD4-CD8"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule annotation_CD4_CD8:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name.capitalize()}.Rmd",
        scGate_done = "scGate.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name.capitalize()}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "annotation"
config_id = "T-cell-subtypes"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule annotation_T_cell_subtypes:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name.capitalize()}.Rmd",
        annotation_CD4_CD8_done = "Annotation_CD4-CD8.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name.capitalize()}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "set_model_variables"

rule set_model_variables:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name}.Rmd",
        annotation_T_cell_subtypes_done = "Annotation_T-cell-subtypes.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "cluster_assessment_umaps"

rule cluster_assessment_umaps:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name}.Rmd",
        set_model_variables_done = "set_model_variables.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "fit_DE_cond_nebula"
config_id = "CD4-CD8"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule fit_DE_cond_nebula_CD4_CD8:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name}.Rmd",
        set_model_variables_done = "set_model_variables.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "stats_DE_cond_nebula"
config_id = "CD4-CD8"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule stats_DE_cond_nebula_CD4_CD8:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name}.Rmd",
        fit_DE_cond_nebula_CD4_CD8_done = "fit_DE_cond_nebula_CD4-CD8.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "GSEA_cond_nebula"
config_id = "CD4-CD8_A"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule GSEA_cond_nebula_CD4_CD8_A:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name}.Rmd",
        stats_DE_cond_nebula_CD4_CD8_done = "stats_DE_cond_nebula_CD4-CD8.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "GSEA_analysis_cond"
config_id = "CD4-CD8_A"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule GSEA_analysis_cond_CD4_CD8_A:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name}.Rmd",
        GSEA_cond_nebula_CD4_CD8_A_done = "GSEA_cond_nebula_CD4-CD8_A.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_base_name = "GSEA_cond_barplots"
config_id = "CD4-CD8_A"
rule_name = rule_base_name + "_" + config_id.replace('-', '_')

rule GSEA_cond_barplots_CD4_CD8_A:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_base_name}.Rmd",
        GSEA_cond_nebula_CD4_CD8_A_done = "GSEA_cond_nebula_CD4-CD8_A.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
        config_id=config_id
    threads: 6
    output:
        f"{rule_base_name}_{config_id}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}', params=list(config = '{params.config}', config_id = '{params.config_id}'), output_file = '{output}')\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "get_cell_counts"

rule get_cell_counts:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name}.Rmd",
        set_model_variables_done = "set_model_variables.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "DA_glm_nb"

rule DA_glm_nb:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name}.Rmd",
        get_cell_counts_done = "get_cell_counts.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

rule_name = "umap_plots_manuscript"

rule umap_plots_manuscript:
    input:
        config_check = f"{SNAKEDIR}.run_{rule_name}",
        script =  f"{rule_name}.Rmd",
        set_model_variables_done = "set_model_variables.html"
    params:
        config = SCRIPTS_DIR + "config.yaml",
        walltime="01:00",
        mem="10GB",
        queue="verylong",
        jobname=rule_name,
    threads: 6
    output:
        f"{rule_name}.html"
    run:
        print("Running " + str(input.script) + " with " + str(threads) + " cores.")
        shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
        write_rule_config(flat_new, rule_name, f".workflow_call/.{rule_name}.yaml")

if flat_new['run_liana']:
  rule liana:
      input:
          config_check = SNAKEDIR + '.run_liana',
          script =  "Liana_interactions.Rmd",
          set_model_variables_done = "set_model_variables.html"
      params:
          config = SCRIPTS_DIR + "config.yaml",
          walltime="01:00",
          mem="10GB",
          queue="verylong",
          jobname="liana"
      threads: 6
      output:
          "Liana_interactions.html"
      run:
          print("Running Liana_interactions.Rmd with " + str(threads) + " cores.")
          shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
          write_rule_config(flat_new, 'liana', '.workflow_call/.liana.yaml')
else:
  Path("Liana_interactions.html").touch()
        
if flat_new['run_shinycell']:
  rule shinycell:
      input:
          config_check = SNAKEDIR + '.run_shinycell',
          script =  "ShinyCell.Rmd",
          set_model_variables_done = "set_model_variables.html"
      params:
          config = SCRIPTS_DIR + "config.yaml",
          walltime="01:00",
          mem="10GB",
          queue="verylong",
          jobname="shinycell"
      threads: 6
      output:
          "ShinyCell.html"
      run:
          print("Running ShinyCell.Rmd with " + str(threads) + " cores.")
          shell(container_command + " Rscript -e \"rmarkdown::render('{input.script}',params=list(config = '{params.config}'))\"")
          write_rule_config(flat_new, 'shinycell', '.workflow_call/.shinycell.yaml')
else:
  Path("ShinyCell.html").touch()

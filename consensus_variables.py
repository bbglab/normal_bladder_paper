import json
import matplotlib as mpl

mpl.rcParams.update({
    'font.family': 'Arial',            # Enforce Arial
    'pdf.fonttype': 42,                # TrueType for PDF
    'ps.fonttype': 42,                 # TrueType for PS/EPS
    'svg.fonttype': 'none',            # Keep text as editable text
})



def load_samples_info(rundir):

    with open(f'{rundir}/table2group/all_groups.json', 'r') as f:
        groups_n_samples_definition = json.load(f)
    all_unit_names = list(groups_n_samples_definition.keys())

    with open(f'{rundir}/table2group/samples.json', 'r') as f:
        samples_definition = json.load(f)
    all_samples_names = list(samples_definition.keys())

    try:
        with open(f'{rundir}/table2group/groups.json', 'r') as f:
            groups_definition = json.load(f)
        all_group_names = list(groups_definition.keys())
    except:
        all_group_names = ["all_samples"]

    return all_samples_names, all_group_names, all_unit_names, groups_n_samples_definition

# Sample names
all_sample_names_dirty = ['P19_0001_BDO_01', 'P19_0001_BTR_01', 'P19_0002_BDO_01', 'P19_0002_BTR_01', 'P19_0003_BDO_01', 'P19_0004_BDO_01', 'P19_0004_BTR_01', 'P19_0005_BDO_01', 'P19_0005_BTR_01', 'P19_0006_BDO_01', 'P19_0007_BDO_01', 'P19_0007_BTR_01', 'P19_0008_BDO_01', 'P19_0008_BTR_01', 'P19_0009_BDO_01', 'P19_0009_BTR_01',
                          'P19_0011_BDO_01', 'P19_0011_BTR_01', 'P19_0012_BDO_01', 'P19_0012_BTR_01', 'P19_0013_BDO_01', 'P19_0013_BTR_01', 'P19_0014_BDO_01', 'P19_0014_BTR_01', 'P19_0015_BDO_01', 'P19_0015_BTR_01', 'P19_0016_BDO_01', 'P19_0016_BTR_01', 'P19_0018_BDO_01', 'P19_0018_BTR_01', 'P19_0019_BDO_01', 'P19_0019_BTR_01', 'P19_0020_BDO_01', 'P19_0020_BTR_01', 'P19_0023_BDO_01', 'P19_0023_BTR_01', 'P19_0024_BDO_01', 'P19_0024_BTR_01', 'P19_0025_BDO_01', 'P19_0025_BTR_01', 'P19_0026_BDO_01', 'P19_0026_BTR_01', 'P19_0027_BTR_01', 'P19_0028_BDO_01', 'P19_0028_BTR_01', 'P19_0029_BDO_01', 'P19_0029_BTR_01', 'P19_0030_BTR_01', 'P19_0031_BTR_01', 'P19_0033_BDO_01', 'P19_0033_BTR_01', 'P19_0034_BDO_01', 'P19_0034_BTR_01', 'P19_0035_BDO_01', 'P19_0036_BTR_01', 'P19_0038_BDO_01', 'P19_0039_BTR_01', 'P19_0040_BDO_01', 'P19_0040_BTR_01', 'P19_0041_BDO_01', 'P19_0041_BTR_01', 'P19_0042_BDO_01', 'P19_0042_BTR_01', 'P19_0043_BTR_01', 'P19_0045_BDO_01', 'P19_0045_BTR_01', 'P19_0046_BDO_01', 'P19_0046_BTR_01', 'P19_0047_BDO_01', 'P19_0047_BTR_01', 'P19_0048_BTR_01',
# new samples
'P19_0050_BDO_01', 'P19_0050_BTR_01', 'P19_0051_BDO_01','P19_0051_BTR_01',
'P19_0052_BDO_01', 'P19_0052_BTR_01', 'P19_0053_BDO_01', 'P19_0053_BTR_01']

repeated_donors = ['0001', '0002', '0004', '0005', '0007', '0008', '0009', '0010', '0011', '0012', '0013', '0014', '0015', '0016', '0018', '0019', '0020', '0023', '0024', '0025', '0026', '0028', '0029', '0033', '0034', '0040', '0041', '0042', '0045', '0046', '0047']

repeated_samples = ['P19_0001_BDO_01', 'P19_0001_BTR_01', 'P19_0002_BDO_01', 'P19_0002_BTR_01', 'P19_0004_BDO_01', 'P19_0004_BTR_01', 'P19_0005_BDO_01', 'P19_0005_BTR_01', 'P19_0007_BDO_01', 'P19_0007_BTR_01', 'P19_0008_BDO_01', 'P19_0008_BTR_01', 'P19_0009_BDO_01', 'P19_0009_BTR_01', 'P19_0010_BDO_01', 'P19_0010_BTR_01', 'P19_0011_BDO_01', 'P19_0011_BTR_01', 'P19_0012_BDO_01', 'P19_0012_BTR_01', 'P19_0013_BDO_01', 'P19_0013_BTR_01', 'P19_0014_BDO_01', 'P19_0014_BTR_01', 'P19_0015_BDO_01', 'P19_0015_BTR_01', 'P19_0016_BDO_01', 'P19_0016_BTR_01', 'P19_0018_BDO_01', 'P19_0018_BTR_01', 'P19_0019_BDO_01', 'P19_0019_BTR_01', 'P19_0020_BDO_01', 'P19_0020_BTR_01', 'P19_0023_BDO_01', 'P19_0023_BTR_01', 'P19_0024_BDO_01', 'P19_0024_BTR_01', 'P19_0025_BDO_01', 'P19_0025_BTR_01', 'P19_0026_BDO_01', 'P19_0026_BTR_01', 'P19_0028_BDO_01', 'P19_0028_BTR_01', 'P19_0029_BDO_01', 'P19_0029_BTR_01', 'P19_0033_BDO_01', 'P19_0033_BTR_01', 'P19_0034_BDO_01', 'P19_0034_BTR_01', 'P19_0040_BDO_01', 'P19_0040_BTR_01', 'P19_0041_BDO_01', 'P19_0041_BTR_01', 'P19_0042_BDO_01', 'P19_0042_BTR_01', 'P19_0045_BDO_01', 'P19_0045_BTR_01', 'P19_0046_BDO_01', 'P19_0046_BTR_01', 'P19_0047_BDO_01', 'P19_0047_BTR_01']

updated_sample_names = ['01_DO', '01_TR', '02_DO', '02_TR', '03_DO', '04_DO', '04_TR', '05_DO', '05_TR', '06_DO', '07_DO', '07_TR', '08_DO', '08_TR', '09_DO', '09_TR', '10_DO', '10_TR', '11_DO', '11_TR', '12_DO', '12_TR', '13_DO', '13_TR', '14_DO', '14_TR', '15_DO', '15_TR', '16_DO', '16_TR', '18_DO', '18_TR', '19_DO', '19_TR', '20_DO', '20_TR', '23_DO', '23_TR', '24_DO', '24_TR', '25_DO', '25_TR', '26_DO', '26_TR', '27_TR', '28_DO', '28_TR', '29_DO', '29_TR', '30_TR', '31_TR', '33_DO', '33_TR', '34_DO', '34_TR', '35_DO', '36_TR', '38_DO', '39_TR', '40_DO', '40_TR', '41_DO', '41_TR', '42_DO', '42_TR', '43_TR', '45_DO', '45_TR', '46_DO', '46_TR', '47_DO', '47_TR', '48_TR']

old2new_sample_names = {'P19_0001_BDO_01': '01_DO', 'P19_0001_BTR_01': '01_TR', 'P19_0002_BDO_01': '02_DO', 'P19_0002_BTR_01': '02_TR', 'P19_0003_BDO_01': '03_DO', 'P19_0004_BDO_01': '04_DO', 'P19_0004_BTR_01': '04_TR', 'P19_0005_BDO_01': '05_DO', 'P19_0005_BTR_01': '05_TR', 'P19_0006_BDO_01': '06_DO', 'P19_0007_BDO_01': '07_DO', 'P19_0007_BTR_01': '07_TR', 'P19_0008_BDO_01': '08_DO', 'P19_0008_BTR_01': '08_TR', 'P19_0009_BDO_01': '09_DO', 'P19_0009_BTR_01': '09_TR', 
'P19_0011_BDO_01': '11_DO', 'P19_0011_BTR_01': '11_TR', 'P19_0012_BDO_01': '12_DO', 'P19_0012_BTR_01': '12_TR', 'P19_0013_BDO_01': '13_DO', 'P19_0013_BTR_01': '13_TR', 'P19_0014_BDO_01': '14_DO', 'P19_0014_BTR_01': '14_TR', 'P19_0015_BDO_01': '15_DO', 'P19_0015_BTR_01': '15_TR', 'P19_0016_BDO_01': '16_DO', 'P19_0016_BTR_01': '16_TR', 'P19_0018_BDO_01': '18_DO', 'P19_0018_BTR_01': '18_TR', 'P19_0019_BDO_01': '19_DO', 'P19_0019_BTR_01': '19_TR', 'P19_0020_BDO_01': '20_DO', 'P19_0020_BTR_01': '20_TR', 'P19_0023_BDO_01': '23_DO', 'P19_0023_BTR_01': '23_TR', 'P19_0024_BDO_01': '24_DO', 'P19_0024_BTR_01': '24_TR', 'P19_0025_BDO_01': '25_DO', 'P19_0025_BTR_01': '25_TR', 'P19_0026_BDO_01': '26_DO', 'P19_0026_BTR_01': '26_TR', 'P19_0027_BTR_01': '27_TR', 'P19_0028_BDO_01': '28_DO', 'P19_0028_BTR_01': '28_TR', 'P19_0029_BDO_01': '29_DO', 'P19_0029_BTR_01': '29_TR', 'P19_0030_BTR_01': '30_TR', 'P19_0031_BTR_01': '31_TR', 'P19_0033_BDO_01': '33_DO', 'P19_0033_BTR_01': '33_TR', 'P19_0034_BDO_01': '34_DO', 'P19_0034_BTR_01': '34_TR', 'P19_0035_BDO_01': '35_DO', 'P19_0036_BTR_01': '36_TR', 'P19_0038_BDO_01': '38_DO', 'P19_0039_BTR_01': '39_TR', 'P19_0040_BDO_01': '40_DO', 'P19_0040_BTR_01': '40_TR', 'P19_0041_BDO_01': '41_DO', 'P19_0041_BTR_01': '41_TR', 'P19_0042_BDO_01': '42_DO', 'P19_0042_BTR_01': '42_TR', 'P19_0043_BTR_01': '43_TR', 'P19_0045_BDO_01': '45_DO', 'P19_0045_BTR_01': '45_TR', 'P19_0046_BDO_01': '46_DO', 'P19_0046_BTR_01': '46_TR', 'P19_0047_BDO_01': '47_DO', 'P19_0047_BTR_01': '47_TR', 'P19_0048_BTR_01': '48_TR',
# new samples
'P19_0050_BDO_01': '50_DO', 'P19_0050_BTR_01': '50_TR', 'P19_0051_BDO_01': '51_DO','P19_0051_BTR_01': '51_TR',
'P19_0052_BDO_01': '52_DO', 'P19_0052_BTR_01': '52_TR', 'P19_0053_BDO_01': '53_DO', 'P19_0053_BTR_01': '53_TR'
}


single_sample_per_donor = ['01_TR', '02_TR', '03_DO', '04_TR', '05_TR', '06_DO', '07_TR', '08_TR', '09_TR',
'11_TR', '12_TR', '13_TR', '14_TR', '15_TR', '16_TR', '18_TR', '19_TR', '20_TR', '23_TR', '24_TR', '25_TR', '26_TR', '27_TR', '28_TR', '29_TR', '30_TR', '31_TR', '33_TR', '34_TR', '35_DO', '36_TR', '38_DO', '39_TR', '40_TR', '41_TR', '42_TR', '43_TR', '45_TR', '46_TR', '47_TR', '48_TR',
# new samples
'50_TR', '51_TR', '52_TR', '53_TR']

females_age_ordered = ['33', '52', '43', '12', '06', '30', '20', '38', '09', '41', '19', '25', '05', '13', '02', '03', '29']
males_age_ordered = ['36', '01', '24', '31', '27', '40', '47', '50', '04', '26', '14', '35', '46', '11', '48', '28', '07', '15', '08', '34', '53', '39', '16', '51', '23', '45', '42', '18']

# paths
origin_dir = '/data/bbg'
bladder_results = f"{origin_dir}/nobackup/bladder_ts/results"
deepcsa_run_dir = f"{bladder_results}/2025-05-14_deepCSA_45_donors"
maf_file = f"{deepcsa_run_dir}/germline_somatic/all_samples.filtered.tsv.gz"
clean_maf_file = f"{deepcsa_run_dir}/clean_germline_somatic/all_samples.clean.mutations.tsv"
somatic_maf_file = f"{deepcsa_run_dir}/clean_somatic/all_samples.somatic.mutations.tsv"

clinvars_file = f"{origin_dir}/projects/bladder_ts/data/complete_cohort/samples_metadata/20250516_metadata_bladder.with_depths.tsv" #TODO: add a metadata.tsv file in the repo with all the metadata compiled

# oncodrive3d datasets
## TODO, make sure that this are absolute paths to the datasets
o3d_alt_datasets = "/data/bbg/nobackup/scratch/oncodrive3d/datasets_240506" # These "alt" are used to rerieve annotations in equivalent Uniprot ID that are missing in the MANE related ones
o3d_datasets = "/data/bbg/nobackup/scratch/oncodrive3d/datasets_mane_240506"
o3d_annotations = "/data/bbg/nobackup/scratch/oncodrive3d/annotations_mane_240506"




####
# genes
####
panel_all_genes = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
                   "RBM10","KDM6A","TP53","FGFR3","CDKN1A","FOXQ1",
                   "PIK3CA","TERTpromoter"
                  ]

gene_order_positive_selection = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
                                 "RBM10","KDM6A","TP53","CDKN1A","FOXQ1",
                                 # "FGFR3",
                                 # "PIK3CA","TERT"
                                ]
gene_order_omega_truncating_decreasing_males = ['RBM10', 'KDM6A', 'STAG2', 'KMT2D', 'ARID1A', 'CDKN1A', 'TP53', 'EP300', 'NOTCH2', 'CREBBP', 'FOXQ1', 'KMT2C', 'RB1', 'FGFR3']

panel_no_hotspots = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
                       "RBM10","KDM6A","TP53","FGFR3","CDKN1A","FOXQ1",
                    ]

gene_order_main = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
                   "RBM10","KDM6A","TP53","FGFR3","CDKN1A","FOXQ1",
                   "PIK3CA","TERT"
                  ]
gene_order = sorted(gene_order_main)

gene_order_main_with_promoter = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
                                 "RBM10","KDM6A","TP53","FGFR3","CDKN1A","FOXQ1",
                                 "PIK3CA","TERTpromoter"
                                ]

gene_order_mutation_rate_normal = ["CDKN1A","FOXQ1", "KDM6A", "RBM10", "KMT2D","STAG2",
                                   "ARID1A","TP53","EP300","CREBBP","NOTCH2","RB1", "FGFR3"]

gene_list_fig5 = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","RB1","RBM10","TP53","CDKN1A","FOXQ1", "KDM6A", "STAG2", "FGFR3"
             # "PIK3CA","TERT", "KMT2C", "FGFR3", 
             ]
genes_mut_epithelium_order = ['KMT2D', 'KDM6A', 'ARID1A', 'RBM10', 'FOXQ1',
                              'CDKN1A', 'STAG2', 'EP300', 'NOTCH2', 'CREBBP',
                              'KMT2C', 'RB1', 'TP53']

gene_order_agebias = [
'RBM10',
'KDM6A',
'KMT2D',
'TP53',
'STAG2',
'CDKN1A',
'ARID1A',
'CREBBP',
'EP300',
'NOTCH2',
'KMT2C',
'RB1',
'FOXQ1']

gene_order_sexbias = [
'RBM10',
'CDKN1A',
'ARID1A',
'STAG2',
'KDM6A',
'KMT2D',
'TP53',
'CREBBP',
'EP300',
'NOTCH2',
'KMT2C',
'RB1',
'FOXQ1']

genes_regressions = ["KMT2D","EP300","ARID1A","CREBBP","NOTCH2","KMT2C","STAG2","RB1",
                   "RBM10","KDM6A","TP53",
                  #  "FGFR3",
                   "CDKN1A","FOXQ1",
                  #  "PIK3CA",
                   "TERTpromoter",
                   "ALL_GENES"
                  ]

####
# clinvars
####
clinvar_id2name = {
   "BLADDER_LOCATION": "Bladder location",
   "dome": "Dome",
   "trigone": "Trigone",
   "AGE": "Age",
   "BMI": "BMI",
   "SEX": "Sex",
   "M": "Male",
   "F": "Female",
   "SMOKING_STATUS": "Smoking",
   "current": "Current",
   "former": "Former",
   "never": "Never",
   "HISTORY_OF_DRINKING": "Drinking",
   "CHEMOTHERAPY_HISTORY_SIMPLE": "Chemotherapy",
   "had_prior_chemotherapy": "Chemotherapy",
   "yes": "Yes",
   "no": "No"
}

clinvar2dummy = {
   "BLADDER_LOCATION": "is_dome",
   "SEX": "is_male",
   "HISTORY_OF_SMOKING": "history_smoking",
   "SMOKING_STATUS": "smoking_status",
   "HISTORY_OF_DRINKING": "history_drinking",
   "CHEMOTHERAPY_HISTORY_SIMPLE": "had_prior_chemotherapy"
}

# metrics
metric_id2name_ylabel = {
   "N_MUTS": "# Somatic\n mutations",
   "MUTRATE_MB": "muts / Mb"
}

metric_id2name_title = {
   "": ""
}

# color dictionaries
gene2color = {
 'KMT2D': '#fa6c8b',
 'ARID1A': '#6b57e3',
 'KDM6A': '#fffcad',
 'RBM10': '#b3e38c',
 'EP300': '#faac66',
 'STAG2': '#e66be4',
 'CREBBP': '#fdd76a',
 'NOTCH2': '#8fcff5',
 'FOXQ1': '#bababa',
 'CDKN1A': '#d3afe5',
 'KMT2C': '#66b966',
 'TP53': '#e48c7a',
 'RB1': '#c5d3f8',
 'FGFR3': '#c6e5ce',
 'PIK3CA': '#d9608f',
 'TERT': '#a65ea2',
 'pTERT': '#a65ea2',
 'TERTpromoter': '#a65ea2',
 'sample': '#D3D3D3',
 'total': '#D3D3D3',
 'ALL_GENES': '#D3D3D3'}

clinvar2color = {
   "BLADDER_LOCATION": {
      "trigone": "#E0E0E0",
      "dome": "#B5B5B5",
   },
   "AGE": {
      "min": "#C5FA8F",
      "max": "#61A856"
   },
   "AGE_above55": {
      "<55": "#C5FA8F",
      ">55": "#61A856"
   },
   "BMI": {
      "min": "#FACC9C",
      "max": "#FFA82E"
   },
   "SEX": {
      "F": "#E5BBFA",
      "M": "#9EC9F9"
      
   },
   "SMOKING_STATUS": {
      "never": "#DACBFF",
      "former": "#9767CF",
      "current": "#694A94"
   },
   "HISTORY_OF_SMOKING": {
      "no": "#DACBFF",
      "yes": "#694A94"   
   },
   "AGE_above55_smoking": {
      "<55 y.o.": "#C5FA8F", 
      ">55 y.o.\nNever\nsmoker": "#DACBFF",
      ">55 y.o.\nEver\nsmoker": "#694A94"   
   },
   "HISTORY_OF_DRINKING": {
      "no": "#F8D4E6",
      "yes": "#D8437F"
   },
   "CHEMOTHERAPY_HISTORY_SIMPLE": {
      "no": "#FAECAD",
      "yes": "#FFEC1D"  
   }
}

metrics_colors_dictionary = {"ofml"        : "viridis_r",
                             "ofml_score"  : "#6A33E0",
                             "omega_trunc" : "#FA5E32",
                             "omega_synon" : "#89E4A2",
                             "omega_miss"  : "#FABE4A",
                             "o3d_score"   : "#6DBDCC",
                             "o3d_cluster" : "skyblue",
                             "o3d_prob"    : "darkgray",
                             "frameshift"  : "#E4ACF4",
                             "inframe"     : "C5",
                             "hv_lines"    : "lightgray", # horizontal and vertical lines,
                             "hv_lines_needle" : "gray",
                             "needle_obs"  : "#003366",
                             "omega_miss_tert" : "#f5840c",
                             "omega_synon_tert": "#378c12",
                             "nonsense" : "#FA5E32",
                             "synonymous" : "#89E4A2",
                             "missense"  : "#FABE4A",
                             #"nonsense"    : "#FB8E6F",  
                             # "synonymous"  : "#ACECBD", 
                             #"missense"    : "#FBD180", 
                             "indel"       : "#ECC4F7", 
                             "splicing"    : "#A1C5DF",
                            }
# plot configs
plots_general_config = {

                        # fonsizes
                        "ylabel_fontsize": 6,
                        "xlabel_fontsize": 6,
                        "xylabel_fontsize": 6,
                        "title_fontsize": 7,
                        "xyticks_fontsize": 5,
                        "xticks_fontsize": 5,
                        "yticks_fontsize": 5,
                        "legend_fontsize": 5,
                        "annots_fontsize": 5,

                        "dot_size_scplot": 15,
                        "dot_size_coeffplot": 5,
                        "dot_sizebelow_coeffplot": 40,
                        "dot_color_coeffplot": "#D3D3D3",
                        "dot_colorabove_coeffplot": "#D62728",
                        "dot_colorbelow_coeffplot": "#f29c9e",
                        "dot_edgethres_coeffplot": 0.2,
                        "dot_edgewidth_coeffplot": 0.5
                        }


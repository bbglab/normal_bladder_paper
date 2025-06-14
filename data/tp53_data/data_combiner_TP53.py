import os
import pandas as pd
import argparse
import re
import sys

from utils import vartype, pyrimidine_alt, pyrimidine_ref


from bgreference import hg38

def update_pyrimidine_based(chrom, pos, base, genome = hg38):
    ref = genome(chrom, pos, size = 1)
    return (ref, base) if ref in ("T", "C") else (pyrimidine_alt[ref], pyrimidine_alt[base])


def fetch_intogen_mutations(file):
    return pd.read_csv(file,sep='\t')

def fetch_boostdm(file):
    boostdm = pd.read_csv(file, sep='\t')
    boostdm['chr'] = 'chr' + boostdm['chr'].astype('str')            
    return boostdm


def fetch_genie_mutations(file, pyrimidine_ref,pyrimidine_alt):
    genie = pd.read_csv(file,sep='\t')

    genie_transformed = []

    for k,r in genie.iterrows():
        r['ref'],r['alt'] = convert_pyrimidines(r['ref'], r['alt'], pyrimidine_ref, pyrimidine_alt)
        genie_transformed.append(r)

    genie_df = pd.DataFrame(genie_transformed)
    genie_df['oncogenic'] = genie_df['Oncogenicity databases'].apply(lambda x: 1 if 'cgi' in x or
                                                                                    'clinvar' in x or
                                                                                    'oncokb' in x else 0)
    return genie_df



def process_normal_TP53_mutations(normal_path, pyrimidine_ref,pyrimidine_alt):
    normal_muts = pd.read_csv(os.path.join(normal_path + '/clean_somatic',
                                           'all_samples.somatic.mutations.tsv'),
                              sep='\t')
    normal_muts = normal_muts[(normal_muts['canonical_SYMBOL']=='TP53') 
                                & (normal_muts['TYPE']=='SNV')]
    normal_muts = normal_muts[['CHROM','POS','REF','ALT','ALT_DEPTH','SAMPLE_ID',
                               'SYMBOL','canonical_Consequence_broader','Protein_position']].drop_duplicates()
    normal_muts =  normal_muts.rename(columns={'CHROM':'chr',
                                       'POS':'pos',
                                       'REF':'ref',
                                       'ALT':'alt',
                                       'SAMPLE_ID':'sample',
                                       'Protein_position':'prot_pos',
                                       'canonical_Consequence_broader':'consequence'})

    normal_muts_format = []

    for k,r in normal_muts.iterrows():
        ref,alt = convert_pyrimidines(r['ref'], r['alt'], pyrimidine_ref, pyrimidine_alt)
        normal_muts_format.append([r['chr'],int(r['pos']),ref,alt,r['sample'],r['prot_pos'],r['consequence']])

    normal_muts_format_df = pd.DataFrame(normal_muts_format,columns=['chr','pos','ref','alt','sample','prot_pos','consequence'])
    return normal_muts_format_df


def mutation_format(hg38_name, pyrimidine_ref, pyrimidine_alt):
    chr_base, location = hg38_name.split(':')
    chr='chr17'
    m = re.search(r'(\d+)([ACGT])>([ACGT])', location)
    if m is not None:
        pos = m.group(1)
        # ref,alt = m.group(2), m.group(3)
        ref,alt = convert_pyrimidines(m.group(2),m.group(3),pyrimidine_ref,pyrimidine_alt)
    else:
        pos,ref,alt = None,None,None
    return chr,pos,ref,alt


def process_saturation_data(saturation_file, pyrmidine_ref, pyrimidine_alt):
    saturation_data = pd.read_excel(saturation_file,sheet_name=2,skiprows=2)
    saturation_format = []

    for k,r in saturation_data.iterrows():
        chr,pos,ref,alt = mutation_format(r['hg38_genomic'], pyrmidine_ref, pyrimidine_alt)

        # this filters out insertions and deletions
        if pos is not None:
            saturation_format.append([chr,int(pos),ref,alt,r['transformed_score_enrich2'],r['p_adjusted (one-sided, all)']])

    saturation_format_df = pd.DataFrame(saturation_format,columns=['chr','pos','ref','alt','RFS','p'])
    return saturation_format_df


def count_mut_occurrence(mutations_df):
    return pd.DataFrame({'count':mutations_df.groupby(by=['chr','pos','ref','alt',
                                                            'prot_pos','consequence'])['sample'].count()}).reset_index()


def fix_chromosome(r):
    return 'chr' + str(r['chr'])

def convert_pyrimidines(ref, alt, pyrimidine_ref, pyrimidine_alt):
    converted_ref = ''
    converted_alt = ''

    if len(ref) == len(alt):
        for i in range(len(ref)):
            new_ref = pyrimidine_ref[ref[i]]
            if new_ref != ref[i]:
                converted_ref += new_ref
                converted_alt += pyrimidine_alt[alt[i]]
            else:
                converted_ref += ref[i]
                converted_alt += alt[i]

    else:
        converted_ref = ref
        converted_alt = alt

    return converted_ref, converted_alt


def reformat_site_selection_data(data_point, pyrimidine_ref, pyrimidine_alt):
    converted_ref, converted_alt = convert_pyrimidines(data_point['ref'], data_point['alt'], pyrimidine_ref, pyrimidine_alt)
    return [data_point['chr'],data_point['pos'],converted_ref,converted_alt,data_point['site_selection'],data_point['p_value']]

def reformat_normal_data(data_point, pyrimidine_ref, pyrimidine_alt):
    converted_ref, converted_alt = convert_pyrimidines(data_point['ref'], data_point['alt'], pyrimidine_ref, pyrimidine_alt)
    return [data_point['chr'],data_point['pos'],converted_ref,converted_alt,data_point['sample']]

def categorize_normal_mutations_counts(r, threshold=1):
    if r['count'] > threshold:
        return '>' + str(threshold)
    elif r['count'] == threshold:
        return str(threshold)
    else:
        return str(0)


def process_site_selection_TP53(normal_path, pyrimidine_ref, pyrimidine_alt):
    site_selection = pd.read_csv(os.path.join(normal_path + '/sitecomparisongloballoc',
                                            'all_samples.global_loc.site.comparison.tsv.gz'),
                                    sep='\t')
    print(len(site_selection[['CHROM','POS','ALT']].drop_duplicates()))
    site_selection = site_selection[(site_selection['GENE'] == 'TP53')&
                                    (site_selection['OBSERVED_MUTS'] > 0)]

    site_selection = site_selection[['CHROM','POS','REF','ALT','OBS/EXP','p_value']].rename(columns={'CHROM':'chr',
                                                                                'POS':'pos',
                                                                                'REF':'ref',
                                                                                'ALT':'alt',
                                                                                'OBS/EXP':'site_selection'})
    site_selection_format = []
    for k,r in site_selection.iterrows():
        ref,alt = convert_pyrimidines(r['ref'], r['alt'], pyrimidine_ref, pyrimidine_alt)
        site_selection_format.append([r['chr'],int(r['pos']),ref,alt,r['site_selection'],r['p_value']])

    site_selection_format_df = pd.DataFrame(site_selection_format,columns=['chr','pos','ref','alt','site_selection','p_value'])
    return site_selection_format_df

def convert_relative_position(r, TERT_gene_TSS):
    return float(r['pos']) - TERT_gene_TSS


def fetch_saturation_data(saturation_file):
    saturation_data = pd.read_csv(saturation_file, sep='\t')
    saturation_data['chr'] = 'chr' + saturation_data['chr'].astype(str)
    return saturation_data


def main(normal_version):

    #Directory with deepCSA output (normal mutations)
    normal_path = f'/data/bbg/nobackup/bladder_ts/results/{normal_version}'

    #Directory with processed data files needed to build the plots
    output_data_dir = '.'
    input_data_dir = 'inputs/'


    total_intogen_cancer = 33218
    total_intogen_bladder = 867
    total_normal = 79

    ##################################TP53 data###################
    #Processed_files
    data_dir = '.'
    outfile_TP53 = os.path.join(output_data_dir, 'TP53_merged_data.tsv')
    intogen_cancer_file = os.path.join(input_data_dir, 'TP53_mutations_intogen.tsv')
    intogen_bladder_file = os.path.join(input_data_dir, 'TP53_mutations_intogen_bladder.tsv')
    genie_file = os.path.join(input_data_dir,'TP53_mut_Genie_annotations.tsv')
    boostdm_cancer_file = os.path.join(input_data_dir,'TP53.model.CANCER.features.CANCER.prediction.tsv.gz')
    boostdm_bladder_file = os.path.join(input_data_dir,'TP53.model.BLCA.features.BLCA.prediction.tsv.gz')
    TP53_saturation_file = os.path.join(input_data_dir,'41588_2024_2039_MOESM4_ESM.xlsx')


    ######## TUMORS
    #Fetch intogen data
    intogen_cancer_muts = fetch_intogen_mutations(intogen_cancer_file)
    intogen_cancer_muts = intogen_cancer_muts.rename(columns={'total':'count_intogen_cancer'})

    print(str(len(intogen_cancer_muts)) + ' TP53 mutations in all tumor types in intogen')

    intogen_bladder_muts = fetch_intogen_mutations(intogen_bladder_file)
    intogen_bladder_muts = intogen_bladder_muts.rename(columns={'total':'count_intogen_bladder'})

    print(str(len(intogen_bladder_muts)) + ' TP53 mutations in bladder in intogen')

    #Merge intogen data
    intogen_merged_muts = pd.merge(intogen_cancer_muts,intogen_bladder_muts,on=['chr','pos','ref','alt','prot_pos','consequence'],how='outer')
    intogen_merged_muts = intogen_merged_muts.rename(columns={'consequence':'consequence_intogen'})
    intogen_merged_muts = intogen_merged_muts.fillna(0)

    print(str(len(intogen_merged_muts)) + ' TP53 mutations after merging cancer and bladder in intogen')

    #Fetch genie data
    genie_TP53_muts = fetch_genie_mutations(genie_file, pyrimidine_ref, pyrimidine_alt)
    genie_TP53_muts = genie_TP53_muts.rename(columns={'hg38_chromosome':'chr',
                                                        'hg38_pos':'pos',
                                                        'Consequence':'consequence_genie'})
    print(str(len(genie_TP53_muts)) + ' TP53 mutations in genie')

    #Merge tumor data
    tumor_data = pd.merge(intogen_merged_muts, genie_TP53_muts,
                            on=['chr','pos','ref','alt'],
                            how='outer')
    print(str(len(tumor_data)) + ' TP53 mutations after merging genie and intogen')

    tumor_data = tumor_data.rename(columns={'ref':'REF', 'alt':'ALT'})
    tumor_data["VARTYPE"] = tumor_data[["REF", "ALT"]].apply(vartype, axis = 1)
    tumor_data = tumor_data.rename(columns={'REF':'ref','ALT':'alt'})
    tumor_data = tumor_data[ tumor_data["VARTYPE"] == 'SNV'].reset_index(drop = True)
    tumor_data = tumor_data.drop([# 'HGVSc', 'HGVSp',
                                  # 'protein_change',
                                  # 'Transcript_ID', 'RefSeq',
                                  # 'Protein_position', 'consequence_genie', 'Variant_Classification',
                                  'Variant_Type'], axis = 'columns')
    print(str(len(tumor_data)) + ' TP53 mutations after filtering by SNVs')
    print()

    tumor_data.to_csv(f"{input_data_dir}/TP53.tumor_data.tsv", sep = '\t', index = False)



    ######## BoostDM
    #Fetch boostDM data
    boostdm_cols = ['chr','pos','alt','boostDM_score','boostDM_class']
    boostdm_cancer = fetch_boostdm(boostdm_cancer_file)
    boostdm_cancer = boostdm_cancer[boostdm_cols]
    boostdm_cancer = boostdm_cancer.rename(columns={'boostDM_score':'boostdm_score_cancer',
                                                        'boostDM_class':'boostdm_class_cancer'})

    print(str(len(boostdm_cancer)) + ' TP53 mutations in boostdm cancer')

    boostdm_bladder = fetch_boostdm(boostdm_bladder_file)
    boostdm_bladder = boostdm_bladder[boostdm_cols]
    boostdm_bladder = boostdm_bladder.rename(columns={'boostDM_score':'boostdm_score_bladder',
                                                        'boostDM_class':'boostdm_class_bladder'})
    print(str(len(boostdm_bladder)) + ' TP53 mutations in boostdm bladder')

    boostdm_data = pd.merge(boostdm_cancer, boostdm_bladder, on=['chr','pos','alt'], how='outer')
    print(str(len(boostdm_data)) + ' TP53 mutations in boostdm in total')

    boostdm_data[["ref","alt"]] = boostdm_data.apply(lambda x: update_pyrimidine_based(x['chr'], x['pos'], x['alt']), axis = 1).to_list()
    boostdm_data.to_csv(f"{input_data_dir}/TP53.boostdm_data.tsv", sep = '\t', index = False)


    #Merge with table data
    table_data = pd.merge(tumor_data, boostdm_data, on=['chr','pos','ref', 'alt'], how='outer')


    #Filter for variants considered in boostdm
    table_data = table_data.dropna(subset=['boostdm_class_cancer'])
    print(str(len(table_data)) + ' TP53 mutations after filtering only for variants included in boostdm')


    #Fetch experimental saturation mutagenesis data
    exp_saturation_data = process_saturation_data(TP53_saturation_file, pyrimidine_ref, pyrimidine_alt)
    print(str(len(exp_saturation_data)) + ' TP53 mutations with experimental functional impact data')


    #Merge with table data
    table_data = pd.merge(table_data, exp_saturation_data, on=['chr', 'pos', 'ref', 'alt'], how='outer')
    print(str(len(table_data)) + ' TP53 mutations after merging with experimental functional impact data')


    #Fetch normal data
    normal_TP53_muts = process_normal_TP53_mutations(normal_path, pyrimidine_ref, pyrimidine_alt)
    normal_TP53_muts = count_mut_occurrence(normal_TP53_muts)
    normal_TP53_muts = normal_TP53_muts.rename(columns={'consequence':'consequence_normal',
                                                        'count':'count_normal'})

    print(str(len(normal_TP53_muts)) + ' TP53 different mutated sites in normal samples')


    #Merge tumor and normal data
    table_data_normal = pd.merge(table_data, normal_TP53_muts, on=['chr','pos','ref','alt'],how='left')
    print(str(len(table_data_normal)) + ' TP53 different mutated sites after merging tumor and normal')

    #Extract meaningful columns
    cols2keep = ['chr','pos','ref','alt',
                    'consequence_intogen','count_intogen_cancer','count_intogen_bladder',
                    'consequence_genie','mut_patients','mut_freq','mut_bladder','mut_bladder_freq',
                    'boostdm_score_cancer','boostdm_class_cancer','boostdm_score_bladder','boostdm_class_bladder',
                    'oncogenic','consequence_normal','count_normal', 'RFS', 'p']
    table_data_normal = table_data_normal[cols2keep]
    table_data_normal['freq_intogen_cancer'] = table_data_normal['count_intogen_cancer'] / total_intogen_cancer
    table_data_normal['freq_intogen_bladder'] = table_data_normal['count_intogen_bladder'] / total_intogen_bladder

    table_data_normal['freq_normal'] = table_data_normal['count_normal'] / total_normal
    table_data_normal = table_data_normal.rename(columns={'mut_patients':'count_genie',
                                                                'mut_freq':'freq_genie',
                                                                'mut_bladder':'count_bladder_genie',
                                                                'mut_bladder_freq':'freq_bladder_genie',
                                                                'RFS' : 'experimental_score',
                                                                'p' : 'p_value_experimental'
                                                                })

    print(str(len(table_data_normal)) + ' TP53 different mutated sites after extracting meaningful columns')

    #TP53 site selection file
    TP53_site_selection = process_site_selection_TP53(normal_path, pyrimidine_ref,pyrimidine_alt)
    print(str(len(TP53_site_selection)) + ' TP53 different mutated sites with site selection value')

    #Merge with table data
    table_data_complete = pd.merge(table_data_normal,TP53_site_selection,on=['chr','pos','ref','alt'],how='left')
    print(str(len(table_data_complete)) + ' TP53 different mutated sites after merging with site selection data')


    #Fill count and freq columns
    table_data_complete[['count_genie','freq_genie','count_bladder_genie','freq_bladder_genie',
                        'count_intogen_cancer','count_intogen_bladder','count_normal',
                        'freq_intogen_cancer','freq_intogen_bladder','freq_normal']] = table_data_complete[['count_genie','freq_genie','count_bladder_genie','freq_bladder_genie',
                                                                                                    'count_intogen_cancer','count_intogen_bladder','count_normal',
                                                                                                    'freq_intogen_cancer','freq_intogen_bladder','freq_normal']].fillna(0)

    table_data_complete.to_csv(outfile_TP53,sep='\t',index=0)



if __name__ == '__main__':
    # parser = argparse.ArgumentParser()

    # parser.add_argument('-N',
    #                     '--normal_version',
    #                     required=True,
    #                     help='Name of version of deepCSA run')


    # args = parser.parse_args()

    # main(args.normal_version)

    main("2025-05-14_deepCSA_45_donors")
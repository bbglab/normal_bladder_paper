import os
import pandas as pd
import argparse
import re

from utils import vartype, pyrimidine_alt, pyrimidine_ref


def fetch_genie_TERT(file, TERT_info, pyrimidine_ref, pyrimidine_alt):
    genie_muts = pd.read_csv(file,sep='\t')
    TERT_genie_muts = genie_muts[(genie_muts['hg38_pos']>=TERT_info['start'])&
                                    (genie_muts['hg38_pos']<=TERT_info['end'])]

    TERT_genie_muts_transformed = []

    for k,r in TERT_genie_muts.iterrows():
        r['ref'],r['alt'] = convert_pyrimidines(r['ref'], r['alt'], pyrimidine_ref, pyrimidine_alt)
        TERT_genie_muts_transformed.append(r)

    TERT_genie_df = pd.DataFrame(TERT_genie_muts_transformed)
    return TERT_genie_df



def mutation_format(hg38_name, pyrimidine_ref, pyrimidine_alt):
    chr_base,location = hg38_name.split(':')
    chr='chr17'
    m = re.search(r'(\d+)([ACGT])>([ACGT])', location)
    if m is not None:
        pos = m.group(1)
        ref,alt = convert_pyrimidines(m.group(2),m.group(3),pyrimidine_ref,pyrimidine_alt)
    else:
        pos,ref,alt = None,None,None
    return chr,pos,ref,alt


def process_saturation_data(saturation_file, pyrmidine_ref, pyrimidine_alt):
    saturation_data = pd.read_excel(saturation_file,sheet_name=2,skiprows=2)
    saturation_format = []

    for k,r in saturation_data.iterrows():
        chr,pos,ref,alt = mutation_format(r['hg38_genomic'], pyrmidine_ref, pyrimidine_alt)

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


def process_normal_TERT_mutations(normal_path, TERT_info_38):
    normal_muts = pd.read_csv(os.path.join(normal_path + '/clean_somatic',
                                            'all_samples.somatic.mutations.tsv'),
                                sep='\t')
    normal_muts = normal_muts[['CHROM','POS','REF','ALT','ALT_DEPTH','SAMPLE_ID', "TYPE"]]
    normal_TERT_muts = normal_muts[(normal_muts['CHROM']==('chr' + TERT_info_38['chr'])) &
                                    (normal_muts['POS'] >= int(TERT_info_38['start'])) &
                                    (normal_muts['POS'] <= int(TERT_info_38['end'])) &
                                    (normal_muts['TYPE'] == "SNV")
                                    ]

    return normal_TERT_muts.rename(columns={'CHROM':'chr',
                                            'POS':'pos',
                                            'REF':'ref',
                                            'ALT':'alt',
                                            'SAMPLE_ID':'sample'})


def process_site_selection_TERT(normal_path, normal_version):
    site_selection = pd.read_csv(os.path.join(normal_path  + '/sitecomparisongloballoc',
                                            'all_samples.global_loc.site.comparison.tsv.gz'),
                                    sep='\t')
    site_selection = site_selection[site_selection['GENE']=='TERTpromoter']
    site_selection = site_selection[['CHROM','POS','REF','ALT','OBS/EXP','p_value']].rename(columns={'CHROM':'chr',
                                                                                            'POS':'pos',
                                                                                            'REF':'ref',
                                                                                            'ALT':'alt',
                                                                                            'OBS/EXP':'site_selection'})
    return site_selection


def convert_relative_position(r, TERT_gene_TSS):
    return float(r['pos']) - TERT_gene_TSS

def fetch_normal_TERT_mutations(normal_file):
    return pd.read_csv(normal_file, sep='\t')

def fetch_pcaw_TERT_mutations(pcawg_file):
    return pd.read_csv(pcawg_file, sep='\t')


def fetch_hartwig_TERT_mutations(hartwig_file):
    return pd.read_csv(hartwig_file, sep='\t')

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

    # total_intogen_cancer = 33218
    # total_intogen_bladder = 867
    # total_normal = 79

    ##################################TERT data###################
    #Processed data files
    outfile_TERT = os.path.join(output_data_dir,'TERT_merged_data.tsv')
    # pcawg_file = os.path.join(processed_data_dir,'TERT_mutations_pcawg.tsv')
    # hartwig_file = os.path.join(processed_data_dir,'TERT_mutations_hartwig.tsv')
    # TERT_genie_file = os.path.join(external_data_dir,'TERT_mut_Genie.tsv')
    TERT_tumors_file = os.path.join(input_data_dir,'TERT.tumor_data.tsv')
    saturation_file = os.path.join(input_data_dir,'TERT_saturation_mutagenesis_data.tsv')

    #TERT genomic coordinates
    TERT_info_38 = {'chr':'5',
                    'start':1294942,
                    'end':1295289}

    TERT_gene_TSS = 1295068

    # Regions in consensus panel
    # chr5	1294823	1295093
    # chr5	1295112	1295289



    #### TUMOR data
    total_hartwig = 5582
    total_pcawg = 2554
    total_tumors = total_hartwig + total_pcawg

    #Fetch tumor data
    # pcawg_TERT_muts = fetch_pcaw_TERT_mutations(pcawg_file)
    # hartwig_TERT_muts = fetch_pcaw_TERT_mutations(hartwig_file)
    # tumor_TERT_muts = pd.concat([pcawg_TERT_muts,hartwig_TERT_muts])
    # tumor_TERT_muts['chromosome'] = tumor_TERT_muts.apply(fix_chromosome,axis=1)
    # tumor_TERT_muts = tumor_TERT_muts.drop('chr',axis=1)
    # tumor_TERT_muts = tumor_TERT_muts.rename(columns={'chromosome':'chr'})

    # #Process tumor data
    # unique_tumor_TERT_muts = tumor_TERT_muts.groupby(['chr','pos','ref','alt'])['sample'].count().reset_index().rename(columns={'sample':'count_tumor'})
    # unique_tumor_TERT_muts['freq_tumor'] = unique_tumor_TERT_muts['count_tumor'] / total_tumors
    # print(str(len(unique_tumor_TERT_muts)) + ' TERT promoter mutations in tumor whole-genomes')

    # #Fetch genie data
    # genie_TERT_muts = fetch_genie_TERT(TERT_genie_file, TERT_info_38, pyrimidine_ref, pyrimidine_alt)
    # genie_TERT_muts = genie_TERT_muts[['hg38_chromosome','hg38_pos','ref','alt','Consequence',
    #                                     'mut_patients', 'mut_freq', 'mut_bladder', 'mut_bladder_freq']]
    # genie_TERT_muts = genie_TERT_muts.rename(columns={'hg38_chromosome':'chr',
    #                                                     'hg38_pos':'pos'})
    # print(str(len(genie_TERT_muts)) + ' TERT promoter mutations in GENIE samples')


    # #Merge tumor data
    # tumor_data = pd.merge(unique_tumor_TERT_muts,genie_TERT_muts,on=['chr','pos','ref','alt'],how='outer')
    # print(str(len(tumor_data)) + ' TERT promoter mutations across tumor samples')

    # tumor_data = tumor_data.rename(columns={'ref':'REF', 'alt':'ALT'})
    # tumor_data["VARTYPE"] = tumor_data[["REF", "ALT"]].apply(vartype, axis = 1)
    # tumor_data = tumor_data.rename(columns={'REF':'ref','ALT':'alt'})
    # tumor_data = tumor_data[ tumor_data["VARTYPE"] == 'SNV'].reset_index(drop = True)
    # print(str(len(tumor_data)) + 'SNVs in the TERT promoter across tumor samples')

    # tumor_data = tumor_data.rename({"count_tumor" : "count_WGS_tumor",
    #                                 "freq_tumor" : "freq_WGS_tumor",
    #                                 'mut_patients' : 'genie_mut_patients',
    #                                 'mut_freq' : 'genie_mut_freq',
    #                                 'mut_bladder' : 'genie_mut_bladder',
    #                                 'mut_bladder_freq' : 'genie_mut_bladder_freq'}, axis = 1)
    # tumor_data[['count_WGS_tumor', 'freq_WGS_tumor',
    #             'genie_mut_patients', 'genie_mut_freq',
    #             'genie_mut_bladder', 'genie_mut_bladder_freq']] = tumor_data[['count_WGS_tumor', 'freq_WGS_tumor',
    #                                                                             'genie_mut_patients', 'genie_mut_freq',
    #                                                                             'genie_mut_bladder', 'genie_mut_bladder_freq']].fillna(0)

    tumor_data = pd.read_table(TERT_tumors_file, sep = '\t')
    print(str(len(tumor_data)) + ' SNVs in the TERT promoter across tumor samples')


    #Fetch normal data
    normal_TERT_muts = process_normal_TERT_mutations(normal_path, TERT_info_38)
    normal_data_list = [reformat_normal_data(r, pyrimidine_ref, pyrimidine_alt) for k, r in normal_TERT_muts.iterrows()]
    normal_data_format = pd.DataFrame(normal_data_list,columns=['chr','pos','ref','alt','sample'])

    #Count unique mutations
    unique_normal_TERT_muts = pd.DataFrame({'count_normal':normal_data_format.groupby(['chr','pos','ref','alt'])['sample'].count()}).reset_index()
    print(str(len(unique_normal_TERT_muts)) + ' TERT promoter mutations in normal samples')


    #Fetch site_selection data
    site_selection_data = process_site_selection_TERT(normal_path, normal_version)
    site_selection_data_list = [reformat_site_selection_data(r, pyrimidine_ref, pyrimidine_alt) for k,r in site_selection_data.iterrows()]
    site_selection_data_format = pd.DataFrame(site_selection_data_list,columns=['chr','pos','ref','alt','site_selection','p_value'])

    print(str(len(site_selection_data_format)) + ' TERT promoter mutations with site selection data')

    #Merge normal data and site selection
    normal_site_selection = pd.merge(site_selection_data_format, unique_normal_TERT_muts, on=['chr','pos','ref','alt'],how='left')
    normal_site_selection["count_normal"] = normal_site_selection["count_normal"].fillna(0)
    normal_site_selection["in_consensus"] = 1
    print(str(len(normal_site_selection)) + ' possible TERT promoter mutations defined by site selection data')

    tumor_normal_selection = pd.merge(normal_site_selection,tumor_data,on=['chr','pos','ref','alt'],how='outer')
    print(str(len(tumor_normal_selection)) + ' TERT promoter mutations after merging normal and tumor data')

    #Fetch saturation data
    saturation_data = fetch_saturation_data(saturation_file)
    print(str(len(saturation_data)) + ' TERT promoter mutations from the experimental saturation data before filtering')

    saturation_data = saturation_data.rename(columns={'ref':'REF', 'alt':'ALT'})
    saturation_data["VARTYPE"] = saturation_data[["REF", "ALT"]].apply(vartype, axis = 1)
    saturation_data = saturation_data.rename(columns={'REF':'ref','ALT':'alt'})
    saturation_data = saturation_data[ saturation_data["VARTYPE"] == 'SNV'].drop("VARTYPE", axis = 1).reset_index(drop = True)

    print(str(len(saturation_data)) + ' TERT promoter mutations from the experimental saturation data')

    saturation_data = saturation_data[(saturation_data['pos']>=TERT_info_38['start']) &
                                        (saturation_data['pos']<=TERT_info_38['end']) &
                                        (saturation_data["system"] == "TERT-GBM")
                                        ].reset_index(drop = True)
    saturation_data["in_saturation"] = 1

    print(str(len(saturation_data)) + ' TERT promoter mutations from the experimental saturation data')


    #Merge into final data table
    data_table = pd.merge(tumor_normal_selection, saturation_data,on=['chr','pos','ref','alt'],how='outer')
    print(str(len(data_table)) + ' TERT promoter mutations after merging normal, tumor and mutagenesis data')

    data_table["in_consensus"] = data_table["in_consensus"].fillna(0)
    data_table["in_saturation"] = data_table["in_saturation"].fillna(0)

    data_table = data_table[(data_table['pos']>=TERT_info_38['start'])&
                                                    (data_table['pos']<=TERT_info_38['end'])]
    print(str(len(data_table)) + ' TERT promoter mutations after filtering normal, tumor and mutagenesis data')


    data_table.to_csv(outfile_TERT,sep='\t',index=0)
    print('\n\n\n\n')



if __name__ == '__main__':
    # parser = argparse.ArgumentParser()

    # parser.add_argument('-N',
    #                     '--normal_version',
    #                     required=True,
    #                     help='Name of version of deepCSA run')


    # args = parser.parse_args()

    # main(args.normal_version)
    main("2025-05-14_deepCSA_45_donors")

# %%
import sys
import pandas as pd
import customeDic

help = '''
    Usage: python3 AS_APA_detect.py [integration_summary] [feature_table]
'''
version = '''
'''

if len(sys.argv) == 1:
    print(help)
    sys.exit()

for i in range(1, len(sys.argv)):
    if sys.argv[i] == '-v':
        print(version)
        sys.exit()
    elif sys.argv[i] == '-h':
        print(help)
        sys.exit()

integration_s = sys.argv[1]
feature_table = sys.argv[2]

# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)

# pd.reset_option('display.max_rows')
# pd.reset_option('display.max_columns')

# for debug
# integration_s = '/home/user/volumes/data/as_apa/all.integration.s'
# feature_table = '/home/user/volumes/data/genome/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_feature_table.txt'

print('Loading integration summary file')
df_integration = pd.read_table(integration_s)
print('Extract transcript id ...')
df_integration = df_integration[~ df_integration['feature'].isin(['Intergenic', 'intron', 'flnc'])]
print('Grab all transcript id ...')
integration_transcript_list = df_integration['PHAS_gene'].drop_duplicates().tolist()

print('Construction Locus-transcript hash table')
df_feature_table = pd.read_table(feature_table, low_memory=False)
df_feature_table = df_feature_table[['product_accession', 'related_accession', 'GeneID']]
df_feature_table.fillna('NA', inplace=True)

Locus_dic = customeDic.OneSetDic()
Transcript_dic = customeDic.nestedDic()

for r_index, row in df_feature_table.iterrows():
    Locus = 'LOC' + str(row['GeneID'])
    product_accession = row['product_accession']
    related_accession = row['related_accession']
    if product_accession != 'NA' and 'P' not in product_accession:
        Locus_dic[Locus].add(product_accession)
        Transcript_dic[product_accession] = Locus
    if related_accession != 'NA' and 'P' not in related_accession:
        Locus_dic[Locus].add(related_accession)
        Transcript_dic[related_accession] = Locus

print('Statics ...')
transcript_can_generating_phasiRNA_dic = customeDic.OneSetDic()

for transcript in integration_transcript_list:
    Locus = Transcript_dic[transcript]
    transcript_can_generating_phasiRNA_dic[Locus].add(transcript)
    Locus_dic[Locus].remove(transcript)

dic = {'Locus':[], 'Transcript_can_generating_phasiRNA': [], 'Transcript_can_not_generating_phasiRNA': [], 'Total': [], 'percentage_generating_phasiRNA_transcript': [], 'phasiRNA_transcript': [], 'non_phasiRNA_transcript': []}
# print(f'Locus\tTranscript_can_generating_phasiRNA\tTranscript_can_not_generating_phasiRNA\ttotal_number\tpercentage_generating_phasiRNA_transcript')
for Locus in transcript_can_generating_phasiRNA_dic:
    Transcript_can_generating_phasiRNA = len(transcript_can_generating_phasiRNA_dic[Locus])
    Transcript_can_not_generating_phasiRNA = len(Locus_dic[Locus])
    Total = Transcript_can_generating_phasiRNA + Transcript_can_not_generating_phasiRNA
    Percentage = Transcript_can_generating_phasiRNA / Total
    dic['Locus'].append(Locus)
    dic['Transcript_can_generating_phasiRNA'].append(Transcript_can_generating_phasiRNA)
    dic['Transcript_can_not_generating_phasiRNA'].append(Transcript_can_not_generating_phasiRNA)
    dic['Total'].append(Total)
    dic['percentage_generating_phasiRNA_transcript'].append(Percentage)
    dic['phasiRNA_transcript'].append(';'.join(transcript_can_generating_phasiRNA_dic[Locus]))
    dic['non_phasiRNA_transcript'].append(';'.join(Locus_dic[Locus]))
    # print(f"{Locus}\t{Transcript_can_generating_phasiRNA}\t{Transcript_can_not_generating_phasiRNA}\t{Total}\t{Percentage}")

df = pd.DataFrame(dic)
df = df[df['percentage_generating_phasiRNA_transcript'] < 1]
dv = df.to_csv('/home/user/volumes/data/as_apa/as_apa.txt', sep='\t', header=True, index=False)

for i in df['phasiRNA_transcript']:
    print('\n'.join(i.split(';')))
for i in df['non_phasiRNA_transcript']:
    print('\n'.join(i.split(';')))
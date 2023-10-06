# %%
import sys
import pandas as pd
import customeDic
import os

help = '''
    Usage: python3 generation_capacity.py [integration_summary] [gff3] [output]
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
output = sys.argv[3]

# pd.set_option('display.max_rows', None)
# pd.set_option('display.max_columns', None)

# pd.reset_option('display.max_rows')
# pd.reset_option('display.max_columns')

# for debug
# feature_table = '/home/user/volumes/data/test_osa/oryza_sativa_gdna.gff3'
# integration_s = '/home/user/volumes/data/test_osa/integration_s.txt'
# output = 'tmp'

print('Loading integration summary file')
df_integration = pd.read_table(integration_s)
print('Extract transcript id ...')
df_integration = df_integration[~ df_integration['feature'].isin(['Intergenic', 'intron', 'flnc'])]
print('Grab all transcript id ...')
integration_transcript_list = df_integration['PHAS_gene'].drop_duplicates().tolist()
print('Construction Locus-transcript hash table')
df_feature_table = pd.read_csv(feature_table, sep='\t', comment='#', header=None)
print('Parse gff3 file')
gene_index = ~df_feature_table[2].str.contains('Gene', case=False) & ~df_feature_table[2].str.contains('region', case=False) & ~df_feature_table[2].str.contains('cds', case=False) & ~df_feature_table[2].str.contains('exon', case=False)
df_feature_table = df_feature_table[gene_index]
# catupre gene and transcript_id 
df_feature_table = df_feature_table[8].str.extract(r'.*(gene=\w*).*(transcript_id=[^;]*)').dropna().drop_duplicates()
df_feature_table[0] = df_feature_table[0].str.replace('gene=', '', case=False).astype(str)
df_feature_table[1] = df_feature_table[1].str.replace('transcript_id=', '', case=False).astype(str)

Locus_dic = customeDic.OneSetDic()
Transcript_dic = customeDic.nestedDic()

for ri, r in df_feature_table.iterrows():
    Locus = r[0]
    transcript = r[1]
    Locus_dic[Locus].add(transcript)
    Transcript_dic[transcript] = Locus


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
df.to_csv(output, sep='\t', header=True, index=False)
cmd = f"sed -i '1i # percentage_generating_phasiRNA_transcript < 1 means only partial transcripts from this locus can generate phasiRNAs.\\n# This phenomenon may be caused by alternative splicing or alternative polyadenylation event' {output}"
os.system(cmd)


# for i in df['phasiRNA_transcript']:
#     print('\n'.join(i.split(';')))
# for i in df['non_phasiRNA_transcript']:
#     print('\n'.join(i.split(';')))
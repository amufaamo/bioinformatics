# データベースの作成
## データベースのリストを表示する
```bash
update_blastdb.pl --showall
```
Output
```bash
16S_ribosomal_RNA
28S_fungal_sequences
18S_fungal_sequences
Betacoronavirus
ITS_RefSeq_Fungi
LSU_eukaryote_rRNA
ITS_eukaryote_sequences
LSU_prokaryote_rRNA
SSU_eukaryote_rRNA
env_nt
env_nr
human_genome
landmark
mito
mouse_genome
nr
nt
nt_others
nt_euk
nt_viruses
nt_prok
pataa
patnt
pdbaa
pdbnt
ref_euk_rep_genomes
ref_prok_rep_genomes
ref_viroids_rep_genomes
ref_viruses_rep_genomes
refseq_select_rna
refseq_select_prot
refseq_protein
refseq_rna
swissprot
tsa_nr
tsa_nt
taxdb
core_nt
```

## ダウンロード
```bash
update_blastdb.pl --decompress swissprot
```

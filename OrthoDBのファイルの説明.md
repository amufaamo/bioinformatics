
### 概要

まず、OrthoDBのデータは大きく分けて2種類あります。

* **タブ区切りファイル (`.tab.gz`)**: 生物種、遺伝子、オルソロググループなどの情報が、表形式でまとめられています。Excelなどで開くこともできますが、プログラムで扱うのが一般的です。**各ファイルに列名（ヘッダー）は含まれていない**ので、この説明がとても重要になりますね！
* **FASTA形式ファイル (`.fasta.gz`, `.fasta.tgz`)**: 遺伝子やタンパク質の配列データが格納されています。

それでは、一つずつ見ていきましょう！

### タブ区切りファイル (.tab) の詳細

#### `odb12v1_species.tab.gz`
生物種に関する基本的な情報がまとまっています。
1.  `NCBI tax id`: NCBIのタクソノミーID（生物の分類ID）です。
2.  `OrthoDB organism id`: OrthoDBが独自に付けている生物のIDです。
3.  `scientific name`: 学名です。（例: *Homo sapiens*）
4.  `genome assembly id`: ゲノムアセンブリのIDです。
5.  `clustered genes count`: この生物種でクラスタリングされた遺伝子の総数です。
6.  `OGs count`: この生物種が属しているオルソロググループ（OG）の総数です。
7.  `mapping type`: マッピングのタイプです。`C`はClustered（クラスタリングされた）、`M`はMapped（マッピングされた）を意味します。

#### 例
```
9       9_0     Buchnera aphidicola     GCF_900128725.1                 
17      417_0   Methylophilus methylotrophus    GCF_003854875.1                 
23      23_0    Shewanella colwelliana  GCF_001735525.1                 
25      25_0    Shewanella hanedai      GCF_007197645.1                 
33      33_0    Myxococcus fulvus       GCF_007991095.1                 
41      41_0    Stigmatella aurantiaca  GCF_900109545.1                 
48      48_0    Archangium gephyra      GCF_001027285.1
```

#### `odb12v1_levels.tab.gz`
オルソロググループ（OG）が計算された分類階層（クレード）に関するファイルです。例えば、「哺乳類」や「真核生物」といった大きなくくりですね。
1.  `level NCBI tax id`: その分類階層のNCBIタクソノミーIDです。
2.  `scientific name`: その分類階層の学名です。
3.  `genes count`: その階層に属する全ての生物の遺伝子数（重複なし）です。
4.  `OGs count`: その階層で構築されたオルソロググループ（OG）の総数です。
5.  `species count`: その階層に属する生物種の数（重複なし）です。

#### 例
```
2       Bacteria        66010659        644134  17551
18      Pelobacter      12014   2745    2
22      Shewanella      209206  10373   54
32      Myxococcus      35831   7715    5
68      Lysobacter      124862  8596    39
75      Caulobacter     129918  8858    31
81      Hyphomicrobium  39398   5059    12
85      Hyphomonas      41951   4819    13
126     Planctomycetaceae       56971   7034    34
136     Spirochaetales  141419  10583   68
137     Spirochaetaceae 119117  9915    15
138     Borrelia        8223    1155    9
157     Treponema       63438   6715    27
170     Leptospiraceae  243681  9227    67
171     Leptospira      238356  9206    65
191     Azospirillum    127786  9749    23
```

#### `odb12v1_level2species.tab.gz`
最上位の分類階層と、そこに属する各生物種との関係を示しています。
1.  `top-most level NCBI tax id`: 最上位の分類階層のIDです。{2(細菌), 2157(古細菌), 2759(真核生物), 10239(ウイルス)}のいずれかです。
2.  `OrthoDB organism id`: OrthoDBの生物IDです。
3.  `hops`: 最上位からその生物の階層までの、分類階層をいくつ経由したかを示す数です。
4.  `intermediate levels`: 最上位からその生物までの、中間の分類階層のリストです。
#### 例
```
2       1052212_0       1       {2,1052212}
2       1052212_1       1       {2,1052212}
2       1069642_0       1       {2,1069642}
2       1078050_0       1       {2,1078050}
2       1104602_0       1       {2,1104602}
2       1104602_1       1       {2,1104602}
2       1120746_0       1       {2,1120746}
2       1121390_0       1       {2,1121390}
2       1121391_0       1       {2,1121391}
2       1121399_0       1       {2,1121399}
2       1121402_0       1       {2,1121402}
2       1121403_0       1       {2,1121403}
2       1121404_0       1       {2,1121404}
2       1121409_0       1       {2,1121409}
```


#### `odb12v1_genes.tab.gz`
遺伝子一つ一つに関する詳細な情報です。
1.  `OrthoDB unique gene id`: OrthoDBが付けたユニークな遺伝子IDです（**注意：バージョンが変わるとIDも変わる可能性があります**）。
2.  `OrthoDB organism id`: この遺伝子を持つ生物のIDです。
3.  `protein original sequence id`: 元々のタンパク質配列のIDです。
4.  `synonyms`: 遺伝子名の同義語リスト（セミコロン区切り）です。
5.  `Uniprot id`: UniprotのIDです。
6.  `Ensembl ids`: EnsemblのIDリスト（セミコロン区切り）です。
7.  `NCBI gid or gene name`: NCBIの遺伝子IDまたは遺伝子名です。
8.  `description`: 遺伝子の機能などの説明です。
9.  `genomic coordinates`: ゲノムDNA上での位置（座標）です。
10. `genomic DNA id`: ゲノムDNAのIDです。
11. `chromosome`: 染色体名です。


#### 例
```
100_0:000000    100_0   WP_131833652.1  EV667_RS01940   A0A4R1I8Y6                      ATP-binding cassette domain-containing protein  [406224:407052](-)      NZ_SMFY01000001.1  Unknown
100_0:000001    100_0   WP_131834425.1  EV667_RS06425   A0A4R1IBT8                      extracellular solute-binding protein    [1412675:1413725](+)    NZ_SMFY01000001.1       Unknown
100_0:000002    100_0   WP_131835073.1  EV667_RS07960                           GNAT family N-acetyltransferase [1718702:1719233](+)    NZ_SMFY01000001.1       Unknown
100_0:000003    100_0   WP_131834175.1  EV667_RS05015   A0A4V2PK40                      ABC transporter permease        [1083734:1084541](-)    NZ_SMFY01000001.1       Unknown
100_0:000004    100_0   WP_131834578.1  EV667_RS07250   A0A4R1IAQ2                      Biopolymer transport protein ExbB       [1572487:1573651](+)    NZ_SMFY01000001.1       Unknown
100_0:000005    100_0   WP_131834150.1  EV667_RS04890   A0A4R1IB33                      sigma-70 family RNA polymerase sigma factor     [1055373:1055874](+)    NZ_SMFY01000001.1  Unknown
100_0:000006    100_0   WP_131834172.1  EV667_RS05000                           ribulose 1,5-bisphosphate carboxylase   [1079260:1080520](-)    NZ_SMFY01000001.1       Unknown
100_0:000007    100_0   WP_131834248.1  EV667_RS05440   A0A4R1I9N8                      DUF2188 domain-containing protein       [1181606:1181873](-)    NZ_SMFY01000001.1       Unknown
```

#### `odb12v1_gene_xrefs.tab.gz`
OrthoDBの遺伝子IDと、他のデータベースのIDとを関連付ける「相互参照」ファイルです。
1.  `OrthoDB gene id`: OrthoDBの遺伝子IDです。
2.  `external gene identifier`: 外部データベースでの遺伝子IDです。
3.  `external DB name`: 外部データベースの名前です。{GOterm, InterPro, NCBIproteinGI, UniProt, ENSEMBL, NCBIgid, NCBIgenename}のいずれかが入ります。

### 例
```
100_0:000000    GO:0005886      GOterm
100_0:000000    GO:0000166      GOterm
100_0:000000    GO:0016887      GOterm
100_0:000000    GO:0005524      GOterm
100_0:000000    IPR017871       InterPro
100_0:000000    IPR003439       InterPro
100_0:000000    IPR003593       InterPro
100_0:000000    IPR027417       InterPro
100_0:000000    WP_131833652.1  NCBIproteinAcc
100_0:000000    A0A4R1I8Y6      UniProt
100_0:000001    GO:0030288      GOterm
100_0:000001    GO:0046872      GOterm
100_0:000001    IPR026045       InterPro
100_0:000001    WP_131834425.1  NCBIproteinAcc
100_0:000001    A0A4R1IBT8      UniProt
```

#### `odb12v1_OGs.tab.gz`
オルソロググループ（OG）自体の情報です。
1.  `OG unique id`: OGのユニークIDです（**注意：バージョンが変わるとIDも変わる・再利用される可能性があります**）。
2.  `level tax_id`: このOGがどの分類階層で定義されたかを示すタクソノミーIDです。
3.  `OG name`: OGの名前です。多くの場合、グループ内で最も一般的な遺伝子名が付けられます。

#### `odb12v1_OG2genes.tab.gz`
どのオルソロググループ（OG）に、どの遺伝子が含まれているかの対応リストです。とても重要なファイルですね！
1.  `OG unique id`: OGのユニークIDです。
2.  `OrthoDB gene id`: OrthoDBの遺伝子IDです。

#### `odb12v1_OG_pairs.tab.gz`
オルソロググループ（OG）間の階層関係を示します。「哺乳類のOG」は「脊椎動物のOG」の子孫である、といった関係が分かります。
1.  `descendant OG id`: 子孫にあたるOGのIDです。
2.  `antecedent OG id`: 祖先にあたるOGのIDです。

#### `odb12v1_OG_xrefs.tab.gz`
オルソロググループ（OG）と、GO（遺伝子オントロジー）やInterPro（タンパク質ドメイン）などの外部情報を関連付けます。
1.  `OG unique id`: OGのユニークIDです。
2.  `external DB or DB section`: 外部データベース名などです。
3.  `external identifier`: 外部データベースでのID（例: GOタームのID）です。
4.  `number of genes`: その外部IDを持つ遺伝子が、このOG内にいくつあるかを示します。

---

### FASTA形式ファイルの詳細

FASTA形式は、`>`から始まる行（ヘッダー行）と、次の行からの配列データで構成されます。ヘッダーにはOrthoDBの内部遺伝子IDと公開IDが含まれているので、どの遺伝子の配列かすぐに分かりますね！

* **`odb12v1_aa_fasta.gz`**: 全ての生物の、全ての遺伝子の**アミノ酸配列**が入っています。
* **`odb12v1_og_aa_fasta.gz`**: オルソロググループ（OG）に含まれている遺伝子の**アミノ酸配列**のみが入っています。系統解析などでよく使います。
* **`odb12v1_cds_fasta.gz`**: 全ての生物の、全ての遺伝子の**CDS（翻訳領域）の塩基配列**が入っています。
* **`odb12v1_dna_fasta.tgz`**: 全ての生物の、全ての遺伝子の**ゲノムDNA配列**が入っています。

---

### その他のファイル

* **`README.txt`**: これらのファイルに関する説明が書かれたメインのファイルです。困ったらまずここを読むのが良いですね！

この情報がお役に立てたら嬉しいです！解析、がんばってくださいね！応援しています💖

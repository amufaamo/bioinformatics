
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

#### `odb12v1_levels.tab.gz`
オルソロググループ（OG）が計算された分類階層（クレード）に関するファイルです。例えば、「哺乳類」や「真核生物」といった大きなくくりですね。
1.  `level NCBI tax id`: その分類階層のNCBIタクソノミーIDです。
2.  `scientific name`: その分類階層の学名です。
3.  `genes count`: その階層に属する全ての生物の遺伝子数（重複なし）です。
4.  `OGs count`: その階層で構築されたオルソロググループ（OG）の総数です。
5.  `species count`: その階層に属する生物種の数（重複なし）です。

#### `odb12v1_level2species.tab.gz`
最上位の分類階層と、そこに属する各生物種との関係を示しています。
1.  `top-most level NCBI tax id`: 最上位の分類階層のIDです。{2(細菌), 2157(古細菌), 2759(真核生物), 10239(ウイルス)}のいずれかです。
2.  `OrthoDB organism id`: OrthoDBの生物IDです。
3.  `hops`: 最上位からその生物の階層までの、分類階層をいくつ経由したかを示す数です。
4.  `intermediate levels`: 最上位からその生物までの、中間の分類階層のリストです。

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

#### `odb12v1_gene_xrefs.tab.gz`
OrthoDBの遺伝子IDと、他のデータベースのIDとを関連付ける「相互参照」ファイルです。
1.  `OrthoDB gene id`: OrthoDBの遺伝子IDです。
2.  `external gene identifier`: 外部データベースでの遺伝子IDです。
3.  `external DB name`: 外部データベースの名前です。{GOterm, InterPro, NCBIproteinGI, UniProt, ENSEMBL, NCBIgid, NCBIgenename}のいずれかが入ります。

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

---
title: 遺伝子発現差解析のツールedgeRの使い方
emoji: "🧬"
type: tech
topics: ["R", "edgeR", "バイオインフォマティクス", "RNA-seq"]
published: true
---
# はじめに

RNA-Seqデータ解析、特に発現変動遺伝子（Differentially Expressed Genes, DEGs）の探索は、生命科学研究において非常に重要なステップです。数あるツールの中でも、Rの`edgeR`パッケージは、その信頼性と柔軟性から広く利用されています。

この記事では、`edgeR`パッケージの`exactTest`を用いて、RNA-Seqのカウントデータから発現変動遺伝子を検出する基本的な流れを、サンプルデータを使いながらステップバイステップで解説します。RやRNA-Seq解析が初めての方でも理解しやすいように、各処理の意味や出力例も交えて説明していきます。

# 準備するもの

本格的な解析を始める前に、以下の準備が必要です。

1.  **R と RStudio のインストール**
    * Rは統計解析のためのプログラミング言語です。 [The Comprehensive R Archive Network (CRAN)](https://cran.r-project.org/) からダウンロードしてインストールしてください。
    * RStudioはRをより使いやすくするための統合開発環境（IDE）です。[RStudio Desktop](https://posit.co/download/rstudio-desktop/) からダウンロードしてインストールすることをおすすめします。

2.  **edgeR パッケージのインストール**
    * RまたはRStudioのコンソールで以下のコマンドを実行し、`edgeR`パッケージをインストールします。

    ```r
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("edgeR")
    ```

    * もし`dplyr`パッケージも使いたい場合は、同様にインストールしておきましょう（今回のメイン関数では直接使用しませんが、データ操作に便利です）。

    ```r
    install.packages("dplyr")
    ```

# サンプルデータの準備

まずは解析に使うためのサンプルデータを作成しましょう。実際のRNA-Seqデータはリードカウントと呼ばれる数値の羅列で、遺伝子ごと、サンプルごとの発現量を示します。

ここでは、説明のために小さなダミーデータセットをRで生成します。
コントロール群3サンプル、トリートメント群3サンプルの計6サンプル、遺伝子数は100個とします。

```r
# 再現性のためのシード設定
set.seed(123)

# counts_df: 行が遺伝子、列がサンプル のカウントデータフレームを作成
counts_example_df <- data.frame(
  Sample1_Ctrl = rnbinom(100, mu=50, size=10), # Control群サンプル1
  Sample2_Ctrl = rnbinom(100, mu=60, size=10), # Control群サンプル2
  Sample3_Ctrl = rnbinom(100, mu=55, size=10), # Control群サンプル3
  Sample4_Treat = rnbinom(100, mu=150, size=10),# Treatment群サンプル1 (発現量多め)
  Sample5_Treat = rnbinom(100, mu=160, size=10),# Treatment群サンプル2 (発現量多め)
  Sample6_Treat = rnbinom(100, mu=155, size=10) # Treatment群サンプル3 (発現量多め)
)
rownames(counts_example_df) <- paste0("Gene", 1:100) # 遺伝子名を付与

# 作成したデータの一部を表示して確認
print(head(counts_example_df))
```

出力例 (一部):

```text
        Sample1_Ctrl Sample2_Ctrl Sample3_Ctrl Sample4_Treat Sample5_Treat Sample6_Treat
Gene1             43           51           61           163           171           151
Gene2             50           68           46           141           151           154
Gene3             50           63           51           141           150           141
Gene4             63           56           50           164           170           156
Gene5             44           54           53           147           156           155
Gene6             54           56           52           158           149           164
```

この`counts_example_df`が、解析の入力となるカウントデータです。行が遺伝子、列がサンプルに対応しています。
次に、各サンプルがどのグループに属しているかを示す情報と、どのグループ間比較をしたいかを定義します。

```R
# group_vector: 各サンプルが属するグループ ("Control" または "Treatment")
group_info_vector <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
```

```R
# comparison_pair: 比較するグループのペア (Treatment vs Control の変化を見たい)
# 最初の要素が対照群 (分母)、2番目の要素が処理群 (分子) になります。
# つまり、logFC は "Treatment" / "Control" の対数を意味します。
comparison_groups <- c("Control", "Treatment")
```

これでデータの準備は完了です！

# edgeRを使った解析関数

今回は、以下の`run_edgeR_exactTest`関数を使って解析を進めます。この関数は、カウントデータ、グループ情報、比較ペアを引数に取り、差次発現遺伝子のリストを返します。

```R
# edgeR を用いた差次発現遺伝子解析 (exactTest)
run_edgeR_exactTest <- function(counts_df, group_vector, comparison_pair) {

  # 引数のチェック (簡易版)
  if (!is.data.frame(counts_df)) {
    stop("counts_df はデータフレームである必要があります。")
  }
  if (ncol(counts_df) != length(group_vector)) {
    stop("counts_df の列数と group_vector の長さが一致しません。")
  }
  if (length(comparison_pair) != 2) {
    stop("comparison_pair は2つのグループ名を含むベクトルである必要があります。")
  }
  if (!all(comparison_pair %in% unique(group_vector))) {
    stop("comparison_pair で指定されたグループ名が group_vector に存在しません。")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("パッケージ 'edgeR' がインストールされていません。インストールしてください。")
  }

  # 1. カウントデータの行列への変換
  count_data_matrix <- as.matrix(counts_df)

  # 2. グループ情報の定義
  group <- factor(group_vector, levels = comparison_pair)

  # 3. DGEListオブジェクトの作成
  d <- edgeR::DGEList(counts = count_data_matrix, group = group)
  cat("DGEListオブジェクト作成完了。\n")
  cat("サンプル数:", ncol(d), "初期遺伝子数:", nrow(d), "\n")

  # 4. 低発現遺伝子のフィルタリング
  keep <- edgeR::filterByExpr(d, group = group)
  d <- d[keep, , keep.lib.sizes = FALSE]
  cat("フィルタリング後の遺伝子数:", nrow(d), "\n")

  if (nrow(d) == 0) {
    stop("フィルタリングの結果、解析対象の遺伝子が0になりました。フィルタリング基準を確認してください。")
  }

  # 5. 正規化 (TMM)
  d <- edgeR::calcNormFactors(d)
  cat("TMM正規化完了。\n")

  # 6. 分散の推定
  d <- edgeR::estimateCommonDisp(d)
  cat("コモン分散推定完了。\n")
  d <- edgeR::estimateTrendedDisp(d)
  cat("トレンド分散推定完了。\n")
  d <- edgeR::estimateTagwiseDisp(d)
  cat("タグワイズ分散推定完了。\n")

  # 7. 差次発現検定 (exactTest)
  result_test <- edgeR::exactTest(d) # pair は DGEList の group レベルに基づき自動設定
  cat("exactTest実行完了。\n")

  # 8. 結果の取得とデータフレームへの変換
  result_table <- edgeR::topTags(result_test, n = nrow(d$counts))
  result_df <- as.data.frame(result_table)
  cat("結果テーブル作成完了。\n")

  return(result_df)
}
```

# 解析ステップの解説

それでは、`run_edgeR_exactTest`関数の中で行われている各ステップを詳しく見ていきましょう。

## カウントデータの行列への変換

```R
count_data_matrix <- as.matrix(counts_df)
```

edgeRの多くの関数は、入力として数値行列（matrix）を期待します。そのため、まずデータフレーム形式のカウントデータを数値行列に変換します。

## グループ情報の定義

```R
group <- factor(group_vector, levels = comparison_pair)
```

サンプルがどの実験群に属するかを示す`group_vector`を因子（factor）型に変換します。ここで重要なのは`levels`引数です。`comparison_pair`（例: `c("Control", "Treatment")`）を`levels`に指定することで、後の比較検定で「Treatment vs Control」という意図した比較が正しく行われるように基準群を設定します。

## DGEListオブジェクトの作成

```R
d <- edgeR::DGEList(counts = count_data_matrix, group = group)
```

DGEListオブジェクトは、edgeRで解析を行うための基本的なデータ構造です。カウントデータ、グループ情報、ライブラリサイズ（各サンプルの総カウント数、自動計算）などが格納されます。

DGEListオブジェクトの中身は以下のとおりです。
```R
names(d)
```
```text
'counts''samples'
```
d$countsの中身は以下のとおりです。カウントデータが入ってます。

```R
head(d$counts)
```
```text
Sample1_Ctrl Sample2_Ctrl Sample3_Ctrl Sample4_Treat	Sample5_Treat	Sample6_Treat
Gene1	46	28	40	139	104	159
Gene2	25	34	31	133	136	181
Gene3	81	89	54	190	195	97
Gene4	20	84	62	134	181	186
Gene5	71	61	51	180	129	162
Gene6	54	51	45	68	169	106
```

d$samplesの中身は以下のとおりです。
サンプルの名前、グループの他に、総カウント数と、normalization factorが格納されています。normalization factorは計算する前はすべて１になっています。
```R
d$samples
```

```text
	group	lib.size	norm.factors
<fct>	<dbl>	<dbl>
Sample1_Ctrl	Control	5064	1
Sample2_Ctrl	Control	5881	1
Sample3_Ctrl	Control	5388	1
Sample4_Treat	Treatment	15389	1
Sample5_Treat	Treatment	14880	1
Sample6_Treat	Treatment	15464	1
```

## 低発現遺伝子のフィルタリング

```R
keep <- edgeR::filterByExpr(d, group = group)
d <- d[keep, , keep.lib.sizes = FALSE]
```

発現量が極端に低い遺伝子は、統計的な検出力が低く、解析のノイズになることがあります。`filterByExpr`関数は、実験デザイン（グループ情報）を考慮して、発現が十分でない遺伝子を自動的に判断し、除外するための論理ベクトル（`keep`）を返します。
`d[keep, , keep.lib.sizes = FALSE]`で、実際にDGEListオブジェクトからこれらの遺伝子を除外します。`keep.lib.sizes = FALSE`は、ライブラリサイズを再計算しないようにするオプションです（通常はフィルタリング後に正規化を行うため）。

コンソール出力例:

```text
フィルタリング後の遺伝子数: 95
```

※この数値はデータの特性によって変わります。ここでは例として5遺伝子が除去されたとします。
もしフィルタリングで全ての遺伝子が除去されてしまう場合は、データが小さすぎるか、`filterByExpr`のデフォルト閾値が厳しすぎる可能性があります。

keepの中身は以下です。
```R
head(keep)
```
```text
Gene1: TRUE Gene2: TRUE Gene3: TRUE Gene4: TRUE Gene5: TRUE Gene6: TRUE
```

## 発現変動遺伝子の計算
分散の計算などで、３つの方法に分かれます。
### [推奨]準尤度パイプライン（edgeR v3, v4）
もっとも新しく、現在推奨されている計算方法です。

```R
design <- model.matrix(~group)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
```
### 負の二項GLMフレームワーク(edgeR v2)

```R
design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
ltr <- glmFLRT(fit, coef = 2)
topTags(ltr)
```
### クラシックパイプライン(edgeR v1)
```R
y <- estimateDisp(y, design)
et <- exactText(y)
topTags(ex)
```

## 正規化 (TMM)

```R
d <- edgeR::calcNormFactors(d)
cat("TMM正規化完了。\n")
```

RNA-Seqデータでは、サンプル間の総リード数（ライブラリサイズ）の違いや、一部の非常に高発現な遺伝子の影響で、単純な比較が難しい場合があります。正規化はこれらのバイアスを補正し、サンプル間の発現量をより公平に比較できるようにする処理です。
`calcNormFactors`関数は、TMM (Trimmed Mean of M-values) という方法で正規化係数を計算します。

コンソール出力例:

```text
TMM正規化完了。
```
これにより、norm.factorが計算されます。計算される前は、norm.factorはすべて1です。

```R
d$samples
```
```R
	group	lib.size	norm.factors
<fct>	<dbl>	<dbl>
Sample1_Ctrl	Control	5064	1.0107421
Sample2_Ctrl	Control	5881	0.9833187
Sample3_Ctrl	Control	5388	0.9996209
Sample4_Treat	Treatment	15389	0.9846202
Sample5_Treat	Treatment	14880	1.0164202
Sample6_Treat	Treatment	15464	1.0057452
```

## 分散の推定

```R
d <- edgeR::estimateCommonDisp(d)
cat("コモン分散推定完了。\n")
d <- edgeR::estimateTrendedDisp(d)
cat("トレンド分散推定完了。\n")
d <- edgeR::estimateTagwiseDisp(d)
cat("タグワイズ分散推定完了。\n")
```

edgeR (特に`exactTest`) は、負の二項分布モデルに基づいて統計検定を行います。このモデルでは「分散」の推定が非常に重要です。
* `estimateCommonDisp`: 全遺伝子で共通の分散（コモン分散）を推定します。
* `estimateTrendedDisp`: 遺伝子の平均発現量と分散の関係（トレンド分散）を推定します。
* `estimateTagwiseDisp`: 個々の遺伝子ごとの分散（タグワイズ分散）を、コモン分散とトレンド分散の情報を利用しつつ、経験的ベイズ法を用いて安定的に推定します。

これらのステップを経ることで、特に反復サンプル数が少ない場合でも、より信頼性の高い分散推定値を得ることができます。

コンソール出力例:

```text
コモン分散推定完了。
トレンド分散推定完了。
タグワイズ分散推定完了。
```

edgeRにおける負の二項分布についてはこちらをご覧ください。

https://qiita.com/amufaamo/items/a286bc132ab63746e39d

分散推定について、詳しい記事は以下を御覧ください。
[edgeRの分散推定に関する詳しい説明](https://qiita.com/amufaamo/items/425d08112a5ff2ee3992)

また、edgeRの分散推定後のDGEListオブジェクトの中身については以下の記事を御覧ください。
[edgeRの分散推定後のDGEListオブジェクトの中身](https://qiita.com/amufaamo/items/b52fa2dc715dda9b538b)


## 差次発現検定 (exactTest)

```R
result_test <- edgeR::exactTest(d)
cat("exactTest実行完了。\n")
```

いよいよ発現変動遺伝子の検定です。`exactTest`関数は、Fisherの正確検定（に類似した手法）を使い、指定された2群間で各遺伝子の発現量に統計的に有意な差があるかを評価します。
DGEList作成時に`group`ファクターの`levels`を`comparison_pair`（例: `c("Control", "Treatment")`）で指定したため、`exactTest`は自動的に `levels`の2番目のグループ vs 1番目のグループ（Treatment vs Control）の比較を行います。

コンソール出力例:

```text
exactTest実行完了。
```

result_testの中身は以下です。
```R
names(result_test)
```
```
'table''comparison''genes'
```
```R
head(result_test$table)
```
```R
	logFC	logCPM	PValue
<dbl>	<dbl>	<dbl>
Gene1	0.3092994	12.99407	0.43313480
Gene2	0.8292998	12.96356	0.02615609
Gene3	-0.3809196	13.58259	0.37585946
Gene4	0.1256244	13.38271	0.76724013
Gene5	-0.1351120	13.41798	0.69883189
Gene6	-0.3015343	13.05316	0.46444746
```
```R
result_test$comparison
```
```R
'Control''Treatment'
```
```R
result_test$genes
```
```R
NULL
```


## 結果の取得とデータフレームへの変換

```R
result_table <- edgeR::topTags(result_test, n = nrow(d$counts))
result_df <- as.data.frame(result_table)
cat("結果テーブル作成完了。\n")
```

`exactTest`の実行結果は`DGEExact`オブジェクトとして得られます。ここから、発現変動の大きさ（logFC）、p値、FDR（False Discovery Rate、多重検定補正後のp値）などを含む結果テーブルを抽出するのが`topTags`関数です。
`n = nrow(d$counts)`とすることで、フィルタリング後の全遺伝子について結果を取得します。
最後に、扱いやすいようにデータフレーム形式に変換します。

コンソール出力例:

```text
結果テーブル作成完了。
```

`result_df`の中身は以下です。
```R
head(result_df)
```
```R
logFC	logCPM	PValue	FDR
<dbl>	<dbl>	<dbl>	<dbl>
Gene87	-0.9992087	13.35877	0.005826470	0.3130362
Gene46	-0.8812362	13.22718	0.006900652	0.3130362
Gene89	-1.1238634	13.20669	0.009391085	0.3130362
Gene12	0.9247378	13.15863	0.022065282	0.5231217
Gene2	0.8292998	12.96356	0.026156087	0.5231217
Gene29	-0.7561143	13.30031	0.039860245	0.6509024
```

この結果の意味を説明します。


1.  **`logFC` (log2 Fold Change)**
    * **意味・定義**:
        比較している2つの実験群間での遺伝子発現量の変化率を、2を底とする対数（log2スケール）で示した値です。例えば、実験群A（例：コントロール群）と実験群B（例：薬剤処理群）を比較し、B群の変化をA群に対して見ている場合、`logFC` は「B群での発現量 / A群での発現量」のlog2値に相当します。
    * **解釈**:
        * **正の値**: 2番目に指定した群（例：薬剤処理群）で遺伝子の発現量が1番目に指定した群（例：コントロール群）に比べて高い（アップレギュレートされている）ことを意味します。
            * `logFC = 1` は、発現量が約2倍になったことを示します。
            * `logFC = 2` は、発現量が約4倍になったことを示します。
        * **負の値**: 2番目に指定した群で遺伝子の発現量が1番目に指定した群に比べて低い（ダウンレギュレートされている）ことを意味します。
            * `logFC = -1` は、発現量が約1/2倍（半分）になったことを示します。
            * `logFC = -2` は、発現量が約1/4倍になったことを示します。
        * **0に近い値**: 2群間で発現量にほとんど差がないことを示します。
    * **例のデータ (`Gene87`)**: `logFC = -0.9992087` は、この遺伝子の発現量が比較対象の群に対して約半分 ( $2^{-0.999} \approx 0.5$ ) に減少していることを示唆します。

---

2.  **`logCPM` (log2 Counts Per Million)**
    * **意味・定義**:
        各遺伝子の平均的な発現量を、log2スケールで示したものです。CPM (Counts Per Million) は、ライブラリサイズ（サンプルごとの総リード数）の違いを補正するために、各遺伝子のカウント数をそのサンプルの総カウント数で割り、100万を掛けた値です。`logCPM` は、このCPM値をさらにlog2変換し、通常は全サンプルにわたる平均値や、モデルに基づいた平均的な値として計算されます。
    * **解釈**:
        * 値が大きいほど、その遺伝子の平均的な発現量が高いことを示します。
        * 値が小さいほど、その遺伝子の平均的な発現量が低いことを示します。
        この値は、遺伝子の絶対的な発現レベルの目安となります。
    * **例のデータ (`Gene87`)**: `logCPM = 13.35877` は、この遺伝子がデータセット全体で比較的高めに発現していることを示唆しています（具体的な発現量はデータのスケールに依存します）。

---

3.  **`PValue` (P-value)**
    * **意味・定義**:
        観測された遺伝子の発現量の差（またはそれ以上に極端な差）が、実際には2群間で発現量に差がないという帰無仮説のもとで、偶然生じる確率を示します。`edgeR`では、`exactTest`の場合はフィッシャーの正確確率検定に類似した手法、`glmLRT`や`glmQLFTest`の場合は尤度比検定やF検定などを用いて計算されます。
    * **解釈**:
        * P値が小さいほど、観測された発現量の差が偶然であるとは考えにくく、統計的に有意な差があると判断されます。
        * 一般的に、P値が0.05や0.01を下回る場合に「統計的に有意」とされますが、これはあくまで目安であり、多重検定の問題を考慮する必要があります。
    * **例のデータ (`Gene87`)**: `PValue = 0.005826470` は、この遺伝子の発現変動が偶然生じた確率が約0.58%であることを示しており、このP値だけを見れば比較的小さいと言えます。

---

4.  **`FDR` (False Discovery Rate)**
    * **意味・定義**:
        偽発見率。RNA-Seqのように多数の遺伝子（数千～数万）に対して同時に統計検定を行うと、実際には差がないのに偶然P値が小さくなってしまう「偽陽性」の数が増加します（多重検定の問題）。FDRは、この多重検定の問題を補正したP値の一種で、あるFDRの閾値（例：0.05）で有意と判定された遺伝子群のうち、実際に偽陽性である遺伝子の割合の期待値をコントロールするための指標です。`edgeR`では、Benjamini-Hochberg (BH) 法などがデフォルトで用いられてFDRが計算されます。
    * **解釈**:
        * 例えば、`FDR < 0.05` を有意の基準とした場合、その基準を満たして選ばれた遺伝子リストのうち、平均して5%程度が偽陽性（実際には変動していないのに変動していると判定されたもの）であると許容することを意味します。
        * P値と同様に、FDRの値が小さいほど、その遺伝子の発現変動の統計的な信頼性が高いと考えられます。通常、発現変動遺伝子を選定する際には、生のP値ではなく、このFDRを用います。
    * **例のデータ (`Gene87`)**: `FDR = 0.3130362` は、この遺伝子の発現変動を「有意」と判定した場合、その判定が誤りである（偽陽性である）確率が約31.3%と推定されることを意味します。一般的に`FDR < 0.05` や `FDR < 0.1` を有意の閾値とすることが多いため、この遺伝子はこのFDR値だけを見ると「統計的に有意な変動とは言えない」と判断される可能性が高いです。

---

**まとめると、`result_df` の各行は一つの遺伝子に対応し、以下の情報を提供します。**

* その遺伝子が2群間でどれだけ発現量が変わったか (`logFC`)
* その遺伝子がデータ全体で平均的にどれくらい発現しているか (`logCPM`)
* 観測された発現量の差が偶然である確率 (`PValue`)
* 多数の遺伝子を検定したことを考慮した上で、その発現量の差が統計的に信頼できるか (`FDR`)

研究者はこれらの情報を総合的に評価し、生物学的に意味のある発現変動遺伝子を絞り込んでいきます。特に `logFC` の大きさと `FDR` の小ささが重要な判断基準となります。

# 実際に解析を実行してみよう

それでは、準備したサンプルデータと`run_edgeR_exactTest`関数を使って、実際に解析を実行してみましょう。

```R
# 必要なパッケージがロードされているか確認 (edgeRは関数内でチェックされます)
if (requireNamespace("edgeR", quietly = TRUE)) {
  cat("edgeR パッケージがロードされています。\n")

  # tryCatchでエラーが発生した場合の処理を記述
  tryCatch({
    deg_results <- run_edgeR_exactTest(
      counts_df = counts_example_df,
      group_vector = group_info_vector,
      comparison_pair = comparison_groups
    )

    # 結果の最初の数行を表示
    cat("\n--- 解析結果 (上位数件) ---\n")
    print(head(deg_results))

    # logFC の向きを確認 (Treatment群で発現が高い遺伝子はlogFCが正になるはず)
    # 例として、Gene1 の平均カウントを比較
    cat("\n--- logFCの向き確認 (Gene1) ---\n")
    gene1_control_mean <- mean(as.numeric(counts_example_df["Gene1", group_info_vector == "Control"]))
    gene1_treatment_mean <- mean(as.numeric(counts_example_df["Gene1", group_info_vector == "Treatment"]))
    cat("Gene1 の Control 群平均カウント:", gene1_control_mean, "\n")
    cat("Gene1 の Treatment 群平均カウント:", gene1_treatment_mean, "\n")
    cat("Gene1 の deg_results での logFC:", deg_results["Gene1", "logFC"], "\n")

  }, error = function(e) {
    cat("エラーが発生しました:", e$message, "\n")
  })

} else {
  cat("edgeR パッケージがインストールされていません。サンプルコードを実行できません。\n")
}
```

実行時のコンソール出力例 (関数内の`cat`文も含む):

```text
edgeR パッケージがロードされています。
DGEListオブジェクト作成完了。
サンプル数: 6 初期遺伝子数: 100 
フィルタリング後の遺伝子数: (実際のフィルタリング結果による数値)
TMM正規化完了。
コモン分散推定完了。
トレンド分散推定完了。
タグワイズ分散推定完了。
exactTest実行完了。
結果テーブル作成完了。

--- 解析結果 (上位数件) ---
            logFC   logCPM       PValue          FDR
GeneX  X.XXXXX XX.XXXXX X.XXXXXXXXe-XX X.XXXXXXXXe-XX
GeneY  X.XXXXX XX.XXXXX X.XXXXXXXXe-XX X.XXXXXXXXe-XX
GeneZ -X.XXXXX XX.XXXXX X.XXXXXXXXe-XX X.XXXXXXXXe-XX
... (以下略) ...

--- logFCの向き確認 (Gene1) ---
Gene1 の Control 群平均カウント: (例: 51.66) 
Gene1 の Treatment 群平均カウント: (例: 161.66) 
Gene1 の deg_results での logFC: (例: 1.64) 
```

※ GeneX, X.XXXXX の部分は実際の遺伝子名と数値に置き換わります。`set.seed(123)` を使っているので、同じ環境であれば同じような傾向の結果が得られるはずです。


# おわりに

この記事では、`edgeR`パッケージの`exactTest`を用いてRNA-Seqデータから発現変動遺伝子を解析する基本的な流れを解説しました。
`edgeR`には、より複雑な実験デザインに対応できるGLM (Generalized Linear Model) アプローチ（`glmFit`, `glmLRT`など）も用意されています。また、得られた結果を可視化するMAプロットやボルケーノプロットを作成することで、結果の全体像を把握しやすくなります。
RNA-Seqデータ解析は奥が深いですが、`edgeR`のような強力なツールを使いこなすことで、生命現象の理解に繋がる多くの知見を得ることができます。ぜひ、ご自身のデータで試してみてください！

# 付録：今回使用したRコード全体

```R
# --- ライブラリのロード ---
# edgeRは関数内で存在チェックとロード推奨メッセージが出ます
# BiocManager::install("edgeR") # 未インストールの場合
# install.packages("dplyr") # 未インストールの場合

# --- サンプルデータの準備 ---
set.seed(123) # 再現性のためのシード設定

# counts_df: 行が遺伝子、列がサンプル
counts_example_df <- data.frame(
  Sample1_Ctrl = rnbinom(100, mu=50, size=10),
  Sample2_Ctrl = rnbinom(100, mu=60, size=10),
  Sample3_Ctrl = rnbinom(100, mu=55, size=10),
  Sample4_Treat = rnbinom(100, mu=150, size=10),
  Sample5_Treat = rnbinom(100, mu=160, size=10),
  Sample6_Treat = rnbinom(100, mu=155, size=10)
)
rownames(counts_example_df) <- paste0("Gene", 1:100)

# group_vector: 各サンプルが属するグループ
group_info_vector <- c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")

# comparison_pair: 比較するグループのペア (Treatment vs Control を見たい場合)
comparison_groups <- c("Control", "Treatment") # 最初の要素が対照

# --- 解析関数の定義 ---
run_edgeR_exactTest <- function(counts_df, group_vector, comparison_pair) {

  # 引数のチェック (簡易版)
  if (!is.data.frame(counts_df)) {
    stop("counts_df はデータフレームである必要があります。")
  }
  if (ncol(counts_df) != length(group_vector)) {
    stop("counts_df の列数と group_vector の長さが一致しません。")
  }
  if (length(comparison_pair) != 2) {
    stop("comparison_pair は2つのグループ名を含むベクトルである必要があります。")
  }
  if (!all(comparison_pair %in% unique(group_vector))) {
    stop("comparison_pair で指定されたグループ名が group_vector に存在しません。")
  }
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("パッケージ 'edgeR' がインストールされていません。インストールしてください。")
  }

  # 1. カウントデータの行列への変換
  count_data_matrix <- as.matrix(counts_df)

  # 2. グループ情報の定義
  group <- factor(group_vector, levels = comparison_pair)

  # 3. DGEListオブジェクトの作成
  d <- edgeR::DGEList(counts = count_data_matrix, group = group)
  cat("DGEListオブジェクト作成完了。\n")
  cat("サンプル数:", ncol(d), "初期遺伝子数:", nrow(d), "\n")

  # 4. 低発現遺伝子のフィルタリング
  keep <- edgeR::filterByExpr(d, group = group)
  d <- d[keep, , keep.lib.sizes = FALSE]
  cat("フィルタリング後の遺伝子数:", nrow(d), "\n")

  if (nrow(d) == 0) {
    stop("フィルタリングの結果、解析対象の遺伝子が0になりました。フィルタリング基準を確認してください。")
  }

  # 5. 正規化 (TMM)
  d <- edgeR::calcNormFactors(d)
  cat("TMM正規化完了。\n")

  # 6. 分散の推定
  d <- edgeR::estimateCommonDisp(d)
  cat("コモン分散推定完了。\n")
  d <- edgeR::estimateTrendedDisp(d)
  cat("トレンド分散推定完了。\n")
  d <- edgeR::estimateTagwiseDisp(d)
  cat("タグワイズ分散推定完了。\n")

  # 7. 差次発現検定 (exactTest)
  result_test <- edgeR::exactTest(d)
  cat("exactTest実行完了。\n")

  # 8. 結果の取得とデータフレームへの変換
  result_table <- edgeR::topTags(result_test, n = nrow(d$counts))
  result_df <- as.data.frame(result_table)
  cat("結果テーブル作成完了。\n")

  return(result_df)
}

# --- 関数の実行と結果表示 ---
if (requireNamespace("edgeR", quietly = TRUE)) {
  cat("edgeR パッケージがロードされています。\n")
  tryCatch({
    deg_results <- run_edgeR_exactTest(
      counts_df = counts_example_df,
      group_vector = group_info_vector,
      comparison_pair = comparison_groups
    )
    cat("\n--- 解析結果 (上位数件) ---\n")
    print(head(deg_results))

    cat("\n--- logFCの向き確認 (Gene1) ---\n")
    gene1_control_mean <- mean(as.numeric(counts_example_df["Gene1", group_info_vector == "Control"]))
    gene1_treatment_mean <- mean(as.numeric(counts_example_df["Gene1", group_info_vector == "Treatment"]))
    cat("Gene1 の Control 群平均カウント:", gene1_control_mean, "\n")
    cat("Gene1 の Treatment 群平均カウント:", gene1_treatment_mean, "\n")
    # Gene1がフィルタリングで除去されていないか確認してからアクセス
    if ("Gene1" %in% rownames(deg_results)) {
        cat("Gene1 の deg_results での logFC:", deg_results["Gene1", "logFC"], "\n")
    } else {
        cat("Gene1 はフィルタリングで除去されたか、結果テーブルに存在しません。\n")
    }

  }, error = function(e) {
    cat("エラーが発生しました:", e$message, "\n")
  })
} else {
  cat("edgeR または dplyr パッケージがインストールされていません。サンプルコードを実行できません。\n")
}

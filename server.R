library(shiny)
library(shinyAce)
library(sciplot)
library(ggplot2)



shinyServer(function(input, output) {


    
    options(warn=-1)
    
    
    anovakun <- reactive({
        
        dat <- read.csv(text=input$text, sep="\t")
    
    #--------------------------------------------------------------
    # ここからANOVA君の関数貼付け
    # anovakun_451.txt（2014-0303）
    
        # 【ANOVA君：要因計画のタイプと水準数を入力することにより，分散分析を行う関数】
        # 1)フリーの統計ソフトウェア「R」で動作する関数
        # 2)被験者間要因（独立測度），被験者内要因（反復測度）のいずれか，または，両方を含む各タイプの分散分析を扱う
        # 3)引数としては，最初にデータフレーム名，次に計画のタイプ（""で囲むこと）を入力し，その後，各要因の水準数を順に入力する
        # （作成：井関龍太）
        #
        # 【使用法】
        # anovakunに読み込むためのデータフレームは，以下のような形式で作っておく
        # 1)被験者間要因はタテ，被験者内要因はヨコに並べる
        # 2)被験者内要因を表す行はデータとして読み込まない（下の例では，点線より上の行は読み込まず，点線より下のデータのみ読み込む）
        # 3)被験者間要因を表す列はデータとして読み込む（下の例では，a1などの文字を含む列；被験者間要因数が増えるたびにラベル用の列を増やす）
        # 4)被験者間要因を表すラベルは，水準ごとに別の文字または数字を用いる（同じラベル＝同じ水準と見なされる）
        # 5)被験者内，被験者間要因ともに，前の方の要因から順に各水準のデータを入れ子状に整理して並べること（例を参照）
        # 6)下の被験者１〜６の列は例の説明のためにつけたものなので，実際のデータフレームには必要ない
        #
        # ［AsBC計画の例］（被験者間要因でデータ数が同じ）
        # 		b1	b1	b2	b2
        # 		c1	c2	c1	c2
        # 	------------------------------------
        # 	a1	12	9	14	13	---被験者１
        # 	a1	13	10	14	12	---被験者２
        # 	a1	11	10	13	15	---被験者３
        # 	a2	18	12	16	15	---被験者４
        # 	a2	17	14	15	14	---被験者５
        # 	a2	15	13	18	15	---被験者６
        #
        # データフレームを代入した変数名をxとすると，
        #
        # > anovakun(x, "AsBC", 2, 2, 2)
        #
        # のようにして関数を呼び出す
        #
        # ［ABsC計画の例］（被験者間要因でデータ数がふぞろい）
        # 			c1	c2
        # 	-------------------------------
        # 	a1	b1	3.5	4.2	---被験者１
        # 	a1	b1	2.7	3.2	---被験者２
        # 	a1	b2	2.5	3.8	---被験者３
        # 	a1	b3	4.0	3.9	---被験者４
        # 	a2	b1	3.3	4.0	---被験者５
        # 	a2	b1	1.4	2.5	---被験者６
        # 	a2	b2	3.7	4.2	---被験者７
        # 	a2	b2	2.2	4.2	---被験者８
        # 	a2	b3	1.3	2.1	---被験者９
        # 	a2	b3	3.4	3.9	---被験者10
        #
        # データフレームを代入した変数名をxとすると，
        #
        # > anovakun(x, "ABsC", 2, 3, 2)
        #
        # のようにして関数を呼び出す
        #
        # 【オプション】
        # 1)long……long = Tとすると，ロング形式のデータを読み込んで処理する
        # 2)type2……type2 = Tとすると，平方和の計算方法をタイプⅡに切り替える
        # 3)nopost……nopost = Tとすると，下位検定を実行しない
        # 4)tech……テクニカルアウトプット；tech = Tとすると，データフレームをリストでつないだ形式で結果を出力する
        # 5)data.frame……data.frame = Tとすると，計算に使用したデータフレームを出力する（関数中でdatと表現されているデータフレーム）
        # 6)copy……copy = Tとすると，出力結果をクリップボードにもコピーする（クリップボードの内容は上書きされる）
        # 7)holm……holm = Tとすると，多重比較の方法がHolmの方法になる
        # 8)hc……hc = Tとすると，多重比較の方法がHolland-Copenhaverの方法になる
        # 9)s2r……s2r = Tとすると，Shafferの方法のための仮説数の計算法を具体的な棄却のパターンに基づく方法に変更する（Rasmussenのアルゴリズムに基づく）
        # 10)s2d……s2d = Tとすると，Shafferの方法のための仮説数の計算法を具体的な棄却のパターンに基づく方法に変更する（Donoghueのアルゴリズムに基づく）
        # 11)fs1……fs1 = Tとすると，ステップ１の基準をステップ２の基準に置き換えた方法でShafferの方法を行う
        # 12)fs2r……fs2r = Tとすると，Shaffer2の方法とF-Shaffer1の方法を組み合わせた方法でShafferの方法を行う（Rasmussenのアルゴリズムに基づく）
        # 13)fs2d……fs2d = Tとすると，Shaffer2の方法とF-Shaffer1の方法を組み合わせた方法でShafferの方法を行う（Donoghueのアルゴリズムに基づく）
        # 14)hc, s2r……hc = Tかつs2r = T（または，fs2r = T，s2d = T，fs2d = T）とすると，Shaffer2の基準を用いてHolland-Copenhaverの方法を行う
        # 15)holm, hc……holm = Tかつhc = Tとすると，Holmの調整基準にSidakの不等式を用いた多重比較（Holm-Sidak法）を行う
        # 16)criteria……criteria = Tとすると，多重比較の出力において，調整済みｐ値の変わりに調整済みの有意水準を表示する
        # 17)lb……lb = Tとすると，すべての被験者内効果に対してイプシロンの下限値（Lower Bound）による自由度の調整を行う（Geisser-Greenhouseの保守的検定）
        # 18)gg……gg = Tとすると，すべての被験者内効果に対してGreehnouse-Geisserのイプシロンによる自由度の調整を行う
        # 19)hf……hf = Tとすると，すべての被験者内効果に対してHuynh-Feldtのイプシロンによる自由度の調整を行う
        # 20)auto……auto = Tとすると，球面性検定が有意であった被験者内効果に対してGreehnouse-Geisserのイプシロンによる自由度の調整を行う
        # 21)mau……mau = Tとすると，球面性検定の方法がMauchlyの球面性検定になる
        # 22)har……har = Tとすると，球面性検定の方法がHarrisの多標本球面性検定になる
        # 23)iga……iga = Tとすると，改良版一般近似検定を行う；イプシロンの代わりに各種の推定値を算出し，分散分析の際に適用する
        # 24)ciga……ciga = Tとすると，修正改良版一般近似検定を行う；イプシロンの代わりに各種の推定値を算出し，分散分析の際に適用する
        # 25)eta……eta = Tとすると，分散分析表にイータ二乗を追加する
        # 26)peta……peta = Tとすると，分散分析表に偏イータ二乗を追加する
        # 27)geta……geta = Tとすると，分散分析表に一般化イータ二乗を追加する；geta = "要因ラベル（A, B, C...）"とすると，指定した要因を測定要因と見なして一般化イータ二乗を計算する；
        # 複数の要因を測定要因として指定するには，例えば，geta = c("A", "C")のようにする（計画に含まれない要因の指定は無効になる）
        # 28)eps……eps = Tとすると，分散分析表にイプシロン二乗を追加する
        # 29)peps……peps = Tとすると，分散分析表に偏イプシロン二乗を追加する
        # 30)geps……geps = Tとすると，分散分析表に一般化イプシロン二乗を追加する；geps = "要因ラベル（A, B, C...）"とすると，指定した要因を測定要因と見なして一般化イプシロン二乗を計算する；
        # 複数の要因を測定要因として指定するには，例えば，geps = c("A", "C")のようにする（計画に含まれない要因の指定は無効になる）
        # 31)omega……omega = Tとすると，分散分析表にオメガ二乗（加算モデル）を追加する
        # 32)omegana……omegana = Tとすると，分散分析表にオメガ二乗（非加算モデル）を追加する
        # 33)pomega……pomega = Tとすると，分散分析表に偏オメガ二乗を追加する
        # 34)gomega……gomega = Tとすると，分散分析表に一般化オメガ二乗（加算モデル）を追加する；gomega = "要因ラベル（A, B, C...）"とすると，指定した要因を測定要因と見なして一般化オメガ二乗を計算する；
        # 複数の要因を測定要因として指定するには，例えば，gomega = c("A", "C")のようにする（計画に含まれない要因の指定は無効になる）
        # 35)gomegana……gomegana = Tとすると，分散分析表に一般化オメガ二乗（非加算モデル）を追加する；gomegana = "要因ラベル（A, B, C...）"とすると，指定した要因を測定要因と見なして一般化オメガ二乗を計算する；
        # 複数の要因を測定要因として指定するには，例えば，gomegana = c("A", "C")のようにする（計画に含まれない要因の指定は無効になる）
        # 36)prep……prep = Tとすると，分散分析表にp_repを追加する
        # 37)cind……cind = Tとすると，記述統計量の表に差分調整型の信頼区間を追加する
        # 38)cin……cin = Tとすると，記述統計量の表に信頼区間を追加する
        #
        # ［オプション使用の例］（テクニカルアウトプットによる出力とHolmの方法による多重比較を指定）
        #
        # > anovakun(x, "AsB", 2, 2, tech = T, holm = T)
        #
        # 【技術情報】
        # 1)anovakunを構成する関数は，仕様上は，最大で26要因までの計画に対応できる；この上限は，要因を表すラベルとしてアルファベット26文字（LETTERSとletters）を使用していることによる
        # 2)ci.calc関数による信頼区間の計算法は，被験者内要因を含むデータについてはCousineau-Moreyの正規化データに基づく方法を用いている
        # 3)ss.calc関数は，デフォルトでは，タイプⅢ平方和の計算法に基づいて分散分析を行う
        # 4)epsilon.calc関数は，デフォルトでは，被験者内要因を含むデータに対してMendozaの多標本球面性検定を行う（近似カイ二乗による）
        # 5)epsilon.calc関数は，オプション指定により，被験者内要因を含むデータに対してMauchlyの球面性検定を行う（近似カイ二乗による）
        # 6)epsilon.calc関数は，オプション指定により，被験者内要因を含むデータに対してHarrisの多標本球面性検定を行う（近似カイ二乗による）
        # 7)epsilon.calc関数によるHuynh-Feldtのイプシロンの計算法は，Lecoutre（1991）の修正に基づく
        # 8)epsilon.calc関数によるIGAとCIGAは，非加重平均を想定した計算法に基づく（Algina, 1997; Algina & Oshima, 1995）
        # 9)anova.modeler関数は，IGAとCIGAを用いた際には，b_hat，c_hatによって調整した後のF値をF値の列に表示する
        # 10)anova.modeler関数による加算モデルに基づくオメガ二乗，一般化オメガ二乗の計算式は，被験者と各要因の間のすべての交互作用が存在しないことを仮定するモデルによる（Dodd & Schultz, 1973を参照）
        # 11)mod.Bon関数は，デフォルトでは，Shafferの方法による多重比較を行う（任意の棄却パターンにおける可能な真の帰無仮説の最大数に基づく方法）
        # 12)mod.Bon関数におけるShafferの方法，Holland-Copenhaverの方法の有意水準の計算は，Rasmussen（1993），Donoghue（2004）のアルゴリズムに基づく（オプションによる）
        # 13)mod.Bon関数による多重比較では，p値の低い対から順に基準値よりも値が低いか否かを判定し，いったん有意でない対が見つかったら，
        # 以降の対についてはすべて判断を保留する（p値が基準値を下回っても*マークを表示しない）；調整済みｐ値の表示の際には，調整済みｐ値が
        # 既定の有意水準（５％）を下回った対にのみ*マークを表示する
        # 14)pro.fraction関数は，単純主効果の検定において誤差項をプールしない（水準別誤差項を使用；サブセットに分散分析を再適用するのと同じ）
        #
        # 【このファイルの含む関数】
        # hmean……調和平均を計算する
        # read.clip……クリップボードの情報を読み込む
        # anovakun……プロセス全体の制御を行う
        # uni.long……データフレームをロング形式に変形する
        # ci.calc……信頼区間を計算する
        # elematch……文字列中のすべての要素を含む文字列をマッチングする
        # ss.calc……平方和を計算する
        # sig.sign……有意水準に合わせて記号を表示する
        # epsilon.calc……Greenhouse-GeisserとHuynh-Feldtのイプシロンを計算する
        # anova.modeler……分散分析表を作る
        # mod.Bon……修正Bonferroniの方法による多重比較を行う
        # post.analyses……下位検定を行う
        # pro.fraction……効果のタイプに適した下位検定を割り当てる
        # info.out……基本情報を出力する
        # bstat.out……記述統計量を出力する
        # anovatab.out……分散分析表を出力する
        # mod.Bon.out……修正Bonferroniの方法による多重比較の結果を出力する
        # simeffects.out……単純主効果の検定の結果を出力する
        # post.out……下位検定の結果を出力する
        # each.out……効果のタイプによって異なる出力方式を割り当てる
        # anovatan……指定した要因についてデータを分割して単純効果の検定を行う
        #
        # 【バージョン情報】
        # 1)anovakun version 1.0.0（R 2.5.1作成；2007/09/03公開）
        # ・分散分析，単純主効果の検定，多重比較
        # 2)anovakun version 2.0.0（R 2.5.1作成；2007/10/01公開）
        # ・球面性検定とイプシロン調整；epsilon.calc関数の追加とそれに伴う変更
        # ・平方和の計算方法を修正；ss.calc関数の修正とelematch関数の追加，それに伴う変更
        # ・多重比較を行う際に，他の要因の効果を考慮した上でのMSeを用いてｔ統計量を計算するように修正；mod.Bon関数と関連する部分の改修
        # 3)anovakun version 2.1.0（R 2.5.1作成；2007/11/01公開）
        # ・平方和の計算方法を変更し，高速化を試みる；ss.calc関数の変更；その他，最適化
        # 4)anovakun version 2.2.0（R 2.5.1，R 2.6.0作成；2007/12/03公開）
        # ・QR分解の方法をLAPACKに変更し，高速化を試みる；その他の修正
        # 5)anovakun version 3.0.0（R 2.5.1，R 2.6.0作成；2008/01/04公開）
        # ・タイプⅢ平方和を計算する機能を追加（タイプⅢをデフォルトに設定）；ss.calc関数の変更とそれに伴う変更
        # ・Shafferの方法のための可能な真の帰無仮説の数の計算アルゴリズムを追加；mod.Bon関数の変更とshaffer2，fshaffer1，fshaffer2オプションの追加
        # ・多重比較の出力において調整済みｐ値を表示する機能を追加（デフォルトに設定）；mod.Bon関数の変更とそれに伴う変更
        # 6)anovakun version 3.1.0（R 2.5.1，R 2.6.0作成；2008/02/01公開）
        # ・aov関数のアルゴリズムを利用して誤差平方和の計算の高速化を試みる；ss.calc関数の変更；その他の修正
        # 7)anovakun version 3.2.0（R 2.5.1，R 2.6.0作成；2008/04/01公開）
        # ・Shafferの方法のための別バージョンのアルゴリズムを追加；mod.Bon関数の変更とs2d，fs2dオプションの追加
        # ・オプション名の変更；“shaffer2，fshaffer1，fshaffer2”を“s2r，fs1，fs2r”に変更
        # ・効果量とp_repを計算する機能を追加；anova.modeler関数，anovatab.out関数ほかの変更
        # 8)anovakun version 4.0.0（R 2.5.1，R 2.6.0，R 2.7.0作成；2008/06/02公開）
        # ・Mendozaの多標本球面性検定とHarrisの多標本球面性検定を行う機能を追加（Mendozaをデフォルトに設定）；epsilon.calc関数と関連する部分の変更
        # ・Huynhの改良版一般近似検定とAlgina-Lecoutreの修正改良版一般近似検定を行う機能を追加；epsilon.calc関数及び関連する箇所の変更
        # ・取り扱い可能な要因の数を拡張；anovakun関数，anova.modeler関数，epsilon.calc関数，pro.fraction関数ほかの変更
        # ・複数の効果量をオプション指定したときに，すべての指標が同時に出力されるように変更；anova.modeler関数，pro.fraction関数ほかの変更
        # 9)anovakun version 4.1.0（R 2.9.0，R 2.9.1，R 2.9.2作成；2009/09/01公開）
        # ・欠損値があるケース（数値以外のデータを含む行）を除外して分析するように変更
        # ・混合要因計画で被験者間要因のみを取り出して下位検定を行う際に警告メッセージが表示される点を修正
        # 10)anovakun version 4.1.1（R 2.9.2，R 2.10.0，R 2.10.1作成；2010/01/04公開）
        # ・Rのバージョン情報の表示を修正
        # ・出力関数の表示エラーを修正；anovatab.out関数，simeffects.out関数の修正
        # 11)anovakun version 4.2.0（R 2.10.0，R 2.10.1，R 2.11.1作成；2010/07/01公開）
        # ・効果量の指標にイータ二乗，一般化イータ二乗，オメガ二乗，一般化オメガ二乗を追加；anova.modeler関数，anovatab.out関数ほかの変更
        # ・データフレームの作成方法を変更；被験者間ラベルに文字列を使用しても指定した水準とデータがずれないようにする
        # ・Lecoutre（1991）に基づいてHyunh-Feldtのイプシロンの計算式を変更；epsilon.calc関数の修正
        # 12)anovakun version 4.3.0（R 2.14.1，R 2.14.2，R 2.15.0作成；2012/07/02公開）
        # ・効果量の指標にイプシロン二乗，偏イプシロン二乗，一般化イプシロン二乗を追加；anova.modeler関数，pro.fraction関数ほかの変更
        # ・オメガ二乗，偏オメガ二乗，一般化オメガ二乗の計算式を修正；anova.modeler関数，pro.fraction関数ほかの修正
        # ・非加算モデルに基づくオメガ二乗，一般化オメガ二乗のオプションを追加；anova.modeler関数，pro.fraction関数ほかの修正
        # ・結果をクリップボードにも出力するオプションcopyを追加（Linuxでは無効）；anovakun関数の変更
        # 13)anovakun version 4.3.1（R 2.15.0，R 2.15.1作成；2012/09/03公開）
        # ・copy機能がMacで働かないエラーを修正；anovakun関数の修正
        # 14)anovakun version 4.3.2（R 2.15.2作成；2013/01/04公開）
        # ・read.clip関数を追加
        # ・nopostオプションを追加；anovakun関数，post.analysis関数の変更
        # ・copy機能がLinuxで働くように修正（xclipをインストールしている場合のみ機能）；anovakun関数の変更
        # ・出力の桁数に合わせて表の長さを調節するように変更；bstat.out関数，anovatab.out関数，mod.Bon.out関数，simeffects.out関数，post.out関数，each.out関数の変更
        # 15)anovakun version 4.3.3（R 2.15.3，R 3.0.0作成；2013/05/07公開）
        # ・R 3.0.0で動作するように修正；ss.calc関数の修正
        # ・３要因以上の分析で一次の交互作用が有意だった際の多重比較の独立・非独立測度への割り当てが誤っていたのを修正；pro.fraction関数の修正
        # ・出力桁数の調整；anovatab.out関数ほかの修正
        # 16)anovakun version 4.4.0（R 3.0.1，R 3.0.2作成；2013/11/01公開）
        # ・データフレームの変形手順を変更し，ロング形式のデータも扱えるように修正；anovakun関数の変更
        # ・出力時の要因名，水準名を任意に指定できるように変更；anovakun関数，anova.modeler関数，post.analyses関数，pro.fraction関数ほかの変更
        # ・出力表示の調整；bstat.out関数，anovatab.out関数，mod.Bon.out関数，simeffects.out関数，post.out関数，each.out関数の変更
        # ・一般化効果量における複数の測定要因を指定する方法を変更（cまたはlist形式を使用する）；anova.modeler関数，pro.fraction関数の変更
        # 17)anovakun version 4.5.0（R 3.0.2作成；2014/02/03公開）
        # ・信頼区間を計算する機能を追加；ci.calc関数の追加，anovakun関数，bstat.out関数の変更
        # ・データフレーム変形部分をuni.long関数として分離；uni.long関数の追加，anovakun関数の変更
        # ・３要因以上の計画について単純効果を扱うための関数を追加；anovatan関数の追加
        #
        # 【この関数の使用に関して】
        # 1)anovakunとこれを構成する関数（コード群）は，自由に使用，改変，再配布していただいて結構です。
        # 2)ただし，改変を加えたものを公開する際には，改変を加えたことを明記し，メインの関数の名前をanovakun以外の名前に変更してください。
        # 3)anovakunとこれを構成する関数（コード群）の使用によって生じるいかなる結果に関しても作成者は責任を負いかねますのでご了承ください。
        
        # 調和平均を計算する関数
        hmean <- function(datvector){
            return(length(datvector)/sum(1/datvector))
        }
        
        
        # クリップボードの情報を読み込む関数
        # read.tableのラッパー関数；read.tableの出力先以外のオプションをすべて指定できる
        # OSにかかわらず同じ書式で機能するようにしてある；ただし，LinuxはUbuntuでのみテストしており，xclipをインストールしていることを前提とする
        read.clip <- function(file, header = FALSE, sep = "", quote = "\"'", dec = ".", row.names, col.names, as.is = !stringsAsFactors,
        na.strings = "NA", colClasses = NA, nrows = -1, skip = 0, check.names = TRUE, fill = !blank.lines.skip,
        strip.white = FALSE, blank.lines.skip = TRUE, comment.char = "#", allowEscapes = FALSE, flush = FALSE,
        stringsAsFactors = default.stringsAsFactors(), fileEncoding = "", encoding = "unknown", text){
            
            # OSごとのクリップボードを出力先に指定
            plat.info <- .Platform
            if(sum(grep("windows", plat.info)) !=0){# Windowsの場合
                outboard <- "clipboard"
            }else if(sum(grep("mac", plat.info)) != 0){# Macの場合
                outboard <- pipe("pbpaste")
            }else if(sum(grep("linux", R.version$system)) != 0){# Linxの場合（xclipをインストールしている必要がある）
                system("xclip -o | xclip -sel primary")
                outboard <- "clipboard"
            }else{# いずれのOSでもない場合
                stop(message = "Unknown operating system!!")
            }
            
            # read.table関数を実行
            read.table(file = outboard, header = header, sep = sep, quote = quote, dec = dec, row.names = row.names, col.names =
            col.names, as.is = as.is, na.strings = na.strings, colClasses = colClasses, nrows = nrows, skip = skip,
            check.names = check.names, fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip,
            comment.char = comment.char, allowEscapes = allowEscapes, flush = flush, stringsAsFactors = stringsAsFactors,
            fileEncoding = fileEncoding, encoding = encoding, text = text)
        }
        
        
        # ANOVA君本体：プロセス全体の制御を行う関数
        anovakun <- function(dataset, design, ..., long = FALSE, type2 = FALSE, nopost = FALSE, tech = FALSE, data.frame = FALSE, copy = FALSE,
        holm = FALSE, hc = FALSE, s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE,
        lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE,
        eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE,
        gomega = NA, gomegana = NA, prep = FALSE, cind = FALSE, cin = FALSE){
            maxfact <- nchar(design) - 1# 実験計画全体における要因数
            
            # データフレームの変形
            datform <- uni.long(dataset = dataset, design = design, ... = ..., long = long)
            dat <- datform$dat
            factnames <- datform$factnames
            eachlev <- datform$eachlev
            miscase <- datform$miscase
            
            # 記述統計量を計算する
            baseresults <- ci.calc(dat = dat, design = design, factnames = factnames, cind = cind, cin = cin)
            
            # 記述統計量の表を出力する際の改行位置を指定する
            if(maxfact < 3) margin <- prod(eachlev)
            else margin <- prod(eachlev[-(1:(maxfact-2))])
            
            # anova.modelerにデータフレームを送り，分散分析の結果を得る
            mainresults <- anova.modeler(dat = dat, design = design, factnames = factnames, type2 = type2,
            lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
            eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
            gomega = gomega, gomegana = gomegana, prep = prep)
            
            # post.analysesにデータフレームと分散分析の結果を送り，下位検定の結果を得る
            postresults <- post.analyses(dat = dat, design = design, factnames = factnames, mainresults = mainresults, type2 = type2,
            nopost = nopost, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria,
            lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
            eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
            gomega = gomega, gomegana = gomegana, prep = prep)
            
            # 基本情報の取得
            info1 <- paste0("[ ", design, "-Type Design ]")# 要因計画の型
            info2 <- paste0("This output was generated via anovakun 4.5.0 at ", strsplit(R.version$version.string, " \\(")[[1]][1], ".")# バージョン情報など
            info3 <- paste0("It was executed on ", date(), ".")# 実行日時
            
            # Unbalancedデザイン（データ数ふぞろい）の場合，プロンプトを追加
            if(length(unique(baseresults$bstatist$n)) != 1){
                if(sum(is.na(mainresults$ano.info1)) == 1){
                    if(type2 == TRUE) mainresults$ano.info1 <- c("== This data is UNBALANCED!! ==", "== Type II SS is applied. ==")
                    else mainresults$ano.info1 <- c("== This data is UNBALANCED!! ==", "== Type III SS is applied. ==")
                }else{
                    if(type2 == TRUE) mainresults$ano.info1 <- append(mainresults$ano.info1, c("== This data is UNBALANCED!! ==", "== Type II SS is applied. =="))
                    else mainresults$ano.info1 <- append(mainresults$ano.info1, c("== This data is UNBALANCED!! ==", "== Type III SS is applied. =="))
                }
            }
            
            # 除外したケースの報告
            if(miscase != 0){
                if(sum(is.na(baseresults$bstat.info1)) == 1){
                    baseresults$bstat.info1 <- paste0("== The number of removed case is ", miscase, ". ==")
                }
                else{
                    baseresults$bstat.info1 <- append(baseresults$bstat.info1, paste0("== The number of removed case is ", miscase, ". =="))
                }
            }
            
            # 結果を表示する
            if(copy == TRUE){# 指定があった場合，出力をクリップボードにコピー
                plat.info <- .Platform
                if(sum(grep("windows", plat.info)) != 0){# Windowsの場合
                    sink("clipboard", split = TRUE)
                }else if(sum(grep("mac", plat.info)) != 0){# Macの場合
                    tclip <- pipe("pbcopy", "w")
                    sink(tclip, split = TRUE)
                }else if(sum(grep("linux", R.version$system)) != 0){# Linxの場合（xclipをインストールしている必要がある）
                    tclip <- pipe("xclip -selection clipboard")
                    sink(tclip, split = TRUE)
                }
            }
            if(tech == TRUE){# データフレーム形式での出力の場合
                if(data.frame == TRUE){
                    names(dat) <- c("s", factnames, "y")
                    return(list("INFORMATION" = rbind(info1, info2, info3),
                    "DESCRIPTIVE STATISTICS" = baseresults,
                    "SPHERICITY INDICES" = list(mainresults$epsi.info1, mainresults$epsitab),
                    "ANOVA TABLE" = list(mainresults$ano.info1, mainresults$anovatab),
                    "POST ANALYSES" = postresults,
                    "DATA.FRAME" = dat))# 計算に使用したデータフレームを付加
                }else{
                    return(list("INFORMATION" = rbind(info1, info2, info3),
                    "DESCRIPTIVE STATISTICS" = baseresults,
                    "SPHERICITY INDICES" = list(mainresults$epsi.info1, mainresults$epsitab),
                    "ANOVA TABLE" = list(mainresults$ano.info1, mainresults$anovatab),
                    "POST ANALYSES" = postresults))
                }
            }else{# 表形式での出力の場合
                if(data.frame == TRUE){
                    names(dat) <- c("s", factnames, "y")
                    info.out(info1, info2, info3)
                    bstat.out(baseresults, maxfact, margin)
                    anovatab.out(mainresults)
                    post.out(postresults, design)
                    return(list("DATA.FRAME" = dat))# 計算に使用したデータフレームを付加
                }else{
                    info.out(info1, info2, info3)
                    bstat.out(baseresults, maxfact, margin)
                    anovatab.out(mainresults)
                    post.out(postresults, design)
                }
            }
            if(copy == TRUE){
                sink()
                if(plat.info$OS.type != "windows"){# Mac，Linuxの場合
                    close(tclip)
                }
            }
        }
        
        
        # データフレームをロング形式に変形する関数
        uni.long <- function(dataset, design, ..., long = FALSE){
            bet.with <- strsplit(design, "s")
            betN <- nchar(bet.with[[1]])[1]# 被験者間要因がないときは“０”
            withN <- nchar(bet.with[[1]])[2]# 被験者内要因がないときは“NA”
            maxfact <- nchar(design) - 1# 実験計画全体における要因数
            levlist <- list(...)
            
            # データフレームの変形
            if(long == TRUE){# ロング形式のデータフレームのヘッダをANOVA君に適合したものに変更
                # 指定されたデザインとデータフレームの列数が合致しない場合にはメッセージを表示して終了
                if((maxfact + 2) != ncol(dataset)){
                    stop(message = "\"anovakun\" has stopped working...\nThe entered design does not match the data.")
                }
                
                # ヘッダの差し替え
                misv <- suppressWarnings(as.numeric(sapply(dataset[, ncol(dataset)], function(x) as.vector(x))))# 数値以外の値をNAに強制変換
                dat <- cbind(data.frame(lapply(dataset[, 1:(maxfact+1), drop = FALSE], function(x) factor(x, levels = unique(x)))), misv)# ラベルをfactor型に変換
                if(setequal(names(dataset), paste0("V", 1:ncol(dataset))) == TRUE | is.null(names(dataset)) == TRUE){# 入力時に指定されていない場合
                    factnames <- LETTERS[1:maxfact]# 要因名のラベルとしてアルファベットを使用
                    levlist <- lapply(1:maxfact, function(x) paste0(letters[x], 1:length(unique(dataset[, x+1]))))
                    for(i in 1:maxfact){# 被験者間要因の水準名を小文字アルファベットにする
                        levels(dat[, i+1]) <- levlist[[i]]# 水準名を上書き
                    }
                }else{# 入力時に指定されている場合
                    factnames <- names(dataset)[-c(1, maxfact+2)]# もとのヘッダを保存；要因名のラベルとして使用
                }
                
                names(dat) <- c("s", LETTERS[1:maxfact], "y")
                eachlev <- sapply(dat, nlevels)[-c(1, maxfact+2)]# 各要因の水準数
                attributes(eachlev) <- NULL# attributesを消去
                
                # 欠損値の除去
                if(sum(is.na(misv)) == 0){# 欠損がなかった場合
                    miscase <- 0# 欠損ケースの数
                }else{# 欠損があった場合
                    misid <- dataset[is.na(misv), 1]# 欠損のあったID
                    dat <- dat[rowSums(sapply(misid, function(x) x == dat[, 1])) == 0, ]# datからNAを含むIDの行（ケース）を除く
                    miscase <- length(misid)# 欠損ケースの数
                }
            }else{# ワイド形式のデータフレームをロング形式に変換
                # 各要因の水準数と要因ラベルの設定
                if(is.null(names(levlist)) == TRUE){# 入力時に指定されていない場合
                    eachlev <- unlist(levlist)# 各要因の水準数
                    factnames <- LETTERS[1:maxfact]# 要因名のラベル
                    levlist <- lapply(1:length(eachlev), function(x) paste0(letters[x], 1:eachlev[x]))
                }else{# 入力時に指定されている場合
                    eachlev <- sapply(levlist, length)# 各要因の水準数
                    attributes(eachlev) <- NULL# attributesを消去
                    factnames <- names(levlist)# 要因名のラベル
                }
                
                # 水準数１が指定されていたらメッセージを表示して終了；存在しないオプション名が指定されると水準数１と解釈される
                if(min(eachlev) == 1){
                    stop(message = "\"anovakun\" has stopped working...\nSome factor specifies only one level.\nMaybe non-existent option-names are requested.")
                }
                
                # 指定されたデザインとデータフレームの列数が合致しない場合にはメッセージを表示して終了
                if((betN + ifelse(is.na(prod(eachlev[(betN+1):length(eachlev)])), 1, prod(eachlev[(betN+1):length(eachlev)]))) != ncol(dataset)){
                    stop(message = "\"anovakun\" has stopped working...\nThe entered design does not match the data.")
                }
                
                # 欠損値の除去
                misv <- suppressWarnings(as.numeric(sapply(dataset[, (betN+1):ncol(dataset)], function(x) as.vector(x))))# 数値以外の値をNAに強制変換
                cdata <- array(misv, c(nrow(dataset), ncol(dataset) - betN))# NAを含む行列
                compcase <- complete.cases(cdata)# 完全ケース
                dataset <- dataset[compcase,]# datasetからNAを含む行（ケース）を除く
                miscase <- sum(!compcase)# 欠損ケースの数
                depv <- as.vector(cdata[compcase, ])# 従属変数をベクトル化
                slab <- rep(paste0("s", 1:nrow(dataset)), ncol(dataset) - betN)# 被験者のラベル
                dat <- data.frame(factor(slab, levels = unique(slab)))# 被験者のラベルをデータフレームにする
                
                # 被験者間要因を含む場合
                if(betN > 0){
                    betlab <- data.frame(lapply(dataset[, 1:betN, drop = FALSE], function(x) factor(x, levels = unique(x))))# 被験者間要因のラベルをfactor型に変換；水準名をもとの順序にそろえる
                    for(i in 1:betN){# 被験者間要因の水準名を小文字アルファベットにする
                        levels(betlab[, i]) <- levlist[[i]]# 水準名を上書き
                    }
                    betlab <- do.call("rbind", replicate(ncol(dataset) - betN, betlab, simplify = FALSE))# くりかえし数のぶんだけ増やす
                    dat <- cbind(dat, betlab)# データフレームに追加
                }
                # 被験者内要因を含む場合
                if(is.na(withN) == FALSE){
                    eachcue <- c(sapply(1:withN, function(x) prod(eachlev[(betN+x):maxfact])), 1)[-1]
                    timescue <- prod(eachlev[(betN+1):maxfact]) / (eachcue * eachlev[(betN+1):maxfact])
                    withlab <- data.frame(lapply(1:withN, function(x) factor(rep(levlist[[betN+x]], each = nrow(dataset) * eachcue[x],
                    times = timescue[x]), levels = levlist[[betN+x]])))# 各被験者内要因の各水準を表すデータフレーム
                    dat <- cbind(dat, withlab)# データフレームに追加
                }
                dat <- cbind(dat, depv)# 従属変数をデータフレームに追加
                names(dat) <- c("s", LETTERS[1:maxfact], "y")
            }
            return(list("dat" = dat, "factnames" = factnames, "eachlev" = eachlev, "miscase" = miscase))
        }
        
        
        # 信頼区間を計算する関数
        ci.calc <- function(dat, design, factnames = NA, conf.level = 0.95, cind = FALSE, cin = FALSE){
            bet.with <- strsplit(design, "s")
            betN <- nchar(bet.with[[1]])[1]# 被験者間要因がないときは“０”
            withN <- nchar(bet.with[[1]])[2]# 被験者内要因がないときは“NA”
            maxfact <- nchar(design) - 1
            cellN <- length(unique(dat$s))# 被験者間要因をつぶしての全被験者の数
            flev <- sapply(2:(maxfact+1), function(x) nlevels(dat[, x]))# 各要因の水準数
            comlev <- prod(flev[(betN+1):maxfact])# すべての被験者内要因の組み合わせ水準数
            
            # factnamesが省略されているときは名前をつける
            if(sum(is.na(factnames)) != 0) factnames <- LETTERS[1:maxfact]
            
            # 各条件ごとの平均と標準偏差を計算する
            tabbase <- tapply(dat$y, dat[, (maxfact+1):2], function(x) x)
            sncol <- sapply(tabbase, function(x) length(x))# セルごとのデータ数を計算
            mncol <- sapply(tabbase, function(x) mean(x, na.rm = TRUE))# セルごとの平均を計算
            sdcol <- sapply(tabbase, function(x) sd(x, na.rm = TRUE))# セルごとの標準偏差を計算
            
            # 記述統計量の表において各要因の各水準を表すためのラベル（数字列）を作成
            maincols <- expand.grid(lapply((maxfact+1):2, function(x) levels(dat[, x])))
            maincols <- maincols[, order(maxfact:1)]# アルファベット順に列を並べ替え
            
            # 記述統計量をデータフレームにまとめる
            bstatist <- data.frame(maincols, sncol, mncol, sdcol)
            names(bstatist) <- c(factnames, "n", "Mean", "S.D.")# 要因のラベルほかを列名に設定
            bstat.info1 <- NA
            
            # 信頼区間の計算
            if(cind == TRUE || cin == TRUE){
                # 信頼区間の調整用の因数
                if(cind == TRUE){# CIの重なりと有意性検定の結果が対応するように調整
                    diff.factor <- sqrt(2)/2
                    id.part <- "== Difference-Adjusted Confidence Intervals for Independent Means =="
                    rm.part <- "== Cousineau-Morey-Baguley's Difference-Adjusted Confidence Intervals =="
                }else if(cin == TRUE){# 調整なし
                    diff.factor <- 1
                    id.part <- "== Confidence Intervals for Independent Means =="
                    rm.part <- "== Cousineau-Morey's Confidence Intervals =="
                }
                
                # 計画に合った計算法を適用
                if(is.na(withN) == TRUE){# 被験者間計画の場合
                    bstat.info1 <- id.part
                    cdf <-  cellN - prod(flev[1:betN])
                    vcom <- sqrt(sum(tapply(dat$y, dat[, (betN+1):2], function(x) sum((x - mean(x))^2))) / cdf * mean(1/sncol))
                    ci.tval <- diff.factor * qt((1 - conf.level)/2, cdf, lower.tail = FALSE) * vcom
                }else{# 反復測定要因を含む場合
                    bstat.info1 <- rm.part
                    if(betN == 0){# 被験者間要因がないとき（被験者内計画の場合）
                        caldat <- list(dat)
                    }else{# １つ以上の被験者間要因があるとき（混合要因計画の場合）
                        caldat <- split(dat, dat[, (maxfact+1-withN):(betN+1)])
                    }
                    wdat <- lapply(caldat, function(x) data.frame(matrix(x$y[order(eval(parse(text = paste0("x[, ", c(1, (maxfact+1):(betN+2)), "]"))))], length(x$y)/comlev, comlev)))# データをワイド形式のデータフレームに変換；tapplyによるアクセス順に対応させるための並べ替え
                    normdat <- lapply(wdat, function(x) x - matrix(rowMeans(x, na.rm = TRUE), nrow(x), 1) + mean(colMeans(x, na.rm = TRUE), na.rm = TRUE))# 個人ごとのデータによって正規化したデータ
                    nsd <- unlist(lapply(normdat, function(x) sapply(x, function(y) sd(y, na.rm = TRUE))), use.names = FALSE)# 正規化したデータを使って標準偏差を計算
                    ci.tval <- diff.factor * qt((1 - conf.level)/2, sncol - 1, lower.tail = FALSE) * sqrt(comlev / (comlev - 1)) * nsd/sqrt(sncol)
                }
                
                lwlab <- paste0(100 * conf.level, "%CI-L")
                uplab <- paste0(100 * conf.level, "%CI-U")
                citab <- data.frame(mncol - ci.tval, mncol + ci.tval)
                names(citab) <- c(lwlab, uplab)
                bstatist <- cbind(bstatist, citab)
            }
            return(list("bstat.info1" = bstat.info1, "bstatist" = bstatist))
        }
        
        
        # 文字列中のすべての要素を含む文字列をマッチングする関数
        # grepとの違いは“A:C”などの文字列を照合パターンとした場合に“A:B:C”のように間に別の文字を挟んだ文字列もマッチと判定する点
        # 照合パターンが１文字の場合はgrepと同じ結果を返す
        elematch <- function(Mstrings, stex){
            # マッチングする文字列を分解して，それぞれgrep関数を適用
            matchlist <- lapply(strsplit(Mstrings, "")[[1]], function(x) grep(x, stex))
            
            # 文字列の各要素とマッチした値の共通部分のみ取り出す
            buffer <- matchlist[[1]]
            
            # 文字列が１文字のときはgrepの結果をそのまま返す
            if(length(matchlist) != 1){
                for(i in 2:length(matchlist)){
                    buffer <- buffer[is.element(buffer, matchlist[[i]])]
                }
            }
            return(buffer)
        }
        
        
        # 平方和を計算する関数
        ss.calc <- function(full.elem, dat, type2 = FALSE){
            # full.elemの並べ替え
            elem.num <- seq(1, max(nchar(full.elem)), by = 2)
            full.elem <- unlist(sapply(elem.num, function(x) full.elem[nchar(full.elem) == x]))
            
            # 計画行列を作る
            if(sum(grep("s", full.elem)) == 0) eff.elem <- full.elem
            else eff.elem <- full.elem[-grep("s", full.elem)]# 誤差項との交互作用効果を除いたもの
            
            er.elem <- grep("s", full.elem, value = TRUE)# 誤差項のみ
            eff.modeleq <- paste0("~ ", gsub(",", " +", toString(eff.elem)))
            er.modeleq <- paste0("~ ", gsub(",", " +", toString(er.elem)))
            
            def.contr <- options()$contrasts# contrastsのデフォルト設定を保存
            options(contrasts = c("contr.sum", "contr.sum"))# 設定を変更
            dmat <- model.matrix(as.formula(eff.modeleq), dat)
            
            # 計画行列とデータを統合する
            exmat <- as.matrix(cbind(dmat, dat$y))# 拡大行列を作る
            promat <- crossprod(exmat)# 拡大行列の積和行列
            endline <- nrow(promat)# 積和行列の列数
            
            # 各効果に対応する計画行列の列（行）番号を得る
            sepcol <- attr(dmat, "assign")# 計画行列からピボットを表すベクトルを取り出す
            pivot.col <- lapply(1:max(sepcol), function(x) (1:length(sepcol))[sepcol == x])# 各効果を表現する列の番号
            names(pivot.col) <- eff.elem
            
            if(type2 == TRUE){# 線形モデルを用いてタイプⅡ平方和を計算する
                # 各効果の平方和の計算
                # 各モデルのための部分行列を選択
                ss.line1 <- lapply(eff.elem, function(x) c(1, unlist(pivot.col[match(x, names(pivot.col))]), unlist(pivot.col[-elematch(x, names(pivot.col))])))
                ss.line2 <- lapply(eff.elem, function(x) c(1, unlist(pivot.col[-elematch(x, names(pivot.col))])))
                
                # 各モデルの最小二乗解を得てもとのベクトルにかけたものの和を取る
                ss.base1 <- sapply(ss.line1, function(x) sum(qr.coef(qr(promat[x, x], LAPACK = TRUE), promat[x, endline]) * promat[x, endline]))
                ss.base2 <- sapply(ss.line2, function(x) sum(qr.coef(qr(promat[x, x], LAPACK = TRUE), promat[x, endline]) * promat[x, endline]))
                
                # 各効果を含むモデルと含まないモデルの差分を取る
                ss.all <- ss.base1 - ss.base2
                ss.all <- lapply(ss.all, function(x) x)# リスト化
                names(ss.all) <- eff.elem
                
            }else{# 線形モデルを用いてタイプⅢ平方和を計算する
                # 各効果の平方和の計算
                eff.line <- c(1, unlist(pivot.col))
                
                # 各モデルのための部分行列を選択
                ss.line <- lapply(eff.elem, function(x) eff.line[-match(pivot.col[[match(x, names(pivot.col))]], eff.line)])
                
                # 各モデルの最小二乗解を得てもとのベクトルにかけたものを合計
                ss.eff <- sum(qr.coef(qr(promat[eff.line, eff.line], LAPACK = TRUE), promat[eff.line, endline]) * promat[eff.line, endline])
                ss.base <- sapply(ss.line, function(x) sum(qr.coef(qr(promat[x, x], LAPACK = TRUE), promat[x, endline]) * promat[x, endline]))
                
                # 各効果を含むモデルと含まないモデルの差分を取る
                ss.all <- ss.eff - ss.base
                ss.all <- lapply(ss.all, function(x) x)# リスト化
                names(ss.all) <- eff.elem
            }
            
            # 全体平方和の計算
            ss.T <- promat[endline, endline] - sum(qr.coef(qr(promat[1, 1], LAPACK = TRUE), promat[1, endline]) * promat[1, endline])
            
            # 誤差平方和の計算
            if(sum(grep("s", full.elem)) == 0){# 誤差項を分解する必要がない場合
                ss.Er <- promat[endline, endline] - sum(qr.coef(qr(promat[1:(endline-1), 1:(endline-1)], LAPACK = TRUE), promat[1:(endline-1), endline]) * promat[1:(endline-1), endline])
            }else{# 誤差項を分解する必要がある場合
                ss.Er <- NA
                emat <- model.matrix(as.formula(er.modeleq), dat)# 切片と誤差項のみの計画行列
                er.num <- length(er.elem)# 誤差項モデルの項の数
                
                qr.er <- qr(emat, LAPACK = TRUE)# 誤差項モデルのQR分解
                er.col <- attr(emat, "assign")[qr.er$pivot[1:qr.er$rank]]# 有効ピボット
                tq.er <- t(qr.Q(qr.er))
                qty <- as.matrix(tq.er %*% dat$y)
                qtx <- tq.er %*% dmat
                
                ss.er <- rep(list(NA), er.num)# 誤差平方和格納用のリストを宣言
                names(ss.er) <- er.elem
                
                for(i in 1:er.num){
                    select <- er.col == i
                    xi <- qtx[select, , drop = FALSE]
                    cols <- colSums(xi^2) > 1e-05
                    ss.er[[i]] <- sum(qr.resid(qr(xi[, cols, drop = FALSE]), qty[select, , drop = FALSE]) * qty[select, , drop = FALSE])
                }
                ss.all <- c(ss.all, ss.er)
            }
            
            ss.results <- c(ss.all, "ss.T" = ss.T, "ss.Er" = ss.Er)
            options(contrasts = def.contr)# contrastsの設定をもどす
            
            return(ss.results)
        }
        
        
        # 有意水準に合わせて記号を表示する関数
        sig.sign <- function(pvalue){
            ifelse(is.na(pvalue), "",
            ifelse(pvalue < 0.001, "***",
            ifelse(pvalue < 0.01, "**",
            ifelse(pvalue < 0.05, "*",
            ifelse(pvalue < 0.10, "+", "ns")))))
        }
        
        
        # Greenhouse-GeisserとHuynh-Feldtのイプシロンを計算する関数
        # 被験者内要因を含まない計画を投入すると適切に動作しないので注意
        epsilon.calc <- function(dat, design, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE){
            # 要因計画の型から被験者内要因を特定
            bet.with <- strsplit(design, "s")
            
            # 被験者間要因の特定
            othlabel <- match(strsplit(bet.with[[1]][1], "")[[1]], names(dat))
            othnum <- othlabel[order(othlabel, decreasing = TRUE)]# 後の方の要因のラベルを先に並べる
            othmat <- sapply(othnum, function(x) nlevels(dat[,x]))# 各要因の水準数をベクトル化
            ol <- ifelse(length(othmat) == 0, 1, prod(othmat))# 全被験者間要因の組み合わせ水準数を取得；被験者間要因がなければ，１を代入
            
            # 被験者内要因の特定
            replabel <- match(strsplit(bet.with[[1]][2], "")[[1]], names(dat))
            repnum <- replabel[order(replabel, decreasing = TRUE)]# 後の方の要因のラベルを先に並べる
            repmat <- sapply(repnum, function(x) nlevels(dat[,x]))# 各要因の水準数をベクトル化
            rl <- prod(repmat)# 全被験者内要因の組み合わせ水準数を取得
            
            cellN <- length(unique(dat$s))# 被験者間要因をつぶしての全被験者の数を取得
            
            if(length(othmat) == 0) othN <- cellN# 被験者間要因がないときはcellNと同じ値を代入
            else othN <- as.vector(table(dat[names(dat)[othnum]]) / rl)# 被験者間要因の各組み合わせにおける被験者数をベクトル化
            
            # データフレームを分割し，共分散行列を作る
            if(length(othnum) == 0){# 被験者間要因がないときはデータフレームを分割しない
                covmatrices <- cov(as.data.frame(split(dat$y, dat[names(dat)[repnum]])))
            }else{# 被験者間要因の組み合わせ水準ごとにデータフレームを分割
                covmatrices <- lapply(split(dat, dat[names(dat)[othnum]]), function(x) cov(as.data.frame(split(x$y, x[names(x)[repnum]]))))
            }
            
            # データフレームのリストを三次元配列に変換
            covbuffer <- array(unlist(covmatrices), dim = c(rl, rl, ol))
            
            # 複数の共分散行列をプール
            tm <- apply(covbuffer, c(1, 2), function(x) sum((othN-1) * x)) / (cellN - ol)
            
            # 正規直交対比行列を作る；被験者内要因の数に合わせて異なるパターンを得る
            replev <- sapply(replabel, function(x) nlevels(dat[, x]))# 計画タイプの順に各被験者内要因の水準数を得る
            repleng <- length(replabel)# 反復測定要因の数
            
            def.contr <- options()$contrasts# contrastsのデフォルト設定を保存
            options(contrasts = c("contr.helmert", "contr.helmert"))# 設定を変更
            
            ortho.model <- expand.grid(lapply(repnum, function(x) levels(dat[, x])))
            names(ortho.model) <- paste0("R", repleng:1)
            ortho.helm <- model.matrix(as.formula(paste0("~ ", gsub(",", " *", toString(paste0("R", 1:repleng))))), ortho.model)
            matcue <- attr(ortho.helm, "assign")
            matdivider <- sapply(1:max(matcue), function(x) length(matcue[matcue == x]))# 正規直交行列を分割する際の行数
            
            options(contrasts = def.contr)# 設定をもどす
            effect.name <- unlist(sapply(1:repleng, function(y) combn(names(dat)[replabel], y, function(x) gsub(", ", "x", toString(x)))))# 効果のラベルを作る
            
            ortho.coef <- t(ortho.helm[, 2:rl, drop = FALSE])# 直交対比のパターンのみを取り出す；行が１のときにベクトルに変換されないようにdrop = FALSEを使う
            ortho.denomi <- rowSums(ortho.coef^2)^(1/2)
            
            # パターンを直交対比行列にする
            orthoM <- ortho.coef / ortho.denomi
            
            # 共分散行列と正規直交対比行列をかける
            otoM <- orthoM %*% tm %*% t(orthoM)
            
            # 行列全体を使っての球面性検定
            if(mau == TRUE){
                # プロンプトの準備
                epsi.info1 <- paste0("== Mauchly's Sphericity Test and Epsilons ==")
                lamlab <- "W"
                
                # Mauchlyの球面性検定
                eps.Lambda <- det(otoM) / (sum(diag(otoM)) / (rl - 1))^(rl - 1)
                eps.m <- 1 - (2 * rl^2 - 3 * rl + 3) / (6 * (cellN - ol) * (rl - 1))
                
                if(cellN <= rl){# 被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
                    epsChi <- NA
                    epsi.info1 <- paste0(epsi.info1, "\n",
                    "*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n",
                    "*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***")
                }else{
                    epsChi <- -(cellN - ol) * eps.m * log(eps.Lambda)
                }
                
                eps.df <- (((rl - 1) * rl) / 2) - 1
                eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
                
            }else if(har == TRUE){
                # プロンプトの準備
                epsi.info1 <- paste0("== Harris's Multisample Sphericity Test and Epsilons ==")
                lamlab <- "h_hat"
                
                # Harrisの多標本球面性検定
                proA <- array(apply(covbuffer, 3, function(x) orthoM %*% x %*% t(orthoM)), dim = c(rl-1, rl-1, ol))
                
                if(cellN <= rl){# 被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
                    epsChi <- NA
                    eps.Lambda <- NA
                    epsi.info1 <- paste0(epsi.info1, "\n",
                    "*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n",
                    "*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***")
                }else{
                    harTr <- apply(proA, 3, function(x) sum(diag(x)))
                    harSq <- apply(proA, 3, function(x) sum(diag(x %*% x)))
                    eps.Lambda <- sum((othN - 1) * harTr)^2 / sum((othN - 1) * harSq)
                    epsChi <- pmax(0, ((cellN - ol) * (rl - 1) / 2) * ((cellN - ol) * (rl - 1) / eps.Lambda - 1))# 負の値は０にそろえる
                }
                
                eps.df <- (ol * (rl - 1) * rl) / 2 - 1
                eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
                
            }else{
                # プロンプトの準備
                epsi.info1 <- paste0("== Mendoza's Multisample Sphericity Test and Epsilons ==")
                lamlab <- "Lambda"
                
                # Mendozaの多標本球面性検定
                ptm <- covbuffer * rep(othN, each = rl^2)
                proA <- array(apply(ptm, 3, function(x) orthoM %*% x %*% t(orthoM)), dim = c(rl-1, rl-1, ol))
                
                if(cellN <= rl){# 被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
                    epsChi <- NA
                    eps.Lambda <- NA
                    epsi.info1 <- paste0(epsi.info1, "\n",
                    "*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n",
                    "*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***")
                }else{
                    menL1 <- log((cellN-ol)^((cellN-ol)*(rl-1)/2) / prod((othN-1)^((othN-1)*(rl-1)/2)))
                    menL2 <- apply(proA, 3, function(x) det(x))
                    menL2 <- sum(log(ifelse(menL2 < 0, NA, menL2)) * (othN-1)/2)# 行列式が負になったときは対数にできないのでNAを代入
                    menL3 <- log(sum(diag(apply(proA, c(1, 2), function(x) sum(x))/(rl-1)))) * ((cellN-ol)*(rl-1)/2)
                    menL <- menL1 + menL2 - menL3
                    eps.m <- 1 - ((((cellN-ol) * (rl-1)^2 * rl * (2*(rl-1)+1) - (2*(cellN-ol)*(rl-1)^2)) * sum(1/(othN-1)) - 4) / (6*(cellN-ol)*(rl-1) * (ol*(rl-1)*rl-2)))
                    eps.m[is.nan(eps.m)] <- 0# NaNが出たところには０を代入
                    epsChi <- - 2 * eps.m * menL
                    eps.Lambda <- exp(menL)
                }
                
                eps.df <- (ol * (rl - 1) * rl) / 2 - 1
                eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
                
            }
            
            # イプシロンを計算する
            LB.ep <- 1 / (rl - 1)
            GG.ep <- sum(diag(otoM))^2 / (nrow(otoM) * sum(otoM^2))
            HF.ep <- ((cellN - ol + 1)* (rl - 1) * GG.ep - 2) / ((rl - 1) * (cellN - ol - (rl - 1) * GG.ep))# Lecoutre（1991）の修正
            
            # 被験者内要因の数によって処理を変更する
            if(repleng == 1){# 被験者内要因が１つのときは有意性判定のマークを用意するのみ
                sig.mark <- sig.sign(eps.p)
                seportM <- list(orthoM)
            }else{# 被験者内要因が複数あるときは各要因ごとに，検定統計量ととイプシロンを計算
                # ラベルの追加
                effect.name <- c("Global", effect.name)
                
                # 正規直交対比行列を被験者内要因の水準数によって分割
                divpoint <- rbind(cumsum(matdivider) - (matdivider - 1), cumsum(matdivider))
                seportM <- unlist(apply(divpoint, 2, function(x) list(orthoM[x[1]:x[2], , drop = FALSE])), recursive = FALSE)
                sepM <- lapply(seportM, function(x) x %*% tm %*% t(x))
                
                if(mau == TRUE){
                    # Mauchlyの球面性検定
                    sep.Lambda <- sapply(sepM, function(x) det(x) / (sum(diag(x)) / nrow(x))^nrow(x))
                    sep.m <- 1 - (2 * matdivider^2 + matdivider + 2) / (6 * (cellN - ol) * matdivider)
                    sepChi <- -(cellN - ol) * sep.m * log(sep.Lambda)
                    sep.df <- ((matdivider * (matdivider + 1)) / 2) - 1
                    sep.p <- pchisq(sepChi, ifelse(sep.df == 0, NA, sep.df), lower.tail = FALSE)
                }else if(har == TRUE){
                    # Harrisの多標本球面性検定
                    sepA <- lapply(seportM, function(x) array(apply(covbuffer, 3, function(y) x %*% y %*% t(x)), dim = c(nrow(x), nrow(x), ol)))
                    
                    sepTr <- lapply(sepA, function(y) apply(y, 3, function(x) sum(diag(x))))
                    sepSq <- lapply(sepA, function(y) apply(y, 3, function(x) sum(diag(x %*% x))))
                    sep.Lambda <- sapply(1:ncol(divpoint), function(x) sum((othN - 1) * sepTr[[x]])^2 / sum((othN - 1) * sepSq[[x]]))
                    sepChi <- pmax(0, ((cellN - ol) * matdivider / 2) * ((cellN - ol) * matdivider / sep.Lambda - 1))
                    sep.df <- ((ol * matdivider * (matdivider + 1)) / 2) - 1
                    sep.p <- pchisq(sepChi, ifelse(sep.df == 0, NA, sep.df), lower.tail = FALSE)
                }else{
                    # Mendozaの多標本球面性検定
                    sepA <- lapply(seportM, function(x) array(apply(ptm, 3, function(y) x %*% y %*% t(x)), dim = c(nrow(x), nrow(x), ol)))
                    
                    sepL1 <- log((cellN-ol)^((cellN-ol)*(matdivider)/2) / sapply(matdivider, function(x) prod((othN-1)^((othN-1) * x / 2))))
                    sepL2 <- sapply(sepA, function(y) sum(log(apply(y, 3, function(x) det(x))) * (othN-1)/2))
                    sepL3 <- log(sapply(sepA, function(y) sum(diag(apply(y, c(1, 2), function(x) sum(x))/nrow(y))))) * ((cellN-ol) * matdivider / 2)
                    sepL <- sepL1 + sepL2 - sepL3
                    sep.m <- 1 - ((((cellN-ol) * matdivider^2 * (matdivider + 1) * (2 * matdivider + 1) - (2*(cellN-ol) * matdivider^2)) *
                    sum(1/(othN-1)) - 4) / (6 * (cellN-ol) * matdivider * (ol * matdivider * (matdivider + 1) - 2)))
                    sep.m[is.nan(sep.m)] <- 0# NaNが出たところには０を代入
                    sepChi <- - 2 * sep.m * sepL
                    sep.Lambda <- exp(sepL)
                    sep.df <- ((ol * matdivider * (matdivider + 1)) / 2) - 1
                    sep.p <- pchisq(sepChi, ifelse(sep.df == 0, NA, sep.df), lower.tail = FALSE)
                }
                
                # イプシロンを計算する
                sepLB.ep <- 1 / matdivider
                sepGG.ep <- sapply(sepM, function(x) sum(diag(x))^2 / (nrow(x) * sum(x^2)))
                sepHF.ep <- ((cellN - ol + 1) * matdivider * sepGG.ep - 2) / (matdivider * (cellN - ol - matdivider * sepGG.ep))# Lecoutre（1991）の修正
                
                # 行列全体での計算結果に追加する
                eps.Lambda <- append(eps.Lambda, sep.Lambda)
                epsChi <- append(epsChi, sepChi)
                eps.df <- append(eps.df, sep.df)
                eps.p <- append(eps.p, sep.p)
                sig.mark <- sig.sign(eps.p)
                LB.ep <- append(LB.ep, sepLB.ep)
                GG.ep <- append(GG.ep, sepGG.ep)
                HF.ep <- append(HF.ep, sepHF.ep)
            }
            
            # 結果をデータフレームにまとめる
            epsitab <- data.frame("Effect" = effect.name, "Dummy" = eps.Lambda, "approx.Chi" = epsChi, "df" = eps.df,
            "p" = eps.p, "sig.mark" = sig.mark, "LB" = LB.ep, "GG" = GG.ep, "HF" = HF.ep)
            names(epsitab)[2] <- lamlab# ラベルを検定方法に応じたものに変更
            
            # オプション；IGAのための統計量を計算
            if(iga == TRUE | ciga == TRUE){
                # HuynhのImproved General Approximate Test
                proSj <- lapply(seportM, function(y) array(apply(covbuffer, 3, function(x) y %*% x %*% t(y)), dim = c(nrow(y), nrow(y), ol)))
                proSh <- lapply(seportM, function(y) y %*% (apply(covbuffer, c(1, 2), function(x) sum((1/othN) * x)) / sum(1/othN)) %*% t(y))
                
                trDt <- sapply(proSh, function(x) sum(diag(x)))
                trDj <- lapply(proSj, function(y) apply(y, 3, function(x) sum(diag(x))))
                trDj2 <- lapply(proSj, function(y) apply(y, 3, function(x) sum(diag(x %*% x))))
                iga.D <- lapply(seportM, function(x) t(x) %*% x)
                
                iga.h0b <- trDt^2 / sapply(proSh, function(x) sum(diag(x %*% x)))
                iga.b <- ((cellN - ol) * trDt) / sapply(trDj, function(x) sum((othN - 1) * x))
                
                iga.cue <- array(sapply(1:ol, function(x) ifelse(1:(ol^2) == ol * (x-1) + x, 1, 0)), dim = c(ol, ol, ol))
                iga.Sigstr <- array(rowSums(sapply(1:ol, function(x) iga.cue[,,x] %x% (covbuffer[,,x]/othN[x]))), dim = c(ol*rl, ol*rl))
                iga.g1 <- lapply(iga.D, function(y) array(rowSums(sapply(1:ol, function(x) iga.cue[,,x] %x% (othN[x] * (1-othN[x]/cellN) * y))), dim = c(ol*rl, ol*rl)))
                if(ol == 1){
                    iga.h0c <- rep(1, length(matdivider))
                    iga.c <- rep(1, length(matdivider))
                }else{
                    iga.g2 <- lapply(iga.D, function(z) apply(combn(ol, 2, function(x) array(ifelse(1:(ol^2) == prod(x) |
                    1:(ol^2) == prod(x) + (ol - 1) * abs(diff(x)), 1, 0), dim = c(ol, ol)) %x% (-prod(othN[x]) * z / cellN)), c(1, 2), function(y) sum(y)))
                    iga.G <- mapply(function(x, y) x + y, iga.g1, iga.g2, SIMPLIFY = FALSE)
                    iga.GS <- lapply(iga.G, function(x) x %*% iga.Sigstr)
                    iga.h0c <- sapply(iga.GS, function(x) sum(diag(x))^2) / sapply(iga.GS, function(x) sum(diag(x %*% x)))
                    iga.c <- sapply(iga.GS, function(x) ((cellN - ol) * sum(diag(x)))) / sapply(trDj, function(x) ((ol - 1) * sum((othN - 1) * x)))
                }
                
                iga.h1 <- ifelse(iga.h0b == 1, 1, (cellN * iga.h0b - 2) / (cellN - ol - iga.h0b))
                iga.h2 <- ifelse(iga.h0c == 1, 1, ((ol - 1) * (cellN * iga.h0c - 2 * (ol - 1))) / ((cellN - ol) * (ol - 1) - iga.h0c))
                if(ol == 1){
                    iga.eta <- mapply(function(y, z) sum((othN - 1)^3 / ((othN + 1) * (othN - 2)) * (othN * y^2 - 2 * z)), trDj, trDj2)
                }else{
                    iga.eta <- mapply(function(y, z) sum((othN - 1)^3 / ((othN + 1) * (othN - 2)) * (othN * y^2 - 2 * z)) + 2 * sum(combn(ol, 2, function(x) prod((othN[x] - 1) * y[x]))), trDj, trDj2)
                }
                iga.sigma <- mapply(function(x, y) sum((othN - 1)^2 / ((othN + 1) * (othN - 2)) * ((othN - 1) * y - x^2)), trDj, trDj2)
                iga.e <- iga.eta / iga.sigma
                
                # Algina-LecoutreのCorrected Improved General Approximation Testのための指標
                iga.al1 <- ifelse(iga.h0b == 1, 1, ((cellN - ol + 1) * iga.h0b - 2) / (cellN - ol - iga.h0b))
                iga.al2 <- ifelse(iga.h0c == 1, 1, ((ol - 1) * ((cellN - ol + 1) * iga.h0c - 2 * (ol - 1))) / ((cellN - ol) * (ol - 1) - iga.h0c))
                
                # 結果をデータフレームにまとめる
                if(length(effect.name) == 1) iga.name <- effect.name
                else iga.name <- effect.name[-1]
                
                if(iga == TRUE){
                    epsi.info1 <- sub("Epsilons", "Estimates for IGA", epsi.info1)
                    if(repleng == 1) epsitab <- cbind(epsitab[,1:6], "b_hat" = iga.b, "c_hat" = iga.c, "h_d" = iga.h1, "h_dd" = iga.h2, "h" = iga.e)
                    else epsitab <- cbind(epsitab[,1:6], "b_hat" = c(NA, iga.b), "c_hat" = c(NA, iga.c), "h_d" = c(NA, iga.h1), "h_dd" = c(NA, iga.h2), "h" = c(NA, iga.e))
                }else{
                    epsi.info1 <- sub("Epsilons", "Estimates for CIGA", epsi.info1)
                    if(repleng == 1) epsitab <- cbind(epsitab[,1:6], "b_hat" = iga.b, "c_hat" = iga.c, "h_d" = iga.al1, "h_dd" = iga.al2, "h" = iga.e)
                    else epsitab <- cbind(epsitab[,1:6], "b_hat" = c(NA, iga.b), "c_hat" = c(NA, iga.c), "h_d" = c(NA, iga.al1), "h_dd" = c(NA, iga.al2), "h" = c(NA, iga.e))
                }
            }
            return(list("epsi.info1" = epsi.info1, "epsitab" = epsitab))
        }
        
        
        # 分散分析表を作る関数
        anova.modeler <- function(dat, design, factnames = NA, type2 = FALSE, lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE,
        mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA,
        omega = FALSE, omegana = FALSE, pomega = FALSE, gomega = NA, gomegana = NA, prep = FALSE, inter = NA){
            bet.with <- strsplit(design, "s")
            betN <- nchar(bet.with[[1]])[1]# 被験者間要因がないときは“０”
            withN <- nchar(bet.with[[1]])[2]# 被験者内要因がないときは“NA”
            maxfact <- nchar(design)-1
            
            # factnamesが省略されているときは名前をつける（デバッグ用）
            if(sum(is.na(factnames)) != 0) factnames <- LETTERS[1:maxfact]
            
            # 要因計画の型に合わせてss.calc関数を適用し分散分析表を作る
            full.elem <- unlist(lapply(1:(maxfact+1), function(y) combn(c("s", LETTERS[1:maxfact]), y, function(x) gsub(", ", ":", toString(x)))))
            fact.elem <- unlist(lapply(1:(maxfact+1), function(y) combn(c("s", factnames), y, function(x) gsub(", ", " x ", toString(x)))))
            epsi.effect <- c("Global", unlist(sapply((betN+1):maxfact, function(y) combn(factnames[(betN+1):maxfact], y - betN, function(x) gsub(", ", " x ", toString(x))))))
            
            # 効果の自由度
            cellN <- length(unique(dat$s))
            flev <- sapply(2:nchar(design), function(x) length(unique(dat[,x])))# 各要因の水準数
            eff.df <- unlist(sapply(1:maxfact, function(y) combn(maxfact, y, function(x) prod(flev[x] - 1))))# 効果の自由度
            
            # 要因計画のタイプ別の処理
            if(betN != 0 && is.na(withN)){# 被験者間計画の場合
                # full.elemの整理
                er.cue <- gsub(", ", ":", toString(c("s", strsplit(bet.with[[1]][1], "")[[1]])))
                cut.er <- setdiff(grep("s", full.elem), grep(er.cue, full.elem))
                full.elem <- full.elem[-c(cut.er, length(full.elem))]
                fact.elem <- fact.elem[-c(cut.er, length(fact.elem))]
                
                # epsilon.calcの出力に対応する情報
                epsi.info1 <- NA
                epsitab <- NA
                ano.info1 <- NA
                
                # 誤差項の自由度
                er.df <- cellN - prod(flev[1:betN])
                
            }else if(betN != 0){# 混合要因計画の場合
                # full.elemの整理
                er.cue <- gsub(", ", ":", toString(c("s", strsplit(bet.with[[1]][1], "")[[1]])))
                cut.er <- setdiff(grep("s", full.elem), grep(er.cue, full.elem))
                full.elem <- full.elem[-cut.er ]
                fact.elem <- fact.elem[-cut.er]
                
                # epsilon.calcの適用
                epsiresults <- epsilon.calc(dat = dat, design = design, mau = mau, har = har, iga = iga, ciga = ciga)
                epsi.info1 <- epsiresults$epsi.info1
                epsitab <- epsiresults$epsitab
                
                # 誤差項の自由度
                er.pri <- cellN - prod(flev[1:nchar(bet.with[[1]][1])])
                er.df <- er.pri * c(1, unlist(sapply(1:nchar(bet.with[[1]][2]), function(y)
                combn(nchar(bet.with[[1]][2]), y, function(x) prod(flev[-(1:nchar(bet.with[[1]][1]))][x]-1)))))
                
            }else{# 被験者内計画の場合
                # epsilon.calcの適用
                epsiresults <- epsilon.calc(dat = dat, design = design, mau = mau, har = har, iga = iga, ciga = ciga)
                epsi.info1 <- epsiresults$epsi.info1
                epsitab <- epsiresults$epsitab
                
                # 誤差項の自由度
                er.df <- (cellN - 1) * c(1, eff.df)
            }
            
            # ss.calcの適用
            ss.results <- ss.calc(full.elem = full.elem, dat = dat, type2 = type2)# ss.calc関数にモデルの要素を投入し，平方和を得る
            
            # 並べ替えのための指標と自由度の計算
            ieff.elem <- full.elem[-grep("s", full.elem)]# 効果の項のみ
            ier.elem <- grep("s", full.elem, value = TRUE)# 誤差項のみ
            eff.elem <- fact.elem[-grep("s", full.elem)]# 効果の項のみ
            er.elem <- fact.elem[grep("s", full.elem)]# 誤差項のみ
            
            if(length(er.elem) == 0){# 被験者間計画の場合
                internal.cue <- c(full.elem, "ss.Er", "ss.T")
                source.cue <- c(fact.elem, "ss.Er", "ss.T")
                df.col <- c(eff.df, er.df, length(dat$y)-1)
                mse.row <- length(source.cue) - 1
            }else{# 被験者内要因を含む計画の場合
                source.ord <- sapply(ieff.elem, function(x) min(elematch(x, ier.elem)))
                internal.cue <- c(unlist(sapply(1:length(ier.elem), function(x) c(ieff.elem[source.ord == x], ier.elem[x]))), "ss.T")
                source.cue <- c(unlist(sapply(1:length(er.elem), function(x) c(eff.elem[source.ord == x], er.elem[x]))), "ss.T")
                
                df.col <- unlist(sapply(1:length(er.elem), function(x) c(eff.df[source.ord == x], er.df[x])))
                df.col <- c(df.col, length(dat$y)-1)
                mse.row <- match(er.elem, source.cue)
                
                # 各効果の自由度を調整するための係数をベクトル化する
                if(iga == TRUE){
                    mdf <- 1
                    ano.info1 <- "== Huynh's Improved General Approximation Test =="
                }else if(ciga == TRUE){
                    mdf <- 1
                    ano.info1 <- "== Algina-Lecoutre's Corrected Improved General Approximation Test =="
                }else if(lb == TRUE){
                    mdf <- pmin(1, rep(c(1, epsitab$LB[epsitab$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
                    ano.info1 <- "== Geisser-Greenhouse's Conservative Test =="
                }else if(gg == TRUE){
                    mdf <- pmin(1, rep(c(1, epsitab$GG[epsitab$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
                    ano.info1 <- "== Adjusted by Greenhouse-Geisser's Epsilon =="
                }else if(hf == TRUE){
                    mdf <- pmin(1, rep(c(1, epsitab$HF[epsitab$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
                    ano.info1 <- "== Adjusted by Huynh-Feldt's Epsilon =="
                }else if(auto == TRUE){
                    sigepsi <- epsitab
                    sigepsi$GG[((sigepsi$sig.mark == "") | (sigepsi$sig.mark == "ns"))] <- 1
                    mdf <- pmin(1, rep(c(1, sigepsi$GG[sigepsi$Effect != "Global"]), c(mse.row[1], diff(mse.row))))
                    ano.info1 <- "== Adjusted by Greenhouse-Geisser's Epsilon for Suggested Violation =="
                }else{
                    mdf <- 1
                    ano.info1 <- NA
                }
                
                # 自由度を調整する
                df.col <- c(mdf, 1) * df.col
            }
            
            f.denomi <- rep(mse.row, c(mse.row[1], diff(mse.row)))
            f.denomi[mse.row] <- NA
            f.denomi <- append(f.denomi, NA)
            
            # 分散分析表を作る
            internal.lab <- gsub("ss.Er", "Error", gsub("ss.T", "Total", internal.cue))
            source.col <- gsub("ss.Er", "Error", gsub("ss.T", "Total", source.cue))
            ss.col <- unlist(sapply(internal.cue, function(x) ss.results[names(ss.results) == x]))
            attributes(ss.col) <- NULL# もとの変数名を消去
            ms.col <- ss.col / df.col# MSを計算する
            f.col <- ms.col[1:length(ms.col)] / ms.col[f.denomi]# F値を計算する
            
            # IGA，CIGAの適用
            if((iga == TRUE && !is.na(epsi.info1[1])) | (ciga == TRUE && !is.na(epsi.info1[1]))){
                mfv <- c(rep(1, mse.row[1]), rep(epsitab$c_hat[epsitab$Effect != "Global"], c(diff(mse.row))))
                mfv[mse.row[-length(mse.row)] + 1] <- epsitab$b_hat[epsitab$Effect != "Global"]
                mfv <- c(mfv, 1)
                
                iga.df <- c(df.col[1:mse.row[1]], rep(epsitab$h_dd[epsitab$Effect != "Global"], c(diff(mse.row))))
                iga.df[mse.row[-length(mse.row)] + 1] <- epsitab$h_d[epsitab$Effect != "Global"]
                iga.df[mse.row[-1]] <- epsitab$h[epsitab$Effect != "Global"]
                iga.df <- append(iga.df, length(dat$y)-1)
                
                df.col <- ifelse(iga.df >= df.col, df.col, iga.df)# 推定した自由度がもとの自由度よりも高くなった場合は調整を行わない
            }else{
                mfv <- 1# 適用しないときは１（調整なし）
            }
            
            f.col <- f.col / mfv# IGA，CIGAのための推定値によってF値を調整
            p.col <- pf(f.col, df.col, df.col[f.denomi], lower.tail = FALSE)# p値を算出する
            sig.col <- sig.sign(p.col)# p値が有意かどうかを判定して記号を表示する
            anovatab <- data.frame(source.col, ss.col, df.col, ms.col, f.col, p.col, sig.col)# 分散分析表をまとめたデータフレーム
            
            # 効果量の計算と追加
            if(eta == TRUE){# イータ二乗
                eta.col <- ss.col / sum(ss.col[-length(ss.col)])
                eta.col[is.na(f.denomi)] <- NA
                anovatab <- cbind(anovatab, "eta^2" = eta.col)
            }
            if(peta == TRUE){# 偏イータ二乗
                peta.col <- ss.col / (ss.col + ss.col[f.denomi])
                anovatab <- cbind(anovatab, "p.eta^2" = peta.col)
            }
            if(sum(is.na(geta)) == 0){# 一般化イータ二乗
                if(is.logical(geta) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                    source2int <- sapply(geta, function(x) LETTERS[match(x, factnames)])
                    measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                    ss.meas <- sum(ss.col[setdiff(measfact, mse.row)])# 誤差平方和を除く
                    measvec <- rep(1, length(ss.col))
                    measvec[measfact] <- 0# measfactに含まれる効果は分母に加える必要がないので０を代入
                    geta.col <- ss.col / (measvec * ss.col + ss.meas + sum(ss.col[mse.row]))
                }else{
                    geta.col <- ss.col / (ss.col + sum(ss.col[mse.row]))
                }
                geta.col[c(mse.row, length(geta.col))] <- NA
                anovatab <- cbind(anovatab, "G.eta^2" = geta.col)
            }
            if(eps == TRUE){# イプシロン二乗；値が負になったときは０に直す
                eps.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) /  sum(ss.col[-length(ss.col)]))
                eps.col[is.na(f.denomi)] <- NA
                anovatab <- cbind(anovatab, "epsilon^2" = eps.col)
            }
            if(peps == TRUE){# 偏イプシロン二乗；値が負になったときは０に直す
                peps.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (ss.col + ss.col[f.denomi]))
                anovatab <- cbind(anovatab, "p.epsilon^2" = peps.col)
            }
            if(sum(is.na(geps)) == 0){# 一般化イプシロン二乗；値が負になったときは０に直す
                if(is.logical(geps) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                    source2int <- sapply(geps, function(x) LETTERS[match(x, factnames)])
                    measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                    ss.meas <- sum(ss.col[setdiff(measfact, mse.row)])# 誤差平方和を除く
                    measvec <- rep(1, length(ss.col))
                    measvec[measfact] <- 0# measfactに含まれる効果は分母に加える必要がないので０を代入
                    geps.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (measvec * ss.col + ss.meas + sum(ss.col[mse.row])))
                }else{
                    geps.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (ss.col + sum(ss.col[mse.row])))
                }
                geps.col[c(mse.row, length(geps.col))] <- NA
                anovatab <- cbind(anovatab, "G.epsilon^2" = geps.col)
            }
            if(omega == TRUE){# オメガ二乗（加算モデル）；値が負になったときは０に直す；Dodd & Schultz（1973）の加算モデルの場合の計算式に基づく
                omega.denomi <- sum((ss.col - df.col * ms.col[f.denomi])[!is.na(f.denomi)]) + nrow(dat) * (sum(ss.col[mse.row]) / sum(df.col[mse.row]))
                omega.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / omega.denomi)
                anovatab <- cbind(anovatab, "omega^2" = omega.col)
            }
            if(omegana == TRUE){# オメガ二乗（非加算モデル）；値が負になったときは０に直す；Dodd & Schultz（1973）の非加算モデルの場合の計算式に基づく
                if(length(mse.row) == 1){# 被験者間計画の場合
                    omega.dummy <- cellN
                }else{# その他の計画の場合
                    dflev <- flev[(nchar(bet.with[[1]][1]) + 1):length(flev)]
                    omega.dummy <- cellN * c(1, unlist(sapply(1:length(dflev), function(y) combn(1:length(dflev), y, function(x) prod(dflev[x])))))
                }
                omegana.denomi <- sum((ss.col - df.col * ms.col[f.denomi])[!is.na(f.denomi)]) + sum(omega.dummy * ms.col[mse.row])
                omegana.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / omegana.denomi)
                anovatab <- cbind(anovatab, "omega^2_NA" = omegana.col)
            }
            if(pomega == TRUE){# 偏オメガ二乗；値が負になったときは０に直す
                pomega.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (ss.col - df.col * ms.col[f.denomi] + nrow(dat) * ms.col[f.denomi]))
                anovatab <- cbind(anovatab, "p.omega^2" = pomega.col)
            }
            if(sum(is.na(gomega)) == 0){# 一般化オメガ二乗（加算モデル）；値が負になったときは０に直す
                if(is.logical(gomega) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                    source2int <- sapply(gomega, function(x) LETTERS[match(x, factnames)])
                    measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                    ms.copy <- ms.col[f.denomi]
                    ss.meas <- sum(ss.col[setdiff(measfact, mse.row)] - df.col[setdiff(measfact, mse.row)] * ms.copy[setdiff(measfact, mse.row)])
                    measvec <- rep(1, length(ss.col))
                    measvec[measfact] <- 0
                }else{
                    ss.meas <- 0
                    measvec <- rep(1, length(ss.col))
                }
                gomega.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (measvec * (ss.col - df.col * ms.col[f.denomi]) + ss.meas + nrow(dat) * (sum(ss.col[mse.row])/sum(df.col[mse.row]))))
                anovatab <- cbind(anovatab, "G.omega^2" = gomega.col)
            }
            if(sum(is.na(gomegana)) == 0){# 一般化オメガ二乗（非加算モデル）；値が負になったときは０に直す
                if(length(mse.row) == 1){# 被験者間計画の場合
                    omega.dummy <- cellN
                }else{# その他の計画の場合
                    dflev <- flev[(nchar(bet.with[[1]][1]) + 1):length(flev)]
                    omega.dummy <- cellN * c(1, unlist(sapply(1:length(dflev), function(y) combn(1:length(dflev), y, function(x) prod(dflev[x])))))
                }
                if(is.logical(gomegana) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                    source2int <- sapply(gomegana, function(x) LETTERS[match(x, factnames)])
                    measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                    ms.copy <- ms.col[f.denomi]
                    ss.meas <- sum(ss.col[setdiff(measfact, mse.row)] - df.col[setdiff(measfact, mse.row)] * ms.copy[setdiff(measfact, mse.row)])
                    measvec <- rep(1, length(ss.col))
                    measvec[measfact] <- 0
                }else{
                    ss.meas <- 0
                    measvec <- rep(1, length(ss.col))
                }
                gomega.denomi <- ss.meas + sum(omega.dummy * ms.col[mse.row])
                gomegana.col <- pmax(0, (ss.col - df.col * ms.col[f.denomi]) / (measvec * (ss.col - df.col * ms.col[f.denomi]) + gomega.denomi))
                anovatab <- cbind(anovatab, "G.omega^2_NA" = gomegana.col)
            }
            if(prep == TRUE){# p_rep（両側）
                prep.col <- pmin(0.9999, pnorm(qnorm(1 - p.col/2) / sqrt(2)))
                anovatab <- cbind(anovatab, "p_rep" = prep.col)
            }
            
            # 結果を返す；interに入力がない場合はanovatabとmse.row，入力がある場合はintertabを返す
            if(is.na(inter) == TRUE){
                if(is.na(withN) == FALSE) epsitab$Effect <- epsi.effect[(2 * sum(withN == 1)):length(epsi.effect)]# 一要因の場合のみ“Global”を省略
                return(list("epsi.info1" = epsi.info1, "epsitab" = epsitab, "ano.info1" = ano.info1, "anovatab" = anovatab, "mse.row" = mse.row, "internal.lab" = internal.lab))
            }else{
                # 単純主効果の検定用の出力を用意する
                sim.row <- charmatch(inter, internal.lab)
                intertab <- rbind(anovatab[sim.row,], anovatab[f.denomi[sim.row],])
                
                # 球面性検定の結果からinterに関する部分のみ取り出す
                if(charmatch(inter, strsplit(design, "")[[1]]) < charmatch("s", strsplit(design, "")[[1]])){
                    # 被験者間要因ならNAの行を返す
                    if(iga == TRUE | ciga == TRUE) interepsi <- rep(NA, 11)# IGA，CIGAを使ったときは列数が多い
                    else interepsi <- rep(NA, 9)
                }else{
                    interepsi <- epsitab[charmatch(inter, epsitab$Effect),]
                }
                return(list("intertab" = intertab, "interepsi" = interepsi, "internal.lab" = internal.lab))
            }
        }
        
        
        # 修正Bonferroniの方法（Holmの方法，Shafferの方法，Holland-Copenhaverの方法）による多重比較を行う関数
        # デフォルトはShafferの方法を適用し，holm = TとするとHolmの方法，hc = TとするとHolland-Copenhaverの方法を適用する
        # s2r = T，s2d = Tとすると，具体的な棄却のパターンを反映した有意水準の調整によるShafferの方法を適用する
        mod.Bon <- function(dat, design, taref, bet.mse, factlabel = NA, factnames = NA, type2 = FALSE, holm = FALSE, hc = FALSE,
        s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, alpha = 0.05, criteria = FALSE){
            # 対象となる要因のラベルを得る
            if(is.na(factlabel)) factlabel <- taref
            
            bonl <- nlevels(dat[, taref])# 水準数の取得
            h0size <- bonl * (bonl-1)/2# 帰無仮説の個数
            bon.name <- combn(bonl, 2, function(ij) paste0(levels(dat[, taref])[ij][1], "-", levels(dat[, taref])[ij][2]))# 可能な対の組み合わせのラベル
            bon.num <- charmatch(taref, names(dat))-1# 分析対象となる要因が何番目の要因かを特定
            
            # 周辺平均を計算する
            cont.means <- tapply(dat$y, dat[, 2:nchar(design)], mean)# 各セルの平均を求める
            bon.means <- apply(cont.means, bon.num, mean)# 分析対象となる周辺平均を求める
            
            factlevels <- sapply(names(dat), function(x) nlevels(dat[, x]))# 各要因の水準数
            factlevels[charmatch(taref, names(dat))] <- 2# 多重比較の対象となる効果の水準数を２に固定
            factlevels <- factlevels[!(factlevels == factlevels[1] | factlevels == factlevels[length(factlevels)])]# 最初と最後（sとyの列）を除く
            
            cont.N <- table(dat[, 2:nchar(design)])# 各セルのデータ数
            bon.denomi <- apply(1/cont.N, bon.num, mean) / (prod(factlevels)/2)# セルごとに重み付けしたデータ数の平均を残りの条件数で割ったもの
            
            bon.delta <- combn(bonl, 2, function(ij) (bon.means[ij][1] - bon.means[ij][2]))# 平均偏差
            
            # 分析対象が被験者間要因か被験者内要因かによって標準誤差を得る方法を変える
            if(length(bet.mse) != 1){
                # 被験者間要因の場合；上位の分析のMSeを適用する
                bon.df <- bet.mse$df.col# 自由度
                bon.Ve <- bet.mse$ms.col# 平均平方
            }else{
                # 被験者内要因の場合；比較する２水準ごとにMSeを再計算する
                bon.lev <- combn(levels(dat[, taref]), 2, function(x) x)# 各水準の組み合わせ
                subdat <- apply(bon.lev, 2, function(x) dat[dat[, taref] == x[1] | dat[, taref] == x[2], ])# データを多重比較の対象となる効果について２水準ずつのサブデータにリスト化
                for(i in 1:length(subdat)){
                    subdat[[i]][, taref] <- subdat[[i]][, taref][, drop = TRUE]
                }
                
                bon.anova <- lapply(subdat, function(x) anova.modeler(dat = x, design = design, factnames = factnames, type2 = type2,
                inter = taref)$intertab[2, ])
                bon.df <- sapply(bon.anova, function(x) x$df.col)# 自由度
                bon.Ve <- sapply(bon.anova, function(x) x$ms.col)
            }
            
            # 検定統計量とｐ値を得る
            bon.SE <- sqrt(combn(bonl, 2, function(ij) sum(bon.denomi[ij])) * bon.Ve)
            bon.t <- abs(bon.delta / bon.SE)# ｔ値は絶対値を取る
            bon.p <- pt(bon.t, bon.df, lower.tail = FALSE) * 2# 両側確率
            
            # 結果をデータフレームにまとめる
            bontab <- data.frame("pair" = bon.name, "difference" = bon.delta, "t" = bon.t, "df" = bon.df, "p.value" = bon.p)
            bontab <- bontab[order(bontab$p.value),]# p値の小さい順に並べ替え
            
            # 調整した有意水準を設定する
            if(holm == TRUE){
                p.criteria <- h0size:1# Holmの方法用の調整値
                bon.info1 <- paste0("== Holm's Sequentially Rejective Bonferroni Procedure ==")
                
            }else if(s2d | fs2d == TRUE){# Donoguhe（2004）のアルゴリズムをベースとするShafferの多重比較のための論理ステップの計算
                # Donoghue, J. R. (2004). Implementing Shaffer's multiple comparison procedure for a large number of groups.
                # Recent developments in multiple comparison procedures (Institute of mathematical statistics-Lecture Notes-Monograph Series, 47), pp. 1-23.
                # Donoghueと完全に同じ手順ではないことに注意
                
                # 隣接行列を作る
                bon.comb <- combn(bonl, 2, function(x) x)# 帰無仮説を表す行列
                bon.comb <- bon.comb[,order(bon.p)]# ｐ値の順に並べ替え
                hvec <- 1:bonl
                a.mat <- diag(bonl)
                diag(a.mat) <- 0
                
                shaf.value <- c(h0size, rep(NA, h0size - 1))# ステップごとの仮説数を代入するためのベクトル
                
                allcomb <- lapply((bonl-1):1, function(y) combn(bonl, y, function(x) x, simplify = FALSE))# すべての帰無仮説の組み合わせを示すリスト
                allcomb <- unlist(allcomb, recursive = FALSE)# リストの階層をなくす
                
                for(j in 1:(h0size-1)){
                    # 隣接行列に棄却された仮説を書き込む
                    a.mat[bon.comb[1, j], bon.comb[2, j]] <- 1# 棄却された帰無仮説の部分に１を代入
                    a.mat[bon.comb[2, j], bon.comb[1, j]] <- 1# 対角線を通して反対の側にも代入
                    
                    # 未分化クラスを作る
                    # 隣接行列の下位行列の中から０のみで構成される正方行列を探す
                    # 未分化クラスを表す行列：各行が各水準に相当；各列が帰無仮説を表す（互いに差がない水準に１を代入）
                    undiff <- array(rep(0, bonl), c(bonl, 1))# ダミー
                    cnt <- 1
                    while(cnt <= length(allcomb)){
                        hnum <- allcomb[[cnt]]
                        if(max(colSums(undiff[hnum, , drop = FALSE])) == length(hnum)){
                            # 上位の仮説に包含される仮説は含めない；カットはしない
                            cnt <- cnt + 1
                        }else if(sum(a.mat[hnum, hnum]) == 0){
                            # 正方行列の場合は成立する帰無仮説を表す列を追加
                            undiff <- cbind(undiff, 1 - 0^match(hvec, hnum, nomatch = 0))
                            cnt <- cnt + 1
                        }else{
                            # その他の場合；このパターンは後に支持されることはないので，allcombからカット
                            allcomb <- allcomb[-cnt]
                        }
                    }
                    
                    undiff <- undiff[,-1]# 一列目のダミーを除く
                    gsize <- colSums(undiff)# 各グループの要素数を示すベクトル
                    
                    # sig.minを決定する
                    sig.min <- max(gsize)^2# 最大クラスの要素数を二乗した値
                    nxcand <- undiff# 未分化クラスのコピー
                    gi <- 1
                    while(ncol(nxcand) > 1 && nrow(nxcand) > 1){
                        nxcand <- nxcand[(1:nrow(nxcand))[nxcand[, gi] == 0], , drop = FALSE]
                        lengvec <- colSums(nxcand)# 各クラスの要素数
                        sig.min <- sig.min + max(lengvec)^2# 最大クラスの要素数の二乗値を足す
                        gi <- which.max(lengvec)# 最大クラスの番号
                    }
                    
                    # don.maxを決定する
                    don.smax <- sig.min
                    
                    for(i in 2:min(ncol(undiff)-1, bonl)){
                        don.sig <- gsize[i]^2
                        nxcand <- undiff
                        gi <- i
                        while(ncol(nxcand) > 1 && nrow(nxcand) > 1){
                            nxcand <- nxcand[(1:nrow(nxcand))[nxcand[,gi] == 0], , drop = FALSE]
                            lengvec <- colSums(nxcand)# 各クラスの要素数
                            don.sig <- don.sig + max(lengvec)^2# 最大クラスの要素数の二乗値を足す
                            gi <- which.max(lengvec)# 最大クラスの番号
                        }
                        don.smax <- max(don.sig, don.smax)# より大きい値を残す
                    }
                    shaf.value[j+1] <- (don.smax - bonl) / 2
                }
                
                if(s2d == T){
                    p.criteria <- shaf.value# Shafferの方法用の調整値
                    bon.info1 <- paste0("== Shaffer's Modified Sequentially Rejective Bonferroni Procedure [SPECIFIC] ==", "\n",
                    "== This computation is based on the algorithm by Donoghue (2004). ==")
                    shaf.meth <- paste0(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Donoghue (2004). ==")
                }else{
                    p.criteria <- c(shaf.value[2], shaf.value[2:length(shaf.value)])# F-Shafferの方法用の調整値
                    bon.info1 <- paste0("== Shaffer's F-Modified Sequentially Rejective Bonferroni Procedure [SPECIFIC] ==", "\n",
                    "== This computation is based on the algorithm by Donoghue (2004). ==")
                    shaf.meth <- paste0(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Donoghue (2004). ==")
                }
                
            }else{# Rasmussen（1993）のアルゴリズムによるShafferの多重比較のための論理ステップの計算
                # Rasmussen, J. L. (1993). Algorithm for Shaffer's multiple comparison tests. Educational and Psychological Measurement, 53, 329-335.
                
                # 平均間の異同パターンを表す行列を作る
                hpattern <- 2^(bonl-1)# 可能な真偽の仮説のパターン数
                nbuffer <- 2^((bonl-2):0)
                c.mat <- cbind(rep(0, hpattern), sapply(nbuffer, function(x) rep(c(rep(0, x), rep(1, x)), nbuffer[1]/x)))
                c.mat <- t(apply(c.mat, 1, function(x) cumsum(x)))
                
                f.mat <- combn(bonl, 2, function(x) c.mat[, x[1]] - c.mat[, x[2]])# 各水準の組み合わせを表現する行列
                f.mat[f.mat != 0] <- 1# 帰無仮説が真のときに０，偽のときに１となるようにする
                rebon.p <- bon.p[order(combn(rank(bon.means), 2, function(x) prod(x)))]# ｐ値の順序を平均値の大きさにそって並べ替え
                f.mat <- f.mat[, order(rebon.p)]# ｐ値の小さい順に列を並べ替え
                i.vector <- rowSums(f.mat)# 棄却される帰無仮説の数
                t.vector <- h0size - i.vector# 成立しうる真の帰無仮説の数
                
                if(s2r | fs2r == TRUE){# 各比較までの特定の仮説が偽であったときの可能な真の帰無仮説の最大数
                    shaf.value <- c(max(t.vector), max(t.vector[i.vector >= (2 - 1)][(f.mat[i.vector >= (2 - 1), 1:(2-1)]) == (2 - 1)]))
                    shaf.value <- c(shaf.value, sapply(3:h0size, function(x) max(t.vector[i.vector >= (x - 1)][rowSums(f.mat[i.vector >= (x - 1), 1:(x-1)]) == (x - 1)])))
                    shaf.meth <- paste0(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Rasmussen (1993). ==")
                }else{# 各比較までの任意の仮説が偽であったときの可能な真の帰無仮説の最大数
                    shaf.value <- sapply(1:h0size, function(x) max(t.vector[i.vector >= (x - 1)]))
                    shaf.meth <- " =="
                }
                
                if(fs1 | fs2r == TRUE){
                    p.criteria <- c(shaf.value[2], shaf.value[2:length(shaf.value)])# F-Shafferの方法用の調整値
                    bon.info1 <- paste0("== Shaffer's F-Modified Sequentially Rejective Bonferroni Procedure", shaf.meth)
                }else{
                    p.criteria <- shaf.value# Shafferの方法用の調整値
                    bon.info1 <- paste0("== Shaffer's Modified Sequentially Rejective Bonferroni Procedure", shaf.meth)
                }
            }
            
            # 平均値の差の方向を調べ，不等号のベクトルを作る
            bon.differ <- ifelse(bontab$difference <= 0, sub("-", " < ", bontab$pair), sub("-", " > ", bontab$pair))
            # 差が見られなかった場合の等号のベクトルを作る
            bon.equal <- sub("-", " = ", bontab$"pair")
            
            if(criteria == TRUE){# データフレームに調整済み有意水準の列を加える
                if(hc == TRUE){
                    bontab <- transform(bontab, "criteria" = 1 - (1 - alpha) ^ (1/p.criteria))# Sidakの不等式による有意水準の調整
                    if(holm == TRUE) bon.info1 <- paste0("== Holm's Sequentially Rejective Sidak Procedure ==")
                    else bon.info1 <- paste0("== Holland-Copenhaver's Improved Sequentially Rejective Sidak Procedure", shaf.meth)
                    
                    if(length(bet.mse) == 1) bon.info1 <- append(bon.info1, "*** CAUTION! This procedure might be inappropriate for dependent means. ***")
                }else{
                    bontab <- transform(bontab, "criteria" = alpha/p.criteria)# Bonferroniの不等式による有意水準の調整
                }
                # 有意であった行は不等号，そうでない行は等号を表示する
                bon.sign <- ifelse(cummin(bontab$p.value < bontab$criteria), paste(bon.differ, "*", sep = " "), paste(bon.equal, " ", sep = " "))
            }else{# データフレームに調整済みｐ値の列を加える
                if(hc == TRUE){
                    bontab <- transform(bontab, "adj.p" = pmin(1, cummax((1-(1-bontab$p.value)^p.criteria))))# Sidakの不等式による調整済みｐ値
                    if(holm == TRUE) bon.info1 <- paste0("== Holm's Sequentially Rejective Sidak Procedure ==")
                    else bon.info1 <- paste0("== Holland-Copenhaver's Improved Sequentially Rejective Sidak Procedure", shaf.meth)
                    
                    if(length(bet.mse) == 1) bon.info1 <- append(bon.info1, "*** CAUTION! This procedure might be inappropriate for dependent means. ***")
                }else{
                    bontab <- transform(bontab, "adj.p" = pmin(1, cummax(bontab$p.value * p.criteria)))# Bonferroniの不等式による調整済みｐ値
                }
                # 有意であった行は不等号，そうでない行は等号を表示する
                bon.sign <- ifelse(bontab$adj.p < alpha, paste(bon.differ, "*", sep = " "), paste(bon.equal, " ", sep = " "))
            }
            
            # 判定結果をデータフレームに反映する
            bontab <- transform(bontab, "significance" = bon.sign)
            
            # 記述統計量の計算
            b.sncol <- tapply(dat$y, dat[,taref], length)# セルごとのデータ数を計算
            b.sdcol <- tapply(dat$y, dat[,taref], sd)# セルごとの標準偏差を計算
            
            bonstat <- data.frame("Dummy" = levels(dat[, taref]), "n" = b.sncol, "Mean" = bon.means, "S.D." = b.sdcol)
            
            names(bonstat)[1] <- factlabel# 水準を表すラベルをfactlabelとして入力した値に置き換える
            
            # その他の出力を準備する
            bon.info2 <- if(length(bet.mse) != 1) paste0("== The factor < ", factlabel, " > is analysed as independent means. ==")
            else paste0("== The factor < ", factlabel, " > is analysed as dependent means. ==")
            bon.info3 <- paste0("== Alpha level is ", alpha, ". ==")
            
            bonresults <- list(factlabel, bon.info1, bon.info2, bon.info3, bonstat, bontab)
            names(bonresults) <- c(taref, "bon.info1", "bon.info2", "bon.info3", "bonstat", "bontab")
            return(bonresults)
        }
        
        
        # 下位検定を行う関数
        post.analyses <- function(dat, design, factnames = NA, mainresults, type2 = FALSE, nopost = FALSE, holm = FALSE, hc = FALSE,
        s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE,
        lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE,
        eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE,
        gomega = NA, gomegana = NA, prep = FALSE){
            anovatab <- mainresults$anovatab
            internal.lab <- mainresults$internal.lab
            
            # 要因計画の型から被験者間要因と被験者内要因の情報を得る
            bet.with <- strsplit(design, "s")
            # 効果が有意であった行のソースラベルを得る
            sig.source <- internal.lab[!((anovatab$"sig.col" == "") | (anovatab$"sig.col" == "ns"))]
            
            if(length(sig.source) == 0 | nopost == TRUE) {
                return(NA)# 有意な行が存在しないか，nopostオプションが指定されている場合はここで終了
            }else{
                # 下位検定の結果を格納するための空のデータフレームを宣言
                postresults <- lapply(paste0("post", 1:length(sig.source)), function(x) x)
                
                # pro.fraction関数を反復適用
                for (i in 1:length(sig.source)) postresults[[i]] <- pro.fraction(dat = dat, design = design, factnames = factnames,
                postplan = sig.source[i], bet.with = bet.with, mainresults = mainresults, type2 = type2, holm = holm, hc = hc,
                s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria,
                lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, eta = eta, peta = peta, geta = geta,
                eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega, gomega = gomega, gomegana = gomegana,
                prep = prep)
                return(postresults)
            }
        }
        
        
        # 効果のタイプに適した下位検定を割り当てる関数
        pro.fraction <- function(dat, design, factnames = NA, postplan, bet.with, mainresults, type2 = FALSE, holm = FALSE, hc = FALSE,
        s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE,
        lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE,
        eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE,
        gomega = NA, gomegana = NA, prep = FALSE){
            # 情報の展開
            anovatab <- mainresults$anovatab
            mse.row <- mainresults$mse.row
            internal.lab <- mainresults$internal.lab
            
            # 内部ラベルを文字列に変換し，その文字数を得る
            sig.term <- as.character(postplan)
            sig.num <- nchar(sig.term)
            sig.source <- as.character(anovatab$source.col[match(sig.term, internal.lab)])
            
            # 効果の種類によって違った処理を割り当てる
            if(sig.num > 3){# 高次の交互作用：効果のラベルを返す
                highx <- sig.source
                names(highx) <- sig.term
                return(highx)
            }else if(sig.num == 3){# １次の交互作用については，単純主効果の検討を行う
                each.term <- strsplit(sig.term, ":")# 交互作用の各項を分離する
                interfact <- as.numeric(sapply(each.term[[1]], function(x) grep(x, bet.with[[1]][1])))# 被験者間要因なら１，被験者内要因なら０を返す
                
                # 交互作用を構成する第一項の項の列番号，水準数を取得
                col.num1 <- charmatch(each.term[[1]][1], names(dat))# 列番号
                level.num1 <- nlevels(dat[, col.num1])# 水準数
                each.fact1 <- factnames[grep(each.term[[1]][1], LETTERS)]# 要因ラベル
                
                # 交互作用を構成する第二項の項の列番号，水準数を取得
                col.num2 <- charmatch(each.term[[1]][2], names(dat))# 列番号
                level.num2 <- nlevels(dat[, col.num2])# 水準数
                each.fact2 <- factnames[grep(each.term[[1]][2], LETTERS)]# 要因ラベル
                
                # 単純主効果の検定用の要因計画を作る
                resdesign1 <- sub(each.term[[1]][2], "", design)# つぶす要因のラベルを消去
                bet.with1 <- strsplit(resdesign1, "s")# 残りの要因を被験者間と被験者内に分離
                bet.num1 <- if(nchar(bet.with1[[1]][1]) == 0) 0 else 1:nchar(bet.with1[[1]][1])# 被験者間要因の数を数える
                with.num1 <- if(is.na(bet.with1[[1]][2])) 0 else (max(bet.num1)+1):(nchar(resdesign1)-1)# 被験者内要因の数を数える
                bet1 <- gsub(", ", "", toString(LETTERS[bet.num1]))# 被験者間要因のラベルを新たに作成
                with1 <- gsub(", ", "", toString(LETTERS[with.num1]))# 被験者内要因のラベルを新たに作成
                postdesign1 <- paste0(bet1, "s", with1)# 両要因を合成
                
                resdesign2 <- sub(each.term[[1]][1], "", design)# つぶす要因のラベルを消去
                bet.with2 <- strsplit(resdesign2, "s")# 残りの要因を被験者間と被験者内に分離
                bet.num2 <- if(nchar(bet.with2[[1]][1]) == 0) 0 else 1:nchar(bet.with2[[1]][1])# 被験者間要因の数を数える
                with.num2 <- if(is.na(bet.with2[[1]][2])) 0 else (max(bet.num2)+1):(nchar(resdesign2)-1)# 被験者内要因の数を数える
                bet2 <- gsub(", ", "", toString(LETTERS[bet.num2]))# 被験者間要因のラベルを新たに作成
                with2 <- gsub(", ", "", toString(LETTERS[with.num2]))# 被験者内要因のラベルを新たに作成
                postdesign2 <- paste0(bet2, "s", with2)# 両要因を合成
                
                # 分析対象となる効果が被験者間要因か被験者内要因かによってソース列のラベルを作成
                if(is.na(interfact[1]) == FALSE) err.label1 <- "Er"
                else err.label1 <- paste0("s x ", each.fact1)
                
                if(is.na(interfact[2]) == FALSE) err.label2 <- "Er"
                else err.label2 <- paste0("s x ", each.fact2)
                
                # データフレームを各要因の各水準ごとに分離
                subdat1 <- dat# データフレームをコピー
                rename.vector1 <- c("s", LETTERS[1:(nchar(design)-2)], "y")
                rename.vector1 <- append(rename.vector1, "X", after = col.num2-1)
                names(subdat1) <- rename.vector1# 列名を変更
                target.effect1 <- names(subdat1)[col.num1]# 検討したい効果の変更後の列名を取得
                subdat1 <- split(subdat1, subdat1[,col.num2])# 分析対象でない要因の水準ごとにデータフレームを分割
                subdat1 <- lapply(subdat1, function(x) x[-col.num2])# X列を消去
                
                subdat2 <- dat# データフレームをコピー
                rename.vector2 <- c("s", LETTERS[1:(nchar(design)-2)], "y")
                rename.vector2 <- append(rename.vector2, "X", after = col.num1-1)
                names(subdat2) <- rename.vector2# 列名を変更
                target.effect2 <- names(subdat2)[col.num2]# 検討したい効果の変更後の列名を取得
                subdat2 <- split(subdat2, subdat2[,col.num1])# 分析対象でない要因の水準ごとにデータフレームを分割
                subdat2 <- lapply(subdat2, function(x) x[-col.num1])# X列を消去
                
                # 各効果について要因を１つ落とした分散分析
                sim.effects1 <- rep(list(NA), col.num2)# 結果格納用リスト
                sim.effects2 <- rep(list(NA), col.num1)# 結果格納用リスト
                
                sim.effects1 <- lapply(subdat1, function(x) anova.modeler(dat = x, design = postdesign1, factnames = factnames[-(col.num2-1)], type2 = type2,
                lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
                eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
                gomega = gomega, gomegana = gomegana, prep = prep, inter = target.effect1))
                sim.effects2 <- lapply(subdat2, function(x) anova.modeler(dat = x, design = postdesign2, factnames = factnames[-(col.num1-1)], type2 = type2,
                lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
                eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
                gomega = gomega, gomegana = gomegana, prep = prep, inter = target.effect2))
                
                # 結果を１つのデータフレームにまとめる
                simtab <- data.frame()# 結果格納用データフレーム
                simepsi <- data.frame()# 結果格納用データフレーム
                sim.internal.lab <- c()# ソースラベル格納用ベクトル
                subsource <- c()# ソースラベル格納用ベクトル
                
                # 計画のタイプによって出力を整理し直す
                if(substr(design, nchar(design), nchar(design)) == "s"){
                    # 被験者間計画なら，anovatabから主分析の誤差平方和等を得る
                    between.mse <- anovatab[charmatch("Error", internal.lab),]
                    bet.mse1 <- between.mse
                    bet.mse2 <- between.mse
                    
                    # 主分析の誤差平方和を用いて単純主効果の結果を再計算
                    for (i in 1:level.num2){
                        simtab <- rbind(simtab, sim.effects1[[i]]$intertab[1,])
                        sim.internal.lab <- append(sim.internal.lab, c(paste0(each.term[[1]][1], " at ", names(sim.effects1)[i])))
                        subsource <- append(subsource, c(paste0(each.fact1, " at ", names(sim.effects1)[i])))
                    }
                    for (i in 1:level.num1){
                        simtab <- rbind(simtab, sim.effects2[[i]]$intertab[1,])
                        sim.internal.lab <- append(sim.internal.lab, c(paste0(each.term[[1]][2], " at ", names(sim.effects2)[i])))
                        subsource <- append(subsource, c(paste0(each.fact2, " at ", names(sim.effects2)[i])))
                    }
                    simtab <- rbind(simtab, between.mse)
                    sim.f.col <- simtab$ms.col[1:(nrow(simtab)-1)] / simtab$ms.col[nrow(simtab)]
                    sim.p.col <- pf(sim.f.col[1:(nrow(simtab)-1)], simtab$df.col[1:(nrow(simtab)-1)], simtab$df.col[nrow(simtab)], lower.tail = FALSE)
                    sim.sig.col <- sig.sign(sim.p.col)
                    
                    # 計算結果をデータフレームに反映
                    simtab$sig.col <- ""
                    simtab <- cbind(simtab[,1:4], "f.col" = c(sim.f.col, NA), "p.col" = c(sim.p.col, NA), "sig.col" = c(sim.sig.col, ""))
                    
                    # 効果量の再計算
                    if(eta == TRUE){# イータ二乗
                        eta.col <- simtab$ss.col[1:(nrow(simtab)-1)] / sum(anovatab$ss.col[-nrow(anovatab)])# 主分析の総平方和を使用
                        simtab <- cbind(simtab, "eta^2" = c(eta.col, NA))
                    }
                    if(peta == TRUE){# 偏イータ二乗
                        peta.col <- simtab$ss.col[1:(nrow(simtab)-1)] / (simtab$ss.col[1:(nrow(simtab)-1)]  + simtab$ss.col[nrow(simtab)])
                        simtab <- cbind(simtab, "p.eta^2" = c(peta.col, NA))
                    }
                    if(sum(is.na(geta)) == 0){# 一般化イータ二乗
                        if(is.logical(geta) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                            source2int <- sapply(geta, function(x) LETTERS[match(x, factnames)])
                            measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)])# 誤差平方和を除く
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(source2int, function(x) grep(x, sim.internal.lab))))# 個人差変数に含まれる単純主効果の取り出し
                            measvec[pmeasfact] <- 0
                            geta.col <- simtab$ss.col[1:(nrow(simtab)-1)] / (measvec * simtab$ss.col[1:(nrow(simtab)-1)] + ss.meas + simtab$ss.col[nrow(simtab)])
                        }else{
                            geta.col <- simtab$ss.col[1:(nrow(simtab)-1)] / (simtab$ss.col[1:(nrow(simtab)-1)] + simtab$ss.col[nrow(simtab)])
                        }
                        simtab <- cbind(simtab, "G.eta^2" = c(geta.col, NA))
                    }
                    if(eps == TRUE){# イプシロン二乗；値が負になったときは０に直す
                        eps.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) / sum(anovatab$ss.col[-nrow(anovatab)]))# 主分析の総平方和を使用
                        simtab <- cbind(simtab, "epsilon^2" = c(eps.col, NA))
                    }
                    if(peps == TRUE){# 偏イプシロン二乗；値が負になったときは０に直す
                        peps.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) / (simtab$ss.col[1:(nrow(simtab)-1)]  + simtab$ss.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "p.epsilon^2" = c(peps.col, NA))
                    }
                    if(sum(is.na(geps)) == 0){# 一般化イプシロン二乗；値が負になったときは０に直す
                        if(is.logical(geps) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                            source2int <- sapply(geps, function(x) LETTERS[match(x, factnames)])
                            measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)])# 誤差平方和を除く
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(source2int, function(x) grep(x, sim.internal.lab))))# 個人差変数に含まれる単純主効果の取り出し
                            measvec[pmeasfact] <- 0
                            geps.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) / (measvec * simtab$ss.col[1:(nrow(simtab)-1)] + ss.meas + simtab$ss.col[nrow(simtab)]))
                        }else{
                            geps.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) / (simtab$ss.col[1:(nrow(simtab)-1)] + simtab$ss.col[nrow(simtab)]))
                        }
                        simtab <- cbind(simtab, "G.epsilon^2" = c(geps.col, NA))
                    }
                    if(omega == TRUE){# オメガ二乗（加算モデル）；値が負になったときは０に直す
                        omega.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) /
                        (sum(anovatab$ss.col[-nrow(anovatab)]) + simtab$ms.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "omega^2" = c(omega.col, NA))
                    }
                    if(omegana == TRUE){# オメガ二乗（非加算モデル）；値が負になったときは０に直す
                        omega.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) /
                        (sum(anovatab$ss.col[-nrow(anovatab)]) + simtab$ms.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "omega^2_NA" = c(omegana.col, NA))
                    }
                    if(pomega == TRUE){# 偏オメガ二乗；値が負になったときは０に直す
                        cellN <- length(unique(dat$s))
                        pomega.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * simtab$ms.col[nrow(simtab)]) /
                        (simtab$ss.col[1:(nrow(simtab)-1)] + (cellN - simtab$df.col[1:(nrow(simtab)-1)]) * simtab$ms.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "p.omega^2" = c(pomega.col, NA))
                    }
                    if(sum(is.na(gomega)) == 0){# 一般化オメガ二乗（加算モデル）；値が負になったときは０に直す
                        ms.copy <- anovatab$ms.col[nrow(anovatab)-1]
                        if(is.logical(gomega) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                            source2int <- sapply(gomega, function(x) LETTERS[match(x, factnames)])
                            measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)] - anovatab$df.col[setdiff(measfact, mse.row)] * ms.copy)
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(source2int, function(x) grep(x, sim.internal.lab))))# 個人差変数に含まれる単純主効果の取り出し
                            measvec[pmeasfact] <- 0
                        }else{
                            ss.meas <- 0
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                        }
                        gomega.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * ms.copy) /
                        (measvec * (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * ms.copy) + ss.meas + nrow(dat) * simtab$ms.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "G.omega^2" = c(gomega.col, NA))
                    }
                    if(sum(is.na(gomegana)) == 0){# 一般化オメガ二乗（非加算モデル）；値が負になったときは０に直す
                        ms.copy <- anovatab$ms.col[nrow(anovatab)-1]
                        if(is.logical(gomegana) == FALSE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            eff.internal <- gsub("Error", NA, gsub("Total", NA, internal.lab))
                            source2int <- sapply(gomegana, function(x) LETTERS[match(x, factnames)])
                            measfact <- unique(unlist(lapply(source2int, function(x) grep(x, eff.internal))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)] - anovatab$df.col[setdiff(measfact, mse.row)] * ms.copy)
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(source2int, function(x) grep(x, sim.internal.lab))))# 個人差変数に含まれる単純主効果の取り出し
                            measvec[pmeasfact] <- 0
                        }else{
                            ss.meas <- 0
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                        }
                        gomegana.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * ms.copy) /
                        (measvec * (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * ms.copy) + ss.meas + nrow(dat) * simtab$ms.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "G.omega^2_NA" = c(gomegana.col, NA))
                    }
                    if(prep == TRUE){# p_rep（両側）
                        prep.col <- pmin(0.9999, pnorm(qnorm(1 - sim.p.col/2) / sqrt(2)))
                        simtab <- cbind(simtab, "p_rep" = c(prep.col, NA))
                    }
                    sim.internal.lab <- append(sim.internal.lab, "Error")
                    subsource <- append(subsource, "Error")
                    simepsi <- NA# 球面性検定の結果はなし
                    
                    # ソース列，行番号のラベルを張り替える
                    simtab$source.col <- subsource
                    row.names(simtab) <- c(1:nrow(simtab))
                }else{
                    # 被験者内計画，混合要因計画の場合，結果を１行ずつ取り出してrbindで結合する
                    for (i in 1:level.num2){
                        simtab <- rbind(simtab, sim.effects1[[i]]$intertab)
                        simepsi <- rbind(simepsi, sim.effects1[[i]]$interepsi)
                        sim.internal.lab <- append(sim.internal.lab, c(paste0(each.term[[1]][1], " at ", names(sim.effects1)[i]),
                        paste0(err.label1, " at ", names(sim.effects1)[i])))
                        subsource <- append(subsource, c(paste0(each.fact1, " at ", names(sim.effects1)[i]),
                        paste0(err.label1, " at ", names(sim.effects1)[i])))
                    }
                    names(simepsi) <- names(sim.effects2[[1]]$interepsi)# 名前が異なるとrbindできないので統一；混合要因計画の場合，必ず被験者間要因がsim.effects1に入ることを前提とする
                    for (i in 1:level.num1){
                        simtab <- rbind(simtab, sim.effects2[[i]]$intertab)
                        simepsi <- rbind(simepsi, sim.effects2[[i]]$interepsi)
                        sim.internal.lab <- append(sim.internal.lab, c(paste0(each.term[[1]][2], " at ", names(sim.effects2)[i]),
                        paste0(err.label2, " at ", names(sim.effects2)[i])))
                        subsource <- append(subsource, c(paste0(each.fact2, " at ", names(sim.effects2)[i]),
                        paste0(err.label2, " at ", names(sim.effects2)[i])))
                    }
                    if(match("df", names(simepsi), nomatch = 0) == 0){# 被験者間要因だけを取り出して下位検定を行った場合（データフレームの名前が""になっている）
                        simepsi <- NA
                    }else{
                        simepsi$Effect <- subsource[(1:length(subsource)) %% 2 == 1]# subsourceの奇数番の値をsimepsiの効果ラベルに貼り付ける
                        simepsi <- simepsi[!is.na(simepsi$df),]# dfがNAの行（被験者間効果の行）を除く
                        if(nrow(simepsi) == 0) simepsi <- NA# 数値のある行がなければNAを代入
                    }
                    
                    # ソース列，行番号のラベルを張り替える
                    simtab$source.col <- subsource
                    row.names(simtab) <- c(1:nrow(simtab))
                    
                    # 混合計画における被験者間要因に対して主分析の平均平方を取り出す
                    if(is.na(interfact[1]) == FALSE){
                        bet.mse1 <- simtab[1:level.num2 * 2,]# 被験者間要因が分析対象の場合
                    }else{
                        bet.mse1 <- NA# 被験者内要因が分析対象の場合
                    }
                    
                    if(is.na(interfact[2]) == FALSE){
                        bet.mse2 <- simtab[(level.num2 + 1):(level.num1 + level.num2) * 2,]# 被験者間要因が分析対象の場合
                    }else{
                        bet.mse2 <- NA# 被験者内要因が分析対象の場合
                    }
                }
                
                # 記述統計量の計算
                sim.sncol <- as.vector(tapply(dat$y, list(dat[, col.num2], dat[, col.num1]), length))# セルごとのデータ数を計算
                sim.mncol <- as.vector(tapply(dat$y, list(dat[, col.num2], dat[, col.num1]), mean))# セルごとの平均を計算
                sim.sdcol <- as.vector(tapply(dat$y, list(dat[, col.num2], dat[, col.num1]), sd))# セルごとの標準偏差を計算
                
                sim.stat <- data.frame("Term1" = rep(levels(dat[, col.num1]), each = level.num2),
                "Term2" = rep(levels(dat[, col.num2]), level.num1), "n" = sim.sncol, "Mean" = sim.mncol, "S.D." = sim.sdcol)
                
                # 行ラベルの張り替え
                names(sim.stat)[1] <- each.fact1
                names(sim.stat)[2] <- each.fact2
                
                # 有意であった行のソースをチェック
                sim.sig.col <- sim.internal.lab[!((simtab$sig.col == "") | (simtab$sig.col == "ns"))]
                sim.sig.lab <- simtab$source.col[!((simtab$sig.col == "") | (simtab$sig.col == "ns"))]
                
                # 多重比較を実行
                if(length(sim.sig.col) == 0){
                    sim.multresults <- NA# 有意な行がない場合はNAを代入
                }else{
                    # 多重比較の結果を格納するための空のデータフレームを宣言
                    sim.multresults <- sapply(paste0("mult", 1:length(sim.sig.col)), function(x) list(x))
                    
                    # 多重比較の関数を反復適用
                    for (i in 1:length(sim.sig.col)){
                        # ソースラベルを分解
                        source.term <- strsplit(sim.sig.col[i], " at ")
                        
                        if((source.term[[1]][1] == each.term[[1]][1]) && (level.num1 >= 3) == TRUE){
                            # 第一項の単純主効果が有意；subdat1を使用
                            multdat <- subdat1[[source.term[[1]][2]]]# 対象となる水準のデータのみ抽出
                            # 混合要因計画の被験者間要因の場合は複数の行から適切なMSeを選択；それ以外の場合は単一の行をコピー
                            if(!is.na(bet.mse1) && nrow(bet.mse1) != 1){
                                bet.mse <- bet.mse1[match(paste("Er at", source.term[[1]][2]), bet.mse1$source.col), ]
                            }else{
                                bet.mse <- bet.mse1
                            }
                            sim.multresults[[i]] <- mod.Bon(dat = multdat, design = postdesign1, taref = target.effect1, bet.mse = bet.mse, factlabel = sim.sig.lab[i],
                            factnames = factnames[-(col.num2-1)], type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
                        }else if((source.term[[1]][1] == each.term[[1]][2]) && (level.num2 >= 3) == TRUE){
                            # 第二項の単純主効果が有意；subdat2を使用
                            multdat <- subdat2[[source.term[[1]][2]]]# 対象となる水準のデータのみ抽出
                            # 混合要因計画の被験者間要因の場合は複数の行から適切なMSeを選択；それ以外の場合は単一の行をコピー
                            if(!is.na(bet.mse2) && nrow(bet.mse2) != 1){
                                bet.mse <- bet.mse2[match(paste("Er at", source.term[[1]][2]), bet.mse2$source.col), ]
                            }else{
                                bet.mse <- bet.mse2
                            }
                            sim.multresults[[i]] <- mod.Bon(dat = multdat, design = postdesign2, taref = target.effect2, bet.mse = bet.mse, factlabel = sim.sig.lab[i],
                            factnames = factnames[-(col.num1-1)], type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
                        }else{sim.multresults[[i]] <- NA}
                    }
                }
                simresults <- list(sig.source, sim.stat, simepsi, simtab, sim.multresults)
                names(simresults) <- c(sig.term, "sim.stat", "simepsi", "simtab", "sim.multresults")
                return(simresults)
            }else if(sig.num == 1){
                # 単純主効果が有意で，水準数が３以上であれば多重比較を行う
                # 分析対象となるの項の列番号，水準数を取得
                col.num <- charmatch(sig.term, names(dat))# 列番号
                level.num <- nlevels(dat[, col.num])# 水準数
                
                # 水準数が３以上の場合にのみ，mod.Bon関数を適用する
                if(level.num >= 3){
                    multfact <- grep(sig.term, bet.with[[1]][1])# 被験者間要因なら１，被験者内要因なら０を返す
                    
                    # 被験者間要因か被験者内要因かを判定して，mod.Bon関数にデータを送る
                    if(is.na(multfact[1]) == FALSE){
                        if(substr(design, nchar(design), nchar(design)) == "s"){
                            # 被験者間計画なら全体の誤差項のMSを得る
                            bet.mse <- anovatab[charmatch("Error", internal.lab),]
                        }else{
                            # 混合要因計画ならその効果の誤差項のMSを得る
                            f.denomi <- rep(mse.row, c(mse.row[1], diff(mse.row)))
                            bet.mse <- anovatab[f.denomi[charmatch(sig.term, internal.lab)],]
                        }
                    }else{
                        bet.mse <- NA
                    }
                    bonout <- mod.Bon(dat = dat, design = design, taref = sig.term, bet.mse = bet.mse, factlabel = sig.source, factnames = factnames,
                    type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
                    return(bonout)
                }else{
                    return(NA)# 水準数が２ならNAを返す
                }
            }
        }
        
        
        # 基本情報を出力する関数
        info.out <- function(info1, info2, info3){
            cat("\n")# １行空ける
            cat(sprintf(info1), "\n", "\n")
            cat(sprintf(info2), "\n")
            cat(sprintf(info3), "\n", "\n", "\n")
        }
        
        
        # 平均と標準偏差の表を出力する関数
        bstat.out <- function(baseresults, maxfact, margin){
            # 情報を展開
            bstat.info1 <- baseresults$bstat.info1
            bstatist <- baseresults$bstatist
            
            tabline <- pmax(3, sapply(1:ncol(bstatist), function(x) ifelse(is.numeric(bstatist[, x]), max(nchar(sprintf("%.4f", bstatist[, x]), type = "width")), max(nchar(as.character(bstatist[, x]), type = "width")))) + 1)
            tabhead <- nchar(names(bstatist)) + 1
            tabline <- ifelse(tabline > tabhead, tabline, tabhead)
            tabline[match("n", names(bstatist))] <- pmax(3, max(nchar(bstatist$n)))# サンプルサイズは小数点を取ることがないので短くする
            mainline <- sum(tabline) + length(tabline)
            headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
            bodymargin <- headmargin
            bodymargin[(maxfact+2):length(tabline)] <- paste0("%", tabline[(maxfact+2):length(tabline)], ".4f")
            
            cat(sprintf("<< DESCRIPTIVE STATISTICS >>"), sep = "", "\n", "\n")# タイトル
            if(sum(is.na(bstat.info1)) != 1){# bstat.info1がNAでないときのみ表示
                cat(sprintf(bstat.info1), sep = "\n")# プロンプトを表示
                cat("\n")# １行空ける
            }
            
            cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")# ラインを引く；要因の数に合わせて長さを調整
            cat(sprintf(headmargin, names(bstatist)), sep = " ", "\n")# データフレームの列名をスペース区切りで表示する
            cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")
            
            # データフレームのデータ部分を一行ずつスペース区切りで表示する
            for (i in 1:nrow(bstatist)){
                if((i %% margin == 0) && !(i == nrow(bstatist))){
                    cat(sprintf(bodymargin[1:maxfact], unlist(bstatist[i, 1:maxfact])), sprintf(bodymargin[maxfact+1], paste(bstatist[i, maxfact+1])), sprintf(bodymargin[(maxfact+2):ncol(bstatist)], bstatist[i, (maxfact+2):ncol(bstatist)]), sep = " ", "\n")
                    cat("\n")
                }else{
                    cat(sprintf(bodymargin[1:maxfact], unlist(bstatist[i, 1:maxfact])), sprintf(bodymargin[maxfact+1], paste(bstatist[i, maxfact+1])), sprintf(bodymargin[(maxfact+2):ncol(bstatist)], bstatist[i, (maxfact+2):ncol(bstatist)]), sep = " ", "\n")
                }
            }
            cat(sprintf(rep("-", mainline + 2)), sep = "", "\n", "\n", "\n")
        }
        
        
        # 分散分析表を出力する関数
        anovatab.out <- function(mainresults){
            # 情報を展開
            epsi.info1 <- mainresults$epsi.info1
            epsitab <- mainresults$epsitab
            ano.info1 <- mainresults$ano.info1
            anovatab <- mainresults$anovatab
            mse.row <- mainresults$mse.row
            
            # 球面性検定の結果を出力
            if(is.na(epsi.info1) == FALSE){# epsi.info1がNAでないときのみ出力
                tabline <- sapply(1:ncol(epsitab), function(x) ifelse(is.numeric(epsitab[, x]), max(nchar(round(epsitab[, x], 4), type = "width")), max(nchar(as.character(epsitab[,x]), type = "width")))) + 1
                tabline <- ifelse(nchar(names(epsitab)) > tabline, nchar(names(epsitab)) + 1, tabline)
                tabline[5] <- 7
                tabline[6] <- 3
                tabline[7:ncol(epsitab)] <- 6
                mainline <- sum(tabline) + length(tabline)
                headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
                bodymargin <- sapply(tabline, function(x) paste0("%", x, ".4f"))
                bodymargin[1] <- paste0("%", tabline[1], "s")
                bodymargin[4] <- paste0("%", tabline[4], "s")
                bodymargin[6] <- "%-3s"
                
                cat(sprintf("<< SPHERICITY INDICES >>"), sep = "", "\n", "\n")# タイトル
                cat(sprintf(epsi.info1), "\n", "\n")# プロンプト
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1:5], names(epsitab[1:5])), sprintf(headmargin[6], paste("")), sprintf(headmargin[7:ncol(epsitab)], names(epsitab[7:ncol(epsitab)])), sep = " ", "\n")# 列名をスペース区切りで出力
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する
                for (i in 1:nrow(epsitab)){
                    cat(sprintf(bodymargin[1], paste(epsitab[i, 1])), sprintf(bodymargin[2], epsitab[i, 2]), sprintf(bodymargin[3], epsitab[i, 3]), 
                    sprintf(bodymargin[4], epsitab[i, 4]), replace(sprintf(bodymargin[5], epsitab[i, 5]), is.na(epsitab[i, 5]), sprintf(paste("%", tabline[5], "s", sep = ""), "")), 
                    replace(sprintf(bodymargin[6], paste(epsitab[i, 6])), epsitab[i, 6] == "", sprintf(bodymargin[6], "")), sprintf(bodymargin[7], epsitab[i, 7:ncol(epsitab)]), sep = " ", "\n")
                }
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                if(ncol(epsitab) == 9){
                    cat(sprintf(paste0("%", mainline - 1, "s"), "LB = lower.bound, GG = Greenhouse-Geisser, HF = Huynh-Feldt"), sep = " ", "\n", "\n", "\n")
                }else{
                    cat(sprintf(paste0(rep("", 3))), sprintf(paste0("%", mainline - 1, "s"), "b_hat = b^, c_hat = c^, h_d = h~', h_dd = h~'', h = h~"), sep = " ", "\n", "\n", "\n")
                }
            }
            
            # 分散分析表を出力
            cat(sprintf("<< ANOVA TABLE >>"), sep = "", "\n", "\n")# タイトル
            if(sum(is.na(ano.info1)) != 1){# ano.info1がNAでないときのみ表示
                cat(sprintf(ano.info1), sep = "\n")# プロンプトを表示
                cat("\n")# １行空ける
            }
            
            tabline <- sapply(1:ncol(anovatab), function(x) ifelse(is.numeric(anovatab[, x]), max(nchar(sprintf("%.4f", anovatab[, x]), type = "width")), max(nchar(as.character(anovatab[, x]), type = "width")))) + 1
            tabline[3] <- pmax(max(nchar(round(anovatab[, 3], 4), type = "chars")), 3)# ss.colとdf.colが近くなりがちなので最小幅を広めに指定
            tabline[6] <- 8# F-ratioとp-valueの間を少し広く指定
            mainline <- sum(tabline) + length(tabline)
            headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
            bodymargin <- sapply(tabline, function(x) paste0("%", x, ".4f"))
            bodymargin[1] <- paste0("%", tabline[1], "s")
            bodymargin[3] <- paste0("%", tabline[3], "s")
            bodymargin[7] <- paste0("%-", tabline[7], "s")
            
            esn <- ncol(anovatab) - 7
            if(esn == 0){# 分散分析表のみの場合
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"), 
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")), sep = " ", "\n")# 列名をスペース区切りで表示する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する；誤差項の行（mse.row）の後にラインを入れる
                for (i in 1:nrow(anovatab)){
                    if(i == nrow(anovatab)){# 最終行に到達したときには以下を実行
                        cat(sprintf(bodymargin[1], format(anovatab[i, 1], justify = "right")), sprintf(bodymargin[2], anovatab[i, 2]), 
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sep = " ", "\n")
                        cat(sprintf(paste0("%", mainline - 2, "s"), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = " ", "\n")
                    }else if(sum(mse.row == i) == 1){# mse.rowの中に一致する値があるときには以下を実行
                        cat(sprintf(bodymargin[1], format(anovatab[i, 1]), justify = "right"), sprintf(bodymargin[2], anovatab[i, 2]), 
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]), 
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""), 
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""), 
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), sep = " ", "\n")
                        cat(sprintf(rep("-", mainline)), sep = "", "\n")
                    }else{# その他のときは以下を実行
                        cat(sprintf(bodymargin[1], format(anovatab[i, 1]), justify = "right"), sprintf(bodymargin[2], anovatab[i, 2]), 
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]), 
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""), 
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""), 
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), sep = " ", "\n")
                    }
                }
            }else{# 効果量の指標を追加した場合
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"), 
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")),	
                sprintf(headmargin[8:length(headmargin)], names(anovatab[8:ncol(anovatab)])), sep = " ", "\n")# 列名をスペース区切りで表示する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する；誤差項の行（mse.row）の後にラインを入れる
                for (i in 1:nrow(anovatab)){
                    if(i == nrow(anovatab)){# 最終行に到達したときには以下を実行
                        cat(sprintf(bodymargin[1], format(anovatab[i, 1], justify = "right")), sprintf(bodymargin[2], anovatab[i, 2]), 
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), "\n", sprintf(paste0("%", mainline - 2, "s"), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = " ", "\n")
                    }else if(sum(mse.row == i) == 1){# mse.rowの中に一致する値があるときには以下を実行
                        cat(sprintf(bodymargin[1], format(anovatab[i, 1], justify = "right")), sprintf(bodymargin[2], anovatab[i, 2]), 
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]), 
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""), 
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""), 
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), sep = " ", "\n")
                        cat(sprintf(rep("-", mainline)), sep = "", "\n")
                    }else{# その他のときは以下を実行
                        cat(sprintf(bodymargin[1], format(anovatab[i, 1], justify = "right")), sprintf(bodymargin[2], anovatab[i, 2]), 
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]), 
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""), 
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""), 
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), 
                        replace(sprintf(bodymargin[8:length(bodymargin)], anovatab[i, 8:ncol(anovatab)]), is.na(anovatab[i, 8:ncol(anovatab)]), ""), sep = " ", "\n")
                    }
                }
            }
            cat("\n", "\n")# ２行空ける
        }
        
        
        # 修正Bonferroniの方法による多重比較の結果を出力する関数
        mod.Bon.out <- function(bon.list, omit = FALSE){
            # 情報を展開
            factlabel <- bon.list[[1]]
            bon.info1 <- bon.list[[2]]
            bon.info2 <- bon.list[[3]]
            bon.info3 <- bon.list[[4]]
            bonstat <- bon.list[[5]]
            bontab <- bon.list[[6]]
            
            cat(sprintf(paste0("< MULTIPLE COMPARISON for \"", factlabel, "\" >")), sep = "", "\n", "\n")# タイトル
            cat(sprintf(bon.info1), sep = "\n")# プロンプト
            cat(sprintf(bon.info2), "\n")
            cat(sprintf(bon.info3), "\n", "\n")
            
            # omitがFなら記述統計量を出力する
            if(omit == FALSE){
                # 記述統計量を出力する
                tabline <- sapply(1:ncol(bonstat), function(x) ifelse(is.numeric(bonstat[, x]), max(nchar(sprintf("%.4f", bonstat[, x]), type = "width")), max(nchar(as.character(bonstat[, x]), type = "width")))) + 1
                tabhead <- nchar(names(bonstat)) + 1
                tabline <- ifelse(tabline > tabhead, tabline, tabhead)
                tabline[2] <- pmax(3, max(nchar(bonstat[, 2])))# サンプルサイズは小数点を取ることがないので短くする
                mainline <- sum(tabline) + length(tabline)
                headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
                bodymargin <- sapply(tabline, function(x) paste0("%", x, "s"))
                bodymargin[(length(tabline)-1):length(tabline)] <- paste0("%", tabline[(length(tabline)-1):length(tabline)], ".4f")
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")# ラインを引く；要因の数に合わせて長さを調整
                cat(sprintf(headmargin, names(bonstat)), sep = " ", "\n")# データフレームの列名をスペース区切りで表示する
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する
                for (i in 1:nrow(bonstat)){
                    cat(sprintf(bodymargin[1], bonstat[i, 1]), sprintf(bodymargin[2], paste(bonstat[i, 2])), sprintf(bodymargin[3:4], bonstat[i, 3:4]), sep = " ", "\n")
                }
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n", "\n")	
            }
            
            # 多重比較の結果を表形式で出力する
            tabline <- sapply(1:ncol(bontab), function(x) ifelse(is.numeric(bontab[, x]), max(nchar(sprintf("%.4f", bontab[, x]), type = "width")), max(nchar(as.character(bontab[, x]), type = "width")))) + 2
            tabline[4] <- max(nchar(bontab[, 4])) + 2# dfは小数点以下の値を取らないので別に調整
            mainline <- sum(tabline) + length(tabline)
            headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
            bodymargin <- sapply(tabline, function(x) paste0("%", x, ".4f"))
            bodymargin[1] <- paste("%", tabline[1], "s", sep = "")
            bodymargin[4] <- paste("%", tabline[4], "s", sep = "")
            bodymargin[7] <- paste("%", tabline[7], "s", sep = "")
            
            cat(sprintf(rep("-", mainline)), sep = "", "\n")
            cat(sprintf(headmargin[1], "Pair"), sprintf(headmargin[2], paste("Diff")), sprintf(headmargin[3], "t-value"), sprintf(headmargin[4], "df"), 
            sprintf(headmargin[5], "p"), sprintf(headmargin[6], names(bontab[6])), sprintf(headmargin[7], paste("")), sep = " ", "\n")# 列名をスペース区切りで表示する
            cat(sprintf(rep("-", mainline)), sep = "", "\n")
            
            # データフレームのデータ部分を一行ずつスペース区切りで表示する
            for (i in 1:nrow(bontab)){
                cat(sprintf(bodymargin[1], format(bontab[i, 1], justify = "right")), sprintf(bodymargin[2], bontab[i, 2]), sprintf(bodymargin[3], bontab[i, 3]), 
                sprintf(bodymargin[4], paste(" ", bontab[i, 4])), sprintf(bodymargin[5], bontab[i, 5]), sprintf(bodymargin[6], bontab[i, 6]), 
                sprintf(bodymargin[7], format(bontab[i, 7], justify = "right")), sep = " ", "\n")
            }
            
            cat(sprintf(rep("-", mainline)), sep = "", "\n", "\n")
        }
        
        
        # 単純主効果の検定の結果を出力する関数
        simeffects.out <- function(partresult, omit = FALSE){
            # 情報を展開
            part.info1 <- partresult[[1]]
            partstat <- partresult$sim.stat
            partepsi <- partresult$simepsi
            parttab <- partresult$simtab
            partmulttab <- partresult$sim.multresults
            
            cat(sprintf("%s", paste0("< SIMPLE EFFECTS for \"", part.info1, "\" INTERACTION >"), sep = ""), sep = "", "\n", "\n")# タイトル
            
            # omitがFALSEなら記述統計量を出力する
            if(omit == FALSE){
                # 記述統計量を出力する
                tabline <- pmax(3, sapply(1:ncol(partstat), function(x) ifelse(is.numeric(partstat[, x]), max(nchar(sprintf("%.4f", partstat[, x]), type = "width")), max(nchar(as.character(partstat[, x]), type = "width")))) + 1)
                tabhead <- nchar(names(partstat)) + 1
                tabline <- ifelse(tabline > tabhead, tabline, tabhead)
                tabline[3] <- pmax(3, max(nchar(partstat[, 3])))# サンプルサイズは小数点を取ることがないので短くする
                mainline <- sum(tabline) + length(tabline)
                headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
                bodymargin <- headmargin
                bodymargin[(length(tabline)-1):length(tabline)] <- paste0("%", tabline[(length(tabline)-1):length(tabline)], ".4f")
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")# ラインを引く；要因の数に合わせて長さを調整
                cat(sprintf(headmargin, names(partstat)), sep = " ", "\n")# データフレームの列名をスペース区切りで表示する
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する
                for (i in 1:nrow(partstat)){
                    cat(sprintf(bodymargin[1:2], unlist(partstat[i, 1:2])), sprintf(bodymargin[3], paste(partstat[i, 3])), sprintf(bodymargin[4:5], partstat[i, 4:5]), sep = " ", "\n")
                }
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n", "\n")	
                
            }
            
            # 球面性検定の結果を出力
            if(is.na(charmatch("Error", parttab$source.col)) && !is.na(partepsi)){
                # 被験者内計画，混合要因計画なら球面性検定の結果を出力する
                tabline <- sapply(1:ncol(partepsi), function(x) ifelse(is.numeric(partepsi[, x]), max(nchar(round(partepsi[, x], 4), type = "width")), max(nchar(as.character(partepsi[, x]), type = "width")))) + 1
                tabline <- ifelse(nchar(names(partepsi)) > tabline, nchar(names(partepsi)) + 1, tabline)
                tabline[5] <- 7
                tabline[6] <- 3
                tabline[7:ncol(partepsi)] <- 6
                mainline <- sum(tabline) + length(tabline)
                headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
                bodymargin <- sapply(tabline, function(x) paste0("%", x, ".4f"))
                bodymargin[1] <- paste0("%", tabline[1], "s")
                bodymargin[4] <- paste0("%", tabline[4], "s")
                bodymargin[6] <- "%-3s"
                
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1:5], names(partepsi[1:5])), sprintf(headmargin[6], paste("")), sprintf(headmargin[7:ncol(partepsi)], names(partepsi[7:ncol(partepsi)])), sep = " ", "\n")# 列名をスペース区切りで出力
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する
                for (i in 1:nrow(partepsi)){
                    cat(sprintf(bodymargin[1], paste(partepsi[i, 1])), sprintf(bodymargin[2], partepsi[i, 2]), sprintf(bodymargin[3], partepsi[i, 3]), 
                    sprintf(bodymargin[4], partepsi[i, 4]), replace(sprintf(bodymargin[5], partepsi[i, 5]), is.na(partepsi[i, 5]), sprintf(paste0("%", tabline[5], "s"), "")), 
                    replace(sprintf(bodymargin[6], paste(partepsi[i, 6])), partepsi[i, 6] == "", sprintf(bodymargin[6], "")), sprintf(bodymargin[7], partepsi[i, 7:ncol(partepsi)]), sep = " ", "\n")
                }
                
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                if(ncol(partepsi) == 9){
                    cat(sprintf(paste0("%", mainline - 1, "s"), "LB = lower.bound, GG = Greenhouse-Geisser, HF = Huynh-Feldt"), sep = " ", "\n", "\n", "\n")
                }else{
                    cat(sprintf(paste0(rep("", 3))), sprintf(paste0("%", mainline - 1, "s"), "b_hat = b^, c_hat = c^, h_d = h~', h_dd = h~'', h = h~"), sep = " ", "\n", "\n", "\n")
                }
                
            }
            
            # 分散分析の結果を出力
            tabline <- sapply(1:ncol(parttab), function(x) ifelse(is.numeric(parttab[, x]), max(nchar(sprintf("%.4f", parttab[, x]), type = "width")), max(nchar(as.character(parttab[, x]), type = "width")))) + 1
            tabline[3] <- pmax(max(nchar(round(parttab[,3], 4), type = "width")), 3)# ss.colとdf.colが近くなりがちなので最小幅を広めに指定
            tabline[6] <- 8# F-ratioとp-valueの間を少し広く指定
            mainline <- sum(tabline) + length(tabline)
            headmargin <- sapply(tabline, function(x) paste0("%", x, "s"))
            bodymargin <- sapply(tabline, function(x) paste0("%", x, ".4f"))
            bodymargin[1] <- paste0("%", tabline[1], "s")
            bodymargin[3] <- paste0("%", tabline[3], "s")
            bodymargin[7] <- paste0("%-", tabline[7], "s")
            
            esn <- ncol(parttab) - 7
            if(esn == 0){# 分散分析表のみの場合
                # 単純主効果の検定の分散分析表を出力する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"), 
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")), sep = " ", "\n")# 列名をスペース区切りで表示する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する
                if(is.na(charmatch("Error", parttab$source.col))){
                    # 被験者内計画，混合要因計画なら２行ごとにラインを入れる
                    for (i in 1:nrow(parttab)){
                        cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                        sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), 
                        replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""), 
                        replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""), 
                        sprintf(bodymargin[7], paste(parttab[i, 7])), sep = " ", "\n")
                        if((i %% 2) == 0){# 行番号が２で割り切れるときはラインを引く
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                        }
                    }
                    # 表の末尾に有意性判定を示す記号を記す
                    cat(sprintf(paste0("%", mainline - 2, "s"), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = " ", "\n", "\n")
                }else{
                    # 被験者間計画なら誤差項の上の行だけにラインを入れる
                    for (i in 1:nrow(parttab)){
                        if(nrow(parttab) == i){# 行番号が最後の行に一致するときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), sep = " ", "\n")
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                            cat(sprintf(paste0("%", mainline - 2, "s"), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = " ", "\n", "\n")
                        }else{# その他のときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), 
                            replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""), 
                            replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""), 
                            sprintf(bodymargin[7], paste(parttab[i, 7])), sep = " ", "\n")
                        }
                    }
                }
            }else{# 効果量の指標を追加した場合
                # 単純主効果の検定の分散分析表を出力する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"), 
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")), sprintf(headmargin[8:ncol(parttab)], names(parttab[8:ncol(parttab)])), sep = " ", "\n")# 列名をスペース区切りで表示す
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつスペース区切りで表示する
                if(is.na(charmatch("Error", parttab$source.col))){
                    # 被験者内計画，混合要因計画なら２行ごとにラインを入れる
                    for (i in 1:nrow(parttab)){
                        cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                        sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), 
                        replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""), 
                        replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""), 
                        sprintf(bodymargin[7], paste(parttab[i, 7])), replace(sprintf(bodymargin[8:ncol(parttab)], 
                        parttab[i, 8:ncol(parttab)]), is.na(parttab[i, 8:ncol(parttab)]), ""), sep = " ", "\n")
                        if((i %% 2) == 0){# 行番号が２で割り切れるときはラインを引く
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                        }
                    }
                    # 表の末尾に有意性判定を示す記号を記す
                    cat(sprintf(paste0("%", mainline - 2, "s"), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = " ", "\n", "\n")
                }else{
                    # 被験者間計画なら誤差項の上の行だけにラインを入れる
                    for (i in 1:nrow(parttab)){
                        if(nrow(parttab) == i){# 行番号が最後の行に一致するときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), sep = " ", "\n")
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                            cat(sprintf(paste0("%", mainline - 2, "s"), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = " ", "\n", "\n")
                        }else{# その他のときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), 
                            replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""), 
                            replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""), 
                            sprintf(bodymargin[7], paste(parttab[i, 7])), replace(sprintf(bodymargin[8:ncol(parttab)], 
                            parttab[i, 8:ncol(parttab)]), is.na(parttab[i, 8:ncol(parttab)]), ""), sep = " ", "\n")
                        }
                    }
                }
            }
            
            # 多重比較の結果を出力する
            for (i in 1:length(partmulttab)){
                if(is.na(partmulttab[[i]][1]) == FALSE){# リストに中身があったら出力する
                    cat("\n")# １行空ける
                    mod.Bon.out(partmulttab[[i]], omit = TRUE)
                }
            }
            
            cat("\n")# １行空ける
        }
        
        
        # 下位検定の結果を出力する関数
        post.out <- function(postresults, design){
            # リストに含まれる要素の数を特定
            postnum <- length(postresults)
            
            # すべてのリストの中身がNAの場合には，タイトルも表示しない
            if(sum(is.na(postresults)) != postnum){
                cat(sprintf("<< POST ANALYSES >>"), sep = "", "\n", "\n")# タイトル
            }
            
            # リストの中身をひとつずつeach.out関数に送る
            for (i in 1:postnum){
                each.out(postresults[[i]], design)
            }
            
            # 出力が終わったことを知らせるプロンプト
            cat(sprintf("output is over "), sprintf(rep("-", 20)), sprintf("///"), sep = "", "\n")
            cat("\n")# １行空ける
        }
        
        
        # 効果のタイプによって異なる出力方式を割り当てる関数
        each.out <- function(partresult, design){
            if(is.null(names(partresult)[1]) == TRUE){# 何もない場合（２水準の主効果が有意で多重比較の必要なしの場合）
            }else if(nchar(names(partresult)[1]) == 1){# 多重比較の結果がある場合
                mod.Bon.out(partresult)
                cat("\n")# １行空ける
            }else if(nchar(names(partresult)[1]) == 3){# 一次の交互作用が見られた場合
                if(nchar(design) == 3) simeffects.out(partresult, omit = TRUE)# もともと２要因の計画のときには，単純主効果の検定の出力時に記述統計量を出力しない（主分析の記述統計量と重複するので）
                else simeffects.out(partresult)
            }else if(nchar(names(partresult)[1]) >= 5){# 高次の交互作用が見られた場合
                cat(sprintf(paste0("< HIGHER-ORDER \"", partresult[[1]], "\" INTERACTION >")),sep = "", "\n")
                cat(sprintf("*** Split the dataset for further analysis. ***"), sep ="", "\n", "\n")# データ分割を促すプロンプト
            }
        }
        
        
        # 指定した要因についてデータを分割して単純効果の検定を行う関数；３要因以上の計画での使用を想定
        # tfactの引数として分割の基準とする要因のラベルを指定する（tfact = "A"など）
        # 同時に複数の要因について分割する場合，tfact = c("A", "C")のような形式で入力する
        anovatan <- function(dataset, design, ..., tfact = NULL, long = FALSE, type2 = FALSE, nopost = FALSE, tech = FALSE, data.frame = FALSE, copy = FALSE, 
        holm = FALSE, hc = FALSE, s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE, 
        lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, 
        eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE, 
        gomega = NA, gomegana = NA, prep = FALSE, cind = FALSE, cin = FALSE){
            # 分割する要因が指定されていない場合は終了
            if(is.null(tfact) == TRUE){
                stop(message = "\"anovatan\" has stopped working...\nPlease specify some factor to split.")
            }
            
            # 要因計画情報の取得
            bet.with <- strsplit(design, "s")
            betN <- nchar(bet.with[[1]])[1]# 被験者間要因がないときは“０”
            withN <- nchar(bet.with[[1]])[2]# 被験者内要因がないときは“NA”
            
            # データフレームの変形
            datform <- uni.long(dataset = dataset, design = design, ... = ..., long = long)
            dat <- datform$dat
            factnames <- datform$factnames
            miscase <- datform$miscase
            names(dat) <- c("s", factnames, "y")
            targetnm <- sapply(tfact, function(x) grep(x, factnames))# 分割する要因が主分析で何番目の要因に当たるかを表すベクトル
            
            # データと要因計画の準備
            subdat <- split(dat, dat[, tfact])# データフレームを分割
            whichfact <- sapply(LETTERS[targetnm], function(x) grep(x, bet.with[[1]]))# 分割する要因が被験者間要因と被験者内要因のどちらに属するのかを表すベクトル
            bet.withN <- c(betN, replace(withN, is.na(withN), 0)) - c(sum(whichfact == 1), sum(whichfact == 2))# 下位検定における被験者間要因・被験者内要因の数を表すベクトル
            subdesign <- paste0(paste0(LETTERS[(1 * (bet.withN[1] != 0)):bet.withN[1]], collapse = ""), "s", paste0(LETTERS[1 * (bet.withN[2] != 0) * (bet.withN[1]+1:bet.withN[2])], collapse = ""))# 下位検定の計画の型
            
            # copyオプションの指定があった場合，出力をクリップボードにコピー
            if(copy == TRUE){
                plat.info <- .Platform
                if(sum(grep("windows", plat.info)) != 0){# Windowsの場合
                    sink("clipboard", split = TRUE)
                }else if(sum(grep("mac", plat.info)) != 0){# Macの場合
                    tclip <- pipe("pbcopy", "w")
                    sink(tclip, split = TRUE)
                }else if(sum(grep("linux", R.version$system)) != 0){# Linxの場合（xclipをインストールしている必要がある）
                    tclip <- pipe("xclip -selection clipboard")
                    sink(tclip, split = TRUE)
                }
            }
            
            # 除外したケースの報告
            if(miscase != 0){
                cat("[[ Information for Whole Analysis ]]", "\n")
                cat(paste0("== The number of removed case is ", miscase, ". =="), "\n", "\n")
            }
            
            # anovakunの実行
            for(i in 1:length(subdat)){
                cat("\n", sprintf(paste0("[[ Simple Effects for ", names(subdat)[i], " ]]")), "\n", sep = "")
                anovakun(dataset = subdat[[i]][, -(targetnm + 1)], design = subdesign, long = TRUE, type2 = type2, nopost = nopost, tech = tech, 
                data.frame = FALSE, copy = FALSE, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, 
                criteria = criteria, lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, 
                eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega, 
                gomega = gomega, gomegana = gomegana, prep = prep, cind = cind, cin = cin)
            }
            
            # copyを実行していた場合，クリップボードを閉じる
            if(copy == TRUE){
                sink()
                if(plat.info$OS.type != "windows"){# Mac，Linuxの場合
                    close(tclip)
                }
            }
            
            # 指定があった場合は，分割後のデータフレームを出力
            if(data.frame == TRUE){
                return(subdat)
            }
        }




    # ANOVA君ここまで
    #--------------------------------------------------------------
    
    # Results

        # One-way ANOVA
        
        if (input$factor == "oneway") {

            type <- switch(input$one.design,
                        Between = "As",
                        Within = "sA")
            
            level1 <- input$factor1.level

            anovakun(dat, type, level1, mau=T, auto=T, holm = T, peta=T)
            
        } else { # Two-way ANOVA
            
            type <- switch(input$two.design,
                        Factor1Between_Factor2Between = "ABs",
                        Factor1Between_Factor2Within = "AsB",
                        Factor1Within_Factor2Within = "sAB")
            
            level1 <- input$factor1.level
            level2 <- input$factor2.level
            
            anovakun(dat, type, level1, level2, mau=T, auto=T, holm = T, peta=T)
            
        }

    })





    makeANOVAPlot1 <- function(){
        if (input$axis == "default") {
            
            if (input$factor == "oneway") {
                
                if (input$one.design == "Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, ci.fun=function(x) c(mean(x)-1.96*se(x), mean(x)+1.96*se(x)),
                                xlab="Level", ylab="", main="Factor 1 (Error bars show 95% CI)")
                    
                } else {  # Within = "sA"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=1:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, ci.fun=function(x) c(mean(x)-1.96*se(x), mean(x)+1.96*se(x)),
                                xlab="Level", ylab="", main="Factor 1 (Error bars show 95% CI)")

                }
                
            } else { # Two-way ANOVA
                
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat[,1] <- as.factor(dat[,1])
                    dat[,2] <- as.factor(dat[,2])
                    names(dat) <- c("fac1", "fac2", "res")
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(),
                                         panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac1, y=x, colour=fac2, group=fac2))
                    gp + geom_line(aes(linetype=fac2), size=.6, position=xgap) +
                    geom_point(aes(shape=fac2), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 1 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 2") +
                    scale_linetype_discrete(name  ="Factor 2") +
                    scale_shape_discrete(name  ="Factor 2")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    res <- dat[,3]
                    names(dat) <- c("fac1", "fac2", "res")
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac1, y=x, colour=fac2, group=fac2))
                    gp + geom_line(aes(linetype=fac2), size=.6, position=xgap) +
                    geom_point(aes(shape=fac2), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 1 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 2") +
                    scale_linetype_discrete(name  ="Factor 2") +
                    scale_shape_discrete(name  ="Factor 2")
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    z <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    dat <- reshape(z, varying=1:ncol(z), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(dat)/nlva, length.out = nrow(dat))
                    fac1 <- anms
                    dat <- cbind(dat, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(z), length.out = nrow(dat))
                    fac2 <- bnms
                    dat <- cbind(dat, fac2)
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac1, y=x, colour=fac2, group=fac2))
                    gp + geom_line(aes(linetype=fac2), size=.6, position=xgap) +
                    geom_point(aes(shape=fac2), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 1 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18), 
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18), 
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 2") +
                    scale_linetype_discrete(name  ="Factor 2") +
                    scale_shape_discrete(name  ="Factor 2")
                }
            }
            
        } else { # y-axis min-max
            
            if (input$factor == "oneway") {
                
                if (input$one.design == "Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, ylim=c(min(res), max(res)),
                                ci.fun=function(x) c(mean(x)-1.96*se(x), mean(x)+1.96*se(x)),
                                xlab="Level", ylab="", main="Factor 1 (Error bars show 95% CI)")
                    
                } else {  # Within = "sA"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=1:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, ylim=c(min(res), max(res)),
                    ci.fun=function(x) c(mean(x)-1.96*se(x), mean(x)+1.96*se(x)),
                    xlab="Level", ylab="", main="Factor 1 (Error bars show 95% CI)")

                }
                
            } else { # Two-way ANOVA
                
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat[,1] <- as.factor(dat[,1])
                    dat[,2] <- as.factor(dat[,2])
                    names(dat) <- c("fac1", "fac2", "res")
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac1, y=x, colour=fac2, group=fac2))
                    gp + geom_line(aes(linetype=fac2), size=.6, position=xgap) +
                    geom_point(aes(shape=fac2), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 1 (Error bars show 95% CI)", x = "Level", y = "") +

                    ylim(min(dat$res), max(dat$res)) +

                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 2") +
                    scale_linetype_discrete(name  ="Factor 2") +
                    scale_shape_discrete(name  ="Factor 2")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    z <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    res <- z[,3]
                    names(z) <- c("fac1", "fac2", "res")
                    df <- with(z, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(z , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac1, y=x, colour=fac2, group=fac2))
                    gp + geom_line(aes(linetype=fac2), size=.6, position=xgap) +
                    geom_point(aes(shape=fac2), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    ylim(min(z$res), max(z$res)) +
                    labs(title = "Factor 1 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 2") +
                    scale_linetype_discrete(name  ="Factor 2") +
                    scale_shape_discrete(name  ="Factor 2")
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    z <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    dat <- reshape(z, varying=1:ncol(z), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(dat)/nlva, length.out = nrow(dat))
                    fac1 <- anms
                    dat <- cbind(dat, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(z), length.out = nrow(dat))
                    fac2 <- bnms
                    dat <- cbind(dat, fac2)
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac1, y=x, colour=fac2, group=fac2))
                    gp + geom_line(aes(linetype=fac2), size=.6, position=xgap) +
                    geom_point(aes(shape=fac2), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 1 (Error bars show 95% CI)", x = "Level", y = "") +
                    ylim(min(dat$res), max(dat$res)) +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 2") +
                    scale_linetype_discrete(name  ="Factor 2") +
                    scale_shape_discrete(name  ="Factor 2")
                }
            }
        }
    }

    output$AnovaPlot1 <- renderPlot({
        print(makeANOVAPlot1())
    })
    
    
    
    
    makeANOVAPlot2 <- function(){
        
        if (input$axis == "default") {
            
            if (input$factor == "twoway") {
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat[,1] <- as.factor(dat[,1])
                    dat[,2] <- as.factor(dat[,2])
                    names(dat) <- c("fac1", "fac2", "res")
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac2, y=x, colour=fac1, group=fac1))
                    gp + geom_line(aes(linetype=fac1), size=.6, position=xgap) +
                    geom_point(aes(shape=fac1), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 2 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 1") +
                    scale_linetype_discrete(name  ="Factor 1") +
                    scale_shape_discrete(name  ="Factor 1")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    res <- dat[,3]
                    names(dat) <- c("fac1", "fac2", "res")
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3] 
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])

                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac2, y=x, colour=fac1, group=fac1))
                    gp + geom_line(aes(linetype=fac1), size=.6, position=xgap) +
                    geom_point(aes(shape=fac1), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 2 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 1") +
                    scale_linetype_discrete(name  ="Factor 1") +
                    scale_shape_discrete(name  ="Factor 1")

                } else { # input$two.design == "Factor1Within_Factor2Within"
                    z <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    dat <- reshape(z, varying=1:ncol(z), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(dat)/nlva, length.out = nrow(dat))
                    fac1 <- anms
                    dat <- cbind(dat, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(z), length.out = nrow(dat))
                    fac2 <- bnms
                    dat <- cbind(dat, fac2)

                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3] 
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac2, y=x, colour=fac1, group=fac1))
                    gp + geom_line(aes(linetype=fac1), size=.6, position=xgap) +
                    geom_point(aes(shape=fac1), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 2 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 1") +
                    scale_linetype_discrete(name  ="Factor 1") +
                    scale_shape_discrete(name  ="Factor 1")
                }
                
            } else {
                
                cat("\n")
                
            }
            
        } else { # y-axis min-max
            
            if (input$factor == "twoway") {
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat[,1] <- as.factor(dat[,1])
                    dat[,2] <- as.factor(dat[,2])
                    names(dat) <- c("fac1", "fac2", "res")
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac2, y=x, colour=fac1, group=fac1))
                    gp + geom_line(aes(linetype=fac1), size=.6, position=xgap) +
                    geom_point(aes(shape=fac1), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    
                    ylim(min(dat$res), max(dat$res)) +

                    labs(title = "Factor 2 (Error bars show 95% CI)", x = "Level", y = "") +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 1") +
                    scale_linetype_discrete(name  ="Factor 1") +
                    scale_shape_discrete(name  ="Factor 1")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    res <- dat[,3]
                    names(dat) <- c("fac1", "fac2", "res")
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac2, y=x, colour=fac1, group=fac1))
                    gp + geom_line(aes(linetype=fac1), size=.6, position=xgap) +
                    geom_point(aes(shape=fac1), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 2 (Error bars show 95% CI)", x = "Level", y = "") +
                    ylim(min(dat$res), max(dat$res)) +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 1") +
                    scale_linetype_discrete(name  ="Factor 1") +
                    scale_shape_discrete(name  ="Factor 1")
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    z <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    dat <- reshape(z, varying=1:ncol(z), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(dat)/nlva, length.out = nrow(dat))
                    fac1 <- anms
                    dat <- cbind(dat, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(z), length.out = nrow(dat))
                    fac2 <- bnms
                    dat <- cbind(dat, fac2)
                    
                    df <- with(dat, aggregate(res, list(fac1=fac1, fac2=fac2), mean))
                    ci <- with(dat , aggregate(res, list(fac1=fac1, fac2=fac2), function(x) qt(0.975,length(x)-1)*sd(x)/sqrt(length(x))))
                    df$ci <- ci[,3]
                    df$fac1 <- factor(df[,1])
                    df$fac2 <- factor(df[,2])
                    
                    opar <- theme_update(panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(colour = "black"))
                    xgap <- position_dodge(0.15)
                    gp <- ggplot(df, aes(x=fac2, y=x, colour=fac1, group=fac1))
                    gp + geom_line(aes(linetype=fac1), size=.6, position=xgap) +
                    geom_point(aes(shape=fac1), size=4, position=xgap) +
                    geom_errorbar(aes(ymax=x+ci, ymin=x-ci), width=.1, position=xgap) +
                    labs(title = "Factor 2 (Error bars show 95% CI)", x = "Level", y = "") +
                    ylim(min(dat$res), max(dat$res)) +
                    theme(plot.title=element_text(size=18), axis.title.x = element_text(size=18), axis.title.y = element_text(size=18),
                    axis.text.x = element_text(size=18), axis.text.y = element_text(size=18),
                    legend.position = c(1, 1), legend.justification = c(1, 1), legend.key = element_rect(fill = "white",colour = "white"),
                    legend.title = element_text(size=16), legend.text = element_text(size = 16)) +
                    scale_colour_discrete(name  ="Factor 1") +
                    scale_linetype_discrete(name  ="Factor 1") +
                    scale_shape_discrete(name  ="Factor 1")
                }
                
            } else {
                
                cat("\n")
                
            }
        }
    }
    
    output$AnovaPlot2 <- renderPlot({
        print(makeANOVAPlot2())
    })




    output$anovakun.out <- renderPrint({
        anovakun()
    })
    
    output$check.out <- renderPrint({
        check()
    })
    
    output$downloadANOVAPlot1 <- downloadHandler(
    filename = function() {
        paste('Plot1-', Sys.Date(), '.pdf', sep='')
    },
    content = function(FILE=NULL) {
        pdf(file=FILE)
		print(makeANOVAPlot1())
		dev.off()
	})
    
    output$downloadANOVAPlot2 <- downloadHandler(
    filename = function() {
        paste('Plot2-', Sys.Date(), '.pdf', sep='')
    },
    content = function(FILE=NULL) {
        pdf(file=FILE)
		print(makeANOVAPlot2())
		dev.off()
	})

})

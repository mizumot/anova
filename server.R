library(shiny)
library(shinyAce)
library(sciplot)




shinyServer(function(input, output) {


    
    options(warn=-1)
    
    
    anovakun <- reactive({
        
        dat <- read.csv(text=input$text, sep="\t")
        
        
        # 【ANOVA君：要因計画のタイプと水準数を入力することにより，分散分析を行う関数】
        # 1)フリーの統計ソフトウェア「R」で動作する関数
        # 2)被験者間要因（独立測度），被験者内要因（反復測度）のいずれか，または，両方を含む各タイプの分散分析を扱う
        # 3)引数としては，最初にデータフレーム名，次に計画のタイプ（""で囲むこと）を入力し，その後，各要因の水準数を順に入力する
        # （作成：井関龍太）
        # http://www11.atpages.jp/~riseki/pukiwikiplus/index.php?ANOVA%B7%AF
        
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
                warning("Unknown operating system !!")
                return(NA)
            }
            
            # read.table関数を実行
            read.table(file = outboard, header = header, sep = sep, quote = quote, dec = dec, row.names = row.names, col.names =
            col.names, as.is = as.is, na.strings = na.strings, colClasses = colClasses, nrows = nrows, skip = skip,
            check.names = check.names, fill = fill, strip.white = strip.white, blank.lines.skip = blank.lines.skip,
            comment.char = comment.char, allowEscapes = allowEscapes, flush = flush, stringsAsFactors = stringsAsFactors,
            fileEncoding = fileEncoding, encoding = encoding, text = text)
        }
        
        
        # ANOVA君本体：データフレームの作成，プロセス全体の制御を行う関数
        anovakun <- function(dataset, design, ..., type2 = FALSE, nopost = FALSE, tech = FALSE, data.frame = FALSE, copy = FALSE, holm = FALSE,	hc = FALSE, s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE, lb = FALSE, gg = FALSE, hf = FALSE,
        auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, eta = FALSE, peta = FALSE, geta = NA,
        eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE, gomega = NA, gomegana = NA, prep = FALSE) {
            bet.with <- strsplit(design, "s")
            betN <- nchar(bet.with[[1]])[1]# 被験者間要因がないときは“０”
            withN <- nchar(bet.with[[1]])[2]# 被験者内要因がないときは“NA”
            maxfact <- nchar(design) - 1# 実験計画全体における要因数
            eachlev <- unlist(list(...))# 各要因の水準数
            
            # 欠損値の除去
            misv <- suppressWarnings(as.numeric(sapply(dataset[,(betN+1):ncol(dataset)], function(x) as.vector(x))))# 数値以外の値をNAに強制変換
            cdata <- array(misv, c(nrow(dataset), ncol(dataset) - betN))# NAを含む行列
            compcase <- complete.cases(cdata)# 完全ケース
            dataset <- dataset[compcase,]# datasetからNAを含む行（ケース）を除く
            miscase <- sum(!compcase)# 欠損ケースの数
            depv <- as.vector(cdata[compcase,])# 従属変数をベクトル化
            
            # 各要因を表すラベルの作成
            if(betN == 0){# 被験者内計画の場合
                if(withN == 1){# １要因の場合
                    eachcue <- 1
                    timescue <- 1
                }else{# その他の場合
                    eachcue <- c(sapply(2:withN, function(x) prod(eachlev[x:withN])), 1)
                    timescue <- prod(eachlev) / (eachcue * eachlev)
                }
                
                # 各被験者内要因の各水準を表すベクトル
                mainlab <- as.data.frame(sapply(1:withN, function(x) rep(paste(letters[x], 1:eachlev[x], sep = ""),
                each = nrow(dataset) * eachcue[x], times = timescue[x])))
                
            }else if(is.na(withN)){# 被験者間計画の場合
                # もとのデータフレームの水準名を変更
                for(i in 1:betN){
                    levels(dataset[,i])[order(match(levels(dataset[,i]), dataset[,i]))] <- paste(letters[i], 1:nlevels(dataset[,i]), sep = "")
                    dataset[,i] <- reorder(dataset[,i], 1:nrow(dataset), mean)# 水準の順序をもとのデータフレームに合わせる
                }
                
                # 各被験者間要因の各水準を表すベクトル
                mainlab <- dataset[,1:betN, drop = FALSE]
                
            }else{# 混合要因計画の場合
                # もとのデータフレームの水準名を変更
                for(i in 1:betN){
                    levels(dataset[,i])[order(match(levels(dataset[,i]), dataset[,i]))] <- paste(letters[i], 1:nlevels(dataset[,i]), sep = "")
                    dataset[,i] <- reorder(dataset[,i], 1:nrow(dataset), mean)# 水準の順序をもとのデータフレームに合わせる
                }
                withlev <- prod(eachlev[-(1:betN)])# 被験者内要因のすべての水準数をかけたもの
                
                # 各被験者間要因の各水準を表すベクトル
                mainlab1 <- dataset[,1:betN, drop = FALSE]
                mainlab1 <- sapply(1:ncol(mainlab1), function(x) rep(mainlab1[,x], withlev))
                
                if(withN == 1){# １要因の場合
                    eachcue <- 1
                    timescue <- 1
                }else{# その他の場合
                    eachcue <- c(sapply(2:withN, function(x) prod(eachlev[(x+betN):(withN+betN)])), 1)
                    timescue <- withlev / (eachcue * eachlev[-(1:betN)])
                }
                
                # 各被験者内要因の各水準を表すベクトル
                mainlab2 <- as.data.frame(sapply(1:withN, function(x) rep(paste(letters[x+betN], 1:eachlev[x+betN], sep = ""),
                each = nrow(dataset) * eachcue[x], times = timescue[x])))
                
                mainlab <- cbind(mainlab1, mainlab2)
            }
            
            names(mainlab) <- LETTERS[1:maxfact]# 要因のラベルを列名に設定
            
            # データフレームにまとめる
            dat <- data.frame("s" = rep(paste("s", 1:nrow(dataset), sep = ""), ncol(dataset) - betN), mainlab, "y" = depv)
            
            # 作成したデータフレームをもとに，各条件ごとの平均と標準偏差を計算する
            sncol <- as.vector(tapply(dat$y, dat[,(maxfact+1):2], function(x) length(x)))# セルごとのデータ数を計算
            mncol <- as.vector(tapply(dat$y, dat[,(maxfact+1):2], function(x) mean(x)))# セルごとの平均を計算
            sdcol <- as.vector(tapply(dat$y, dat[,(maxfact+1):2], function(x) sd(x)))# セルごとの標準偏差を計算
            
            # 記述統計量の表において各要因の各水準を表すためのラベル（数字列）を作成
            maincols <- expand.grid(lapply((maxfact+1):2, function(x) levels(dat[, x])))
            maincols <- maincols[, order(maxfact:1)]# アルファベット順に並べ替え
            
            # 記述統計量をデータフレームにまとめる
            bstatist <- data.frame(maincols, sncol, mncol, sdcol)
            names(bstatist) <- c(LETTERS[1:maxfact], "N", "Mean", "S.D.")# 要因のラベルほかを列名に設定
            
            # 記述統計量の表を出力する際の改行位置を指定する
            if(maxfact < 3) margin <- prod(eachlev)
            else margin <- prod(eachlev[-(1:(maxfact-2))])
            
            # anova.modelerにデータフレームを送り，分散分析の結果を得る
            mainresults <- anova.modeler(dat = dat, design = design, type2 = type2,
            lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
            eta = eta, peta = peta, eps = eps, peps = peps, geps = geps, geta = geta, omega = omega, omegana = omegana, pomega = pomega,
            gomega = gomega, gomegana = gomegana, prep = prep)
            
            # post.analysesにデータフレームと分散分析の結果を送り，下位検定の結果を得る
            postresults <- post.analyses(dat = dat, design = design, mainresults = mainresults, type2 = type2, nopost = nopost,
            holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria,
            lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
            eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
            gomega = gomega, gomegana = gomegana, prep = prep)
            
            # 基本情報の取得
            info1 <- paste("[ ", design, "-Type Design ]", sep = "")# 要因計画の型
            info2 <- paste("This output was generated via anovakun 4.3.3 at ", strsplit(R.version$version.string, " \\(")[[1]][1], ".", sep = "")# バージョン情報など
            info3 <- paste("It was executed on ", date(), ".", sep = "")# 実行日時
            
            # Unbalancedデザイン（データ数ふぞろい）の場合，プロンプトを追加
            if(length(unique(sncol)) != 1){
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
                if(sum(is.na(mainresults$ano.info1)) == 1){
                    mainresults$ano.info1 <- paste("== The number of removed case is ", miscase, ". ==", sep = "")
                }
                else{
                    mainresults$ano.info1 <- append(mainresults$ano.info1, paste("== The number of removed case is ", miscase, ". ==", sep = ""))
                }
            }
            
            # 結果を表示する
            if(copy == TRUE){# 指定があった場合，出力をクリップボードにコピー
                plat.info <- .Platform
                if(sum(grep("windows", plat.info)) !=0){# Windowsの場合
                    sink("clipboard", split =TRUE)
                }else if(sum(grep("mac", plat.info)) != 0){# Macの場合
                    tclip <- pipe("pbcopy", "w")
                    sink(tclip, split = TRUE)
                }else if(sum(grep("linux", R.version$system)) != 0){# Linxの場合（xclipをインストールしている必要がある）
                    tclip <- pipe("xclip -selection clipboard")
                    sink(tclip, split = TRUE)
                }
            }
            if(tech == TRUE){# データフレーム形式での出力の場合
                if(data.frame == TRUE){return(list("INFORMATION" = rbind(info1, info2, info3),
                    "DESCRIPTIVE STATISTICS" = bstatist,
                    "SPHERICITY INDICES" = list(mainresults$epsi.info1, mainresults$epsitab),
                    "ANOVA TABLE" = list(mainresults$ano.info1, mainresults$anovatab),
                    "POST ANALYSES" = postresults,
                    "DATA.FRAME" = dat))# 計算に使用したデータフレームを付加
                }else{return(list("INFORMATION" = rbind(info1, info2, info3),
                    "DESCRIPTIVE STATISTICS" = bstatist,
                    "SPHERICITY INDICES" = list(mainresults$epsi.info1, mainresults$epsitab),
                    "ANOVA TABLE" = list(mainresults$ano.info1, mainresults$anovatab),
                    "POST ANALYSES" = postresults))
                }
            }else{# 表形式での出力の場合
                if(data.frame == TRUE){info.out(info1, info2, info3)
                    bstat.out(bstatist, maxfact, margin, copy = copy)
                    anovatab.out(mainresults, copy = copy)
                    post.out(postresults, design, copy = copy)
                    return(list("DATA.FRAME" = dat))# 計算に使用したデータフレームを付加
                }else{info.out(info1, info2, info3)
                    bstat.out(bstatist, maxfact, margin, copy = copy)
                    anovatab.out(mainresults, copy = copy)
                    post.out(postresults, design, copy = copy)
                }
            }
            if(copy == TRUE){
                sink()
                if(plat.info$OS.type != "windows"){# Mac，Linuxの場合
                    close(tclip)
                }
            }
        }
        
        
        # 文字列中のすべての要素を含む文字列をマッチングする関数
        # grepとの違いは“A:C”などの文字列を照合パターンとした場合に“A:B:C”のように間に別の文字を挟んだ文字列もマッチと判定する点
        # 照合パターンが１文字の場合はgrepと同じ結果を返す
        elematch <- function(Mstrings, stex){
            # マッチングする文字列を分解して，それぞれgrep関数を適用
            matchlist <- lapply(strsplit(Mstrings, "")[[1]], function(x) grep(x, stex))
            
            # 文字列の各要素とマッチした値の共通部分のみ取り出す
            buffer <- matchlist[[1]]
            if(length(matchlist) != 1){
                for(i in 2:length(matchlist)){
                    buffer <- buffer[is.element(buffer, matchlist[[i]])]
                }
            }# 文字列が１文字のときはgrepの結果をそのまま返す
            
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
            eff.modeleq <- paste("~ ", gsub(",", " +", toString(eff.elem)), sep = "")
            er.modeleq <- paste("~ ", gsub(",", " +", toString(er.elem)), sep = "")
            
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
            names(ortho.model) <- paste("R", repleng:1, sep = "")
            ortho.helm <- model.matrix(as.formula(paste("~ ", gsub(",", " *", toString(paste("R", 1:repleng, sep = ""))), sep = "")), ortho.model)
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
                epsi.info1 <- paste("== Mauchly's Sphericity Test and Epsilons ==", sep = "")
                lamlab <- "W"
                
                # Mauchlyの球面性検定
                eps.Lambda <- det(otoM) / (sum(diag(otoM)) / (rl - 1))^(rl - 1)
                eps.m <- 1 - (2 * rl^2 - 3 * rl + 3) / (6 * (cellN - ol) * (rl - 1))
                
                if(cellN <= rl){# 被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
                    epsChi <- NA
                    epsi.info1 <- paste(epsi.info1, "\n",
                    "*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n",
                    "*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***", sep = "")
                }else{
                    epsChi <- -(cellN - ol) * eps.m * log(eps.Lambda)
                }
                
                eps.df <- (((rl - 1) * rl) / 2) - 1
                eps.p <- pchisq(epsChi, ifelse(eps.df == 0, NA, eps.df), lower.tail = FALSE)
                
            }else if(har == TRUE){
                # プロンプトの準備
                epsi.info1 <- paste("== Harris's Multisample Sphericity Test and Epsilons ==", sep = "")
                lamlab <- "h_hat"
                
                # Harrisの多標本球面性検定
                proA <- array(apply(covbuffer, 3, function(x) orthoM %*% x %*% t(orthoM)), dim = c(rl-1, rl-1, ol))
                
                if(cellN <= rl){# 被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
                    epsChi <- NA
                    eps.Lambda <- NA
                    epsi.info1 <- paste(epsi.info1, "\n",
                    "*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n",
                    "*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***", sep = "")
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
                epsi.info1 <- paste("== Mendoza's Multisample Sphericity Test and Epsilons ==", sep = "")
                lamlab <- "Lambda"
                
                # Mendozaの多標本球面性検定
                ptm <- covbuffer * rep(othN, each = rl^2)
                proA <- array(apply(ptm, 3, function(x) orthoM %*% x %*% t(orthoM)), dim = c(rl-1, rl-1, ol))
                
                if(cellN <= rl){# 被験者数が被験者内要因の組み合わせ水準数を下回るときは妥当なカイ二乗値を計算できない
                    epsChi <- NA
                    eps.Lambda <- NA
                    epsi.info1 <- paste(epsi.info1, "\n",
                    "*** CAUTION! The test of GLOBAL SPHERICITY is INVALID because of small sample size. ***", "\n",
                    "*** The minimum sample size for valid computation is N = ", rl + 1, " at each group. ***", sep = "")
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
        anova.modeler <- function(dat, design, type2 = FALSE, lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE,
        mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE, eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA,
        omega = FALSE, omegana = FALSE, pomega = FALSE, gomega = NA, gomegana = NA, prep = FALSE, inter = NA){
            # 要因計画の型に合わせてss.calc関数を適用し分散分析表を作る
            bet.with <- strsplit(design, "s")
            full.elem <- unlist(lapply(1:nchar(design), function(y) combn(c("s", LETTERS[1:(nchar(design)-1)]), y, function(x) gsub(", ", ":", toString(x)))))
            
            # 効果の自由度
            cellN <- length(unique(dat$s))
            flev <- sapply(2:nchar(design), function(x) length(unique(dat[,x])))# 各要因の水準数
            eff.df <- unlist(sapply(1:(nchar(design)-1), function(y) combn(nchar(design)-1, y, function(x) prod(flev[x]-1))))# 効果の自由度
            
            # 要因計画のタイプ別の処理
            if(nchar(bet.with[[1]][1]) != 0 && is.na(bet.with[[1]][2])){# 被験者間計画の場合
                # full.elemの整理
                er.cue <- gsub(", ", ":", toString(c("s", strsplit(strsplit(design, "s")[[1]][1], "")[[1]])))
                full.elem <- full.elem[-setdiff(grep("s", full.elem), grep(er.cue, full.elem))]
                full.elem <- full.elem[-length(full.elem)]
                
                # epsilon.calcの出力に対応する情報
                epsi.info1 <- NA
                epsitab <- NA
                ano.info1 <- NA
                
                # 誤差項の自由度
                er.df <- cellN - prod(flev[1:nchar(bet.with[[1]][1])])
                
            }else if(nchar(bet.with[[1]][1]) != 0){# 混合要因計画の場合
                # full.elemの整理
                er.cue <- gsub(", ", ":", toString(c("s", strsplit(strsplit(design, "s")[[1]][1], "")[[1]])))
                full.elem <- full.elem[-setdiff(grep("s", full.elem), grep(er.cue, full.elem))]
                
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
            eff.elem <- full.elem[-grep("s", full.elem)]# 効果の項のみ
            er.elem <- grep("s", full.elem, value = TRUE)# 誤差項のみ
            if(length(er.elem) == 0){# 被験者間計画の場合
                source.cue <- c(full.elem, "ss.Er", "ss.T")
                df.col <- c(eff.df, er.df, length(dat$y)-1)
                mse.row <- length(source.cue) - 1
            }else{# 被験者内要因を含む計画の場合
                source.ord <- sapply(eff.elem, function(x) min(elematch(x, er.elem)))
                source.cue <- unlist(sapply(1:length(er.elem), function(x) c(eff.elem[source.ord == x], er.elem[x])))
                source.cue <- c(source.cue, "ss.T")
                
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
            source.col <- gsub("ss.Er", "Error", gsub("ss.T", "Total", gsub(":", "x", source.cue)))
            ss.col <- unlist(sapply(source.cue, function(x) ss.results[names(ss.results) == x]))
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
            if(is.na(geta) == FALSE){# 一般化イータ二乗
                if(is.character(geta) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    measfact <- unique(unlist(lapply(strsplit(geta, "")[[1]], function(x) grep(x, source.col))))# 個人差変数を含む効果の取り出し
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
            if(is.na(geps) == FALSE){# 一般化イプシロン二乗；値が負になったときは０に直す
                if(is.character(geps) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    measfact <- unique(unlist(lapply(strsplit(geps, "")[[1]], function(x) grep(x, source.col))))# 個人差変数を含む効果の取り出し
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
            if(is.na(gomega) == FALSE){# 一般化オメガ二乗（加算モデル）；値が負になったときは０に直す
                if(is.character(gomega) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    measfact <- unique(unlist(lapply(strsplit(gomega, "")[[1]], function(x) grep(x, source.col))))# 個人差変数を含む効果の取り出し
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
            if(is.na(gomegana) == FALSE){# 一般化オメガ二乗（非加算モデル）；値が負になったときは０に直す
                if(length(mse.row) == 1){# 被験者間計画の場合
                    omega.dummy <- cellN
                }else{# その他の計画の場合
                    dflev <- flev[(nchar(bet.with[[1]][1]) + 1):length(flev)]
                    omega.dummy <- cellN * c(1, unlist(sapply(1:length(dflev), function(y) combn(1:length(dflev), y, function(x) prod(dflev[x])))))
                }
                if(is.character(gomegana) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                    measfact <- unique(unlist(lapply(strsplit(gomegana, "")[[1]], function(x) grep(x, source.col))))# 個人差変数を含む効果の取り出し
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
            if(is.na(inter)){
                return(list("epsi.info1" = epsi.info1, "epsitab" = epsitab, "ano.info1" = ano.info1, "anovatab" = anovatab, "mse.row" = mse.row))
            }else{
                # 単純主効果の検定用の出力を用意する
                sim.row <- charmatch(inter, anovatab$source.col)
                intertab <- rbind(anovatab[sim.row,], anovatab[f.denomi[sim.row],])
                
                # 球面性検定の結果からinterに関する部分のみ取り出す
                if(charmatch(inter, strsplit(design, "")[[1]]) < charmatch("s", strsplit(design, "")[[1]])){
                    # 被験者間要因ならNAの行を返す
                    if(iga == TRUE | ciga == TRUE) interepsi <- rep(NA, 11)# IGA，CIGAを使ったときは列数が多い
                    else interepsi <- rep(NA, 9)
                }else{
                    interepsi <- epsitab[charmatch(inter, epsitab$Effect),]
                }
                return(list("intertab" = intertab, "interepsi" = interepsi))
            }
        }
        
        
        # 修正Bonferroniの方法（Holmの方法，Shafferの方法，Holland-Copenhaverの方法）による多重比較を行う関数
        # デフォルトはShafferの方法を適用し，holm = TとするとHolmの方法，hc = TとするとHolland-Copenhaverの方法を適用する
        # s2r = T，s2d = Tとすると，具体的な棄却のパターンを反映した有意水準の調整によるShafferの方法を適用する
        mod.Bon <- function(dat, design, taref, bet.mse, factlabel = NA, type2 = FALSE, holm = FALSE, hc = FALSE,
        s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, alpha = 0.05, criteria = FALSE){
            # 対象となる要因のラベルを得る
            if(is.na(factlabel)) factlabel <- taref
            
            bonl <- nlevels(dat[, taref])# 水準数の取得
            h0size <- bonl * (bonl-1)/2# 帰無仮説の個数
            bon.name <- combn(bonl, 2, function(ij) paste(levels(dat[, taref])[ij][1], "-", levels(dat[, taref])[ij][2], sep = ""))# 可能な対の組み合わせのラベル
            bon.num <- charmatch(taref, names(dat))-1# 分析対象となる要因が何番目の要因かを特定
            
            # 周辺平均を計算する
            cont.means <- tapply(dat$y, dat[, 2:nchar(design)], mean)# 各セルの平均を求める
            bon.means <- apply(cont.means, bon.num, mean)# 分析対象となる周辺平均を求める
            
            factlevels <- sapply(names(dat), function(x) nlevels(dat[,x]))# 各要因の水準数
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
                subdat <- apply(bon.lev, 2, function(x) dat[dat[, taref] == x[1] | dat[, taref] == x[2],])# データを多重比較の対象となる効果について２水準ずつのサブデータにリスト化
                for(i in 1:length(subdat)){
                    subdat[[i]][, taref] <- subdat[[i]][, taref][, drop = TRUE]
                }
                
                bon.anova <- lapply(subdat, function(x) anova.modeler(dat = x, design = design, type2 = type2,
                inter = taref)$intertab[2,])
                bon.df <- sapply(bon.anova, function(x) x$df.col)# 自由度
                bon.Ve <- sapply(bon.anova, function(x) x$ms.col)
            }
            
            # 検定統計量とｐ値を得る
            bon.SE <- sqrt(combn(bonl, 2, function(ij) sum(bon.denomi[ij])) * bon.Ve)
            bon.t <- abs(bon.delta / bon.SE)# ｔ値は絶対値を取る
            bon.p <- pt(bon.t, bon.df, lower.tail = FALSE) * 2# 両側確率
            
            # 結果をデータフレームにまとめる
            bontab <- data.frame("pair" = bon.name, "interval" = bon.delta, "t" = bon.t, "df" = bon.df, "p.value" = bon.p)
            bontab <- bontab[order(bontab$p.value),]# p値の小さい順に並べ替え
            
            # 調整した有意水準を設定する
            if(holm == TRUE){
                p.criteria <- h0size:1# Holmの方法用の調整値
                bon.info1 <- paste("== Holm's Sequentially Rejective Bonferroni Procedure ==", sep = "")
                
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
                    bon.info1 <- paste("== Shaffer's Modified Sequentially Rejective Bonferroni Procedure [SPECIFIC] ==", "\n",
                    "== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
                    shaf.meth <- paste(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
                }else{
                    p.criteria <- c(shaf.value[2], shaf.value[2:length(shaf.value)])# F-Shafferの方法用の調整値
                    bon.info1 <- paste("== Shaffer's F-Modified Sequentially Rejective Bonferroni Procedure [SPECIFIC] ==", "\n",
                    "== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
                    shaf.meth <- paste(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Donoghue (2004). ==", sep = "")
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
                    shaf.meth <- paste(" [SPECIFIC] ==", "\n", "== This computation is based on the algorithm by Rasmussen (1993). ==", sep = "")
                }else{# 各比較までの任意の仮説が偽であったときの可能な真の帰無仮説の最大数
                    shaf.value <- sapply(1:h0size, function(x) max(t.vector[i.vector >= (x - 1)]))
                    shaf.meth <- " =="
                }
                
                if(fs1 | fs2r == TRUE){
                    p.criteria <- c(shaf.value[2], shaf.value[2:length(shaf.value)])# F-Shafferの方法用の調整値
                    bon.info1 <- paste("== Shaffer's F-Modified Sequentially Rejective Bonferroni Procedure", shaf.meth, sep = "")
                }else{
                    p.criteria <- shaf.value# Shafferの方法用の調整値
                    bon.info1 <- paste("== Shaffer's Modified Sequentially Rejective Bonferroni Procedure", shaf.meth, sep = "")
                }
            }
            
            # 平均値の差の方向を調べ，不等号のベクトルを作る
            bon.differ <- ifelse(bontab$interval <= 0, sub("-", " < ", bontab$pair), sub("-", " > ", bontab$pair))
            # 差が見られなかった場合の等号のベクトルを作る
            bon.equal <- sub("-", " = ", bontab$"pair")
            
            if(criteria == TRUE){# データフレームに調整済み有意水準の列を加える
                if(hc == TRUE){bontab <- transform(bontab, "criteria" = 1 - (1 - alpha) ^ (1/p.criteria))# Sidakの不等式による有意水準の調整
                    if(holm == TRUE) bon.info1 <- paste("== Holm's Sequentially Rejective Sidak Procedure ==", sep = "")
                    else bon.info1 <- paste("== Holland-Copenhaver's Improved Sequentially Rejective Sidak Procedure", shaf.meth, sep = "")
                    
                    if(length(bet.mse) == 1) bon.info1 <- append(bon.info1, "*** CAUTION! This procedure might be inappropriate for dependent means. ***")
                }else{bontab <- transform(bontab, "criteria" = alpha/p.criteria)# Bonferroniの不等式による有意水準の調整
                }
                # 有意であった行は不等号，そうでない行は等号を表示する
                bon.sign <- ifelse(cummin(bontab$p.value < bontab$criteria), paste(bon.differ, "*", sep = " "), paste(bon.equal, " ", sep = " "))
            }else{# データフレームに調整済みｐ値の列を加える
                if(hc == TRUE){bontab <- transform(bontab, "adj.p" = pmin(1, cummax((1-(1-bontab$p.value)^p.criteria))))# Sidakの不等式による調整済みｐ値
                    if(holm == TRUE) bon.info1 <- paste("== Holm's Sequentially Rejective Sidak Procedure ==", sep = "")
                    else bon.info1 <- paste("== Holland-Copenhaver's Improved Sequentially Rejective Sidak Procedure", shaf.meth, sep = "")
                    
                    if(length(bet.mse) == 1) bon.info1 <- append(bon.info1, "*** CAUTION! This procedure might be inappropriate for dependent means. ***")
                }else{bontab <- transform(bontab, "adj.p" = pmin(1, cummax(bontab$p.value * p.criteria)))# Bonferroniの不等式による調整済みｐ値
                }
                # 有意であった行は不等号，そうでない行は等号を表示する
                bon.sign <- ifelse(bontab$adj.p < alpha, paste(bon.differ, "*", sep = " "), paste(bon.equal, " ", sep = " "))
            }
            
            # 判定結果をデータフレームに反映する
            bontab <- transform(bontab, "significance" = bon.sign)
            
            # 記述統計量の計算
            b.sncol <- tapply(dat$y, dat[,taref], length)# セルごとのデータ数を計算
            b.sdcol <- tapply(dat$y, dat[,taref], sd)# セルごとの標準偏差を計算
            
            bonstat <- data.frame("Dummy" = levels(dat[, taref]), "N" = b.sncol, "Mean" = bon.means, "S.D." = b.sdcol)
            
            names(bonstat)[1] <- factlabel# 水準を表すラベルをfactlabelとして入力した値に置き換える
            
            # その他の出力を準備する
            bon.info2 <- if(length(bet.mse) != 1) paste("== The factor < ", factlabel, " > is analysed as independent means. ==", sep  = "")
            else paste("== The factor < ", factlabel, " > is analysed as dependent means. ==", sep  = "")
            bon.info3 <- paste("== Alpha level is ", alpha, ". ==", sep  = "")
            
            return(list(factlabel, bon.info1, bon.info2, bon.info3, bonstat, bontab))
            
        }
        
        
        # 下位検定を行う関数
        post.analyses <- function(dat, design, mainresults, type2 = FALSE, nopost = FALSE, holm = FALSE, hc = FALSE,
        s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE,
        lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE,
        eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE,
        gomega = NA, gomegana = NA, prep = FALSE){
            anovatab <- mainresults$anovatab
            
            # 要因計画の型から被験者間要因と被験者内要因の情報を得る
            bet.with <- strsplit(design, "s")
            # 効果が有意であった行のソースラベルを得る
            sig.source <- anovatab$"source.col"[!((anovatab$"sig.col" == "") | (anovatab$"sig.col" == "ns"))]
            
            if(length(sig.source) == 0 | nopost == TRUE) {
                return(NA)# 有意な行が存在しないか，nopostオプションが指定されている場合はここで終了
            }else{
                # 下位検定の結果を格納するための空のデータフレームを宣言
                postresults <- lapply(paste("post", 1:length(sig.source), sep = ""), function(x) x)
                
                # pro.fraction関数を反復適用
                for (i in 1:length(sig.source)) postresults[[i]] <- pro.fraction(dat = dat, design = design, postplan = sig.source[i],
                bet.with = bet.with, mainresults = mainresults, type2 = type2, holm = holm, hc = hc,
                s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria,
                lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga, eta = eta, peta = peta, geta = geta,
                eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega, gomega = gomega, gomegana = gomegana,
                prep = prep)
                
                return(postresults)
            }
        }
        
        
        # 効果のタイプに適した下位検定を割り当てる関数
        pro.fraction <- function(dat, design, postplan, bet.with, mainresults, type2 = FALSE, holm = FALSE, hc = FALSE,
        s2r = FALSE, s2d = FALSE, fs1 = FALSE, fs2r = FALSE, fs2d = FALSE, criteria = FALSE,
        lb = FALSE, gg = FALSE, hf = FALSE, auto = FALSE, mau = FALSE, har = FALSE, iga = FALSE, ciga = FALSE,
        eta = FALSE, peta = FALSE, geta = NA, eps = FALSE, peps = FALSE, geps = NA, omega = FALSE, omegana = FALSE, pomega = FALSE,
        gomega = NA, gomegana = NA, prep = FALSE){
            # 情報の展開
            anovatab <- mainresults$anovatab
            mse.row <- mainresults$mse.row
            
            # ソースラベルを文字列に変換し，その文字数を得る
            sig.term <- as.character(postplan)
            sig.num <- nchar(sig.term)
            
            # 効果の種類によって違った処理を割り当てる
            if(sig.num > 3){
                # 高次の交互作用：効果のラベルを返す
                return(sig.term)
                
            }else if(sig.num == 3){
                # １次の交互作用については，単純主効果の検討を行う
                each.term <- strsplit(sig.term, "x")# 交互作用の各項を分離する
                interfact <- as.numeric(sapply(each.term[[1]], function(x) grep(x, bet.with[[1]][1])))# 被験者間要因なら１，被験者内要因なら０を返す
                
                # 交互作用を構成する第一項の項の列番号，水準数を取得
                col.num1 <- charmatch(each.term[[1]][1], names(dat))# 列番号
                level.num1 <- nlevels(dat[, col.num1])# 水準数
                
                # 交互作用を構成する第二項の項の列番号，水準数を取得
                col.num2 <- charmatch(each.term[[1]][2], names(dat))# 列番号
                level.num2 <- nlevels(dat[, col.num2])# 水準数
                
                # 単純主効果の検定用の要因計画を作る
                resdesign1 <- sub(each.term[[1]][2], "", design)# つぶす要因のラベルを消去
                bet.with1 <- strsplit(resdesign1, "s")# 残りの要因を被験者間と被験者内に分離
                bet.num1 <- if(nchar(bet.with1[[1]][1]) == 0) 0 else 1:nchar(bet.with1[[1]][1])# 被験者間要因の数を数える
                with.num1 <- if(is.na(bet.with1[[1]][2])) 0 else (max(bet.num1)+1):(nchar(resdesign1)-1)# 被験者内要因の数を数える
                bet1 <- gsub(", ", "", toString(LETTERS[bet.num1]))# 被験者間要因のラベルを新たに作成
                with1 <- gsub(", ", "", toString(LETTERS[with.num1]))# 被験者内要因のラベルを新たに作成
                postdesign1 <- paste(bet1, "s", with1, sep = "")# 両要因を合成
                
                resdesign2 <- sub(each.term[[1]][1], "", design)# つぶす要因のラベルを消去
                bet.with2 <- strsplit(resdesign2, "s")# 残りの要因を被験者間と被験者内に分離
                bet.num2 <- if(nchar(bet.with2[[1]][1]) == 0) 0 else 1:nchar(bet.with2[[1]][1])# 被験者間要因の数を数える
                with.num2 <- if(is.na(bet.with2[[1]][2])) 0 else (max(bet.num2)+1):(nchar(resdesign2)-1)# 被験者内要因の数を数える
                bet2 <- gsub(", ", "", toString(LETTERS[bet.num2]))# 被験者間要因のラベルを新たに作成
                with2 <- gsub(", ", "", toString(LETTERS[with.num2]))# 被験者内要因のラベルを新たに作成
                postdesign2 <- paste(bet2, "s", with2, sep = "")# 両要因を合成
                
                # 分析対象となる効果が被験者間要因か被験者内要因かによってソース列のラベルを作成
                if(is.na(interfact[1]) == FALSE) err.label1 <- "Er"
                else err.label1 <- paste("sx", each.term[[1]][1], sep = "")
                
                if(is.na(interfact[2]) == FALSE) err.label2 <- "Er"
                else err.label2 <- paste("sx", each.term[[1]][2], sep = "")
                
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
                sim.effects1 <- lapply(subdat1, function(x) anova.modeler(dat = x, design = postdesign1, type2 = type2,
                lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
                eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
                gomega = gomega, gomegana = gomegana, prep = prep, inter = target.effect1))
                sim.effects2 <- lapply(subdat2, function(x) anova.modeler(dat = x, design = postdesign2, type2 = type2,
                lb = lb, gg = gg, hf = hf, auto = auto, mau = mau, har = har, iga = iga, ciga = ciga,
                eta = eta, peta = peta, geta = geta, eps = eps, peps = peps, geps = geps, omega = omega, omegana = omegana, pomega = pomega,
                gomega = gomega, gomegana = gomegana, prep = prep, inter = target.effect2))
                
                # 結果を１つのデータフレームにまとめる
                simtab <- data.frame()# 結果格納用データフレーム
                simepsi <- data.frame()# 結果格納用データフレーム
                subsource <- c()# ソースラベル格納用ベクトル
                
                # 計画のタイプによって出力を整理し直す
                if(substr(design, nchar(design), nchar(design)) == "s"){
                    # 被験者間計画なら，anovatabから主分析の誤差平方和等を得る
                    between.mse <- anovatab[charmatch("Error", anovatab$source.col),]
                    bet.mse1 <- between.mse
                    bet.mse2 <- between.mse
                    
                    # 主分析の誤差平方和を用いて単純主効果の結果を再計算
                    for (i in 1:level.num2){
                        simtab <- rbind(simtab, sim.effects1[[i]]$intertab[1,])
                        subsource <- append(subsource, c(paste(each.term[[1]][1], " at ", names(sim.effects1)[i], sep = "")))
                    }
                    for (i in 1:level.num1){
                        simtab <- rbind(simtab, sim.effects2[[i]]$intertab[1,])
                        subsource <- append(subsource, c(paste(each.term[[1]][2], " at ", names(sim.effects2)[i], sep = "")))
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
                    if(is.na(geta) == FALSE){# 一般化イータ二乗
                        if(is.character(geta) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            measfact <- unique(unlist(lapply(strsplit(geta, "")[[1]], function(x) grep(x, anovatab$source.col))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)])# 誤差平方和を除く
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(strsplit(geta, "")[[1]], function(x) grep(x, subsource))))# 個人差変数に含まれる単純主効果の取り出し
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
                    if(is.na(geps) == FALSE){# 一般化イプシロン二乗；値が負になったときは０に直す
                        if(is.character(geps) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            measfact <- unique(unlist(lapply(strsplit(geps, "")[[1]], function(x) grep(x, anovatab$source.col))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)])# 誤差平方和を除く
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(strsplit(geps, "")[[1]], function(x) grep(x, subsource))))# 個人差変数に含まれる単純主効果の取り出し
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
                    if(is.na(gomega) == FALSE){# 一般化オメガ二乗（加算モデル）；値が負になったときは０に直す
                        ms.copy <- anovatab$ms.col[nrow(anovatab)-1]
                        if(is.character(gomega) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            measfact <- unique(unlist(lapply(strsplit(gomega, "")[[1]], function(x) grep(x, anovatab$source.col))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)] - anovatab$df.col[setdiff(measfact, mse.row)] * ms.copy)
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(strsplit(gomega, "")[[1]], function(x) grep(x, subsource))))# 個人差変数に含まれる単純主効果の取り出し
                            measvec[pmeasfact] <- 0
                        }else{
                            ss.meas <- 0
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                        }
                        gomega.col <- pmax(0, (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * ms.copy) /
                        (measvec * (simtab$ss.col[1:(nrow(simtab)-1)] - simtab$df.col[1:(nrow(simtab)-1)] * ms.copy) + ss.meas + nrow(dat) * simtab$ms.col[nrow(simtab)]))
                        simtab <- cbind(simtab, "G.omega^2" = c(gomega.col, NA))
                    }
                    if(is.na(gomegana) == FALSE){# 一般化オメガ二乗（非加算モデル）；値が負になったときは０に直す
                        ms.copy <- anovatab$ms.col[nrow(anovatab)-1]
                        if(is.character(gomega) == TRUE){# 被験者間要因の中に測定変数（個人差変数）がある場合
                            measfact <- unique(unlist(lapply(strsplit(gomegana, "")[[1]], function(x) grep(x, anovatab$source.col))))# 個人差変数を含む効果の取り出し
                            ss.meas <- sum(anovatab$ss.col[setdiff(measfact, mse.row)] - anovatab$df.col[setdiff(measfact, mse.row)] * ms.copy)
                            measvec <- rep(1, length(simtab$ss.col[1:(nrow(simtab)-1)]))
                            pmeasfact <- unique(unlist(lapply(strsplit(gomegana, "")[[1]], function(x) grep(x, subsource))))# 個人差変数に含まれる単純主効果の取り出し
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
                        subsource <- append(subsource, c(paste(each.term[[1]][1], " at ", names(sim.effects1)[i], sep = ""),
                        paste(err.label1, " at ", names(sim.effects1)[i], sep = "")))
                    }
                    names(simepsi) <- names(sim.effects2[[1]]$interepsi)# 名前が異なるとrbindできないので統一；混合要因計画の場合，必ず被験者間要因がsim.effects1に入ることを前提とする
                    for (i in 1:level.num1){
                        simtab <- rbind(simtab, sim.effects2[[i]]$intertab)
                        simepsi <- rbind(simepsi, sim.effects2[[i]]$interepsi)
                        subsource <- append(subsource, c(paste(each.term[[1]][2], " at ", names(sim.effects2)[i], sep = ""),
                        paste(err.label2, " at ", names(sim.effects2)[i], sep = "")))
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
                sim.sncol <- as.vector(tapply(dat$y, list(dat[,col.num2], dat[,col.num1]), length))# セルごとのデータ数を計算
                sim.mncol <- as.vector(tapply(dat$y, list(dat[,col.num2], dat[,col.num1]), mean))# セルごとの平均を計算
                sim.sdcol <- as.vector(tapply(dat$y, list(dat[,col.num2], dat[,col.num1]), sd))# セルごとの標準偏差を計算
                
                sim.stat <- data.frame("Term1" = rep(levels(dat[,col.num1]), each = level.num2),
                "Term2" = rep(levels(dat[,col.num2]), level.num1), "N" = sim.sncol, "Mean" = sim.mncol, "S.D." = sim.sdcol)
                
                # 行ラベルの張り替え
                names(sim.stat)[1] <- each.term[[1]][1]
                names(sim.stat)[2] <- each.term[[1]][2]
                
                # 有意であった行のソースをチェック
                sim.sig.col <- simtab$source.col[!((simtab$sig.col == "") | (simtab$sig.col == "ns"))]
                
                # 多重比較を実行
                if(length(sim.sig.col) == 0){
                    sim.multresults <- NA# 有意な行がない場合はNAを代入
                }else{
                    # 多重比較の結果を格納するための空のデータフレームを宣言
                    sim.multresults <- sapply(paste("mult", 1:length(sim.sig.col), sep = ""), function(x) list(x))
                    
                    # 多重比較の関数を反復適用
                    for (i in 1:length(sim.sig.col)){
                        # ソースラベルを分解
                        source.term <- strsplit(sim.sig.col[i], " at ")
                        
                        if((source.term[[1]][1] == each.term[[1]][1]) && (level.num1 >= 3) == TRUE){
                            # 第一項の単純主効果が有意；subdat1を使用
                            multdat <- subdat1[[source.term[[1]][2]]]# 対象となる水準のデータのみ抽出
                            # 混合要因計画の被験者間要因の場合は複数の行から適切なMSeを選択；それ以外の場合は単一の行をコピー
                            if(!is.na(bet.mse1) && nrow(bet.mse1) != 1) bet.mse <- bet.mse1[grep(source.term[[1]][2], bet.mse1$source.col),]
                            else bet.mse <- bet.mse1
                            sim.multresults[[i]] <- mod.Bon(dat = multdat, design = postdesign1, taref = target.effect1, bet.mse = bet.mse, factlabel = sim.sig.col[i],
                            type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
                        }else if((source.term[[1]][1] == each.term[[1]][2]) && (level.num2 >= 3) == TRUE){
                            # 第二項の単純主効果が有意；subdat2を使用
                            multdat <- subdat2[[source.term[[1]][2]]]# 対象となる水準のデータのみ抽出
                            # 混合要因計画の被験者間要因の場合は複数の行から適切なMSeを選択；それ以外の場合は単一の行をコピー
                            if(!is.na(bet.mse2) && nrow(bet.mse2) != 1) bet.mse <- bet.mse2[grep(source.term[[1]][2], bet.mse2$source.col),]
                            else bet.mse <- bet.mse2
                            sim.multresults[[i]] <- mod.Bon(dat = multdat, design = postdesign2, taref = target.effect2, bet.mse = bet.mse, factlabel = sim.sig.col[i],
                            type2 = type2, holm = holm, hc = hc, s2r = s2r, s2d = s2d, fs1 = fs1, fs2r = fs2r, fs2d = fs2d, criteria = criteria)
                        }else{sim.multresults[[i]] <- NA}
                    }
                }
                
                return(list("sig.term" = sig.term, "sim.stat" = sim.stat, "simepsi" = simepsi, "simtab" = simtab, "sim.multresults" = sim.multresults))
                
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
                            bet.mse <- anovatab[charmatch("Error", anovatab$source.col),]
                        }else{
                            # 混合要因計画ならその効果の誤差項のMSを得る
                            f.denomi <- rep(mse.row, c(mse.row[1], diff(mse.row)))
                            bet.mse <- anovatab[f.denomi[charmatch(sig.term, anovatab$source.col)],]
                        }
                    }else{
                        bet.mse <- NA
                    }
                    
                    bonout <- mod.Bon(dat = dat, design = design, taref = sig.term, bet.mse = bet.mse, factlabel = NA,
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
        bstat.out <- function(bstatist, maxfact, margin, copy = FALSE){
            #インターフェイスによってラインの長さを調整
            if(.Platform$GUI == "Rgui" | .Platform$GUI == "X11" | copy == TRUE){
                btablength = 4
            }else{
                btablength = 1
            }
            
            tabline <- sapply(1:ncol(bstatist), function(x) ifelse(is.numeric(bstatist[,x]), max(nchar(round(bstatist[,x], 4), type = "width")), max(nchar(as.character(bstatist[,x]), type = "width")))) + 2
            mainline <- sum(tabline) + btablength * length(tabline)
            headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
            bodymargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
            bodymargin[(length(tabline)-1):length(tabline)] <- paste("%", tabline[(length(tabline)-1):length(tabline)], ".4f", sep = "")
            
            cat(sprintf("<< DESCRIPTIVE STATISTICS >>"), sep = "", "\n", "\n")# タイトル
            cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")# ラインを引く；要因の数に合わせて長さを調整
            cat(sprintf(headmargin, names(bstatist)), sep = "\t", "\n")# データフレームの列名をタブ区切りで表示する
            cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")
            
            # データフレームのデータ部分を一行ずつタブ区切りで表示する
            for (i in 1:nrow(bstatist)){
                if((i %% margin == 0) && !(i == nrow(bstatist))){
                    cat(sprintf(bodymargin[1:maxfact], unlist(bstatist[i, 1:maxfact])), sprintf(bodymargin[maxfact+1], paste(bstatist[i, maxfact+1])), sprintf(bodymargin[maxfact+2:3], bstatist[i, maxfact+2:3]), sep = "\t", "\n")
                    cat("\n")
                }else{
                    cat(sprintf(bodymargin[1:maxfact], unlist(bstatist[i, 1:maxfact])), sprintf(bodymargin[maxfact+1], paste(bstatist[i, maxfact+1])), sprintf(bodymargin[maxfact+2:3], bstatist[i, maxfact+2:3]), sep = "\t", "\n")
                }
            }
            cat(sprintf(rep("-", mainline + 2)), sep = "", "\n", "\n", "\n")
            
        }
        
        
        # 分散分析表を出力する関数
        anovatab.out <- function(mainresults, copy = FALSE){
            # 情報を展開
            epsi.info1 <- mainresults$epsi.info1
            epsitab <- mainresults$epsitab
            ano.info1 <- mainresults$ano.info1
            anovatab <- mainresults$anovatab
            mse.row <- mainresults$mse.row
            
            #インターフェイスによってラインの長さを調整
            if(.Platform$GUI == "Rgui" | .Platform$GUI == "X11" | copy == TRUE){
                tablength = 6
                stablength = 3
            }else{
                tablength = 1
                stablength = 1
            }
            
            # 球面性検定の結果を出力
            if(is.na(epsi.info1) == FALSE){# epsi.info1がNAでないときのみ出力
                tabline <- sapply(1:ncol(epsitab), function(x) ifelse(is.numeric(epsitab[,x]), max(nchar(round(epsitab[,x], 4), type = "width")), max(nchar(as.character(epsitab[,x]), type = "width")))) + 1
                tabline <- ifelse(nchar(names(epsitab)) > tabline, nchar(names(epsitab)) + 1, tabline)
                tabline[5] <- 7
                tabline[6] <- 3
                tabline[7:ncol(epsitab)] <- 6
                mainline <- sum(tabline) + stablength * length(tabline)
                headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
                bodymargin <- sapply(tabline, function(x) paste("%", x, ".4f", sep = ""))
                bodymargin[1] <- paste("%", tabline[1], "s", sep = "")
                bodymargin[4] <- paste("%", tabline[4], "s", sep = "")
                bodymargin[6] <- "%-3s"
                
                cat(sprintf("<< SPHERICITY INDICES >>"), sep = "", "\n", "\n")# タイトル
                cat(sprintf(epsi.info1), "\n", "\n")# プロンプト
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1:5], names(epsitab[1:5])), sprintf(headmargin[6], paste("")), sprintf(headmargin[7:ncol(epsitab)], names(epsitab[7:ncol(epsitab)])), sep = "\t", "\n")# 列名をタブ区切りで出力
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する
                for (i in 1:nrow(epsitab)){
                    cat(sprintf(bodymargin[1], paste(epsitab[i, 1])), sprintf(bodymargin[2], epsitab[i, 2]), sprintf(bodymargin[3], epsitab[i, 3]),
                    sprintf(bodymargin[4], epsitab[i, 4]), replace(sprintf(bodymargin[5], epsitab[i, 5]), is.na(epsitab[i, 5]), sprintf(paste("%", tabline[5], "s", sep = ""), "")),
                    replace(sprintf(bodymargin[6], paste(epsitab[i, 6])), epsitab[i, 6] == "", sprintf(bodymargin[6], "")), sprintf(bodymargin[7], epsitab[i, 7:ncol(epsitab)]), sep = "\t", "\n")
                }
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                if(ncol(epsitab) == 9){
                    cat(sprintf(paste("%", mainline - 1, "s", sep = ""), "LB = lower.bound, GG = Greenhouse-Geisser, HF = Huynh-Feldt"), sep = "\t", "\n", "\n", "\n")
                }else{
                    cat(sprintf(paste(rep("", 3), sep = "")), sprintf(paste("%", mainline - 1, "s", sep = ""), "b_hat = b^, c_hat = c^, h_d = h~', h_dd = h~'', h = h~"), sep = "\t", "\n", "\n", "\n")
                }
            }
            
            # 分散分析表を出力
            cat(sprintf("<< ANOVA TABLE >>"), sep = "", "\n", "\n")# タイトル
            if(sum(is.na(ano.info1)) != 1){# ano.info1がNAでないときのみ表示
                cat(sprintf(ano.info1), sep = "\n")# プロンプトを表示
                cat("\n")# １行空ける
            }
            
            tabline <- sapply(1:ncol(anovatab), function(x) ifelse(is.numeric(anovatab[,x]), max(nchar(round(anovatab[,x]), type = "width")) + 5, max(nchar(as.character(anovatab[,x]), type = "width")))) + 2
            tabline[3] <- max(nchar(round(anovatab[,3], 4), type = "width"))
            mainline <- sum(tabline) + tablength * length(tabline)
            headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
            bodymargin <- sapply(tabline, function(x) paste("%", x, ".4f", sep = ""))
            bodymargin[1] <- paste("%", tabline[1], "s", sep = "")
            bodymargin[3] <- paste("%", tabline[3], "s", sep = "")
            bodymargin[7] <- paste("%-", tabline[7], "s", sep = "")
            
            esn <- ncol(anovatab) - 7
            if(esn == 0){# 分散分析表のみの場合
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"),
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")), sep = "\t", "\n")# 列名をタブ区切りで表示する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する；誤差項の行（mse.row）の後にラインを入れる
                for (i in 1:nrow(anovatab)){
                    if(i == nrow(anovatab)){# 最終行に到達したときには以下を実行
                        cat(sprintf(bodymargin[1], paste(anovatab[i, 1])), sprintf(bodymargin[2], anovatab[i, 2]),
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sep = "\t", "\n")
                        cat(sprintf(paste("%", mainline - 2, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n")
                    }else if(sum(mse.row == i) == 1){# mse.rowの中に一致する値があるときには以下を実行
                        cat(sprintf(bodymargin[1], paste(anovatab[i, 1])), sprintf(bodymargin[2], anovatab[i, 2]),
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]),
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""),
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""),
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), sep = "\t", "\n")
                        cat(sprintf(rep("-", mainline)), sep = "", "\n")
                    }else{# その他のときは以下を実行
                        cat(sprintf(bodymargin[1], paste(anovatab[i, 1])), sprintf(bodymargin[2], anovatab[i, 2]),
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]),
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""),
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""),
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), sep = "\t", "\n")
                    }
                }
            }else{# 効果量の指標を追加した場合
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"),
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")),
                sprintf(headmargin[8:length(headmargin)], names(anovatab[8:ncol(anovatab)])), sep = "\t", "\n")# 列名をタブ区切りで表示する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する；誤差項の行（mse.row）の後にラインを入れる
                for (i in 1:nrow(anovatab)){
                    if(i == nrow(anovatab)){# 最終行に到達したときには以下を実行
                        cat(sprintf(bodymargin[1], paste(anovatab[i, 1])), sprintf(bodymargin[2], anovatab[i, 2]),
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), "\n", sprintf(paste("%", mainline - 2, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n")
                    }else if(sum(mse.row == i) == 1){# mse.rowの中に一致する値があるときには以下を実行
                        cat(sprintf(bodymargin[1], paste(anovatab[i, 1])), sprintf(bodymargin[2], anovatab[i, 2]),
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]),
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""),
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""),
                        sprintf(bodymargin[7], paste(anovatab[i, 7])), sep = "\t", "\n")
                        cat(sprintf(rep("-", mainline)), sep = "", "\n")
                    }else{# その他のときは以下を実行
                        cat(sprintf(bodymargin[1], paste(anovatab[i, 1])), sprintf(bodymargin[2], anovatab[i, 2]),
                        sprintf(bodymargin[3], round(anovatab[i, 3], 2)), sprintf(bodymargin[4], anovatab[i, 4]),
                        replace(sprintf(bodymargin[5], anovatab[i, 5]), is.na(anovatab[i, 5]), ""),
                        replace(sprintf(bodymargin[6], anovatab[i, 6]), is.na(anovatab[i, 6]), ""),
                        sprintf(bodymargin[7], paste(anovatab[i, 7])),
                        replace(sprintf(bodymargin[8:length(bodymargin)], anovatab[i, 8:ncol(anovatab)]), is.na(anovatab[i, 8:ncol(anovatab)]), ""), sep = "\t", "\n")
                    }
                }
            }
            
            cat("\n", "\n")# ２行空ける
            
        }
        
        
        # 修正Bonferroniの方法による多重比較の結果を出力する関数
        mod.Bon.out <- function(bon.list, omit = FALSE, copy = FALSE){
            # 情報を展開
            factlabel <- bon.list[[1]]
            bon.info1 <- bon.list[[2]]
            bon.info2 <- bon.list[[3]]
            bon.info3 <- bon.list[[4]]
            bonstat <- bon.list[[5]]
            bontab <- bon.list[[6]]
            
            #インターフェイスによってラインの長さを調整
            if(.Platform$GUI == "Rgui" | .Platform$GUI == "X11" | copy == TRUE){
                tablength = 6
                btablength = 4
            }else{
                tablength = 1
                btablength = 1
            }
            
            cat(sprintf(paste("< MULTIPLE COMPARISON for FACTOR ", factlabel, " >" ,sep = "")), sep = "", "\n", "\n")# タイトル
            cat(sprintf(bon.info1), sep = "\n")# プロンプト
            cat(sprintf(bon.info2), "\n")
            cat(sprintf(bon.info3), "\n", "\n")
            
            # omitがFなら記述統計量を出力する
            if(omit == FALSE){
                # 記述統計量を出力する
                tabline <- sapply(1:ncol(bonstat), function(x) ifelse(is.numeric(bonstat[,x]), max(nchar(round(bonstat[,x], 4), type = "width")), max(nchar(as.character(bonstat[,x]), type = "width")))) + 2
                mainline <- sum(tabline) + btablength * length(tabline)
                headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
                bodymargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
                bodymargin[(length(tabline)-1):length(tabline)] <- paste("%", tabline[(length(tabline)-1):length(tabline)], ".4f", sep = "")
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")# ラインを引く；要因の数に合わせて長さを調整
                cat(sprintf(headmargin, names(bonstat)), sep = "\t", "\n")# データフレームの列名をタブ区切りで表示する
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する
                for (i in 1:nrow(bonstat)){
                    cat(sprintf(bodymargin[1:2], paste(bonstat[i, 1:2])), sprintf(bodymargin[3:4], bonstat[i, 3:4]), sep = "\t", "\n")
                }
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n", "\n")
                
            }
            
            # 多重比較の結果を表形式で出力する
            tabline <- sapply(1:ncol(bontab), function(x) ifelse(is.numeric(bontab[,x]), max(nchar(round(bontab[,x], 4), type = "width")), max(nchar(as.character(bontab[,x]), type = "width")))) + 3
            mainline <- sum(tabline) + tablength * length(tabline)
            headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
            bodymargin <- sapply(tabline, function(x) paste("%", x, ".4f", sep = ""))
            bodymargin[1] <- paste("%", tabline[1], "s", sep = "")
            bodymargin[4] <- paste("%", tabline[4], "s", sep = "")
            bodymargin[7] <- paste("%", tabline[7], "s", sep = "")
            
            cat(sprintf(rep("-", mainline)), sep = "", "\n")
            cat(sprintf(headmargin[1], "Pair"), sprintf(headmargin[2], paste("Interval")), sprintf(headmargin[3], "t-value"), sprintf(headmargin[4], "df"),
            sprintf(headmargin[5], "p"), sprintf(headmargin[6], names(bontab[6])), sprintf(headmargin[7], paste("")), sep = "\t", "\n")# 列名をタブ区切りで表示する
            cat(sprintf(rep("-", mainline)), sep = "", "\n")
            
            # データフレームのデータ部分を一行ずつタブ区切りで表示する
            for (i in 1:nrow(bontab)){
                cat(sprintf(bodymargin[1], paste(bontab[i, 1])), sprintf(bodymargin[2], bontab[i, 2]), sprintf(bodymargin[3], bontab[i, 3]),
                sprintf(bodymargin[4], paste(" ", bontab[i, 4])), sprintf(bodymargin[5], bontab[i, 5]), sprintf(bodymargin[6], bontab[i, 6]),
                sprintf(bodymargin[7], paste(bontab[i, 7])), sep = "\t", "\n")
            }
            
            cat(sprintf(rep("-", mainline)), sep = "", "\n", "\n")
            
        }
        
        
        # 単純主効果の検定の結果を出力する関数
        simeffects.out <- function(partresult, omit = FALSE, copy = FALSE){
            # 情報を展開
            part.info1 <- partresult$sig.term
            partstat <- partresult$sim.stat
            partepsi <- partresult$simepsi
            parttab <- partresult$simtab
            partmulttab <- partresult$sim.multresults
            
            #インターフェイスによってラインの長さを調整
            if(.Platform$GUI == "Rgui" | .Platform$GUI == "X11" | copy == TRUE){
                tablength = 6
                btablength = 4
                stablength = 4
            }else{
                tablength = 1
                btablength = 1
                stablength = 1
            }
            
            cat(sprintf(paste("< SIMPLE EFFECTS for <", part.info1, "> INTERACTION >"), sep = ""), sep = "", "\n", "\n")# タイトル
            
            # omitがFALSEなら記述統計量を出力する
            if(omit == FALSE){
                # 記述統計量を出力する
                tabline <- sapply(1:ncol(partstat), function(x) ifelse(is.numeric(partstat[,x]), max(nchar(round(partstat[,x], 4), type = "width")), max(nchar(as.character(partstat[,x]), type = "width")))) + 2
                mainline <- sum(tabline) + btablength * length(tabline)
                headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
                bodymargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
                bodymargin[(length(tabline)-1):length(tabline)] <- paste("%", tabline[(length(tabline)-1):length(tabline)], ".4f", sep = "")
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")# ラインを引く；要因の数に合わせて長さを調整
                cat(sprintf(headmargin, names(partstat)), sep = "\t", "\n")# データフレームの列名をタブ区切りで表示する
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する
                for (i in 1:nrow(partstat)){
                    cat(sprintf(bodymargin[1:3], paste(partstat[i, 1:3])), sprintf(bodymargin[4:5], partstat[i, 4:5]), sep = "\t", "\n")
                }
                
                cat(sprintf(rep("-", mainline + 2)), sep = "", "\n", "\n")
                
            }
            
            # 球面性検定の結果を出力
            if(is.na(charmatch("Error", parttab$source.col)) && !is.na(partepsi)){
                # 被験者内計画，混合要因計画なら球面性検定の結果を出力する
                tabline <- sapply(1:ncol(partepsi), function(x) ifelse(is.numeric(partepsi[,x]), max(nchar(round(partepsi[,x], 4), type = "width")), max(nchar(as.character(partepsi[,x]), type = "width")))) + 1
                tabline <- ifelse(nchar(names(partepsi)) > tabline, nchar(names(partepsi)) + 1, tabline)
                tabline[5] <- 7
                tabline[6] <- 3
                tabline[7:ncol(partepsi)] <- 6
                mainline <- sum(tabline) + stablength * length(tabline)
                headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
                bodymargin <- sapply(tabline, function(x) paste("%", x, ".4f", sep = ""))
                bodymargin[1] <- paste("%", tabline[1], "s", sep = "")
                bodymargin[4] <- paste("%", tabline[4], "s", sep = "")
                bodymargin[6] <- "%-3s"
                
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1:5], names(partepsi[1:5])), sprintf(headmargin[6], paste("")), sprintf(headmargin[7:ncol(partepsi)], names(partepsi[7:ncol(partepsi)])), sep = "\t", "\n")# 列名をタブ区切りで出力
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する
                for (i in 1:nrow(partepsi)){
                    cat(sprintf(bodymargin[1], paste(partepsi[i, 1])), sprintf(bodymargin[2], partepsi[i, 2]), sprintf(bodymargin[3], partepsi[i, 3]),
                    sprintf(bodymargin[4], partepsi[i, 4]), replace(sprintf(bodymargin[5], partepsi[i, 5]), is.na(partepsi[i, 5]), sprintf(paste("%", tabline[5], "s", sep = ""), "")),
                    replace(sprintf(bodymargin[6], paste(partepsi[i, 6])), partepsi[i, 6] == "", sprintf(bodymargin[6], "")), sprintf(bodymargin[7], partepsi[i, 7:ncol(partepsi)]), sep = "\t", "\n")
                }
                
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                if(ncol(partepsi) == 9){
                    cat(sprintf(paste("%", mainline - 1, "s", sep = ""), "LB = lower.bound, GG = Greenhouse-Geisser, HF = Huynh-Feldt"), sep = "\t", "\n", "\n", "\n")
                }else{
                    cat(sprintf(paste(rep("", 3), sep = "")), sprintf(paste("%", mainline - 1, "s", sep = ""), "b_hat = b^, c_hat = c^, h_d = h~', h_dd = h~'', h = h~"), sep = "\t", "\n", "\n", "\n")
                }
                
            }
            
            # 分散分析の結果を出力
            tabline <- sapply(1:ncol(parttab), function(x) ifelse(is.numeric(parttab[,x]), max(nchar(round(parttab[,x]), type = "width")) + 5, max(nchar(as.character(parttab[,x]), type = "width")))) + 2
            tabline[3] <- max(nchar(round(parttab[,3], 4), type = "width"))
            mainline <- sum(tabline) + tablength * length(tabline)
            headmargin <- sapply(tabline, function(x) paste("%", x, "s", sep = ""))
            bodymargin <- sapply(tabline, function(x) paste("%", x, ".4f", sep = ""))
            bodymargin[1] <- paste("%", tabline[1], "s", sep = "")
            bodymargin[3] <- paste("%", tabline[3], "s", sep = "")
            bodymargin[7] <- paste("%-", tabline[7], "s", sep = "")
            
            esn <- ncol(parttab) - 7
            if(esn == 0){# 分散分析表のみの場合
                # 単純主効果の検定の分散分析表を出力する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"),
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")), sep = "\t", "\n")# 列名をタブ区切りで表示する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する
                if(is.na(charmatch("Error", parttab$source.col))){
                    # 被験者内計画，混合要因計画なら２行ごとにラインを入れる
                    for (i in 1:nrow(parttab)){
                        cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]),
                        sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]),
                        replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""),
                        replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""),
                        sprintf(bodymargin[7], paste(parttab[i, 7])), sep = "\t", "\n")
                        if((i %% 2) == 0){# 行番号が２で割り切れるときはラインを引く
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                        }
                    }
                    # 表の末尾に有意性判定を示す記号を記す
                    cat(sprintf(paste("%", mainline - 2, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
                }else{
                    # 被験者間計画なら誤差項の上の行だけにラインを入れる
                    for (i in 1:nrow(parttab)){
                        if(nrow(parttab) == i){# 行番号が最後の行に一致するときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]),
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), sep = "\t", "\n")
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                            cat(sprintf(paste("%", mainline - 2, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
                        }else{# その他のときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]),
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]),
                            replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""),
                            replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""),
                            sprintf(bodymargin[7], paste(parttab[i, 7])), sep = "\t", "\n")
                        }
                    }
                }
            }else{# 効果量の指標を追加した場合
                # 単純主効果の検定の分散分析表を出力する
                cat(sprintf(rep("-", mainline)), sep = "", "\n")# ラインを引く
                cat(sprintf(headmargin[1], "Source"), sprintf(headmargin[2], paste("SS")), sprintf(headmargin[3], "df"), sprintf(headmargin[4], "MS"),
                sprintf(headmargin[5], "F-ratio"), sprintf(headmargin[6], paste("p-value")), sprintf(headmargin[7], paste("")), sprintf(headmargin[8:ncol(parttab)], names(parttab[8:ncol(parttab)])), sep = "\t", "\n")# 列名をタブ区切りで表示す
                cat(sprintf(rep("-", mainline)), sep = "", "\n")
                
                # データフレームのデータ部分を一行ずつタブ区切りで表示する
                if(is.na(charmatch("Error", parttab$source.col))){
                    # 被験者内計画，混合要因計画なら２行ごとにラインを入れる
                    for (i in 1:nrow(parttab)){
                        cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]),
                        sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]),
                        replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""),
                        replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""),
                        sprintf(bodymargin[7], paste(parttab[i, 7])), replace(sprintf(bodymargin[8:ncol(parttab)],
                        parttab[i, 8:ncol(parttab)]), is.na(parttab[i, 8:ncol(parttab)]), ""), sep = "\t", "\n")
                        if((i %% 2) == 0){# 行番号が２で割り切れるときはラインを引く
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                        }
                    }
                    # 表の末尾に有意性判定を示す記号を記す
                    cat(sprintf(paste("%", mainline - 2, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
                }else{
                    # 被験者間計画なら誤差項の上の行だけにラインを入れる
                    for (i in 1:nrow(parttab)){
                        if(nrow(parttab) == i){# 行番号が最後の行に一致するときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]),
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), sep = "\t", "\n")
                            cat(sprintf(rep("-", mainline)), sep = "", "\n")
                            cat(sprintf(paste("%", mainline - 2, "s", sep = ""), "+p < .10, *p < .05, **p < .01, ***p < .001"), sep = "\t", "\n", "\n")
                        }else{# その他のときは以下を実行
                            cat(sprintf(bodymargin[1], paste(parttab[i, 1])), sprintf(bodymargin[2], parttab[i, 2]), 
                            sprintf(bodymargin[3], round(parttab[i, 3], 2)), sprintf(bodymargin[4], parttab[i, 4]), 
                            replace(sprintf(bodymargin[5], parttab[i, 5]), is.na(parttab[i, 5]), ""), 
                            replace(sprintf(bodymargin[6], parttab[i, 6]), is.na(parttab[i, 6]), ""), 
                            sprintf(bodymargin[7], paste(parttab[i, 7])), replace(sprintf(bodymargin[8:ncol(parttab)], 
                            parttab[i, 8:ncol(parttab)]), is.na(parttab[i, 8:ncol(parttab)]), ""), sep = "\t", "\n")
                        }
                    }
                }
            }
            
            # 多重比較の結果を出力する
            for (i in 1:length(partmulttab)){
                if(is.na(partmulttab[[i]][1]) == FALSE){# リストに中身があったら出力する
                    cat("\n")# １行空ける
                    mod.Bon.out(partmulttab[[i]], omit = TRUE, copy = copy)
                }
            }
            
            cat("\n")# １行空ける
            
        }
        
        
        # 下位検定の結果を出力する関数
        post.out <- function(postresults, design, copy = FALSE){
            # リストに含まれる要素の数を特定
            postnum <- length(postresults)
            
            # すべてのリストの中身がNAの場合には，タイトルも表示しない
            if(sum(is.na(postresults)) != postnum){
                cat(sprintf("<< POST ANALYSES >>"), sep = "", "\n", "\n")# タイトル
            }
            
            # リストの中身をひとつずつeach.out関数に送る
            for (i in 1:postnum){
                each.out(postresults[[i]], design, copy = copy)
            }
            
            # 出力が終わったことを知らせるプロンプト
            cat(sprintf("output is over "), sprintf(rep("-", 20)), sprintf("///"), sep = "", "\n")
            cat("\n")# １行空ける
        }
        
        
        # 効果のタイプによって異なる出力方式を割り当てる関数
        each.out <- function(partresult, design, copy = FALSE){
            if(is.na(partresult[[1]]) == TRUE){
                # 何もない場合（２水準の主効果が有意で多重比較の必要なしの場合）
                
            }else if(nchar(partresult[[1]]) == 1){
                # 多重比較の結果がある場合
                mod.Bon.out(partresult, copy = copy)
                cat("\n")# １行空ける
                
            }else if(nchar(partresult[[1]]) == 3){
                # 一次の交互作用が見られた場合
                if(nchar(design) == 3) simeffects.out(partresult, omit = TRUE, copy)# もともと２要因の計画のときには，単純主効果の検定の出力時に記述統計量を出力しない（主分析の記述統計量と重複するので）
                else simeffects.out(partresult, copy = copy)
                
            }else if(nchar(partresult[[1]]) >= 5){
                # 高次の交互作用が見られた場合
                cat(sprintf(paste("< HIGHER-ORDER < ", partresult[[1]], " > INTERACTION >", sep = "")),sep = "", "\n")
                cat(sprintf("*** Split the dataset for further analysis. ***"), sep ="", "\n", "\n")# データ分割を促すプロンプト
            }
        }



        #--------------------------------------------------------------
        # ここから出力の指示

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
                    
                    lineplot.CI(x.factor = fac1, response = res, xlab="Level (Error bars show 95% CI)", ylab="", main="")
                    
                } else {  # Within = "sA"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=1:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, xlab="Level (Error bars show 95% CI)", ylab="", main="")
                }
                
            } else { # Two-way ANOVA
                
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    fac1 <- as.factor(dat[,1])
                    fac2 <- as.factor(dat[,2])
                    res <- dat[,3]
                    
                    lineplot.CI(fac1, res, group = fac2, cex = 2,
                    xlab = "Factor 1 Level (Error bars show 95% CI)", ylab = "", main="Factor 1")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    x <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(x[,1])
                    fac2 <- as.factor(x[,2])
                    res <- x[,3]
                    
                    lineplot.CI(fac1, res, group = fac2, cex = 2,
                    xlab = "Factor 1 Level (Error bars show 95% CI)", ylab = "", , main="Factor 1")
                    
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    x <- reshape(dat, varying=1:ncol(dat), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(x)/nlva, length.out = nrow(x))
                    fac1 <- anms
                    x <- cbind(x, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(dat), length.out = nrow(x))
                    fac2 <- bnms
                    x <- cbind(x, fac2)
                    
                    fac1 <- as.factor(x$fac1)
                    fac2 <- as.factor(x$fac2)
                    res <- x$res
                    
                    lineplot.CI(fac1, res, group = fac2, cex = 2,
                    xlab = "Factor 1 Level (Error bars show 95% CI)", ylab = "", main="Factor 1")
                    
                }
            }
            
        } else { # y-axis min-max
            
            if (input$factor == "oneway") {
                
                if (input$one.design == "Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, ylim=c(min(res), max(res)), xlab="Level (Error bars show 95% CI)", ylab="", main="")
                    
                } else {  # Within = "sA"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    dat <- reshape(dat, varying=1:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(dat[,1])
                    res <- dat[,2]
                    
                    lineplot.CI(x.factor = fac1, response = res, ylim=c(min(res), max(res)), xlab="Level (Error bars show 95% CI)", ylab="", main="")
                }
                
            } else { # Two-way ANOVA
                
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    fac1 <- as.factor(dat[,1])
                    fac2 <- as.factor(dat[,2])
                    res <- dat[,3]
                    
                    lineplot.CI(fac1, res, group = fac2, cex = 2, ylim=c(min(res), max(res)),
                    xlab = "Factor 1 Level (Error bars show 95% CI)", ylab = "", main="Factor 1")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    x <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(x[,1])
                    fac2 <- as.factor(x[,2])
                    res <- x[,3]
                    
                    lineplot.CI(fac1, res, group = fac2, cex = 2, ylim=c(min(res), max(res)),
                    xlab = "Factor 1 Level (Error bars show 95% CI)", ylab = "", , main="Factor 1")
                    
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    x <- reshape(dat, varying=1:ncol(dat), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(x)/nlva, length.out = nrow(x))
                    fac1 <- anms
                    x <- cbind(x, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(dat), length.out = nrow(x))
                    fac2 <- bnms
                    x <- cbind(x, fac2)
                    
                    fac1 <- as.factor(x$fac1)
                    fac2 <- as.factor(x$fac2)
                    res <- x$res
                    
                    lineplot.CI(fac1, res, group = fac2, cex = 2, ylim=c(min(res), max(res)),
                    xlab = "Factor 1 Level (Error bars show 95% CI)", ylab = "", main="Factor 1")
                    
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
                    
                    fac1 <- as.factor(dat[,1])
                    fac2 <- as.factor(dat[,2])
                    res <- dat[,3]
                    
                    lineplot.CI(fac2, res, group = fac1, cex = 2,
                    xlab = "Factor 2 Level (Error bars show 95% CI)", ylab = "", main="Factor 2")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    x <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(x[,1])
                    fac2 <- as.factor(x[,2])
                    res <- x[,3]
                    
                    lineplot.CI(fac2, res, group = fac1, cex = 2,
                    xlab = "Factor 2 Level (Error bars show 95% CI)", ylab = "", main="Factor 2")
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    x <- reshape(dat, varying=1:ncol(dat), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(x)/nlva, length.out = nrow(x))
                    fac1 <- anms
                    x <- cbind(x, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(dat), length.out = nrow(x))
                    fac2 <- bnms
                    x <- cbind(x, fac2)
                    
                    fac1 <- as.factor(x$fac1)
                    fac2 <- as.factor(x$fac2)
                    res <- x$res
                    
                    lineplot.CI(fac2, res, group = fac1, cex = 2,
                    xlab = "Factor 2 Level (Error bars show 95% CI)", ylab = "", main="Factor 2")
                    
                }
                
            } else {
                
                cat("\n")
                
            }
            
        } else {
            
            if (input$factor == "twoway") {
                if (input$two.design == "Factor1Between_Factor2Between") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    fac1 <- as.factor(dat[,1])
                    fac2 <- as.factor(dat[,2])
                    res <- dat[,3]
                    
                    lineplot.CI(fac2, res, group = fac1, cex = 2, ylim=c(min(res), max(res)),
                    xlab = "Factor 2 Level (Error bars show 95% CI)", ylab = "", main="Factor 2")
                    
                } else if (input$two.design == "Factor1Between_Factor2Within") {
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    x <- reshape(dat, varying=2:length(dat), v.names="score", direction = "long")
                    fac1 <- as.factor(x[,1])
                    fac2 <- as.factor(x[,2])
                    res <- x[,3]
                    
                    lineplot.CI(fac2, res, group = fac1, cex = 2, ylim=c(min(res), max(res)),
                    xlab = "Factor 2 Level (Error bars show 95% CI)", ylab = "", main="Factor 2")
                    
                } else { # input$two.design == "Factor1Within_Factor2Within"
                    dat <- read.csv(text=input$text, sep="\t")
                    
                    nlva <- input$factor1.level # factor A の水準
                    nlvb <- input$factor2.level # factor B の水準
                    
                    x <- reshape(dat, varying=1:ncol(dat), v.names="res", direction = "long")
                    
                    anms <- paste("a", 1:nlva, sep = "")
                    anms <- rep(anms, each = nrow(x)/nlva, length.out = nrow(x))
                    fac1 <- anms
                    x <- cbind(x, fac1)
                    
                    bnms <- paste("b", 1:nlvb, sep = "")
                    bnms <- rep(bnms, each = nrow(dat), length.out = nrow(x))
                    fac2 <- bnms
                    x <- cbind(x, fac2)
                    
                    fac1 <- as.factor(x$fac1)
                    fac2 <- as.factor(x$fac2)
                    res <- x$res
                    
                    lineplot.CI(fac2, res, group = fac1, cex = 2, ylim=c(min(res), max(res)),
                    xlab = "Factor 2 Level (Error bars show 95% CI)", ylab = "", main="Factor 2")
                    
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

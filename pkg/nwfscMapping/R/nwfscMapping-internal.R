.addBubbles.ach <-
function (events, type = c("perceptual", "surface", "volume"), 
    z.max = NULL, max.size = 0.8, symbol.zero = "+", symbol.fg = rgb(0, 
        0, 0, 0.6), symbol.bg = rgb(0, 0, 0, 0.3), legend.pos = "bottomleft", 
    legend.breaks = NULL, show.actual = FALSE, legend.type = c("nested", 
        "horiz", "vert"), legend.title = "Abundance", legend.cex = 0.8,legend.title.cex=legend.cex+0.2, 
    ...) 
{
    #A sligth modification from PBSmapping addBubbles to call a modified legend function (fixes an error)
    events <- .validateEventData(events)
    if (is.character(events)) 
        stop(paste("Invalid EventData 'events'.\n", events, sep = ""))
    if (!is.element("Z", names(events))) 
        stop("EventData is missing required column 'Z'.\n")
    type <- match.arg(type)
    if (!is.null(legend.pos)) 
        legend.type <- match.arg(legend.type)
    if (is.null(z.max) || is.na(z.max)) 
        z.max <- max(events$Z, na.rm = TRUE)
    if (is.null(legend.breaks) || is.na(legend.breaks)) 
        legend.breaks <- pretty(range(events$Z), 3)[-1]
    if (show.actual) 
        legend.breaks <- signif(legend.breaks/max(legend.breaks) * 
            max(events$Z, na.rm = TRUE), 3)
    usr.xdiff <- par("usr")[2] - par("usr")[1]
    usr.ydiff <- par("usr")[4] - par("usr")[3]
    stand.rad <- (max.size/2)/par("pin")[1] * usr.xdiff
    events <- events[order(events$Z, decreasing = TRUE), ]
    type <- match.arg(type)
    switch(type, volume = {
        radii <- ((events$Z/z.max)^(1/3)) * stand.rad
        radii.leg <- ((legend.breaks/z.max)^(1/3)) * stand.rad
    }, surface = {
        radii <- sqrt(events$Z/z.max) * stand.rad
        radii.leg <- sqrt(legend.breaks/z.max) * stand.rad
    }, perceptual = {
        radii <- ((events$Z/z.max)^0.57) * stand.rad
        radii.leg <- ((legend.breaks/z.max)^0.57) * stand.rad
    })
    isZero <- unlist(lapply(events$Z, all.equal, current = 0)) == 
        "TRUE"
    symbols(events$X[!isZero], events$Y[!isZero], circles = radii[!isZero], 
        inches = FALSE, bg = symbol.bg, fg = symbol.fg, add = TRUE)
    if (any(isZero) && (!is.logical(symbol.zero) || symbol.zero)) {
        if (is.logical(symbol.zero)) 
            symbol.zero <- "+"
        dots <- list(...)
        if (!is.null(dots$pch)) 
            stop("Specify 'pch' through 'symbol.zero'")
        points(events$X[isZero], events$Y[isZero], pch = symbol.zero, 
            ...)
    }
    if (!is.null(legend.pos)) {
        if (!any(isZero)) 
            symbol.zero <- FALSE
        .addBubblesLegend.ach(radii.leg, usr.xdiff, usr.ydiff, symbol.zero, 
            symbol.fg, symbol.bg, legend.pos, legend.breaks, 
            legend.type, legend.title, legend.cex,legend.title.cex, ...)
    }
    invisible()
}
.addBubblesLegend.ach <-
function (radii.leg, usr.xdiff, usr.ydiff, symbol.zero, symbol.fg, 
    symbol.bg, legend.pos, legend.breaks, legend.type, legend.title, 
    legend.cex,legend.title.cex=legend.cex+0.2, ...) 
{
    #fixes an error in the vertical legend
    ratio.y.x = (usr.ydiff/par("pin")[2])/(usr.xdiff/par("pin")[1])
    gap.x <- par("cxy")[1] * legend.cex/2
    gap.y <- par("cxy")[2] * legend.cex/2
    radii.leg.y <- radii.leg * ratio.y.x
    leg.tex.w <- strwidth(legend.breaks, units = "user") * legend.cex
    title.w = strwidth(legend.title)
    max.tex.w <- max(leg.tex.w)
    switch(legend.type, nested = {
        legend.height <- 2 * max(radii.leg.y) + 3 * gap.y
        legend.width <- 2 * max(radii.leg) + gap.x + max.tex.w
    }, horiz = {
        legend.height <- 2 * max(radii.leg.y) + 3 * gap.y
        legend.width <- 2 * sum(radii.leg) + (length(legend.breaks) - 
            1) * gap.x
    }, vert = {
        legend.height <- 2 * sum(radii.leg.y) + (length(legend.breaks) - 
            1) * gap.y + 3 * gap.y
        legend.width <- 2 * max(radii.leg) + gap.x + max.tex.w
    })
    if (title.w > legend.width) {
        w.adj <- (title.w - legend.width)/2
    }
    else {
        w.adj <- 0
    }
    if (class(legend.pos) == "numeric") {
        legend.loc <- legend.pos
    }
    else {
        corners <- c("bottomleft", "bottomright", "topleft", 
            "topright")
        if (legend.pos %in% corners) {
            legend.loc <- switch(legend.pos, bottomleft = c(par("usr")[1] + 
                0.025 * usr.xdiff + w.adj, par("usr")[3] + 0.025 * 
                usr.ydiff + legend.height), bottomright = c(par("usr")[2] - 
                (0.025 * usr.xdiff + legend.width + w.adj), par("usr")[3] + 
                0.025 * usr.ydiff + legend.height), topleft = c(par("usr")[1] + 
                0.025 * usr.xdiff + w.adj, par("usr")[4] - 0.025 * 
                usr.ydiff), topright = c(par("usr")[2] - (0.025 * 
                usr.xdiff + legend.width + w.adj), par("usr")[4] - 
                0.025 * usr.ydiff))
        }
    }
    switch(legend.type, nested = {
        legend.loc[1] <- legend.loc[1] + max(radii.leg)
        legend.loc[2] <- legend.loc[2] - legend.height
        r <- rev(radii.leg)
        bb <- rev(legend.breaks)
        x.text.leg <- legend.loc[1] + max(r) + gap.x + max.tex.w
        for (i in 1:length(r)) {
            symbols(legend.loc[1], legend.loc[2] + r[i] * ratio.y.x, 
                circles = r[i], inches = FALSE, add = TRUE, bg = symbol.bg, 
                fg = symbol.fg)
            lines(c(legend.loc[1], legend.loc[1] + r[1] + gap.x), 
                rep(legend.loc[2] + 2 * r[i] * ratio.y.x, 2))
            text(x.text.leg, legend.loc[2] + 2 * r[i] * ratio.y.x, 
                bb[i], adj = c(1, 0.5), cex = legend.cex)
        }
        x.title.leg <- legend.loc[1] - max(radii.leg) + (legend.width/2)
        text(x.title.leg, legend.loc[2] + legend.height, legend.title, 
            adj = c(0.5, 0.5), cex = legend.title.cex, col = "black")
        zlab <- c(x.title.leg, legend.loc[2] + legend.height/4)
    }, horiz = {
        legend.loc[2] <- legend.loc[2] + max(radii.leg.y) - legend.height
        offset <- vector()
        for (i in 1:length(radii.leg)) offset[i] <- 2 * sum(radii.leg[1:i]) - 
            radii.leg[i] + (i - 1) * gap.x
        symbols(legend.loc[1] + offset, rep(legend.loc[2], length(radii.leg)), 
            circles = radii.leg, inches = FALSE, bg = symbol.bg, 
            fg = symbol.fg, add = TRUE)
        text(legend.loc[1] + offset, legend.loc[2] + radii.leg.y + 
            gap.y, legend.breaks, adj = c(0.5, 0.5), cex = legend.cex)
        text(legend.loc[1] + legend.width/2, legend.loc[2] + 
            legend.height - max(radii.leg.y), legend.title, adj = c(0.5, 
            0.5), cex = legend.title.cex, col = "black")
        zlab <- c(legend.loc[1], legend.loc[2] - legend.height/8)
    }, vert = {
        if (any(legend.pos == c("bottomleft", "topleft"))) legend.loc[1] <- legend.loc[1] + 
            0.05 * usr.xdiff
        offset <- vector()
        for (i in 1:length(legend.breaks)) offset[i] <- gap.y + 
            2 * sum(radii.leg.y[1:i]) - radii.leg.y[i] + i * 
            gap.y
        symbols(rep(legend.loc[1], length(legend.breaks)), legend.loc[2] - 
            offset, circles = radii.leg, bg = symbol.bg, fg = symbol.fg, 
            inches = FALSE, add = TRUE)
        x.text.leg <- legend.loc[1] + max(radii.leg) + gap.x + 
            max.tex.w
        text(rep(x.text.leg, length(legend.breaks)), legend.loc[2] - 
            offset, legend.breaks, cex = legend.cex, adj = c(1,0.5), col = "black")
        text(legend.loc[1] + legend.width/2 - max(radii.leg), 
            legend.loc[2], legend.title, adj = c(0.5, 0.5), cex = legend.title.cex, col = "black")
        x.title.leg <- legend.loc[1] - max(radii.leg) + (legend.width/2)  #ach
        zlab <- c(x.title.leg, legend.loc[2])
    })
    if (!is.logical(symbol.zero)) 
        legend(zlab[1], zlab[2], legend = "zero", pch = symbol.zero, 
            xjust = 0, yjust = 1, bty = "n", cex = 0.8, x.intersp = 0.5)
    invisible()
}
.lat.to.km <-
function(lat){
# lat in degrees
    lat.rad <- (lat * pi)/180
    return(111.14 - 0.56 * cos(2 * lat.rad))
}
.long.to.km <-
function(lat){
# lat in degrees
    lat.rad <- (lat * pi)/180
    return(111.41 * cos(lat.rad) - 0.1 * cos(3 * lat.rad))
}


.Random.seed <-
c(403L, 10L, 2060855613L, -839351881L, -1855232418L, -124588136L, 
-748118277L, 461035161L, 788152104L, -857573930L, 354089553L, 
1845801219L, 1789619346L, 733045284L, 820468967L, -148369091L, 
1705887156L, 1215167002L, -196525723L, -313602049L, 930162310L, 
461708816L, 1839031475L, -1626841263L, 485954304L, -2142688802L, 
-553717207L, -775904805L, -262910230L, 2142944396L, -1953520689L, 
40278917L, 927930428L, 919001362L, -726397427L, 1263926599L, 
-2017306834L, 118860424L, 191650699L, 1449737065L, -733530824L, 
-51750938L, 1175891841L, -1329314029L, 1837224130L, -85660812L, 
1446467447L, -1543770003L, -634302716L, 1795083818L, -2010602603L, 
2132977711L, 982491446L, 1056866656L, 1329962915L, 323725505L, 
1474107568L, 1624664014L, 1658239801L, -462613685L, 308283642L, 
-28683780L, -550325697L, 2025723669L, -1094786004L, 1029119234L, 
-1307576867L, 833340631L, -1085919298L, -2052142600L, 1361178075L, 
-1217739719L, 853452168L, -1787121354L, 1507134321L, 82143651L, 
548836594L, 1881737284L, 2003249927L, 1768336093L, -434555948L, 
-919342342L, 1676256837L, 958594719L, -1946418138L, 73920688L, 
179464531L, 1135866993L, 622406048L, 1526916414L, -607214519L, 
-1107414789L, -1853425462L, -1824629204L, -1642520017L, -977069339L, 
-1603879012L, -1678599694L, 2076357741L, -2087053913L, 1795276750L, 
-587277016L, -517554389L, 220374345L, 589220440L, -737148154L, 
2037722849L, 1974758003L, -1023126814L, -186741164L, -2019813545L, 
1340009293L, -1039655196L, 4775626L, -1892704459L, -1013312113L, 
1918309014L, 665627200L, -185533693L, -416675039L, 684780560L, 
-968401938L, 2130332569L, -1561848789L, 26841370L, 503358428L, 
1357317599L, 1853649333L, -568025268L, -1885063774L, -1973873411L, 
1015249399L, -314231010L, -1867167272L, -431545541L, -1226357543L, 
1834485864L, -1587815018L, 1174324497L, -803264829L, 2098570962L, 
465805156L, -1889394009L, 748358653L, 865238516L, -526656806L, 
-2016420827L, 1625313599L, 1236289222L, 1168690000L, 163558259L, 
-780229231L, -942020672L, -1574268514L, 485514089L, 1443057563L, 
798701226L, 697991372L, -125959793L, -441821755L, 985921916L, 
1057409490L, 1435743565L, -926701561L, -165935122L, 556838600L, 
-1929723445L, 1001796649L, -1001381256L, 944836518L, -350877247L, 
-22102445L, 585815298L, -1000483788L, 1982277047L, -1551644243L, 
-1227617596L, -1315270934L, 1087643989L, -1327704337L, -1612421898L, 
-683363296L, 431269347L, -1041557887L, -461086480L, 1115594638L, 
337735929L, 671775371L, 666655418L, 1854182844L, -206274817L, 
1366815317L, -584104596L, 78253506L, -1878659555L, 3110551L, 
-1612726274L, -974804040L, -1578858597L, 51621369L, 94793032L, 
-1953191818L, 1261393329L, 161897443L, 1723558066L, 750060804L, 
-1238038969L, 1003835677L, 346775444L, -1666301126L, -1720021115L, 
1844845663L, 1231181542L, -1595581072L, 1663297171L, 486850353L, 
-336624160L, -1127030146L, -884944631L, 560952123L, 1897351946L, 
-170163476L, -127664273L, -463855451L, -2084140452L, -593810894L, 
-1084522195L, 2115241959L, -1801515250L, -246568620L, -940226798L, 
245409184L, -219182388L, 595805792L, 1691398082L, -1426767064L, 
-43541076L, -779022020L, 1200508786L, -502685392L, 1624391460L, 
-1364880072L, -522159398L, -666544240L, 1146335996L, 1832648596L, 
-1880426430L, -496259840L, 158567484L, -1281450672L, -1594464094L, 
-626882744L, -1119856884L, 489388764L, -420871342L, 280158848L, 
1751976020L, 306568168L, 1263441658L, -2063144624L, 1133372204L, 
885815156L, 140546962L, 1777694304L, 1321787276L, -936443040L, 
128002018L, 1316188072L, -451445556L, -922178052L, -114979918L, 
-1721800208L, -1506776636L, -1405062824L, -1656920710L, 748122544L, 
1888604092L, 514781268L, 623920066L, -527834080L, 364973244L, 
-4593072L, 1601834978L, -21547416L, 1070598732L, -488029956L, 
-279121934L, -151159552L, -815982732L, 1123704840L, 222533050L, 
971365840L, 849871852L, 164405140L, 1615195090L, 1164801888L, 
-618359476L, -803455712L, -810642558L, -1363170584L, -514820884L, 
-117753284L, 672783218L, -134288784L, 848594980L, -1077201096L, 
1138818202L, 1227502288L, -293412676L, -501340524L, 1358579330L, 
1861573824L, 1257017660L, -928522096L, 685398562L, -1277405688L, 
835539212L, 1248707228L, -1132634606L, 423947648L, -59052524L, 
1077569064L, -1413900998L, 943491024L, -1821851028L, -1221957836L, 
-824307054L, -385517216L, -577846964L, 213245920L, 390696866L, 
-100911192L, -1045775476L, -1173137732L, 651418482L, 621422832L, 
-279652284L, 882974040L, 227361850L, 1466114928L, 1199088700L, 
1696717588L, -1348761150L, -2049396512L, 392865660L, 1199107792L, 
-782743902L, -2061310936L, -413791988L, -1996217412L, -549262030L, 
-1858494528L, 365646388L, -160511032L, 1117209786L, 406220688L, 
-219363860L, 616284244L, 784436882L, 1238940448L, -238027316L, 
20229088L, -193998142L, -1186321880L, -889004372L, -312610116L, 
-1817058830L, -1101076176L, 1617206436L, 365187768L, 2145495898L, 
-560403056L, -1031609988L, 1327848724L, 1618594114L, -259882240L, 
1047747516L, -156497072L, -1117183710L, 1536429384L, -1168254964L, 
-1439909796L, 1060623570L, -1327527936L, 1602456404L, 1112257000L, 
-198105094L, -749869360L, 489690412L, -2134976268L, -1959438446L, 
407800928L, 1118211340L, 1463674592L, 425635170L, 918184L, 199273036L, 
579209084L, -1259123918L, -1275720080L, -383623868L, -1650185256L, 
-1295720326L, 1932514864L, -972170564L, -569346220L, -1116702142L, 
575986592L, -1728882500L, -1722519088L, -1328659102L, 1770349544L, 
2144536524L, -2144762116L, 1537617906L, -932938624L, 2098424052L, 
1332130952L, 866170170L, 1481909072L, 87698540L, 985895060L, 
188967250L, 1142980320L, -414465844L, 1954076832L, 1346673026L, 
-1034643224L, -1025242900L, -170555332L, 1940349042L, -2128099216L, 
1573577252L, 1718944312L, 1859131930L, -1464827696L, -1065843140L, 
1165044884L, 1387316354L, -839841856L, -1566409412L, -170922352L, 
224592930L, 1346695176L, -1347232628L, 880726044L, 1945195666L, 
-1337086848L, -1972419564L, 1505665832L, -1464822726L, 762547408L, 
-1491364756L, 141690548L, 700076050L, 614424544L, -555605392L, 
2050672665L, 520143371L, 25594108L, 2059580714L, -1146320033L, 
-1043211543L, 632406526L, -1350234644L, -2031657667L, -66542329L, 
-505109264L, -1817831282L, 1848144571L, -550417763L, -405082550L, 
844026984L, 116280401L, -1732176061L, -1340820924L, -1046109966L, 
1087196359L, -617310575L, 296069142L, 1458573940L, 877257989L, 
1922259791L, -198838968L, 1416718566L, -903922157L, -657505931L, 
-1699752622L, 625604000L, -1022712695L, 1751467195L, -1842825460L, 
-285303142L, -965565425L, 1915183193L, 1477364974L, -157251716L, 
1316991981L, -1999029097L, 1333217120L, -183902402L, 678690827L, 
852953229L, -794865222L, 2004942840L, -597509215L, -432818861L, 
-436163468L, -1219265022L, -1480568681L, -1130945951L, -832226650L, 
-576806236L, 2105730901L, -1405801473L, 595737240L, 1865260982L, 
-1174842301L, 556780933L, -2024221150L, 215010576L, -483418951L, 
698343787L, 1741131868L, -1680368246L, 1061448063L, -1130089271L, 
1025411422L, -1605162484L, -427422243L, -853551897L, -20975344L, 
-1222766546L, -701880613L, 1218305981L, 1397684522L, -180097592L, 
-851873039L, 530046051L, -1628936092L, -1799511150L, -1465767833L, 
-1475053519L, -15005322L, 496760468L, 27393189L, 1282302383L, 
1294858856L, 1420608070L, 1977116659L, 830936917L, 678299954L, 
-898261760L, 1746230761L, -1432950693L, 675297964L, 997074106L, 
-1408965009L, -1359252487L, 1262017166L, -1687381156L, -512272435L, 
1851234103L, 717931968L, 1282468510L, -168066197L, 579164909L, 
-1993781414L, 1926333464L, 164135297L, -1732774093L, -1896403628L, 
1091479394L, 813300727L, -1724349375L, 1177759942L, -2117698044L, 
106248629L, -2068986209L, -1518893320L, -990745898L, -2027520157L, 
1872512805L, 1273806018L, 1200077232L, -1270148647L, 317169483L, 
34854332L, 1905919466L, 631669279L, -1024907095L, 421584958L, 
-2050584276L, 1186225021L, -1541852729L, -1276615504L, 2041082574L, 
-698816517L, -633278371L, 1676285834L, 32806824L, -378036975L, 
-1557257213L, -876343548L, -1912167630L, 402014599L, 1295427537L, 
1859127638L, -781916364L, -1559508411L, -1232193905L, -231022584L, 
341722150L, 803658963L, -51626443L, -239340142L, -719124128L, 
916762825L, 932168699L, -2049137716L, 1188191834L, -1035713713L, 
913627673L, -453404626L, 1356912188L, 1351675763L)

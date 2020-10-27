library(ggplot2)

.th <- theme(
    axis.title.x=element_blank(), axis.title.y=element_blank(), 
    strip.text.x = element_text(size=12, face="bold"),
    strip.text.y = element_text(size=12, face="bold"),
    strip.background = element_rect(colour="red", fill="#CCCCFF"),
    legend.title=element_blank(), legend.position='bottom')

## cap the values
.cp <- function(dat, grp, val='val', cap=0.01, mtd=c('both', 'upper', 'lower'))
{
    grp <- split(dat, dat[, grp])
    mtd <- match.arg(mtd, c('both', 'upper', 'lower'))
    grp <- lapply(grp, function(g)
    {
        v <- g[, val]
        if(mtd == 'upper')
            v <- pmin(v, quantile(v, 1-cap, na.rm=TRUE))
        else if(mtd == 'lower')
            v <- pmax(v, quantile(v, 0+cap, na.rm=TRUE))
        else
        {
            v <- pmin(v, quantile(v, 1-cap/2, na.rm=TRUE))
            v <- pmax(v, quantile(v, 0+cap/2, na.rm=TRUE))
        }
        g[, val] <- v
        g
    })
    dat <- do.call(rbind, grp)
    dat
}

get.rpt <- function(sim, cache=TRUE)
{
    rds=paste0(sim, '.rds')
    if(file.exists(rds) && cache)
        agg <- readRDS(rds)
    else
    {
        agg <- lapply(dir(sim, '^[0-9]+.rds$', full=TRUE), function(f)
        {
            print(f)
            readRDS(f)
        })
        agg <- do.call(rbind, agg)
        saveRDS(agg, rds)
    }
    invisible(agg)
}

get.pow <- function(sim, cache=TRUE)
{
    rds=paste0(sim, '.pow')
    if(file.exists(rds) && cache)
        pow <- readRDS(rds)
    else
    {
        rpt <- get.rpt(sim)
        rpt <- subset(rpt, se=-seed)
        grp <- subset(rpt, se=c(key, tag, N, M, mtd))
        pow <- by(rpt, grp, function(g)
        {
            cfg <- subset(g, se=-c(pow, egv, rep))[1, ]
            pow <- with(g, sum(pow * rep) / sum(rep))
            egv <- with(g, sum(pow * egv) / sum(rep))
            rep <- with(g, sum(rep))
            cbind(cfg, pow=pow, egv=egv, rep=rep)
        })
        pow <- do.call(rbind, pow)
        saveRDS(pow, rds)
    }
    invisible(pow)
}

plt.pow <- function(sim, out=paste0(sim, '.pdf'))
{
    rpt <- get.pow(sim)
    
    ## rpt <- subset(rpt, mtd %in% c("tsq", "dot_ch2", "skt", "mra", "bon"))
    ## g <- ggplot(rpt, aes(x=N, y=pow))
    g <- ggplot(rpt, aes(x=M, y=pow))
    g <- g + geom_line(aes(color=mtd), alpha=.5, size=1)
    g <- g + facet_grid(key ~ tag)
    g <- g + .th

    nfy <- length(unique(rpt$key))
    ufy <- 4 # 10 / nfy
    nfx <- length(unique(rpt$tag))
    ufx <- 5 # 19 / nfx
    ## if(ufx / ufy < 19 / 10)
    ##     ufy <- ufx / 19 * 10
    ## else
    ##     ufx <- ufy / 10 * 19
    options(bitmapType = 'cairo', device = 'pdf')
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy), dpi=300, scale=.8)
    invisible(g)
}


plt.egv <- function(sim, out=paste0(sim, '_egv.pdf'))
{
    rpt <- get.pow(sim)
    
    ## rpt <- subset(rpt, mtd %in% c("tsq", "dot_ch2", "skt", "mra", "bon"))
    g <- ggplot(rpt, aes(x=M, y=egv))
    g <- g + geom_line(aes(color=mtd), alpha=.5, size=1)
    g <- g + facet_grid(key ~ tag)
    g <- g + .th

    nfy <- length(unique(rpt$key))
    ufy <- 4 # 10 / nfy
    nfx <- length(unique(rpt$tag))
    ufx <- 5 # 19 / nfx
    options(bitmapType = 'cairo', device = 'pdf')
    ggsave(out, g, width=min(19, ufx * nfx), height=min(10, ufy * nfy), dpi=300, scale=.8)
    invisible(g)
}

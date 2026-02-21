# Load external GMTs
list_gmt <- list(
    "h"      = "msigdb_v2023.2.Hs_GMTs/h.all.v2023.2.Hs.entrez.gmt",
    "c2_cgp" = "msigdb_v2023.2.Hs_GMTs/c2.cgp.v2023.2.Hs.entrez.gmt",
    "c2_cp"  = "msigdb_v2023.2.Hs_GMTs/c2.cp.v2023.2.Hs.entrez.gmt",
    "GOBP"   = "msigdb_v2023.2.Hs_GMTs/c5.go.bp.v2023.2.Hs.entrez.gmt",
    "GOCC"   = "msigdb_v2023.2.Hs_GMTs/c5.go.cc.v2023.2.Hs.entrez.gmt",
    "GOMF"   = "msigdb_v2023.2.Hs_GMTs/c5.go.mf.v2023.2.Hs.entrez.gmt"
) %>%
    map(\(filename) {
        here::here(
            "input",
            "genesets",
            filename
        )
    }) %>%
    map(\(x) getGmt(con = x))

## Define genes related to NF-kB signalling
geneids_nfkb <- c(
    list_gmt[["h"]][["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]]@geneIds
) %>%
    unique() %>%
    set_names(mapIds(org.Hs.eg.db,
        keys = .,
        column = "SYMBOL",
        keytype = "ENTREZID"
    ))

geneids_tgfb <- c(
    list_gmt[["h"]][["HALLMARK_TGF_BETA_SIGNALING"]]@geneIds
) %>%
    unique() %>%
    set_names(mapIds(org.Hs.eg.db,
        keys = .,
        column = "SYMBOL",
        keytype = "ENTREZID"
    ))

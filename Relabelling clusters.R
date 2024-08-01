##For lavelling clusters
OM_new.cluster.ids <- c("T cells", "PreAd2", "T cells", "MSC2s", "APCs", "PreAd1", "MSC1s", "Myeloid cells", "Myeloid cells", "NK cells", "PreAd4", "ECs", "SMCs", "B cells", "PreAd3", "MSC3s", "?")

names(OM_new.cluster.ids) <- levels(OM_merge)

OM_merge <- RenameIdents(OM_merge, OM_new.cluster.ids)

##Re-order clusters
levels(x = myobject) <- c("10", "3", "1", "4", "7", "11", "5", "8", "6", "9", "2", "12")
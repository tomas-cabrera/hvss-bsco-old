import mywheels.gccatalogs.aggregate as ag

mwcat = ag.GCCatalog(
    ["/hildafs/projects/phy200025p/tcabrera/hvss/holger_baumgardt_clean.txt",
    "/hildafs/projects/phy200025p/tcabrera/hvss/harris2010_II_clean.txt"],
)
print(mwcat.df)
mwcat.match_to_cmc_models(
    "/hildafs/projects/phy200025p/share/catalog_files",
)

expect_no_error(run_app_ABI())
# expect_warning(ABI_Helminth(),"")

# distance == 0
expect_warning(ABI_Helminth(),"• Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker\n• Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation\n• Suggest nematode 16S primer from nematode systematics paper")


# Genus-Species
expect_warning(ABI_Helminth(0.21,"NAS","16S rRNA"),"• Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker \n• Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation \n• Suggest nematode 16S primer from nematode systematics paper")
expect_warning(ABI_Helminth(0.141,"NS","COI"),"• Suggest to use mt 16S rRNA gene or mt COII as an alternative genetic marker \n• Although 18S has the smallest gap, but the low sequence variation at the genus-species level is challenging for species delimitation \n• Suggest nematode 16S primer from nematode systematics paper")
expect_warning(ABI_Helminth(0.07,"NT","18S rRNA"),"• Suggest to use mt 12S")
expect_warning(ABI_Helminth(0.16,"TR","12S rRNA"),"• Suggest to use mt 16S \n •Although 18S has small gap ")
expect_warning(ABI_Helminth(0.111,"CE","28S rRNA"),"• Suggest to use cytB or 12S ")


# Family-Genus
expect_warning(ABI_Helminth(0.047,"NAS","18S rRNA"),"• Suggest to use mt 12S rRNA gene as an alternative genetic marker, with primer from nematode systematics paper")
expect_warning(ABI_Helminth(0.167,"NS","16S rRNA"),"• Suggest to use mt COI or mt 12S or mt 16S\n•	But caution for COI primers ")
expect_warning(ABI_Helminth(0.172,"TR","ITS2"),"• Suggest 28S, but need to caution ")
expect_warning(ABI_Helminth(0.283,"CE","COII"),"• Suggest COI or 16S")

expect_warning(ABI_Helminth(0.1,"NT","28S rRNA"),"Out Of Bounds")

expect_no_warning(ABI_Helminth(0.01,"NAS","28S rRNA"))



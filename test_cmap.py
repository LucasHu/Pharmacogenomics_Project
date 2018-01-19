from cmapPy.pandasGEXpress import parse

my_col_metadata = parse("GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx", col_meta_only=True)
my_col_metadata[0:5]
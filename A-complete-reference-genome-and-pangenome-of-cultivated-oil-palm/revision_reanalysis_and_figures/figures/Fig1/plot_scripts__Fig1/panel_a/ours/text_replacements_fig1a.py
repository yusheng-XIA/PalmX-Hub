# Fig. 1a shares, scheme A (nine major oil crops, symmetric). FAOSTAT QCL World 2023, accessed 24 September 2026.
# (label, find, new) -- find strings taken from fix/text/base_paragraphs.txt, unique in their paragraph.
EDITS_FIG1A = [
    ("INTRO [17] · oil-palm shares (FAOSTAT 2023, nine major oil crops)",
     "supplies approximately 39% of global vegetable oil while occupying only 8.5% of oil-crop land area",
     "supplies approximately 40% of the vegetable oil produced by the nine major oil crops while occupying only "
     "8.6% of their combined harvested area"),
    ("LEGEND 1a [858] · contribution bars",
     "The right-hand bars show the contributions of oil palm to the harvested area of the selected major crops "
     "(8.6%) and the corresponding vegetable-oil output (39.5%) in 2023.",
     "The right-hand bars show the contributions of oil palm to the combined harvested area of nine major oil "
     "crops (8.6%) and to their combined vegetable-oil output (40.1%; palm oil plus palm kernel oil) in 2023, "
     "calculated from FAOSTAT (see Methods)."),
]
# Methods: replace this sentence inside the new text of REPLACE_PARAGRAPHS entry
# "MUST 11b · Methods Fig. 1a paragraph (whole paragraph)" in deliver_build/build_main_text.py
METHODS_OLD = ("The contributions of oil palm to the total harvested area of the selected major crops (8.6%) and to the "
               "corresponding vegetable-oil output (39.5%) in 2023 were calculated from FAOSTAT^{2}.")
METHODS_NEW = ("The contributions of oil palm to the combined harvested area of nine major oil crops (8.6%) and to "
               "their combined vegetable-oil output (40.1%) in 2023 were calculated from the FAOSTAT Crops and "
               "livestock products domain (World; accessed 24 September 2026)^{2}. The nine crops were oil palm, "
               "soybean, rapeseed, sunflower, groundnut, cotton, coconut, olive and sesame; for each crop, harvested "
               "area of the primary crop (oil palm fruit, soya beans, rape or colza seed, sunflower seed, groundnuts, "
               "seed cotton, coconuts, olives and sesame seed) was paired with production of the corresponding crude "
               "oil, and oil palm output was taken as palm oil plus palm kernel oil.")
